/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "iridium_frame_printer_impl.h"
#include <gnuradio/io_signature.h>

#include "bits.pb.h"

#include "iridium.h"

#include <boost/format.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include <zmq.hpp>

#if defined(CPPZMQ_VERSION) && defined(ZMQ_MAKE_VERSION) && \
    CPPZMQ_VERSION >= ZMQ_MAKE_VERSION(4, 3, 1)
#define USE_NEW_CPPZMQ_SEND_RECV 1
#else
#define USE_NEW_CPPZMQ_SEND_RECV 0
#endif

#if defined(CPPZMQ_VERSION) && defined(ZMQ_MAKE_VERSION) && \
    CPPZMQ_VERSION >= ZMQ_MAKE_VERSION(4, 8, 0)
#define USE_NEW_CPPZMQ_SET_GET 1
#else
#define USE_NEW_CPPZMQ_SET_GET 0
#endif

using boost::format;

namespace gr {
namespace iridium {

iridium_frame_printer::sptr iridium_frame_printer::make(std::string file_info)
{
    return gnuradio::make_block_sptr<iridium_frame_printer_impl>(file_info);
}


/*
 * The private constructor
 */
iridium_frame_printer_impl::iridium_frame_printer_impl(std::string file_info)
    : gr::block("iridium_frame_printer",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0)),
      d_file_info(file_info),
      d_t0(0),
      d_context(1),
      d_socket(d_context, ZMQ_PUB)
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    bool bind = true;
    auto address = "tcp://127.0.0.1:2323";

    int time = 0;
#if USE_NEW_CPPZMQ_SET_GET
    d_socket.set(zmq::sockopt::linger, time);
#else
    d_socket.setsockopt(ZMQ_LINGER, &time, sizeof(time));
#endif
    if (bind) {
        d_socket.bind(address);
    } else {
        d_socket.connect(address);
    }

    int hwm = 1000;
    /* Set high watermark */
    if (hwm >= 0) {
#if USE_NEW_CPPZMQ_SET_GET
        d_socket.set(zmq::sockopt::sndhwm, hwm);
#else
#ifdef ZMQ_SNDHWM
        d_socket.setsockopt(ZMQ_SNDHWM, &hwm, sizeof(hwm));
#else // major < 3
        uint64_t tmp = hwm;
        d_socket.setsockopt(ZMQ_HWM, &tmp, sizeof(tmp));
#endif
#endif
    }


    auto port_name = pmt::mp("pdus");
    message_port_register_in(port_name);
    set_msg_handler(port_name, [this](const pmt::pmt_t& msg) { this->handler(msg); });

    // hack: hijack system message port
    set_msg_handler(pmt::mp("system"),
                    [this](const pmt::pmt_t& msg) { this->handle_msg_sys(msg); });
}

/*
 * Our virtual destructor.
 */
iridium_frame_printer_impl::~iridium_frame_printer_impl() {}

void iridium_frame_printer_impl::handle_msg_sys(const pmt::pmt_t& msg)
{
    // ignore system message(s)
}

void iridium_frame_printer_impl::handler(const pmt::pmt_t& msg)
{
    if (msg == pmt::PMT_EOF) {
        // pass EOF to (original) system handler
        system_handler(pmt::cons(pmt::mp("done"), pmt::mp(1)));
        return;
    }

    pmt::pmt_t symbols = pmt::cdr(msg);
    pmt::pmt_t meta = pmt::car(msg);

    uint64_t timestamp =
        pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    double center_frequency =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
    float magnitude =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("magnitude"), pmt::PMT_NIL));
    float noise = pmt::to_float(pmt::dict_ref(meta, pmt::mp("noise"), pmt::PMT_NIL));
    uint64_t id = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("id"), pmt::PMT_NIL));
    int confidence =
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("confidence"), pmt::PMT_NIL));
    float level = pmt::to_float(pmt::dict_ref(meta, pmt::mp("level"), pmt::PMT_NIL));
    int n_symbols = pmt::to_long(pmt::dict_ref(meta, pmt::mp("n_symbols"), pmt::PMT_NIL));

    if (d_file_info == "") {
        d_t0 = (timestamp / 1000000000ULL) * 1000000000ULL;
        d_file_info = "i-" + std::to_string(d_t0 / 1000000000ULL) + "-t1";
    }

    std::cout << "RAW: " << d_file_info << " ";
    std::cout << format("%012.4f ") % ((timestamp - d_t0) / 1000000.);
    std::cout << format("%010d ") % int(center_frequency + 0.5);
    std::cout << format("N:%05.2f%+06.2f ") % magnitude % noise;
    std::cout << format("I:%011d ") % id;
    std::cout << format("%3d%% ") % confidence;
    std::cout << format("%.5f ") % level;
    std::cout << format("%3d ") % (n_symbols - ::iridium::UW_LENGTH);


    auto bits = pmt::u8vector_elements(symbols);
    for (const auto i : bits) {
        std::cout << std::to_string(i);
    }

    std::cout << std::endl;

    Frame frame;
    frame.set_fileinfo(d_file_info);
    frame.set_timestamp(timestamp - d_t0);
    frame.set_frequency(int(center_frequency + 0.5));
    frame.set_magnitude(magnitude);
    frame.set_noise(noise);
    frame.set_id(id);
    frame.set_confidence(confidence);
    frame.set_level(level);
    frame.set_len(n_symbols * 2);

    for (const auto i : bits) {
        std::cout << std::to_string(i);
    }

    char packed[(n_symbols * 2 + 7) / 8];
    memset(packed, 0, sizeof(packed));

    for(int i=0; i<n_symbols * 2; i++) {
        if(bits[i]) {
            packed[i/8] |= 1 << (7 - (i % 8));
        }
    }

    frame.set_data(std::string(packed, sizeof(packed)));

    std::string msg_str;
    frame.SerializeToString(&msg_str);
    zmq::message_t zmsg(msg_str.size());
    memcpy(zmsg.data(), msg_str.c_str(), msg_str.size());

#if USE_NEW_CPPZMQ_SEND_RECV
    d_socket.send(zmsg, zmq::send_flags::none);
#else
    d_socket.send(zmsg);
#endif

}

} /* namespace iridium */
} /* namespace gr */
