/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "frame_sorter_cpp_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace iridium {

frame_sorter_cpp::sptr frame_sorter_cpp::make()
{
    return gnuradio::make_block_sptr<frame_sorter_cpp_impl>();
}


/*
 * The private constructor
 */
frame_sorter_cpp_impl::frame_sorter_cpp_impl()
    : gr::sync_block("frame_sorter_cpp",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0))
{
    message_port_register_out(pmt::mp("pdus"));
    auto port_name = pmt::mp("pdus");
    message_port_register_in(port_name);
    set_msg_handler(port_name, [this](const pmt::pmt_t& msg) { this->handler(msg); });
}

/*
 * Our virtual destructor.
 */
frame_sorter_cpp_impl::~frame_sorter_cpp_impl() {}

void frame_sorter_cpp_impl::handler(const pmt::pmt_t& msg)
{
    pmt::pmt_t symbols = pmt::cdr(msg);
    pmt::pmt_t meta = pmt::car(msg);

    frame f;

    f.timestamp = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    f.center_frequency = pmt::to_float(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));


    auto it = frames.begin();

    while(it != frames.end()) {
        if(it->first.timestamp < f.timestamp - 1000000000) {
            message_port_pub(pmt::mp("pdus"), it->second);
            it = frames.erase(it);
        } else {
            break;
        }
    }


    const auto [it2, success] = frames.insert({f, msg});
    if(!success) {
        int confidence = pmt::to_long(pmt::dict_ref(meta, pmt::mp("confidence"), pmt::PMT_NIL));

        pmt::pmt_t meta_old = pmt::car(it2->second);
        int confidence_old = pmt::to_long(pmt::dict_ref(meta_old, pmt::mp("confidence"), pmt::PMT_NIL));

        if(confidence > confidence_old) {
            // insert_or_assign with it2 as a hint would be nice but is C++17
            frames.erase(it2);
            frames.insert({f, msg});
        }
    }
}

bool frame_sorter_cpp_impl::stop()
{
    // Flush remaining messages
    auto it = frames.begin();
    while(it != frames.end()) {
        message_port_pub(pmt::mp("pdus"), it->second);
        it = frames.erase(it);
    }

    // Signal end of messages
    message_port_pub(pmt::mp("pdus"), pmt::PMT_EOF);
    return true;
}

int frame_sorter_cpp_impl::work(int noutput_items,
                                gr_vector_const_void_star& input_items,
                                gr_vector_void_star& output_items)
{
    return 0;
}

} /* namespace iridium */
} /* namespace gr */
