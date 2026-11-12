/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "frame_sorter_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace iridium {

frame_sorter::sptr frame_sorter::make()
{
    return gnuradio::make_block_sptr<frame_sorter_impl>();
}


/*
 * The private constructor
 */
frame_sorter_impl::frame_sorter_impl()
    : gr::block("frame_sorter",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0))
{
    message_port_register_out(pmt::mp("pdus"));
    d_upstream_done.install(this);
    auto port_name = pmt::mp("pdus");
    message_port_register_in(port_name);
    set_msg_handler(port_name, [this](const pmt::pmt_t& msg) { this->handler(msg); });
}

/*
 * Our virtual destructor.
 */
frame_sorter_impl::~frame_sorter_impl() {}

void frame_sorter_impl::handler(const pmt::pmt_t& msg)
{
    pmt::pmt_t symbols = pmt::cdr(msg);
    pmt::pmt_t meta = pmt::car(msg);

    frame f;

    f.timestamp = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    f.center_frequency =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));

    int confidence =
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("confidence"), pmt::PMT_NIL));

    // A TCH_C frame with a valid CRC proves the burst is TCH_C: legacy frames at the
    // same time and frequency are the legacy demodulator misreading it and are dropped.
    // Otherwise legacy and TCH_C frames are never duplicates of each other.
    // Among duplicates of the same kind a frame with a valid CRC beats one without.
    bool tch_c = pmt::dict_has_key(meta, pmt::mp("tch_type"));
    bool crc_ok = pmt::to_bool(pmt::dict_ref(meta, pmt::mp("crc_ok"), pmt::PMT_F));

    auto it = frames.begin();

    while (it != frames.end()) {
        if (it->first.timestamp < f.timestamp - 2000000000) {
            message_port_pub(pmt::mp("pdus"), it->second);
            it = frames.erase(it);
        } else {
            break;
        }
    }

    // Frames within this window will be de-duplicated
    const int64_t max_dt = 1000000; // 1 ms
    const int64_t max_df = 10000;   // 10 kHz

    // Get first element which is after the start of our search window (time).
    it = frames.upper_bound((frame){ .timestamp = f.timestamp - max_dt });

    // Loop over all elements possibly within our search window (time and frequency)
    while (it != frames.end() && abs(it->first.timestamp - f.timestamp) < max_dt) {
        if (std::abs(it->first.center_frequency - f.center_frequency) < max_df) {
            pmt::pmt_t meta_old = pmt::car(it->second);
            bool tch_c_old = pmt::dict_has_key(meta_old, pmt::mp("tch_type"));
            bool crc_ok_old =
                pmt::to_bool(pmt::dict_ref(meta_old, pmt::mp("crc_ok"), pmt::PMT_F));
            if (tch_c_old != tch_c) {
                if (tch_c && crc_ok) {
                    it = frames.erase(it);
                    continue;
                }
                if (tch_c_old && crc_ok_old) {
                    return;
                }
                it = std::next(it);
                continue;
            }
            int confidence_old = pmt::to_long(
                pmt::dict_ref(meta_old, pmt::mp("confidence"), pmt::PMT_NIL));

            if (std::tie(crc_ok, confidence) > std::tie(crc_ok_old, confidence_old)) {
                // insert_or_assign with it2 as a hint would be nice but is C++17
                frames.erase(it);
                frames.insert({ f, msg });
            }
            return;
        }
        it = std::next(it);
    }
    frames.insert({ f, msg });
}

bool frame_sorter_impl::stop()
{
    // Flush remaining messages
    auto it = frames.begin();
    while (it != frames.end()) {
        message_port_pub(pmt::mp("pdus"), it->second);
        it = frames.erase(it);
    }

    // Signal end of messages
    message_port_pub(pmt::mp("pdus"), pmt::PMT_EOF);
    return true;
}

} /* namespace iridium */
} /* namespace gr */
