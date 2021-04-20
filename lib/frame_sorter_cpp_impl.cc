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
    set_msg_handler(port_name, [this](pmt::pmt_t msg) { this->handler(msg); });
}

/*
 * Our virtual destructor.
 */
frame_sorter_cpp_impl::~frame_sorter_cpp_impl() {}

void frame_sorter_cpp_impl::handler(pmt::pmt_t msg)
{
}

int frame_sorter_cpp_impl::work(int noutput_items,
                                gr_vector_const_void_star& input_items,
                                gr_vector_void_star& output_items)
{
    return 0;
}

} /* namespace iridium */
} /* namespace gr */
