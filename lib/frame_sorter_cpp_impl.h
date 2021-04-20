/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_FRAME_SORTER_CPP_IMPL_H
#define INCLUDED_IRIDIUM_FRAME_SORTER_CPP_IMPL_H

#include <iridium/frame_sorter_cpp.h>

namespace gr {
namespace iridium {

class frame_sorter_cpp_impl : public frame_sorter_cpp
{
private:

    void handler(const pmt::pmt_t& msg);

public:
    frame_sorter_cpp_impl();
    ~frame_sorter_cpp_impl();

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FRAME_SORTER_CPP_IMPL_H */
