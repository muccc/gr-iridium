/* -*- c++ -*- */
/*
 * Copyright 2022 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H
#define INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H

#include <iridium/fft_channelizer.h>

namespace gr {
namespace iridium {

class fft_channelizer_impl : public fft_channelizer
{
private:
    // Nothing to declare in this block.

public:
    fft_channelizer_impl(int fft_size, int decimation);
    ~fft_channelizer_impl();

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H */
