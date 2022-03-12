/* -*- c++ -*- */
/*
 * Copyright 2022 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "fft_channelizer_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace iridium {

#pragma message("set the following appropriately and remove this warning")
using input_type = float;
#pragma message("set the following appropriately and remove this warning")
using output_type = float;
fft_channelizer::sptr fft_channelizer::make(int fft_size, int decimation)
{
    return gnuradio::make_block_sptr<fft_channelizer_impl>(fft_size, decimation);
}


/*
 * The private constructor
 */
fft_channelizer_impl::fft_channelizer_impl(int fft_size, int decimation)
    : gr::sync_decimator("fft_channelizer",
                         gr::io_signature::make(
                             1 /* min inputs */, 1 /* max inputs */, sizeof(input_type)),
                         gr::io_signature::make(1 /* min outputs */,
                                                1 /*max outputs */,
                                                sizeof(output_type)),
                         2 /*<+decimation+>*/)
{
}

/*
 * Our virtual destructor.
 */
fft_channelizer_impl::~fft_channelizer_impl() {}

int fft_channelizer_impl::work(int noutput_items,
                               gr_vector_const_void_star& input_items,
                               gr_vector_void_star& output_items)
{
    const input_type* in = reinterpret_cast<const input_type*>(input_items[0]);
    output_type* out = reinterpret_cast<output_type*>(output_items[0]);

#pragma message("Implement the signal processing in your block and remove this warning")
    // Do <+signal processing+>

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace iridium */
} /* namespace gr */
