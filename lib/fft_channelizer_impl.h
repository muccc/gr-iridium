/* -*- c++ -*- */
/*
 * Copyright 2022 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H
#define INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H

#include <iridium/fft_channelizer.h>
#include <gnuradio/fft/fft.h>
#include <fftw3.h>

namespace gr {
namespace iridium {

class fft_channelizer_impl : public fft_channelizer
{
private:
      const int d_fft_size;
      const int d_ifft_size;
      const int d_inverse_overlap;
      const int d_output_step;
      const int d_decimation;
      fftwf_complex * d_out;
      std::map<fftwf_complex*, fftwf_plan> d_fft_plans;

      gr::fft::fft_complex_fwd          d_fft;
      gr::fft::fft_complex_rev          d_ifft;

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
