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

    struct burst_data {
        uint64_t id;
        double offset;
        float magnitude;
        float relative_frequency;
        double center_frequency;
        float sample_rate;
        uint64_t timestamp;
        float noise;
        size_t len;
        gr_complex * data;
    };


class fft_channelizer_impl : public fft_channelizer
{
private:
      const int d_fft_size;
      const int d_ifft_size;
      const int d_inverse_overlap;
      const int d_output_step;
      const int d_decimation;
      const int d_channels;
      const int d_pdu_ports;
      const int d_max_burst_size;

      uint64_t d_channel_active;
      fftwf_complex * d_out;
      std::map<fftwf_complex*, fftwf_plan> d_fft_plans;

      gr::fft::fft_complex_fwd          d_fft;
      gr::fft::fft_complex_rev          d_ifft;

      float channel_center(int channel);
      float channel_lower_border(int channel);
      float channel_upper_border(int channel);
      uint64_t activated_channels(tag_t &new_burst);

      std::vector<std::map<uint64_t, burst_data>> d_bursts;
      burst_data create_burst(tag_t &new_burst, int channel);
      void append_to_burst(burst_data &burst, const gr_complex * data, size_t n);
      void publish_burst(burst_data &burst);

      int d_outstanding;
      void burst_handled(pmt::pmt_t msg);
public:
    fft_channelizer_impl(int fft_size, int decimation, bool activate_streams, int pdu_ports);
    ~fft_channelizer_impl();

    uint64_t get_n_dropped_bursts();
    int get_output_queue_size();
    int get_output_max_queue_size();


    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FFT_CHANNELIZER_IMPL_H */
