/* -*- c++ -*- */
/*
 * Copyright 2020 Free Software Foundation, Inc.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_IRIDIUM_BURST_DOWNMIX_IMPL_H
#define INCLUDED_IRIDIUM_BURST_DOWNMIX_IMPL_H

#include <gnuradio/blocks/rotator.h>
#include <gnuradio/filter/fir_filter.h>

#include "tch_c_format.h"
#include <iridium/burst_downmix.h>
#include <memory>

namespace gr {
namespace iridium {

class burst_downmix_impl : public burst_downmix
{
private:
    static const int RRC_TAPS = 51;
    static const int RC_TAPS = 51;

    const tch_c_format* d_tch_c; // TCH_C burst format, nullptr for legacy bursts
    int d_symbols_per_second;
    int d_output_sample_rate;
    int d_output_samples_per_symbol;
    float d_samples_per_symbol; // at the actual burst sample rate
    float d_burst_sample_rate;
    size_t d_max_burst_size;
    int d_search_depth;
    int d_pre_start_samples;
    int d_cfo_est_fft_size;
    int d_fft_over_size_facor;
    int d_corr_fft_size;
    int d_sync_search_len;
    int d_hard_max_queue_len;
    uint64_t d_n_dropped_bursts;
    bool d_handle_multiple_frames_per_burst;
    bool d_debug;
    int64_t d_debug_id;

    gr_complex* d_input;
    gr_complex* d_coarse;
    gr_complex* d_frame;
    gr_complex* d_tmp_a;
    gr_complex* d_tmp_b;
    gr_complex* d_dl_preamble_reversed_conj_fft;
    gr_complex* d_ul_preamble_reversed_conj_fft;

    float* d_magnitude_f;
    float* d_magnitude_filtered_f;
    float* d_cfo_est_window_f;

    gr::fft::fft_complex_fwd* d_corr_fft;
    gr::fft::fft_complex_rev* d_corr_dl_ifft;
    gr::fft::fft_complex_rev* d_corr_ul_ifft;

    filter::kernel::fir_filter_ccf d_input_fir;
    filter::kernel::fir_filter_fff d_start_finder_fir;
    filter::kernel::fir_filter_ccf d_rrc_fir;
    filter::kernel::fir_filter_fff d_start_finder_coarse_fir;
    filter::kernel::fir_filter_ccf d_rrc_coarse_fir;

    volk::vector<gr_complex> d_dl_preamble_reversed_conj;
    volk::vector<gr_complex> d_ul_preamble_reversed_conj;

    blocks::rotator d_r;

    // The input filter runs on demand, see filter_input()
    blocks::rotator d_input_r;
    const gr_complex* d_burst; // the burst being handled
    int d_decimation;
    size_t d_frame_size;  // its size after decimation
    size_t d_frame_valid; // samples of d_frame filtered so far
    size_t d_input_valid; // samples of d_input shifted so far
    size_t d_coarse_valid; // samples of d_coarse filtered so far
    int d_coarse_step;
    void shift_input(size_t n);
    void filter_input(size_t from, size_t to);
    size_t filter_input_coarse(size_t n);
    void initialize_coarse_filters(const std::vector<float>& start_finder_taps);

    // Normalized correlation a TCH_C UW needs to be taken for one
    static constexpr float UW_MIN_CORRELATION = 0.25f;
    float d_sync_word_energy;
    gr::fft::fft_complex_fwd d_cfo_est_fft;

    // Coarse UW search over a whole tracked burst window (TCH_C only)
    int d_track_fft_size;
    int d_track_decimation;
    bool d_tracked; // the burst being handled is a window on a tracked channel
    std::unique_ptr<gr::fft::fft_complex_fwd> d_track_fwd;
    std::unique_ptr<gr::fft::fft_complex_rev> d_track_inv;
    std::vector<gr_complex> d_track_dl_fft;
    std::vector<gr_complex> d_track_ul_fft;
    std::vector<gr_complex> d_track_buf;
    int find_uw_coarse(size_t burst_size);

    bool d_wait_for_demod;
    uint64_t d_burst_id; // the burst being handled
    int d_n_frames;      // frames published for it so far
    void burst_done(uint64_t id);

    void handler(pmt::pmt_t msg);
    int find_start(size_t burst_size, uint64_t id);
    int process_next_frame(float sample_rate,
                           double center_frequency,
                           uint64_t timestamp,
                           uint64_t sub_id,
                           size_t burst_size,
                           int start,
                           float noise,
                           float magnitude);

    void update_buffer_sizes(size_t burst_size);
    void initialize_cfo_est_fft(void);
    void initialize_correlation_filter(void);
    void free_correlation_filter(void);
    void set_burst_sample_rate(float sample_rate);
    int pre_start_samples() const;
    int sync_search_len() const;
    double rrc_alpha() const;
    int min_width_class() const;
    volk::vector<gr_complex> generate_sync_word(::iridium::direction direction);
    int fft_shift_index(int index, int fft_size);
    int fft_unshift_index(int index, int fft_size);
    float interpolate(float alpha, float beta, float gamma);

public:
    burst_downmix_impl(int sample_rate,
                       int search_depth,
                       size_t hard_max_queue_len,
                       const std::vector<float>& input_taps,
                       const std::vector<float>& start_finder_taps,
                       bool handle_multiple_frames_per_burst,
                       bool wait_for_demod,
                       int tch_c_type);
    ~burst_downmix_impl();

    size_t get_input_queue_size();
    uint64_t get_n_dropped_bursts();
    void debug_id(uint64_t id);

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_BURST_DOWNMIX_IMPL_H */
