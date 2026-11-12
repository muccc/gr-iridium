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

#ifndef INCLUDED_IRIDIUM_FFT_BURST_TAGGER_IMPL_H
#define INCLUDED_IRIDIUM_FFT_BURST_TAGGER_IMPL_H

#include <gnuradio/fft/fft.h>
#include <iridium/fft_burst_tagger.h>

namespace gr {
namespace iridium {

struct burst {
    uint64_t start;
    uint64_t stop;
    uint64_t last_active;
    int center_bin;
    float magnitude;
    float noise;
    uint64_t id;
    bool wide;
    bool tracked; // opened by a watch: does not mask, ends at hard_stop
    uint64_t hard_stop;
};

/*
 * Timing of the tracking of TCH_C channels in seconds. All TCH_C bursts (TCH_C1, C2 and
 * C8) last 8.27 ms and a channel has one burst per 90 ms TDMA frame.
 */
// TDMA frame length
const double TRACK_FRAME_LENGTH = 0.09;
// A wideband burst active this long looks like a TCH_C burst and starts tracking
const double TRACK_MIN_BURST_LENGTH = 0.006;
const double TRACK_MAX_BURST_LENGTH = 0.012;
// A tracked window opens this long before the expected start of the burst ...
const double TRACK_LEAD = 0.004;
// ... and stays open for the burst plus margin and the same time after it
const double TRACK_SLOT_LENGTH = 0.0085;

// A channel on which a TCH_C like burst was seen: open a burst every TDMA frame
struct watch {
    int bin;
    uint64_t next_trigger;
    uint64_t until;
};

struct peak {
    int bin;
    float relative_magnitude;
};

class fft_burst_tagger_impl : public fft_burst_tagger
{
private:
    bool d_history_primed;
    bool d_debug;
    bool d_offline;

    int d_fft_size;
    int d_burst_pre_len;
    int d_history_size;
    int d_burst_width;
    int d_history_index;
    int d_burst_post_len;
    int d_max_bursts;
    int d_max_burst_len;
    int d_sample_rate;
    int d_squelch_count;
    uint64_t d_index;
    uint64_t d_burst_id;
    uint64_t d_n_tagged_bursts;
    uint64_t d_sample_count;
    uint64_t d_last_rx_time_offset;
    uint64_t d_last_rx_time_timestamp;

    float* d_window_f;
    float* d_magnitude_shifted_f;
    float* d_baseline_sum_f;
    float* d_baseline_history_f;
    float* d_relative_magnitude_f;
    float* d_burst_mask_f;
    float* d_ones_f;
    float d_threshold;
    float d_wide_threshold; // 0: wideband detection disabled
    int d_wide_bins;
    std::vector<float> d_wide_f;
    std::vector<float> d_wide_masked_f;
    std::vector<float> d_wide_tmp_f;
    int d_track_frames;    // 0: tracking disabled
    uint64_t d_frame_len;  // TDMA frame length in samples (90 ms)
    uint64_t d_track_lead; // open the window this many samples before the slot
    std::vector<watch> d_watches;
    double d_center_frequency;
    float d_window_enbw;

    FILE* d_burst_debug_file;

    gr::fft::fft_complex_fwd* d_fft;
    std::vector<peak> d_peaks;
    std::vector<burst> d_bursts;
    std::vector<burst> d_new_bursts;
    std::vector<burst> d_gone_bursts;

    bool update_filters_pre(void);
    void update_filters_post(bool force);
    void extract_peaks(void);
    void save_peaks_to_debug_file(char* filename);
    void remove_peaks_around_bursts(void);
    void update_burst_mask(void);
    void update_bursts(void);
    void delete_gone_bursts(void);
    void create_new_bursts(void);
    void moving_average(const float* in, float* out);
    void create_wide_bursts(void);
    void add_watch(const burst& b);
    void create_tracked_bursts(void);
    void mask_burst(burst& b);
    void tag_new_bursts(void);
    void tag_gone_bursts(int noutput_items);

public:
    fft_burst_tagger_impl(double center_frequency,
                          int fft_size,
                          int sample_rate,
                          int burst_pre_len,
                          int burst_post_len,
                          int burst_width,
                          int max_bursts,
                          int max_burst_len,
                          float threshold,
                          int history_size,
                          bool offline,
                          bool debug,
                          float wideband_threshold,
                          int wideband_width,
                          int track_frames);
    ~fft_burst_tagger_impl();

    uint64_t get_n_tagged_bursts();
    uint64_t get_sample_count();

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FFT_BURST_TAGGER_IMPL_H */
