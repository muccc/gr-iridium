/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_IMPL_H
#define INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_IMPL_H

#include "tch_c_frame_timing.h"
#include "tch_c_slot_survey.h"
#include <iridium/tch_c_slot_prober.h>
#include <deque>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace gr {
namespace iridium {

// Ids of the probes, above those of the tagger (fits the I:%011d of the output)
const uint64_t TCH_C_PROBE_ID_BASE = 50000000000ULL;

// Simplex ring alert channel
const double TCH_C_F_RING = 1626270833.33;

// A probe published and not forgotten yet
struct tch_c_probe_info {
    int type;
    uint64_t t_ns; // first UW symbol
    double f_hz;
    float rho2;
    int track; // -1: none
    int position;
    int slot;
    bool uplink; // probe of an uplink slot
    bool acked;  // tch_c_demod has handled it
};

// A TCH_C frame from tch_c_demod: one of our probes or from the energy path
struct tch_c_hit {
    uint64_t id;
    uint64_t t_ns; // first UW symbol
    double f_hz;   // carrier after the fit to the known symbols
    int type;      // 1, 2, 8
    int bits_per_symbol;
    bool uplink;
    bool crc_ok;
    int confidence;
    float level;
    float magnitude;
    float noise;
    int variant;                   // carrier recovery variant of tch_c_demod
    const tch_c_probe_info* probe; // nullptr for frames of the energy path
};

class tch_c_slot_prober_impl : public tch_c_slot_prober
{
private:
    const double d_center_frequency;
    const double d_sample_rate;
    const std::vector<int> d_types;
    const float d_history_s;
    const float d_threshold[3]; // per type index
    const int d_outstanding_limit;
    const bool d_drop_overflow;
    const bool d_uplink;
    const bool d_debug;
    // Surveyed positions with at least this rho2 are logged (IRIDIUM_PROBE_DEBUG_RHO2),
    // above 1: none
    const float d_debug_rho2;
    // All probe times shifted by this (IRIDIUM_PROBE_OFFSET_US), for the false alarm
    // control
    const int64_t d_offset_ns;
    const float d_variants_rho2;

    // Ring buffer of the input: absolute sample n is at n % d_ring.size()
    std::vector<gr_complex> d_ring;
    uint64_t d_end; // absolute index of the next sample
    std::vector<gr_complex> d_scratch;
    // Input appended per step at most, and how far the probing lags the input
    uint64_t d_chunk;
    uint64_t d_lag;
    bool d_flushed;

    // Time bases from the time_base tags of fft_burst_tagger: from sample `start` on,
    // sample n is at timestamp + (n - offset) / fs
    struct time_base {
        uint64_t start;
        uint64_t offset;
        uint64_t timestamp;
    };
    std::vector<time_base> d_time_bases;

    // Ring alerts accepted by the frames handler, handed to on_ira() in work()
    std::deque<std::pair<uint64_t, double>> d_iras;

    // Satellites and, per track and type, the next DL slot to probe
    tch_c_track_table d_tracks;
    struct cursor {
        int64_t frame;
        int slot;
        bool valid;
    };
    std::map<int, std::vector<cursor>> d_cursors;
    /*
     * Uplink (TCH_C1). A UL slot is searched on the channels whose DL slot of the same
     * number carried a frame recently (the UL burst of a terminal is on its DL
     * channel), widely until a UL frame was found there, then around it. The UL slots
     * lie at v = W0 + UL_SLOT_MS * slot (see tch_c_frame_timing.h), W0 is learned from
     * the UL frames; until then a sweep searches all of the frame on the channels with
     * DL frames (acquisition). All of it lags the DL slots by UL_DELAY_S, so the DL
     * frames of the frames before are known.
     */
    const bool d_probe_uplink;
    struct ul_link {
        int64_t last_dl_ns;
        int64_t last_ul_ns;
        double v_ms;        // v of the newest UL frame
        double f_offset_hz; // its carrier minus nominal * (1 - Doppler / F_RING)
    };
    std::map<std::tuple<int, int, int>, ul_link> d_ul_links; // (track, slot, position)
    struct ul_cursor {
        cursor c;
        int64_t last_ns; // nominal time of the last UL slot probed
        double w0;       // W0 the cursor was set up with
    };
    std::map<int, ul_cursor> d_ul_cursors;
    std::deque<std::pair<int64_t, double>> d_ul_w; // (t, W0) of the recent UL frames
    bool d_ul_w_valid;
    double d_ul_w0;
    struct ul_sweep {
        int64_t frame; // next frame to search
        int window;    // window within the frame
        double v;      // start of the next window of v, ms
        uint64_t n;    // windows swept
        bool valid;
    };
    std::map<int, ul_sweep> d_ul_sweeps;
    uint64_t d_n_ul_sweep_frames;
    void probe_ul_sweep(const tch_c_sat_track& tr,
                        uint64_t ul_limit,
                        uint64_t begin,
                        int64_t begin_ns);

    void ul_link_dl(int track, int slot, int position, int64_t t_ns);
    void probe_ul_slots(const tch_c_sat_track& tr,
                        uint64_t ul_limit,
                        uint64_t begin,
                        int64_t begin_ns);
    void
    probe_ul_slot(const tch_c_sat_track& track, int64_t frame, int slot, int64_t t_ns);
    void ul_slot_window(double p, double s, int64_t* first, int64_t* end) const;
    uint64_t ul_search(const gr_complex* x,
                       int64_t first,
                       const tch_c_sat_track& track,
                       int slot,
                       int position,
                       double p_c,
                       double s,
                       double f_c,
                       double f_s,
                       tch_c_slot_survey::channel_candidate* cc,
                       double* p_found);
    bool ul_match(int64_t t_ns,
                  double f_hz,
                  double tol_hz,
                  int* track,
                  int* position,
                  int* slot) const;
    uint64_t d_n_ul_slots;
    uint64_t d_n_ul_searches;
    uint64_t d_n_ul_probes;
    uint64_t d_n_ul_crc_ok;

    // Slots probed recently, so no satellite is probed twice in a slot even while two
    // tracks of it exist
    struct probed_slot {
        int type;
        int64_t t_ns;
        double doppler_hz;
    };
    std::deque<probed_slot> d_probed;

    tch_c_slot_survey d_survey;
    std::vector<tch_c_slot_survey::candidate> d_candidates;
    std::vector<gr_complex> d_symbols;

    // Probes
    uint64_t d_next_id;
    std::map<uint64_t, tch_c_probe_info> d_probes;
    std::multiset<uint64_t> d_unacked_times;
    int d_outstanding;
    uint64_t d_newest;
    uint64_t d_n_probes;
    uint64_t d_n_dropped;
    uint64_t d_n_iras;
    uint64_t d_n_hits;
    uint64_t d_n_probe_hits;
    uint64_t d_n_probe_crc_ok;
    uint64_t d_n_energy_associated;
    // Statistics per type index
    uint64_t d_n_surveys[3];
    uint64_t d_n_type_probes[3];
    uint64_t d_n_type_crc_ok[3];
    uint64_t d_n_missed_slots;
    uint64_t d_n_shared_slots;
    // Delay of the ring alerts behind the input: sum and maximum in s
    double d_ira_delay_sum;
    double d_ira_delay_max;
    double d_cpu[4]; // s: survey, fine timing and extraction, all of work(), uplink

    // Hooks for the track table
    void on_ira(uint64_t t_ns, double f_hz);
    void on_hit(const tch_c_hit& hit);
    void on_ul_hit(const tch_c_hit& hit);

    void handle_frame(const pmt::pmt_t& msg);
    void handle_hit(const pmt::pmt_t& msg);
    void handle_frame_handled(const pmt::pmt_t& msg);

    bool blocked() const;
    void wait_for_change(int ms);
    void forget_old_probes();

    void append(const gr_complex* in, size_t n);
    void update_time_base(uint64_t start, uint64_t end);

    // Samples [first, end) a probe of a slot of the type at p_uw may read
    void slot_window(int type, double p_uw, int64_t* first, int64_t* end) const;
    // Probe every slot whose window starts before limit and lies in the buffer
    void probe_slots(uint64_t limit, uint64_t ul_limit);
    void probe_slot(
        const tch_c_sat_track& track, int type, int64_t frame, int slot, int64_t t_ns);
    bool already_probed(int type, int64_t t_ns, double doppler_hz) const;

public:
    tch_c_slot_prober_impl(double center_frequency,
                           double sample_rate,
                           const std::vector<int>& types,
                           float history_s,
                           float threshold,
                           int outstanding_limit,
                           bool drop_overflow,
                           bool uplink);
    ~tch_c_slot_prober_impl();

    bool stop() override;

    // Sample index <-> time, see fft_burst_tagger_impl::tag_time_base()
    double sample_index(uint64_t t_ns) const;
    uint64_t time_ns(double sample_index) const;

    // Oldest sample still in the ring buffer and the one after the newest
    uint64_t ring_begin() const;
    uint64_t ring_end() const { return d_end; }
    // Samples [start, start + len), which must be in the ring buffer. Points into it or,
    // if the samples wrap around its end, to a copy which is valid until the next call.
    const gr_complex* samples(uint64_t start, size_t len);

    /*
     * Publish a probe: n samples at one sample per symbol, the first one at the first UW
     * symbol (n = n_symbols + 1 of the type, see tch_c_demod). Returns its id, 0 if it
     * was dropped (drop_overflow and too many outstanding probes).
     */
    uint64_t publish_probe(int type,
                           uint64_t t_first_uw_ns,
                           double freq_hz,
                           const gr_complex* symbols,
                           size_t n,
                           float magnitude,
                           float noise,
                           float rho2,
                           int track = -1,
                           int position = -1,
                           int slot = -1,
                           bool uplink = false);

    uint64_t get_n_probes() override { return d_n_probes; }
    uint64_t get_n_dropped_probes() override { return d_n_dropped; }
    uint64_t get_n_ring_alerts() override { return d_n_iras; }
    int get_n_outstanding() override { return d_outstanding; }

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items) override;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_IMPL_H */
