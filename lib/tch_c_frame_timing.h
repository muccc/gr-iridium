/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_FRAME_TIMING_H
#define INCLUDED_IRIDIUM_TCH_C_FRAME_TIMING_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <vector>

namespace gr {
namespace iridium {

/*
 * TDMA frame timing and Doppler of the satellites, as received, from ring alerts (IRA,
 * simplex ring channel) and from TCH_C downlink hits.
 *
 * Phase: the time of a frame's reference (the start of the IRA of that frame) mod 90 ms.
 * It follows the propagation delay: d(phase)/dt = -Doppler / F_RING (phase in s, Doppler
 * of the ring channel in Hz). Downlink slot k of a frame has its first UW symbol at
 * frame * 90 ms + phase + offset(type) + bias(type) + k * 8.380 ms.
 *
 * Times are int64 ns of the stream time base, phases double ms. Within a track the phase
 * is kept continuous (not wrapped to [0, 90)), so frame numbers stay the same over the
 * lifetime of a track: frame n of a track starts at n * 90 ms + phase.
 */
namespace tch_c_timing {

const double F_RING = 1626270833.3333333; // 1616 MHz + 246.5 * 10 MHz / 240
const double SIMPLEX_CHANNEL_HZ = 10e6 / 240; // spacing of the simplex channels, ~41.7 kHz
const int64_t FRAME_NS = 90000000;
const double FRAME_MS = 90.0;
const double SLOT_MS = 8.380; // spacing of the downlink slots
const int N_SLOTS = 4;        // downlink slots per frame
const double MS_PER_HZ_S =
    1e3 / F_RING; // phase change in ms per Hz of Doppler and second

// Index into per-type arrays: TCH_C1 0, TCH_C2 1, TCH_C8 (QPSK and 16APSK) 2
inline int type_index(int type) { return type == 8 ? 2 : type == 2 ? 1 : 0; }

// Offset of the first UW symbol of DL slot 0 relative to the frame reference, in ms
inline double offset_ms(int type)
{
    static const double offset[3] = { 53.017, 53.008, 53.002 };
    return offset[type_index(type)];
}

/*
 * Uplink slots as received near the terminal: v = (t + phase(t)) mod 90 ms of the first
 * UW symbol is W0 + UL_SLOT_MS * slot, phase the DL frame phase of the satellite. W0 is
 * twice the epoch of the satellite frames in the stream's time base (plus a constant),
 * the same for all satellites, but it moves with the time origin of the stream and
 * drifts with its clock: it has to be learned from UL frames (sat-147 from its start:
 * 23.07 ms; slots at 23.07, 31.41, 39.85, 48.29 ms). Terminals farther from the
 * receiver differ by up to ~0.2 ms.
 */
const double UL_SLOT_MS = 8.407;

// x wrapped to [-FRAME_MS / 2, FRAME_MS / 2)
inline double wrap_ms(double x) { return x - FRAME_MS * std::floor(x / FRAME_MS + 0.5); }

// t mod 90 ms in ms, exact for negative t too
inline double frame_phase_ms(int64_t t_ns)
{
    int64_t r = t_ns % FRAME_NS;
    if (r < 0) {
        r += FRAME_NS;
    }
    return r * 1e-6;
}

// Weighted median of (value, weight) pairs; reorders v
inline double weighted_median(std::vector<std::pair<double, double>>& v)
{
    if (v.empty()) {
        return 0;
    }
    std::sort(v.begin(), v.end());
    double total = 0;
    for (const auto& p : v) {
        total += p.second;
    }
    double sum = 0;
    for (size_t i = 0; i < v.size(); i++) {
        sum += v[i].second;
        if (sum > total / 2) {
            return v[i].first;
        }
        if (sum == total / 2 && i + 1 < v.size()) {
            return (v[i].first + v[i + 1].first) / 2;
        }
    }
    return v.back().first;
}

} // namespace tch_c_timing

struct tch_c_timing_config {
    double window_s = 4.0;        // measurements kept: the last window_s seconds ...
    size_t max_measurements = 64; // ... but at most this many
    size_t min_fit = 3;           // Doppler rate fitted from >= min_fit measurements ...
    double min_fit_span_s = 0.5;  // ... spanning >= min_fit_span_s, else rate 0
    double doppler_outlier_hz =
        100; // Doppler residuals beyond max(this, 4 sigma) are ignored
    double hit_interval_s =
        0.25;               // at most one hit measurement per track in this interval
    float hit_weight = 1.f; // weight of a hit measurement in the Doppler fit, IRAs 1
    // Weight of a hit in the phase median, IRAs 1. Small: a hit's phase depends on the
    // bias of its type, which is learned against the IRA phase, so the IRAs must dominate
    // the phase while there are any. Without IRAs the hits carry the track.
    double hit_phase_weight = 0.05;
    double gate_phase_us = 40; // association gate
    double gate_doppler_hz = 400;
    double gate_doppler_hz_per_s = 400;
    size_t max_tracks = 8;
    double expiry_s = 20; // track removed after this long without a measurement
    double unconfirmed_expiry_s = 20; // the same for tracks with a single IRA
    double stale_s = 3;               // track stale after this long without a measurement
    double bias_gain = 1.0 / 16;      // bias update per fine timing result
    double bias_step_ns =
        10000;                  // timing residuals clamped to +- this for the bias update
    double max_bias_ns = 20000; // bias clamped to +- this
    // New tracks start from the median bias of the tracks with at least this many bias
    // updates of the type (a track's bias has then converged from where it started)
    uint64_t bias_min_updates = 32;
};

struct tch_c_timing_measurement {
    int64_t t_ns;      // stream time
    double phase_ms;   // in [0, 90)
    double doppler_hz; // Doppler at the ring channel frequency
    float weight;
    bool hit; // from a TCH_C hit, else from an IRA
};

struct tch_c_timing_prediction {
    double phase_ms;      // in [0, 90)
    double phase_cont_ms; // continuous phase of the track
    double doppler_hz;    // at the ring channel frequency
    double doppler_rate;  // Hz/s
};

/*
 * One satellite. Measurements may arrive in any time order; the model is refitted after
 * every accepted measurement and predictions are O(1).
 *
 * Model: Doppler d(t) = a + b x, x = (t - t_ref) in s, t_ref = newest measurement.
 * Phase(t) = median_i [phase_i - MS_PER_HZ_S * integral_{t_i}^{t} d(s) ds], unwrapped mod
 * 90 around the previous estimate. With the linear Doppler the integral separates:
 * phase(t) = phase0 - MS_PER_HZ_S * (a x + b x^2 / 2), phase0 = median_i c_i,
 * c_i = phase_i + MS_PER_HZ_S * (a x_i + b x_i^2 / 2).
 */
class tch_c_sat_track
{
public:
    tch_c_sat_track(int id,
                    const tch_c_timing_config& cfg,
                    const double* bias_ns = nullptr)
        : d_id(id), d_cfg(cfg)
    {
        for (int i = 0; i < 3; i++) {
            const double b = bias_ns ? bias_ns[i] : 0;
            d_bias_ns[i] = std::max(-d_cfg.max_bias_ns, std::min(d_cfg.max_bias_ns, b));
        }
    }

    int id() const { return d_id; }
    const std::vector<tch_c_timing_measurement>& measurements() const { return d_meas; }
    bool empty() const { return d_meas.empty(); }
    int64_t first_ns() const { return d_first_ns; } // oldest measurement ever
    int64_t last_ns() const { return d_meas.empty() ? INT64_MIN : d_meas.back().t_ns; }
    int64_t last_ira_ns() const { return d_last_ira_ns; }
    uint64_t n_ira() const { return d_n_ira; } // IRAs accepted over the lifetime
    uint64_t n_hit() const
    {
        return d_n_hit;
    } // hits seen over the lifetime (incl. decimated)
    bool confirmed() const { return d_n_ira + d_n_hit >= 2; }
    bool stale(int64_t now_ns) const { return now_ns - last_ns() > d_cfg.stale_s * 1e9; }
    bool expired(int64_t now_ns) const
    {
        const double e =
            d_n_ira + d_n_hit >= 2 ? d_cfg.expiry_s : d_cfg.unconfirmed_expiry_s;
        return d_meas.empty() || now_ns - last_ns() > e * 1e9;
    }

    double bias_ns(int type) const { return d_bias_ns[tch_c_timing::type_index(type)]; }
    void set_bias_ns(int type, double b)
    {
        d_bias_ns[tch_c_timing::type_index(type)] = b;
    }
    // Bias updates of the type (strong hits) over the lifetime
    uint64_t n_bias_updates(int type) const
    {
        return d_n_bias[tch_c_timing::type_index(type)];
    }

    // Measurement from an IRA (frequency as received). Returns false if outside the
    // window.
    bool add_ira(int64_t t_ns, double freq_hz)
    {
        tch_c_timing_measurement m = { t_ns,
                                       tch_c_timing::frame_phase_ms(t_ns),
                                       freq_hz - tch_c_timing::F_RING,
                                       1.f,
                                       false };
        if (!add(m)) {
            return false;
        }
        d_n_ira++;
        d_last_ira_ns = std::max(d_last_ira_ns, t_ns);
        return true;
    }

    /*
     * Measurement from a CRC:OK downlink hit in DL slot `slot` (0..3) of the raster
     * channel with nominal carrier nominal_hz. With update_bias (strong hit, fine timing)
     * the timing residual against the current prediction updates the bias of the type
     * first. Returns the residual t_uw - predicted slot time in ns (0 for an empty
     * track).
     */
    double add_hit(int64_t t_uw_ns,
                   double freq_hz,
                   int type,
                   int slot,
                   double nominal_hz,
                   bool update_bias = false)
    {
        using namespace tch_c_timing;
        double residual = 0;
        if (!d_meas.empty()) {
            const int64_t frame = frame_of_slot(t_uw_ns, slot, type);
            residual = double(t_uw_ns - slot_time_ns(frame, slot, type));
            if (update_bias) {
                const double r =
                    std::max(-d_cfg.bias_step_ns, std::min(d_cfg.bias_step_ns, residual));
                double& b = d_bias_ns[type_index(type)];
                b = std::max(-d_cfg.max_bias_ns,
                             std::min(d_cfg.max_bias_ns, b + d_cfg.bias_gain * r));
                d_n_bias[type_index(type)]++;
            }
        }
        d_n_hit++;
        // decimate: hits are far more frequent than IRAs
        for (const auto& m : d_meas) {
            if (m.hit && std::llabs(m.t_ns - t_uw_ns) < d_cfg.hit_interval_s * 1e9) {
                return residual;
            }
        }
        const double ref_ms = (t_uw_ns % FRAME_NS) * 1e-6 - offset_ms(type) -
                              bias_ns(type) * 1e-6 - SLOT_MS * slot;
        tch_c_timing_measurement m = { t_uw_ns,
                                       ref_ms - FRAME_MS * std::floor(ref_ms / FRAME_MS),
                                       (freq_hz / nominal_hz - 1) * F_RING,
                                       d_cfg.hit_weight,
                                       true };
        add(m);
        return residual;
    }

    // Generic measurement; returns false if it is older than the window
    bool add(const tch_c_timing_measurement& m)
    {
        const int64_t window = int64_t(d_cfg.window_s * 1e9);
        if (!d_meas.empty() && m.t_ns < d_meas.back().t_ns - window) {
            return false;
        }
        auto it = std::upper_bound(
            d_meas.begin(),
            d_meas.end(),
            m.t_ns,
            [](int64_t t, const tch_c_timing_measurement& x) { return t < x.t_ns; });
        d_meas.insert(it, m);
        d_first_ns = std::min(d_first_ns, m.t_ns);
        const int64_t oldest = d_meas.back().t_ns - window;
        size_t drop = 0;
        while (drop < d_meas.size() && d_meas[drop].t_ns < oldest) {
            drop++;
        }
        drop = std::max(drop,
                        d_meas.size() > d_cfg.max_measurements
                            ? d_meas.size() - d_cfg.max_measurements
                            : size_t(0));
        d_meas.erase(d_meas.begin(), d_meas.begin() + drop);
        fit();
        return true;
    }

    // Move the measurements of another track of the same satellite into this one
    void merge(const tch_c_sat_track& o)
    {
        for (const auto& m : o.d_meas) {
            add(m);
        }
        d_n_ira += o.d_n_ira;
        d_n_hit += o.d_n_hit;
        d_first_ns = std::min(d_first_ns, o.d_first_ns);
        d_last_ira_ns = std::max(d_last_ira_ns, o.d_last_ira_ns);
    }

    // Continuous phase in ms at time t
    double phase_cont_ms(int64_t t_ns) const
    {
        const double x = (t_ns - d_t_ref) * 1e-9;
        return d_phase0 - tch_c_timing::MS_PER_HZ_S * x * (d_a + 0.5 * d_b * x);
    }
    double doppler(int64_t t_ns) const { return d_a + d_b * (t_ns - d_t_ref) * 1e-9; }
    double doppler_rate() const { return d_b; }

    tch_c_timing_prediction predict(int64_t t_ns) const
    {
        const double p = phase_cont_ms(t_ns);
        return { p - tch_c_timing::FRAME_MS * std::floor(p / tch_c_timing::FRAME_MS),
                 p,
                 doppler(t_ns),
                 d_b };
    }

    // Received carrier of a raster channel with nominal carrier nominal_hz at time t
    double carrier_hz(double nominal_hz, int64_t t_ns) const
    {
        return nominal_hz * (1 + doppler(t_ns) / tch_c_timing::F_RING);
    }

    // Frame containing t: frame n covers [n * 90 ms + phase, (n + 1) * 90 ms + phase)
    int64_t frame_at(int64_t t_ns) const
    {
        const double rel = (t_ns - int64_t(std::llround(phase_cont_ms(t_ns) * 1e6))) /
                           double(tch_c_timing::FRAME_NS);
        return int64_t(std::floor(rel));
    }

    // First UW symbol of DL slot `slot` of frame `frame` for a TCH_C type (1, 2, 8)
    int64_t slot_time_ns(int64_t frame, int slot, int type) const
    {
        const int64_t base = frame * tch_c_timing::FRAME_NS;
        const double off_ns =
            (tch_c_timing::offset_ms(type) + tch_c_timing::SLOT_MS * slot) * 1e6 +
            bias_ns(type);
        // T = base + phase(T) + offset: two fixed point iterations (|dphase/dt| < 3e-5)
        int64_t t = base + int64_t(std::llround(phase_cont_ms(base) * 1e6 + off_ns));
        for (int i = 0; i < 2; i++) {
            t = base + int64_t(std::llround(phase_cont_ms(t) * 1e6 + off_ns));
        }
        return t;
    }

    // Uplink (see UL_SLOT_MS): v = (t + phase(t)) mod 90 ms
    double ul_phase_ms(int64_t t_ns) const
    {
        const double v =
            std::fmod(tch_c_timing::frame_phase_ms(t_ns) + phase_cont_ms(t_ns),
                      tch_c_timing::FRAME_MS);
        return v < 0 ? v + tch_c_timing::FRAME_MS : v;
    }

    // The UL time with v = v_ms in frame `frame`: T = frame * 90 ms + v - phase(T)
    int64_t ul_time_ns(int64_t frame, double v_ms) const
    {
        const int64_t base = frame * tch_c_timing::FRAME_NS;
        int64_t t = base + int64_t(std::llround((v_ms - phase_cont_ms(base)) * 1e6));
        for (int i = 0; i < 2; i++) {
            t = base + int64_t(std::llround((v_ms - phase_cont_ms(t)) * 1e6));
        }
        return t;
    }

    // Frame whose UL time with v = v_ms lies nearest to t
    int64_t ul_frame_at(int64_t t_ns, double v_ms) const
    {
        return int64_t(std::llround((t_ns + (phase_cont_ms(t_ns) - v_ms) * 1e6) /
                                    tch_c_timing::FRAME_NS));
    }

    // Frame whose slot `slot` of the type lies nearest to t
    int64_t frame_of_slot(int64_t t_ns, int slot, int type) const
    {
        const double off_ns =
            (tch_c_timing::offset_ms(type) + tch_c_timing::SLOT_MS * slot) * 1e6 +
            bias_ns(type);
        const double rel =
            (t_ns - phase_cont_ms(t_ns) * 1e6 - off_ns) / tch_c_timing::FRAME_NS;
        return int64_t(std::llround(rel));
    }

    // The first DL slot of the type whose UW starts at or after t
    void next_slot(int64_t t_ns, int type, int64_t& frame, int& slot) const
    {
        for (int64_t n = frame_at(t_ns) - 1;; n++) {
            for (int k = 0; k < tch_c_timing::N_SLOTS; k++) {
                if (slot_time_ns(n, k, type) >= t_ns) {
                    frame = n;
                    slot = k;
                    return;
                }
            }
        }
    }

    // Wrapped phase difference prediction - measurement in us and Doppler difference in
    // Hz
    void residual(const tch_c_timing_measurement& m,
                  double& dphase_us,
                  double& ddoppler_hz) const
    {
        dphase_us = tch_c_timing::wrap_ms(phase_cont_ms(m.t_ns) - m.phase_ms) * 1e3;
        ddoppler_hz = doppler(m.t_ns) - m.doppler_hz;
    }

    // Distance in s from t to the nearest stored measurement
    double gap_s(int64_t t_ns) const
    {
        int64_t g = INT64_MAX;
        for (const auto& m : d_meas) {
            g = std::min(g, int64_t(std::llabs(m.t_ns - t_ns)));
        }
        return g * 1e-9;
    }

private:
    void fit()
    {
        using namespace tch_c_timing;
        const size_t n = d_meas.size();
        if (n == 0) {
            return;
        }
        const int64_t t_ref = d_meas.back().t_ns;
        // previous continuous phase at the new reference time: keeps the unwrapping
        // stable
        const double prev = d_have_fit ? phase_cont_ms(t_ref) : d_meas.back().phase_ms;
        std::vector<double> x(n);
        for (size_t i = 0; i < n; i++) {
            x[i] = (d_meas[i].t_ns - t_ref) * 1e-9;
        }

        // Doppler: robust line (Theil-Sen), outlier rejection, weighted least squares
        std::vector<std::pair<double, double>> tmp;
        double a = 0, b = 0;
        const bool rate =
            n >= d_cfg.min_fit && x.back() - x.front() >= d_cfg.min_fit_span_s;
        if (rate) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    if (x[j] - x[i] > 0.1) {
                        tmp.emplace_back((d_meas[j].doppler_hz - d_meas[i].doppler_hz) /
                                             (x[j] - x[i]),
                                         1.0);
                    }
                }
            }
            b = weighted_median(tmp);
        }
        tmp.clear();
        for (size_t i = 0; i < n; i++) {
            tmp.emplace_back(d_meas[i].doppler_hz - b * x[i], d_meas[i].weight);
        }
        a = weighted_median(tmp);
        tmp.clear();
        for (size_t i = 0; i < n; i++) {
            tmp.emplace_back(std::fabs(d_meas[i].doppler_hz - a - b * x[i]), 1.0);
        }
        const double limit =
            std::max(d_cfg.doppler_outlier_hz, 4 * 1.4826 * weighted_median(tmp));
        double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (size_t i = 0; i < n; i++) {
            const double y = d_meas[i].doppler_hz;
            if (std::fabs(y - a - b * x[i]) > limit) {
                continue;
            }
            const double w = d_meas[i].weight;
            sw += w;
            sx += w * x[i];
            sy += w * y;
            sxx += w * x[i] * x[i];
            sxy += w * x[i] * y;
        }
        if (sw > 0) {
            const double det = sw * sxx - sx * sx;
            if (rate && det > 1e-12 * sw * sw) {
                b = (sw * sxy - sx * sy) / det;
                a = (sy - b * sx) / sw;
            } else {
                b = 0;
                a = sy / sw;
            }
        }

        // Phase: median of the measurements projected to t_ref, unwrapped around prev
        double ref = prev;
        double phase0 = prev;
        for (int pass = 0; pass < 2; pass++) {
            tmp.clear();
            for (size_t i = 0; i < n; i++) {
                const double c =
                    d_meas[i].phase_ms + MS_PER_HZ_S * x[i] * (a + 0.5 * b * x[i]);
                const double w =
                    d_meas[i].hit ? d_cfg.hit_phase_weight : d_meas[i].weight;
                tmp.emplace_back(ref + wrap_ms(c - ref), w);
            }
            phase0 = weighted_median(tmp);
            ref = phase0;
        }
        d_t_ref = t_ref;
        d_a = a;
        d_b = b;
        d_phase0 = phase0;
        d_have_fit = true;
    }

    int d_id;
    tch_c_timing_config d_cfg;
    std::vector<tch_c_timing_measurement> d_meas; // sorted by time
    double d_bias_ns[3];
    uint64_t d_n_bias[3] = { 0, 0, 0 };
    int64_t d_first_ns = INT64_MAX;
    int64_t d_last_ira_ns = INT64_MIN;
    uint64_t d_n_ira = 0;
    uint64_t d_n_hit = 0;
    bool d_have_fit = false;
    int64_t d_t_ref = 0;
    double d_a = 0;
    double d_b = 0;
    double d_phase0 = 0;
};

/*
 * All tracked satellites. Track ids are unique over the lifetime of the table; a track
 * removed by expiry or merging never comes back, so a caller can key its state by id.
 */
class tch_c_track_table
{
public:
    struct slot_match {
        int track;
        int64_t frame;
        int slot;
        double dt_ns; // measured - predicted time
        double df_hz; // measured - predicted carrier
    };

    explicit tch_c_track_table(const tch_c_timing_config& cfg = tch_c_timing_config())
        : d_cfg(cfg)
    {
        for (double& b : d_bias_ns) {
            b = 0;
        }
    }

    const tch_c_timing_config& config() const { return d_cfg; }
    std::vector<tch_c_sat_track>& tracks() { return d_tracks; }
    const std::vector<tch_c_sat_track>& tracks() const { return d_tracks; }
    std::vector<tch_c_sat_track>::iterator begin() { return d_tracks.begin(); }
    std::vector<tch_c_sat_track>::iterator end() { return d_tracks.end(); }
    std::vector<tch_c_sat_track>::const_iterator begin() const
    {
        return d_tracks.begin();
    }
    std::vector<tch_c_sat_track>::const_iterator end() const { return d_tracks.end(); }
    size_t size() const { return d_tracks.size(); }

    const tch_c_sat_track* find(int id) const
    {
        return const_cast<tch_c_track_table*>(this)->find(id);
    }
    tch_c_sat_track* find(int id)
    {
        for (auto& t : d_tracks) {
            if (t.id() == id) {
                return &t;
            }
        }
        return nullptr;
    }

    // Track an IRA belongs to by the gate, -1 if none
    int match_ira(int64_t t_ns, double freq_hz) const
    {
        const tch_c_timing_measurement m = { t_ns,
                                             tch_c_timing::frame_phase_ms(t_ns),
                                             freq_hz - tch_c_timing::F_RING,
                                             1.f,
                                             false };
        return match(m);
    }

    /*
     * An IRA: added to the matching track, else a new track is created (the least
     * recently updated one is dropped if there are max_tracks). Returns the track id.
     */
    /*
     * A satellite also transmits on the neighbour simplex channels (messaging, not ring
     * alerts). Those bursts share its frame phase but are a channel width off in
     * frequency, so freq - F_RING looks like a Doppler a channel width away and they
     * each seed a spurious track. The IRA frequency window (+-45 kHz) admits the +-1
     * neighbours. A candidate that matches a confirmed track once its frequency is
     * shifted back by a whole channel width (same phase, Doppler then agrees) is such a
     * neighbour: return that track's id, else -1.
     */
    int neighbor_image(int64_t t_ns, double freq_hz) const
    {
        for (int k = -1; k <= 1; k++) {
            if (k == 0) {
                continue;
            }
            const tch_c_timing_measurement m = {
                t_ns,
                tch_c_timing::frame_phase_ms(t_ns),
                freq_hz - tch_c_timing::F_RING - k * tch_c_timing::SIMPLEX_CHANNEL_HZ,
                1.f,
                false
            };
            const int id = match(m);
            const tch_c_sat_track* tr = id >= 0 ? find(id) : nullptr;
            if (tr && tr->confirmed()) {
                return id;
            }
        }
        return -1;
    }

    int add_ira(int64_t t_ns, double freq_hz)
    {
        int id = match_ira(t_ns, freq_hz);
        if (id < 0) {
            // A messaging burst of a satellite we already track, on a neighbour simplex
            // channel: drop it instead of seeding a spurious track
            if (neighbor_image(t_ns, freq_hz) >= 0) {
                d_n_neighbor++;
                return -1;
            }
            if (d_tracks.size() >= d_cfg.max_tracks) {
                auto victim = std::min_element(
                    d_tracks.begin(),
                    d_tracks.end(),
                    [](const tch_c_sat_track& x, const tch_c_sat_track& y) {
                        // unconfirmed tracks first, then the least recently updated
                        if (x.confirmed() != y.confirmed()) {
                            return !x.confirmed();
                        }
                        return x.last_ns() < y.last_ns();
                    });
                d_tracks.erase(victim);
                d_n_evicted++;
            }
            id = d_next_id++;
            d_tracks.emplace_back(id, d_cfg, d_bias_ns);
            d_n_created++;
        }
        find(id)->add_ira(t_ns, freq_hz);
        merge();
        return id;
    }

    /*
     * The DL slot of any track nearest to a burst with first UW symbol at t_uw and
     * carrier freq on the raster channel nominal_hz, within the gate (gate_us < 0: config
     * gate).
     */
    bool match_hit(int64_t t_uw_ns,
                   double freq_hz,
                   int type,
                   double nominal_hz,
                   slot_match& best,
                   double gate_us = -1,
                   double gate_hz = -1) const
    {
        if (gate_us < 0) {
            gate_us = d_cfg.gate_phase_us;
        }
        // score: normalised distance; confirmed tracks always win over unconfirmed ones
        bool found = false;
        double best_score = 1e300;
        for (const auto& tr : d_tracks) {
            if (tr.empty()) {
                continue;
            }
            const double gap = tr.gap_s(t_uw_ns);
            const double ghz =
                gate_hz >= 0
                    ? gate_hz
                    : nominal_hz / tch_c_timing::F_RING *
                          (d_cfg.gate_doppler_hz + d_cfg.gate_doppler_hz_per_s * gap);
            const double df = freq_hz - tr.carrier_hz(nominal_hz, t_uw_ns);
            if (std::fabs(df) > ghz) {
                continue;
            }
            for (int k = 0; k < tch_c_timing::N_SLOTS; k++) {
                const int64_t n = tr.frame_of_slot(t_uw_ns, k, type);
                const double dt = double(t_uw_ns - tr.slot_time_ns(n, k, type));
                if (std::fabs(dt) >= gate_us * 1e3) {
                    continue;
                }
                const double score = (tr.confirmed() ? 0 : 10) +
                                     std::fabs(dt) / (gate_us * 1e3) +
                                     std::fabs(df) / ghz;
                if (score < best_score) {
                    best_score = score;
                    best = { tr.id(), n, k, dt, df };
                    found = true;
                }
            }
        }
        return found;
    }

    // A hit in a known slot of a track; see tch_c_sat_track::add_hit
    double add_hit(int track,
                   int64_t t_uw_ns,
                   double freq_hz,
                   int type,
                   int slot,
                   double nominal_hz,
                   bool update_bias = false)
    {
        tch_c_sat_track* tr = find(track);
        if (!tr) {
            return 0;
        }
        const double r =
            tr->add_hit(t_uw_ns, freq_hz, type, slot, nominal_hz, update_bias);
        if (update_bias) {
            update_bias_ns(type);
        }
        merge();
        return r;
    }

    /*
     * The bias new tracks start from: the median of the tracks whose bias of the type
     * has converged, else unchanged. Not the sum of the updates of all tracks: tracks
     * which converge at the same time from the same start (after a rebuild of all of
     * them) would each add their whole correction, and a value beyond the timing range
     * of the survey would lock the type out of every new track, as no hit can then
     * correct it.
     */
    void update_bias_ns(int type)
    {
        std::vector<std::pair<double, double>> v;
        for (const auto& t : d_tracks) {
            if (t.n_bias_updates(type) >= d_cfg.bias_min_updates) {
                v.emplace_back(t.bias_ns(type), 1.0);
            }
        }
        if (!v.empty()) {
            d_bias_ns[tch_c_timing::type_index(type)] =
                std::max(-d_cfg.max_bias_ns,
                         std::min(d_cfg.max_bias_ns, tch_c_timing::weighted_median(v)));
        }
    }

    // A hit of unknown slot (energy path): associated by the gate, returns the track or
    // -1
    int associate_hit(int64_t t_uw_ns, double freq_hz, int type, double nominal_hz)
    {
        slot_match m;
        if (!match_hit(t_uw_ns, freq_hz, type, nominal_hz, m)) {
            return -1;
        }
        add_hit(m.track, t_uw_ns, freq_hz, type, m.slot, nominal_hz);
        return m.track;
    }

    // Removes expired tracks (now = current stream time); returns the number removed
    size_t expire(int64_t now_ns)
    {
        const size_t before = d_tracks.size();
        d_tracks.erase(
            std::remove_if(d_tracks.begin(),
                           d_tracks.end(),
                           [&](const tch_c_sat_track& t) { return t.expired(now_ns); }),
            d_tracks.end());
        d_n_expired += before - d_tracks.size();
        return before - d_tracks.size();
    }

    double bias_ns(int type) const { return d_bias_ns[tch_c_timing::type_index(type)]; }
    uint64_t n_created() const { return d_n_created; }
    uint64_t n_neighbor() const { return d_n_neighbor; }
    uint64_t n_merged() const { return d_n_merged; }
    uint64_t n_expired() const { return d_n_expired; }
    uint64_t n_evicted() const { return d_n_evicted; }

private:
    bool in_gate(const tch_c_sat_track& tr,
                 const tch_c_timing_measurement& m,
                 double& dphase_us,
                 double& ddoppler_hz) const
    {
        tr.residual(m, dphase_us, ddoppler_hz);
        return std::fabs(dphase_us) < d_cfg.gate_phase_us &&
               std::fabs(ddoppler_hz) <
                   d_cfg.gate_doppler_hz + d_cfg.gate_doppler_hz_per_s * tr.gap_s(m.t_ns);
    }

    // Best track in the gate by normalised distance, confirmed tracks first; -1 if none
    int match(const tch_c_timing_measurement& m) const
    {
        int id = -1;
        double best = 1e300;
        for (const auto& tr : d_tracks) {
            double dp, dd;
            if (tr.empty() || !in_gate(tr, m, dp, dd)) {
                continue;
            }
            const double score =
                (tr.confirmed() ? 0 : 10) + std::fabs(dp) / d_cfg.gate_phase_us +
                std::fabs(dd) / (d_cfg.gate_doppler_hz +
                                 d_cfg.gate_doppler_hz_per_s * tr.gap_s(m.t_ns));
            if (score < best) {
                best = score;
                id = tr.id();
            }
        }
        return id;
    }

    // Merges tracks whose predictions coincide at the newer one's last measurement
    void merge()
    {
        for (size_t i = 0; i < d_tracks.size(); i++) {
            for (size_t j = i + 1; j < d_tracks.size(); j++) {
                tch_c_sat_track& x = d_tracks[i];
                tch_c_sat_track& y = d_tracks[j];
                if (x.empty() || y.empty()) {
                    continue;
                }
                const int64_t t = std::max(x.last_ns(), y.last_ns());
                const double dp =
                    tch_c_timing::wrap_ms(x.phase_cont_ms(t) - y.phase_cont_ms(t));
                const double dd = x.doppler(t) - y.doppler(t);
                if (std::fabs(dp) * 1e3 >= d_cfg.gate_phase_us ||
                    std::fabs(dd) >= d_cfg.gate_doppler_hz) {
                    continue;
                }
                // keep the older track (its id and frame numbering)
                const bool keep_x = x.first_ns() <= y.first_ns();
                tch_c_sat_track& keep = keep_x ? x : y;
                keep.merge(keep_x ? y : x);
                d_tracks.erase(d_tracks.begin() + (keep_x ? j : i));
                d_n_merged++;
                i = size_t(-1); // start over
                break;
            }
        }
    }

    tch_c_timing_config d_cfg;
    std::vector<tch_c_sat_track> d_tracks;
    double d_bias_ns[3];
    int d_next_id = 0;
    uint64_t d_n_created = 0;
    uint64_t d_n_neighbor = 0;
    uint64_t d_n_merged = 0;
    uint64_t d_n_expired = 0;
    uint64_t d_n_evicted = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_FRAME_TIMING_H */
