/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "tch_c_frame_timing.h"
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

using namespace gr::iridium;
using namespace gr::iridium::tch_c_timing;

namespace {

const double C = 299792458.0;

// Gaussian noise; std::mt19937 output is the same on every platform
struct noise {
    std::mt19937 rng;
    explicit noise(uint32_t seed) : rng(seed) {}
    double uniform() { return (rng() + 0.5) / 4294967296.0; }
    double gauss()
    {
        return std::sqrt(-2 * std::log(uniform())) * std::cos(2 * M_PI * uniform());
    }
};

/*
 * A satellite pass on a straight line: range r(t) = sqrt(r_min^2 + (v (t - t_c))^2).
 * The satellite sends frame n at n * 90 ms + phi0; it is received after r / c. Doppler
 * and delay follow from the same range, as for a real satellite.
 */
struct sat_pass {
    double t_c;    // closest approach, s
    double r_min;  // m
    double v;      // m/s
    double phi0;   // ms
    int64_t epoch; // ns added to all times

    double range(double t) const
    {
        const double x = v * (t - t_c);
        return std::sqrt(r_min * r_min + x * x);
    }
    double range_rate(double t) const { return v * v * (t - t_c) / range(t); }
    double doppler(double t) const { return -range_rate(t) / C * F_RING; }

    // Received time of something sent at tx (s, relative to the epoch)
    double rx(double tx) const
    {
        double t = tx;
        for (int i = 0; i < 4; i++) {
            t = tx + range(t) / C;
        }
        return t;
    }
    // Received IRA of frame n and first UW symbol of DL slot k
    double ira_s(int64_t n) const { return rx(n * 0.09 + phi0 * 1e-3); }
    double slot_s(int64_t n, int k, int type) const
    {
        return rx(n * 0.09 + (phi0 + offset_ms(type) + SLOT_MS * k) * 1e-3);
    }
    int64_t ns(double t) const { return epoch + int64_t(std::llround(t * 1e9)); }
};

struct ira {
    int64_t t_ns;
    double f_hz;
};

// IRAs of frames first..last with gaps of 1 to max_gap frames, with timing / Doppler
// jitter
std::vector<ira> make_iras(const sat_pass& p,
                           int64_t first,
                           int64_t last,
                           int max_gap,
                           noise& nz,
                           double jit_ns = 1300,
                           double jit_hz = 50)
{
    std::vector<ira> v;
    for (int64_t n = first; n <= last; n += 1 + nz.rng() % max_gap) {
        const double t = p.ira_s(n);
        v.push_back({ p.ns(t) + int64_t(std::llround(jit_ns * nz.gauss())),
                      F_RING + p.doppler(t) + jit_hz * nz.gauss() });
    }
    return v;
}

struct stats {
    double n = 0, s = 0, ss = 0, max_abs = 0;
    void add(double x)
    {
        n++;
        s += x;
        ss += x * x;
        max_abs = std::max(max_abs, std::fabs(x));
    }
    double mean() const { return s / n; }
    double std() const { return std::sqrt(std::max(0.0, ss / n - mean() * mean())); }
};

/*
 * Feeds the IRAs causally (only those received at least delay_ms before the prediction
 * time) and predicts every DL slot of every frame. Returns the time and carrier errors of
 * the track the slot matches, checks that the frame numbering does not jump.
 */
struct replay_result {
    stats dt, df;
    size_t unmatched = 0;
};

replay_result replay(tch_c_track_table& table,
                     const std::vector<sat_pass>& passes,
                     std::vector<ira> iras,
                     int64_t first,
                     int64_t last,
                     double delay_ms = 50)
{
    std::sort(iras.begin(), iras.end(), [](const ira& a, const ira& b) {
        return a.t_ns < b.t_ns;
    });
    replay_result res;
    size_t next = 0;
    std::vector<std::vector<int64_t>> frame_offset(passes.size());
    for (int64_t n = first; n <= last; n++) {
        for (size_t s = 0; s < passes.size(); s++) {
            const sat_pass& p = passes[s];
            for (int k = 0; k < N_SLOTS; k++) {
                const int64_t t_true = p.ns(p.slot_s(n, k, 1));
                while (next < iras.size() && iras[next].t_ns <= t_true - delay_ms * 1e6) {
                    table.add_ira(iras[next].t_ns, iras[next].f_hz);
                    next++;
                }
                table.expire(t_true);
                if (next < 4) { // warm up
                    continue;
                }
                const double nominal = 1616e6 + 100.5 * 1e7 / 240;
                const double f_true =
                    nominal * (1 + p.doppler(p.slot_s(n, k, 1)) / F_RING);
                tch_c_track_table::slot_match m;
                if (!table.match_hit(t_true, f_true, 1, nominal, m)) {
                    res.unmatched++;
                    continue;
                }
                BOOST_REQUIRE_EQUAL(m.slot, k);
                const tch_c_sat_track* tr = table.find(m.track);
                const int64_t t_pred = tr->slot_time_ns(m.frame, k, 1);
                res.dt.add(double(t_pred - t_true));
                res.df.add(tr->carrier_hz(nominal, t_pred) - f_true);
                // frame numbering of a track: frame - n constant
                auto& off = frame_offset[s];
                if (off.size() <= size_t(m.track)) {
                    off.resize(m.track + 1, INT64_MIN);
                }
                if (off[m.track] == INT64_MIN) {
                    off[m.track] = m.frame - n;
                }
                BOOST_REQUIRE_EQUAL(m.frame - n, off[m.track]);
            }
        }
    }
    return res;
}

const double V = 7460;      // orbital velocity, m/s
const double R_MIN = 790e3; // pass through the zenith: largest Doppler rate (~390 Hz/s)

} // namespace

BOOST_AUTO_TEST_CASE(test_tch_c_timing_helpers)
{
    BOOST_CHECK_CLOSE(F_RING, 1616e6 + 246.5 * 1e7 / 240, 1e-12);
    BOOST_CHECK_EQUAL(frame_phase_ms(90000000LL * 1000 + 1234567), 1.234567);
    BOOST_CHECK_CLOSE(frame_phase_ms(-1000000), 89.0, 1e-12);
    BOOST_CHECK_CLOSE(wrap_ms(89.0), -1.0, 1e-12);
    BOOST_CHECK_CLOSE(wrap_ms(-46.0), 44.0, 1e-12);
    BOOST_CHECK_EQUAL(type_index(1), 0);
    BOOST_CHECK_EQUAL(type_index(2), 1);
    BOOST_CHECK_EQUAL(type_index(8), 2);
    std::vector<std::pair<double, double>> v = { { 3, 1 }, { 1, 1 }, { 100, 1 } };
    BOOST_CHECK_EQUAL(weighted_median(v), 3);
    v = { { 3, 1 }, { 1, 5 }, { 100, 1 } };
    BOOST_CHECK_EQUAL(weighted_median(v), 1);
}

// Noise free: the model must follow the pass to a few ns
BOOST_AUTO_TEST_CASE(test_tch_c_timing_exact)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(1);
    auto iras = make_iras(p, 0, 1500, 13, nz, 0, 0);
    auto r = replay(table, { p }, iras, 0, 1500);
    BOOST_TEST_MESSAGE("exact: dt max " << r.dt.max_abs << " ns, df max " << r.df.max_abs
                                        << " Hz");
    BOOST_CHECK_EQUAL(r.unmatched, 0u);
    BOOST_CHECK_LT(r.dt.max_abs, 20.0);
    BOOST_CHECK_LT(r.df.max_abs, 20.0);
    BOOST_CHECK_EQUAL(table.n_created(), 1u);
}

// Jitter 1.3 us / 50 Hz, IRA gaps up to 1.2 s, at a unix epoch time base
BOOST_AUTO_TEST_CASE(test_tch_c_timing_jitter)
{
    for (double t_c : { 20.0, 60.0, 120.0 }) {
        const sat_pass p = { t_c, R_MIN, V, 30.0, 1700000000LL * 1000000000LL };
        tch_c_track_table table;
        noise nz(2);
        auto iras = make_iras(p, 0, 1500, 13, nz);
        auto r = replay(table, { p }, iras, 0, 1500);
        BOOST_TEST_MESSAGE("jitter, closest approach "
                           << t_c << " s: dt mean " << r.dt.mean() << " std "
                           << r.dt.std() << " max " << r.dt.max_abs << " ns; df std "
                           << r.df.std() << " Hz");
        BOOST_CHECK_EQUAL(r.unmatched, 0u);
        BOOST_CHECK_LT(r.dt.std(), 700.0);
        BOOST_CHECK_LT(std::fabs(r.dt.mean()), 200.0);
        BOOST_CHECK_LT(r.dt.max_abs, 2500.0);
        BOOST_CHECK_LT(r.df.std(), 50.0);
        BOOST_CHECK_EQUAL(table.n_created(), 1u);
        BOOST_CHECK_EQUAL(table.size(), 1u);
    }
}

// The phase crosses 90 -> 0 ms during the pass: frame numbers must not jump (checked in
// replay)
BOOST_AUTO_TEST_CASE(test_tch_c_timing_wrap)
{
    for (double phi0 : { 84.0, 89.9 }) {
        const sat_pass p = { 30, R_MIN, V, phi0, 0 }; // delay 2.6 .. 9 ms
        tch_c_track_table table;
        noise nz(3);
        auto iras = make_iras(p, 0, 1500, 13, nz);
        auto r = replay(table, { p }, iras, 0, 1500);
        BOOST_CHECK_EQUAL(r.unmatched, 0u);
        BOOST_CHECK_LT(r.dt.std(), 700.0);
        BOOST_CHECK_EQUAL(table.n_created(), 1u);
    }
}

// Two satellites at once, phases 0.6 ms apart, one track each
BOOST_AUTO_TEST_CASE(test_tch_c_timing_two_satellites)
{
    const sat_pass a = { 40, R_MIN, V, 30.0, 0 };
    const sat_pass b = { 200, 1500e3, V, 30.6, 0 };
    tch_c_track_table table;
    noise nz(4);
    auto iras = make_iras(a, 0, 1500, 13, nz);
    auto ib = make_iras(b, 0, 1500, 13, nz);
    iras.insert(iras.end(), ib.begin(), ib.end());
    auto r = replay(table, { a, b }, iras, 0, 1500);
    BOOST_TEST_MESSAGE("two satellites: dt std " << r.dt.std() << " max " << r.dt.max_abs
                                                 << " ns; df std " << r.df.std()
                                                 << " Hz");
    BOOST_CHECK_EQUAL(r.unmatched, 0u);
    BOOST_CHECK_LT(r.dt.std(), 700.0);
    BOOST_CHECK_LT(r.df.std(), 55.0);
    BOOST_CHECK_EQUAL(table.n_created(), 2u);
    BOOST_CHECK_EQUAL(table.size(), 2u);
}

// IRAs arrive out of order (up to 0.3 s late) and up to 1 s behind the prediction time
BOOST_AUTO_TEST_CASE(test_tch_c_timing_out_of_order)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(5);
    auto iras = make_iras(p, 0, 1500, 13, nz);
    // feed in arrival order: arrival = time + U(0, 0.3 s)
    std::vector<std::pair<int64_t, ira>> arrival;
    for (const auto& x : iras) {
        arrival.push_back({ x.t_ns + int64_t(nz.uniform() * 3e8), x });
    }
    std::sort(arrival.begin(), arrival.end(), [](const auto& x, const auto& y) {
        return x.first < y.first;
    });
    size_t next = 0;
    stats dt;
    for (int64_t n = 0; n <= 1500; n++) {
        const int64_t t_true = p.ns(p.slot_s(n, 2, 1));
        while (next < arrival.size() && arrival[next].first <= t_true - 1000000000LL) {
            table.add_ira(arrival[next].second.t_ns, arrival[next].second.f_hz);
            next++;
        }
        if (next < 4) {
            continue;
        }
        const tch_c_sat_track& tr = table.tracks().at(0);
        dt.add(double(tr.slot_time_ns(tr.frame_of_slot(t_true, 2, 1), 2, 1) - t_true));
    }
    BOOST_TEST_MESSAGE("out of order: dt std " << dt.std() << " max " << dt.max_abs
                                               << " ns");
    BOOST_CHECK_LT(dt.std(), 800.0);
    BOOST_CHECK_LT(dt.max_abs, 3000.0);
    BOOST_CHECK_EQUAL(table.n_created(), 1u);
}

/*
 * Outliers: 15% garbage IRAs (random time, Doppler +-40 kHz) and 10% with a timing error
 * of 15..35 us (inside the gate). The true track must stay accurate; garbage tracks stay
 * unconfirmed and expire.
 */
BOOST_AUTO_TEST_CASE(test_tch_c_timing_outliers)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_timing_config cfg;
    cfg.unconfirmed_expiry_s = 3;
    tch_c_track_table table(cfg);
    noise nz(6);
    auto iras = make_iras(p, 0, 1500, 13, nz);
    const size_t n_good = iras.size();
    for (size_t i = 0; i < n_good; i++) {
        const double u = nz.uniform();
        if (u < 0.10) {
            iras[i].t_ns +=
                int64_t((nz.uniform() < 0.5 ? -1 : 1) * (15e3 + 20e3 * nz.uniform()));
        } else if (u < 0.13) {
            iras[i].f_hz += (nz.uniform() < 0.5 ? -1 : 1) * (1500 + 3000 * nz.uniform());
        }
    }
    for (size_t i = 0; i < n_good * 15 / 100; i++) {
        iras.push_back(
            { int64_t(nz.uniform() * 135e9), F_RING + 80e3 * (nz.uniform() - 0.5) });
    }
    auto r = replay(table, { p }, iras, 0, 1500);
    BOOST_TEST_MESSAGE("outliers: dt std "
                       << r.dt.std() << " max " << r.dt.max_abs << " ns; df std "
                       << r.df.std() << " Hz; tracks " << table.n_created()
                       << " created, " << table.n_merged() << " merged");
    BOOST_CHECK_EQUAL(r.unmatched, 0u);
    BOOST_CHECK_LT(r.dt.std(), 900.0);
    BOOST_CHECK_LT(r.dt.max_abs, 4000.0);
    BOOST_CHECK_LT(r.df.std(), 70.0);
    BOOST_CHECK_LE(table.size(), 3u);
    size_t confirmed = 0;
    for (const auto& t : table) {
        confirmed += t.confirmed();
    }
    BOOST_CHECK_EQUAL(confirmed, 1u);
}

/*
 * An IRA with a timing error just outside the gate starts a second track; a later IRA
 * between the two pulls the second onto the first and they are merged.
 */
BOOST_AUTO_TEST_CASE(test_tch_c_timing_merge)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    const auto ira_at = [&](int64_t n, int64_t err_ns) {
        const double t = p.ira_s(n);
        return table.add_ira(p.ns(t) + err_ns, F_RING + p.doppler(t));
    };
    BOOST_CHECK_EQUAL(ira_at(0, 0), 0);
    BOOST_CHECK_EQUAL(ira_at(1, 41000), 1); // outside the 40 us gate
    BOOST_CHECK_EQUAL(table.size(), 2u);
    ira_at(2,
           21000); // nearer to track 1: its median moves to 31 us, inside the gate of 0
    BOOST_CHECK_EQUAL(table.n_created(), 2u);
    BOOST_CHECK_EQUAL(table.n_merged(), 1u);
    BOOST_REQUIRE_EQUAL(table.size(), 1u);
    // the older track (and its frame numbering) is kept
    BOOST_CHECK_EQUAL(table.tracks()[0].id(), 0);
    BOOST_CHECK_EQUAL(table.tracks()[0].measurements().size(), 3u);
    for (int64_t n = 3; n < 30; n++) {
        BOOST_CHECK_EQUAL(ira_at(n, 0), 0);
    }
    BOOST_CHECK_EQUAL(table.size(), 1u);
}

// Expiry, stale, eviction beyond max_tracks
BOOST_AUTO_TEST_CASE(test_tch_c_timing_lifetime)
{
    tch_c_track_table table;
    const int64_t s = 1000000000LL;
    for (int i = 0; i < 10; i++) {
        // IRAs 1 ms apart in phase: ten different satellites
        table.add_ira(i * s + i * 1000000LL, F_RING + 1000 * i);
    }
    BOOST_CHECK_EQUAL(table.size(), 8u);
    BOOST_CHECK_EQUAL(table.n_evicted(), 2u);
    const tch_c_sat_track& t = table.tracks().back();
    BOOST_CHECK(!t.stale(12 * s));
    BOOST_CHECK(t.stale(13 * s));
    BOOST_CHECK_EQUAL(table.expire(25 * s), 3u);
    BOOST_CHECK_EQUAL(table.expire(40 * s), 5u);
    BOOST_CHECK_EQUAL(table.size(), 0u);
}

// Frame numbers and slot iteration
BOOST_AUTO_TEST_CASE(test_tch_c_timing_slots)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(8);
    for (const auto& x : make_iras(p, 0, 100, 5, nz, 0, 0)) {
        table.add_ira(x.t_ns, x.f_hz);
    }
    const tch_c_sat_track& tr = table.tracks().at(0);
    for (int type : { 1, 2, 8 }) {
        for (int64_t n = 60; n < 110; n++) { // inside the 4 s window
            for (int k = 0; k < N_SLOTS; k++) {
                const int64_t t = tr.slot_time_ns(n, k, type);
                BOOST_REQUIRE_LT(std::llabs(t - p.ns(p.slot_s(n, k, type))), 50);
                BOOST_REQUIRE_EQUAL(tr.frame_at(t), n);
                BOOST_REQUIRE_EQUAL(tr.frame_of_slot(t + 2000000, k, type), n);
                int64_t f;
                int sl;
                tr.next_slot(t, type, f, sl);
                BOOST_REQUIRE(f == n && sl == k);
                tr.next_slot(t + 1, type, f, sl);
                BOOST_REQUIRE(k < 3 ? f == n && sl == k + 1 : f == n + 1 && sl == 0);
            }
            // the IRA of the frame starts it
            BOOST_REQUIRE_EQUAL(tr.frame_at(p.ns(p.ira_s(n)) + 1000), n);
        }
    }
    const double t = p.ira_s(100);
    const auto pr = tr.predict(p.ns(t));
    BOOST_CHECK_LT(std::fabs(pr.phase_ms - frame_phase_ms(p.ns(t))) * 1e6, 20);
    BOOST_CHECK_LT(std::fabs(pr.doppler_hz - p.doppler(t)), 5); // curvature over 4 s
    BOOST_CHECK_LT(
        std::fabs(pr.doppler_rate - (p.doppler(t + 0.01) - p.doppler(t - 0.01)) / 0.02),
        10);
}

/*
 * Hits with fine timing feedback: a constant timing error of the type (here +3 us for
 * TCH_C2) is learned against the IRA phase and passed to new tracks.
 */
BOOST_AUTO_TEST_CASE(test_tch_c_timing_bias)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(9);
    auto iras = make_iras(p, 0, 1500, 13, nz);
    size_t next = 0;
    const double nominal = 1616e6 + 51 * 1e7 / 240;
    stats late, dt;
    for (int64_t n = 0; n < 1500; n++) {
        const int k = n % 4;
        const double ts = p.slot_s(n, k, 2);
        const int64_t t = p.ns(ts) + 3000 + int64_t(std::llround(300 * nz.gauss()));
        while (next < iras.size() && iras[next].t_ns <= t - 50000000LL) {
            table.add_ira(iras[next].t_ns, iras[next].f_hz);
            next++;
        }
        if (next < 4) {
            continue;
        }
        BOOST_REQUIRE_EQUAL(table.size(), 1u);
        const double f = nominal * (1 + p.doppler(ts) / F_RING) + 20 * nz.gauss();
        tch_c_track_table::slot_match m;
        BOOST_REQUIRE(table.match_hit(t, f, 2, nominal, m));
        BOOST_REQUIRE_EQUAL(m.slot, k);
        table.add_hit(m.track, t, f, 2, k, nominal, true); // fine timing feedback
        if (n > 1000) {
            const tch_c_sat_track& tr = table.tracks()[0];
            late.add(tr.bias_ns(2));
            // a TCH_C1 slot is unaffected by the TCH_C2 bias
            dt.add(double(tr.slot_time_ns(m.frame, k, 1) - p.ns(p.slot_s(n, k, 1))));
        }
    }
    BOOST_TEST_MESSAGE("bias: " << late.mean() << " ns, std " << late.std()
                                << "; TCH_C1 dt mean " << dt.mean() << " ns");
    BOOST_CHECK_LT(std::fabs(late.mean() - 3000), 300);
    BOOST_CHECK_LT(std::fabs(table.bias_ns(2) - 3000), 500);
    BOOST_CHECK_EQUAL(table.bias_ns(1), 0);
    BOOST_CHECK_LT(std::fabs(dt.mean()), 300);
    // a new track starts with the learned bias
    table.add_ira(p.ns(p.ira_s(1500)) + 10000000LL, F_RING);
    BOOST_REQUIRE_EQUAL(table.size(), 2u);
    BOOST_CHECK_CLOSE(table.tracks()[1].bias_ns(2), table.bias_ns(2), 1e-9);
}

/*
 * The bias new tracks start from, over rebuilds of all tracks (a gap longer than the
 * expiry, then every satellite at once, as after a stall). The TCH_C1 bursts are a few
 * us off the table offset, differently per satellite, and a slot is only found within
 * 17 us (+-0.5 symbol) of its prediction, as by the survey. Summing the bias updates of
 * all tracks, every rebuild multiplied the error: after one or two the start bias was
 * beyond 17 us, no new track found a TCH_C1 burst and none could correct it.
 */
BOOST_AUTO_TEST_CASE(test_tch_c_timing_bias_rebuild)
{
    const double bias_us[] = { 4.0, 3.5, 3.0, 4.0, 2.5, 3.5 };
    const int n_sats = 6;
    std::vector<sat_pass> passes;
    for (int i = 0; i < n_sats; i++) {
        passes.push_back({ 100.0 + 120 * i, R_MIN + 150e3 * i, V, 30.0 + 0.7 * i, 0 });
    }
    tch_c_track_table table;
    noise nz(17);
    const double nominal = 1616e6 + 100.5 * 1e7 / 240;
    const int64_t round_frames = 1500, gap_frames = 800; // gap 72 s: longer than expiry_s
    for (int round = 0; round < 5; round++) {
        const int64_t first = round * (round_frames + gap_frames);
        const int64_t last = first + round_frames;
        std::vector<ira> iras;
        for (const auto& p : passes) {
            auto v = make_iras(p, first, last, 13, nz);
            iras.insert(iras.end(), v.begin(), v.end());
        }
        std::sort(iras.begin(), iras.end(), [](const ira& a, const ira& b) {
            return a.t_ns < b.t_ns;
        });
        size_t next = 0;
        uint64_t tried = 0, found = 0;
        for (int64_t n = first; n <= last; n++) {
            const int k = n % N_SLOTS;
            for (int s = 0; s < n_sats; s++) {
                const sat_pass& p = passes[s];
                const double ts = p.slot_s(n, k, 1);
                const int64_t t_true = p.ns(ts) + int64_t(bias_us[s] * 1e3);
                table.expire(t_true);
                while (next < iras.size() && iras[next].t_ns <= t_true - 50000000LL) {
                    table.add_ira(iras[next].t_ns, iras[next].f_hz);
                    next++;
                }
                if (n < first + 100) { // warm up: every satellite has a track
                    continue;
                }
                const double f = nominal * (1 + p.doppler(ts) / F_RING);
                tch_c_track_table::slot_match m;
                if (!table.match_hit(t_true, f, 1, nominal, m)) {
                    continue;
                }
                tried++;
                const int64_t t_pred = table.find(m.track)->slot_time_ns(m.frame, k, 1);
                if (std::llabs(t_true - t_pred) >= 17000) {
                    continue; // not found by the survey
                }
                found++;
                const int64_t t_meas = t_true + int64_t(std::llround(1000 * nz.gauss()));
                table.add_hit(m.track, t_meas, f, 1, k, nominal, true);
            }
        }
        BOOST_TEST_MESSAGE("rebuild " << round << ": TCH_C1 found " << found << " of "
                                      << tried << ", start bias " << table.bias_ns(1)
                                      << " ns, " << table.n_created()
                                      << " tracks created");
        BOOST_CHECK_GT(tried, 0.9 * n_sats * (round_frames - 100));
        BOOST_CHECK_GT(found, 0.95 * tried);
        BOOST_CHECK_LT(std::fabs(table.bias_ns(1) - 3500), 1500);
    }
    // every round rebuilt all tracks
    BOOST_CHECK_EQUAL(table.n_created(), 5u * n_sats);
}

// A new track never starts beyond the bias limit
BOOST_AUTO_TEST_CASE(test_tch_c_timing_bias_start_clamped)
{
    const double far[3] = { 1e6, -1e6, 0 };
    tch_c_timing_config cfg;
    tch_c_sat_track tr(0, cfg, far);
    BOOST_CHECK_EQUAL(tr.bias_ns(1), cfg.max_bias_ns);
    BOOST_CHECK_EQUAL(tr.bias_ns(2), -cfg.max_bias_ns);
    BOOST_CHECK_EQUAL(tr.bias_ns(8), 0);
}

// Without IRAs (after the first 2 s), energy-path hits associated by the gate carry the
// track
BOOST_AUTO_TEST_CASE(test_tch_c_timing_hits_only)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(10);
    for (const auto& x : make_iras(p, 0, 22, 4, nz)) {
        table.add_ira(x.t_ns, x.f_hz);
    }
    BOOST_REQUIRE_EQUAL(table.size(), 1u);
    const int id = table.tracks()[0].id();
    const double nominal = 1616e6 + 151.5 * 1e7 / 240;
    stats dt, df;
    for (int64_t n = 23; n < 1500; n++) {
        const int k = (n * 7) % 4;
        const double ts = p.slot_s(n, k, 1);
        const int64_t t_true = p.ns(ts);
        const double f_true = nominal * (1 + p.doppler(ts) / F_RING);
        const tch_c_sat_track& tr = table.tracks().at(0);
        const int64_t t_pred = tr.slot_time_ns(tr.frame_of_slot(t_true, k, 1), k, 1);
        dt.add(double(t_pred - t_true));
        df.add(tr.carrier_hz(nominal, t_pred) - f_true);
        const int64_t t = t_true + int64_t(std::llround(500 * nz.gauss()));
        BOOST_REQUIRE_EQUAL(table.associate_hit(t, f_true + 20 * nz.gauss(), 1, nominal),
                            id);
        table.expire(t);
    }
    BOOST_TEST_MESSAGE("hits only: dt std " << dt.std() << " max " << dt.max_abs
                                            << " ns; df std " << df.std() << " Hz");
    BOOST_CHECK_LT(dt.max_abs, 5000.0);
    BOOST_CHECK_LT(df.std(), 30.0);
    BOOST_CHECK_EQUAL(table.size(), 1u);
}

// Uplink: a terminal next to the receiver sends so its bursts arrive at the satellite
// at fixed times of the satellite frame. As received v = (t + phase(t)) mod 90 ms is
// constant over the pass, and ul_time_ns() / ul_frame_at() invert it.
BOOST_AUTO_TEST_CASE(test_tch_c_timing_uplink)
{
    const sat_pass p = { 60, R_MIN, V, 30.0, 0 };
    tch_c_track_table table;
    noise nz(1);
    auto iras = make_iras(p, 0, 1500, 13, nz, 0, 0);
    size_t next = 0;

    const double at_sat_ms = 12.0; // UL slot 0 at the satellite, relative to phi0
    double v0 = -1, max_dv = 0, max_dt = 0;
    for (int64_t n = 100; n < 1400; n += 7) {
        // Sent so it arrives at the satellite at n * 90 ms + phi0 + at_sat_ms
        const double at_sat = n * 0.09 + (p.phi0 + at_sat_ms) * 1e-3;
        double t = at_sat;
        for (int i = 0; i < 4; i++) {
            t = at_sat - p.range(t) / C;
        }
        const int64_t t_ns = p.ns(t);
        // The IRAs up to then
        while (next < iras.size() && iras[next].t_ns <= t_ns) {
            table.add_ira(iras[next].t_ns, iras[next].f_hz);
            next++;
        }
        BOOST_REQUIRE_EQUAL(table.size(), 1u);
        const tch_c_sat_track& tr = *table.begin();
        const double v = tr.ul_phase_ms(t_ns);
        if (v0 < 0) {
            v0 = v;
        }
        max_dv = std::max(max_dv, std::fabs(wrap_ms(v - v0)));
        const int64_t back = tr.ul_time_ns(tr.ul_frame_at(t_ns, v0), v0);
        max_dt = std::max(max_dt, std::fabs((double)(back - t_ns)));
    }
    BOOST_TEST_MESSAGE("uplink: v " << v0 << " ms, max dv " << max_dv * 1e3
                                    << " us, max dt " << max_dt << " ns");
    BOOST_CHECK_LT(max_dv, 0.002);
    BOOST_CHECK_LT(max_dt, 2000.0);
}
