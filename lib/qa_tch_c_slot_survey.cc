/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "qa_tch_c_encoder.h"
#include "tch_c_format.h"
#include "tch_c_slot_survey.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace gr::iridium;

namespace {

const double CENTER = 1622e6;
const double DOPPLER = 23456; // at the ring channel
const double SLOT = 8.38e-3;

// Gaussian numbers by Box-Muller from std::mt19937 (the same on every platform)
struct gauss {
    std::mt19937 rng;
    explicit gauss(uint32_t seed) : rng(seed) {}
    double uniform() { return (rng() + 0.5) / 4294967296.0; }
    gr_complex operator()()
    {
        const double r = std::sqrt(-2 * std::log(uniform()));
        const double a = 2 * M_PI * uniform();
        return gr_complex(r * std::cos(a) * M_SQRT1_2, r * std::sin(a) * M_SQRT1_2);
    }
};

double rrc(double t)
{
    const double a = tch_c_slot_survey::RRC_ALPHA;
    if (std::fabs(t) < 1e-9) {
        return 1 - a + 4 * a / M_PI;
    }
    if (std::fabs(std::fabs(4 * a * t) - 1) < 1e-9) {
        return a / M_SQRT2 *
               ((1 + 2 / M_PI) * std::sin(M_PI / (4 * a)) +
                (1 - 2 / M_PI) * std::cos(M_PI / (4 * a)));
    }
    return (std::sin(M_PI * t * (1 - a)) + 4 * a * t * std::cos(M_PI * t * (1 + a))) /
           (M_PI * t * (1 - 16 * a * a * t * t));
}

const gr_complex S0(M_SQRT1_2, M_SQRT1_2);

// The symbols of a QPSK TCH_C burst with a random payload
std::vector<gr_complex> burst_symbols(const tch_c_format& f,
                                      bool uplink,
                                      uint32_t seed,
                                      std::vector<uint8_t>* payload)
{
    std::mt19937 rng(seed);
    payload->assign(f.n_codewords * f.code.info_bits - 10, 0);
    for (auto& b : *payload) {
        b = rng() >> 31;
    }
    const std::vector<uint8_t> bits = qa_tch_c_encoder(f).encode_burst(*payload);

    std::vector<gr_complex> s(f.n_symbols);
    const int* uw = uplink ? TCH_C_UW_UL : TCH_C_UW_DL;
    for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
        s[k] = uw[k] ? -S0 : S0;
    }
    for (int p = 0; p < f.n_pilots; p++) {
        for (int k = 0; k < TCH_C_PILOT_LENGTH; k++) {
            s[f.pilot_start[p] + k] = TCH_C_PILOT[k] ? -S0 : S0;
        }
    }
    int i = 0;
    for (int g = 0; g < f.n_segments; g++) {
        for (int k = f.segment[g][0]; k < f.segment[g][0] + f.segment[g][1]; k++) {
            s[k] = gr_complex(bits[i] ? -M_SQRT1_2 : M_SQRT1_2,
                              bits[i + 1] ? -M_SQRT1_2 : M_SQRT1_2);
            i += 2;
        }
    }
    return s;
}

struct burst {
    const tch_c_format* f;
    int position;
    double p_uw;      // sample index of the centre of UW symbol 0
    double frequency; // Hz
    std::vector<uint8_t> payload;
};

/*
 * Add a burst: RRC pulses (+-6 symbols) at p_uw + k sps, amplitude a (per pulse peak
 * ~1), carrier frequency, random phase
 */
void add_burst(std::vector<gr_complex>& x,
               double fs,
               const std::vector<gr_complex>& s,
               double sps,
               double p_uw,
               double frequency,
               double a,
               double phase)
{
    const double w = 2 * M_PI * (frequency - CENTER) / fs;
    const long n0 = std::max(0L, (long)std::floor(p_uw - 6 * sps));
    const long n1 =
        std::min((long)x.size() - 1, (long)std::ceil(p_uw + (s.size() + 5) * sps));
    for (long n = n0; n <= n1; n++) {
        const double t = (n - p_uw) / sps;
        const int k0 = std::max(0, (int)std::ceil(t - 6));
        const int k1 = std::min((int)s.size() - 1, (int)std::floor(t + 6));
        gr_complex v = 0;
        for (int k = k0; k <= k1; k++) {
            v += s[k] * (float)rrc(t - k);
        }
        x[n] += v * gr_complex(std::polar(a, w * n + phase));
    }
}

/*
 * Stand-in for the soft demodulation of tch_c_demod (which is internal to that block):
 * frequency and phase from the known symbols (grid search, then the phase at symbol 0),
 * gain and noise from the known symbols, QPSK LLRs, decode, descramble, CRC.
 */
bool decode(const tch_c_format& f,
            const gr_complex* y,
            bool uplink,
            std::vector<uint8_t>& msg,
            float* residual_hz)
{
    const int* uw = uplink ? TCH_C_UW_UL : TCH_C_UW_DL;
    std::vector<int> pos;
    std::vector<gr_complex> z;
    for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
        pos.push_back(k);
        z.push_back(y[k] * std::conj(uw[k] ? -S0 : S0));
    }
    for (int p = 0; p < f.n_pilots; p++) {
        for (int k = 0; k < TCH_C_PILOT_LENGTH; k++) {
            pos.push_back(f.pilot_start[p] + k);
            z.push_back(y[f.pilot_start[p] + k] * std::conj(TCH_C_PILOT[k] ? -S0 : S0));
        }
    }
    double best_w = 0, best = -1;
    for (int g = -300; g <= 300; g++) {
        const double w = 2 * M_PI * g * 5.0 / f.symbols_per_second;
        gr_complex acc = 0;
        for (size_t i = 0; i < z.size(); i++) {
            acc += z[i] * gr_complex(std::polar(1.0, -w * pos[i]));
        }
        if (std::norm(acc) > best) {
            best = std::norm(acc);
            best_w = w;
        }
    }
    *residual_hz = best_w / (2 * M_PI) * f.symbols_per_second;
    gr_complex h = 0;
    for (size_t i = 0; i < z.size(); i++) {
        h += z[i] * gr_complex(std::polar(1.0, -best_w * pos[i]));
    }
    const double phi0 = std::arg(h);
    std::vector<gr_complex> c(f.n_symbols);
    for (int k = 0; k < f.n_symbols; k++) {
        c[k] = y[k] * gr_complex(std::polar(1.0, -(phi0 + best_w * k)));
    }

    float gain = 0;
    for (size_t i = 0; i < pos.size(); i++) {
        gain += (z[i] * gr_complex(std::polar(1.0, -(phi0 + best_w * pos[i])))).real();
    }
    gain /= pos.size();
    float sigma2 = 0;
    for (size_t i = 0; i < pos.size(); i++) {
        const gr_complex zi =
            z[i] * gr_complex(std::polar(1.0, -(phi0 + best_w * pos[i]))) / gain;
        sigma2 += std::norm(zi - gr_complex(1, 0));
    }
    sigma2 = std::max(sigma2 / pos.size(), 1e-3f);

    std::vector<float> llr;
    const float scale = 2 * M_SQRT2 / sigma2 / gain;
    for (int g = 0; g < f.n_segments; g++) {
        for (int k = f.segment[g][0]; k < f.segment[g][0] + f.segment[g][1]; k++) {
            llr.push_back(scale * c[k].real());
            llr.push_back(scale * c[k].imag());
        }
    }

    tch_c_decoder decoder(f.code);
    const int K = f.code.info_bits;
    const int L = f.code.coded_bits;
    msg.assign(f.n_codewords * K, 0);
    for (int j = 0; j < f.n_codewords; j++) {
        decoder.decode_raw(&llr[j * L], &msg[j * K]);
    }
    tch_c_decoder::descramble(&msg[0], f.n_codewords * K, 0);
    return tch_c_decoder::check_crc(&msg[0], f.n_codewords * K);
}

const tch_c_format* const TYPES[] = { &TCH_C1, &TCH_C2, &TCH_C8Q };

// Raster positions used, inside the band at 10 MS/s too
const int POSITIONS[3][5] = { { 30, 100, 101, 150, 230 },
                              { 20, 60, 61, 110, -1 },
                              { 5, 15, 16, 25, -1 } };

// Timing offsets of the bursts of a slot from the slot time, in symbols
const double TIMING[5] = { 0, 0.1, -0.15, 0.2, -0.05 };

/*
 * A wideband signal at fs: in slot i (UW at p_slot[i]) bursts of TYPES[i] on
 * POSITIONS[i] with the Doppler, a timing offset of TIMING[] symbols, Es/N0 esn0_db
 * (no noise if esn0_db >= 100).
 */
std::vector<gr_complex> make_signal(double fs,
                                    double esn0_db,
                                    uint32_t seed,
                                    std::vector<burst>& bursts,
                                    double p_slot[3])
{
    const double first = 0.002 * fs + 0.37;
    for (int i = 0; i < 3; i++) {
        p_slot[i] = first + i * SLOT * fs;
    }
    std::vector<gr_complex> x((size_t)(first + 3 * SLOT * fs + 0.002 * fs), 0);

    std::mt19937 rng(seed);
    const double a = 0.01;
    for (int i = 0; i < 3; i++) {
        const tch_c_format& f = *TYPES[i];
        const double sps = fs / f.symbols_per_second;
        for (int j = 0; j < 5 && POSITIONS[i][j] >= 0; j++) {
            burst b;
            b.f = &f;
            b.position = POSITIONS[i][j];
            b.p_uw = p_slot[i] + TIMING[j] * sps;
            b.frequency = tch_c_slot_survey::nominal_frequency(f.type, b.position) *
                              (1 + DOPPLER / tch_c_slot_survey::F_RING) +
                          ((int)(rng() % 200) - 100);
            const std::vector<gr_complex> s = burst_symbols(f, false, rng(), &b.payload);
            // The same Es for all types, amplitudes within 3 dB
            const double ab = a * std::sqrt(fs / TCH_C1.symbols_per_second / sps) *
                              std::pow(10, (rng() % 300) / 2000.0);
            add_burst(x, fs, s, sps, b.p_uw, b.frequency, ab, (rng() % 6283) / 1000.0);
            bursts.push_back(b);
        }
    }

    if (esn0_db < 100) {
        // Es/N0 of the weakest burst: Es = a^2 sps (the pulses have unit energy per
        // symbol), noise variance per sample sigma^2
        const double sps1 = fs / TCH_C1.symbols_per_second;
        const double sigma = a * std::sqrt(sps1 / std::pow(10, esn0_db / 10));
        gauss g(seed + 1);
        for (auto& v : x) {
            v += (float)sigma * g();
        }
    }
    return x;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_tch_c_slot_survey_geometry)
{
    for (double fs : { 12e6, 10e6 }) {
        tch_c_slot_survey s(fs, CENTER);
        for (const tch_c_format* f : TYPES) {
            BOOST_TEST_CONTEXT("fs " << fs << " TCH_C" << f->type)
            {
                const double sps = fs / f->symbols_per_second;
                BOOST_CHECK_CLOSE(s.samples_per_symbol(f->type), sps, 1e-9);
                BOOST_CHECK_EQUAL(s.uw_offset(f->type), (int)std::ceil(6 * sps));
                BOOST_CHECK_EQUAL(s.window_length(f->type),
                                  (int)std::floor(s.uw_offset(f->type) + 37 * sps) + 1);
                BOOST_CHECK(s.fft_size(f->type) >= s.window_length(f->type));
                BOOST_CHECK(s.fft_size(f->type) < 2 * s.window_length(f->type));
            }
        }
        // The window sizes of ARCH.md at 12 MS/s
        if (fs == 12e6) {
            BOOST_CHECK_EQUAL(s.window_length(1), 17201);
            BOOST_CHECK_EQUAL(s.fft_size(1), 32768);
            BOOST_CHECK_EQUAL(s.window_length(2), 8601);
            BOOST_CHECK_EQUAL(s.window_length(8), 2151);
            BOOST_CHECK_EQUAL(s.fft_size(8), 4096);
        }
    }
    BOOST_CHECK_CLOSE(
        tch_c_slot_survey::nominal_frequency(1, 0), 1616e6 + 20833.33, 1e-6);
    BOOST_CHECK_CLOSE(
        tch_c_slot_survey::nominal_frequency(2, 119), 1616e6 + 239 * 41666.667, 1e-6);
    BOOST_CHECK_CLOSE(
        tch_c_slot_survey::nominal_frequency(8, 29), 1616e6 + 236 * 41666.667, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_tch_c_slot_survey_noise)
{
    /*
     * Noise only. rho2 of one frequency cell is about exponentially distributed with a
     * mean of ~1/37 (the energy window spans ~32.5 symbols); the peak search over the
     * bins (+-1 bin and the interpolation) about doubles that. The tail: rarely above
     * 0.25, never above 0.35.
     */
    for (double fs : { 12e6, 10e6 }) {
        tch_c_slot_survey s(fs, CENTER);
        gauss g(7);
        for (const tch_c_format* f : TYPES) {
            BOOST_TEST_CONTEXT("fs " << fs << " TCH_C" << f->type)
            {
                const int L = s.window_length(f->type);
                std::vector<gr_complex> x(L);
                double sum = 0, max = 0;
                int n = 0, n25 = 0;
                const int rounds = f->type == 1 ? 40 : f->type == 2 ? 80 : 300;
                for (int r = 0; r < rounds; r++) {
                    for (auto& v : x) {
                        v = g();
                    }
                    std::vector<tch_c_slot_survey::candidate> c;
                    s.survey(&x[0], f->type, false, DOPPLER, 0, -1, c);
                    for (const auto& ci : c) {
                        sum += ci.rho2;
                        max = std::max(max, (double)ci.rho2);
                        n25 += ci.rho2 >= 0.25;
                        n++;
                    }
                }
                BOOST_TEST_MESSAGE("noise fs " << fs << " TCH_C" << f->type << ": " << n
                                               << " hypotheses, mean rho2 " << sum / n
                                               << ", >= 0.25: " << n25 << ", max "
                                               << max);
                BOOST_CHECK(n > 5000);
                BOOST_CHECK(sum / n > 1 / 25.);
                BOOST_CHECK(sum / n < 1 / 14.);
                BOOST_CHECK(n25 < n / 1000);
                BOOST_CHECK(max < 0.35);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_slot_survey_clean)
{
    // Clean bursts (Es/N0 40 dB): rho2 > 0.9 at their own time; frequency within a
    // fraction of a bin
    for (double fs : { 12e6, 10e6 }) {
        std::vector<burst> bursts;
        double p_slot[3];
        const std::vector<gr_complex> x = make_signal(fs, 40, 1, bursts, p_slot);
        tch_c_slot_survey s(fs, CENTER);
        for (const burst& b : bursts) {
            const int type = b.f->type;
            BOOST_TEST_CONTEXT("fs " << fs << " TCH_C" << type << " position "
                                     << b.position)
            {
                const long i0 = std::lround(b.p_uw) - s.uw_offset(type);
                std::vector<tch_c_slot_survey::candidate> c;
                s.survey(&x[i0], type, false, DOPPLER, 0, 0.5, c);
                const tch_c_slot_survey::candidate* hit = nullptr;
                for (const auto& ci : c) {
                    if (ci.position == b.position) {
                        hit = &ci;
                    }
                }
                BOOST_REQUIRE(hit);
                BOOST_TEST_MESSAGE("clean fs "
                                   << fs << " TCH_C" << type << " pos " << b.position
                                   << ": rho2 " << hit->rho2 << " df "
                                   << hit->estimated_frequency - b.frequency
                                   << " Hz (bin " << s.bin_width(type) << ")");
                BOOST_CHECK(hit->rho2 > 0.9);
                BOOST_CHECK(std::fabs(hit->estimated_frequency - b.frequency) <
                            0.25 * s.bin_width(type));
                // No candidates on the empty positions of the type
                for (const auto& ci : c) {
                    bool used = false;
                    for (const burst& o : bursts) {
                        used |= o.f->type == type && o.position == ci.position;
                    }
                    BOOST_CHECK(used);
                }
                // The UW correlation of the other direction is low
                std::vector<tch_c_slot_survey::candidate> cu;
                s.survey(&x[i0], type, true, DOPPLER, 0, 0.35, cu);
                for (const auto& ci : cu) {
                    BOOST_TEST_MESSAGE("UL pos " << ci.position << " rho2 " << ci.rho2);
                }
                BOOST_CHECK(cu.empty());

                // Fine timing finds the fractional position
                float rho2_fine;
                const double tau = s.fine_timing(&x[i0], type, false, *hit, &rho2_fine);
                const double err = i0 + s.uw_offset(type) + tau - b.p_uw;
                BOOST_TEST_MESSAGE("clean timing error " << err << " samples, rho2 "
                                                         << rho2_fine);
                // Within 1% of a symbol (the data after the UW shifts the peak a bit)
                BOOST_CHECK(std::fabs(err) <
                            std::max(1., 0.01 * s.samples_per_symbol(type)));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_slot_survey_decode)
{
    // Noisy slots surveyed at the slot time: every burst is found at its position, the
    // fine timing finds its own offset, and the extracted symbols decode
    for (double fs : { 12e6, 10e6 }) {
        std::vector<burst> bursts;
        double p_slot[3];
        const std::vector<gr_complex> x = make_signal(fs, 9, 2, bursts, p_slot);
        tch_c_slot_survey s(fs, CENTER);
        for (int i = 0; i < 3; i++) {
            const int type = TYPES[i]->type;
            const double sps = s.samples_per_symbol(type);
            const long i0 = std::lround(p_slot[i]) - s.uw_offset(type);
            std::vector<tch_c_slot_survey::candidate> c;
            s.survey(&x[i0], type, false, DOPPLER, 0, 0.35, c);

            for (const burst& b : bursts) {
                if (b.f->type != type) {
                    continue;
                }
                BOOST_TEST_CONTEXT("fs " << fs << " TCH_C" << type << " position "
                                         << b.position)
                {
                    const tch_c_slot_survey::candidate* hit = nullptr;
                    for (const auto& ci : c) {
                        if (ci.position == b.position) {
                            hit = &ci;
                        }
                    }
                    BOOST_REQUIRE(hit);
                    float rho2_fine;
                    const double tau =
                        s.fine_timing(&x[i0], type, false, *hit, &rho2_fine);
                    const double p_uw = i0 + s.uw_offset(type) + tau;
                    const double timing_error = p_uw - b.p_uw;

                    std::vector<gr_complex> y(b.f->n_symbols + 1);
                    long first, end;
                    s.extract_span(type, p_uw, &first, &end);
                    BOOST_REQUIRE(first >= 0 && end <= (long)x.size());
                    BOOST_CHECK_EQUAL(
                        s.extract(
                            &x[0], p_uw, type, hit->estimated_frequency, false, &y[0]),
                        b.f->n_symbols + 1);
                    std::vector<uint8_t> msg;
                    float residual;
                    const bool ok = decode(*b.f, &y[0], false, msg, &residual);

                    BOOST_TEST_MESSAGE("decode fs "
                                       << fs << " TCH_C" << type << " pos " << b.position
                                       << ": rho2 " << hit->rho2 << " fine " << rho2_fine
                                       << " timing error " << timing_error << " samples ("
                                       << timing_error / sps << " symbols), df "
                                       << hit->estimated_frequency - b.frequency
                                       << " Hz, fit " << residual << " Hz, CRC " << ok);
                    BOOST_CHECK(std::fabs(hit->estimated_frequency - b.frequency) <
                                0.3 * s.bin_width(type));
                    BOOST_CHECK(std::fabs(timing_error) < 0.05 * sps);
                    BOOST_CHECK(std::fabs(residual + hit->estimated_frequency -
                                          b.frequency) < 15);
                    BOOST_CHECK(ok);
                    BOOST_CHECK(
                        std::equal(b.payload.begin(), b.payload.end(), msg.begin()));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_slot_survey_extract_exact)
{
    // Clean single burst at a fractional position and sps: the extracted samples are
    // the symbols (up to the RRC truncation and the neighbouring pulses), phase 0
    for (double fs : { 12e6, 10e6, 7.68e6 }) {
        tch_c_slot_survey s(fs, CENTER);
        for (const tch_c_format* f : TYPES) {
            BOOST_TEST_CONTEXT("fs " << fs << " TCH_C" << f->type)
            {
                const double sps = fs / f->symbols_per_second;
                std::vector<uint8_t> payload;
                const std::vector<gr_complex> sym = burst_symbols(*f, true, 5, &payload);
                const double p_uw = 20 * sps + 0.731;
                const double freq = CENTER + 1234567.8;
                std::vector<gr_complex> x((size_t)((f->n_symbols + 40) * sps), 0);
                add_burst(x, fs, sym, sps, p_uw, freq, 1, 0.3);
                std::vector<gr_complex> y(f->n_symbols + 1);
                s.extract(&x[0], p_uw, f->type, freq, true, &y[0]);
                // Gain of RRC * RRC with taps summing to 1
                gr_complex g = 0;
                for (int k = 0; k < f->n_symbols; k++) {
                    g += y[k] * std::conj(sym[k]);
                }
                g /= f->n_symbols;
                double err = 0;
                for (int k = 0; k < f->n_symbols; k++) {
                    err += std::norm(y[k] / g - sym[k]);
                }
                err /= f->n_symbols;
                BOOST_TEST_MESSAGE("extract fs " << fs << " TCH_C" << f->type << ": gain "
                                                 << g << " EVM " << 10 * std::log10(err)
                                                 << " dB");
                BOOST_CHECK(std::fabs(std::arg(g)) < 0.01);
                // RRC truncation (+-2.5 symbols)
                BOOST_CHECK(err < 3e-3);
            }
        }
    }
}
