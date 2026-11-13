/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tch_c_demod_impl.h"
#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>
#include <algorithm>
#include <cstring>

namespace gr {
namespace iridium {

tch_c_demod::sptr tch_c_demod::make()
{
    return gnuradio::make_block_sptr<tch_c_demod_impl>();
}

tch_c_demod_impl::tch_c_demod_impl()
    : gr::block("tch_c_demod",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0)),
      d_alpha(1 / 5.),
      d_n_handled_bursts(0),
      d_n_access_ok_bursts(0),
      d_n_access_ok_sub_bursts(0),
      d_channel_id(1, UINT64_MAX),
      d_c1_decoder(TCH_C1.code),
      d_c2_decoder(TCH_C2.code),
      d_c8q_decoder(TCH_C8Q.code),
      d_c8a_decoder(TCH_C8A.code)
{
    message_port_register_out(pmt::mp("pdus"));
    d_upstream_done.install(this);

    auto port_name = pmt::mp("probes");
    message_port_register_in(port_name);
    d_input_ports.push_back(port_name);
    set_msg_handler(port_name, [this](pmt::pmt_t msg) { this->handler(0, msg); });
    d_frame_handled_ports.push_back(pmt::mp("probe_handled"));
    for (const auto& port : d_frame_handled_ports) {
        message_port_register_out(port);
    }
}

tch_c_demod_impl::~tch_c_demod_impl() {}

uint64_t tch_c_demod_impl::get_n_handled_bursts() { return d_n_handled_bursts; }

uint64_t tch_c_demod_impl::get_n_access_ok_bursts() { return d_n_access_ok_bursts; }

uint64_t tch_c_demod_impl::get_n_access_ok_sub_bursts()
{
    return d_n_access_ok_sub_bursts;
}

/*
 * Decision directed first order PLL, as in iridium_qpsk_demod
 * (taken from synchronizer_v4_impl.cc of gr-burst)
 */
float tch_c_demod_impl::qpsk_first_order_pll(const gr_complex* x,
                                             int size,
                                             float alpha,
                                             gr_complex* y)
{
    gr_complex phiHat = gr_complex(1, 0);
    float total_phase = 0;
    gr_complex xHat, er, phiHatT;
    for (int ii = 0; ii < size; ii++) {
        // correct w/ estimated phase
        y[ii] = x[ii] * phiHat;

        // demodulating circuit
        if (y[ii].real() >= 0 && y[ii].imag() >= 0) {
            xHat.real(M_SQRT1_2);
            xHat.imag(M_SQRT1_2);
        } else if (y[ii].real() >= 0 && y[ii].imag() < 0) {
            xHat.real(M_SQRT1_2);
            xHat.imag(-M_SQRT1_2);
        } else if (y[ii].real() < 0 && y[ii].imag() < 0) {
            xHat.real(-M_SQRT1_2);
            xHat.imag(-M_SQRT1_2);
        } else {
            xHat.real(-M_SQRT1_2);
            xHat.imag(M_SQRT1_2);
        }

        // loop filter to update phase estimate
        er = std::conj(xHat) * y[ii];
        phiHatT = er / std::abs(er);
        total_phase += std::arg(std::pow(phiHatT, alpha));
        phiHat = std::conj(std::pow(phiHatT, alpha)) * phiHat;
        phiHat = phiHat / std::abs(phiHat);
    }
    return total_phase;
}

/*
 * Hard decisions: quadrants counted counter-clockwise from (+I,+Q). Level is the mean
 * magnitude, confidence the share of symbols within 22 degrees of the diagonals.
 */
void tch_c_demod_impl::demod_qpsk(
    const gr_complex* burst, size_t n_symbols, int* out, float* level, int* confidence)
{
    float sum = 0;
    int n_ok = 0;

    for (size_t i = 0; i < n_symbols; i++) {
        sum += std::abs(burst[i]);

        if (burst[i].real() >= 0 && burst[i].imag() >= 0) {
            out[i] = 0;
        } else if (burst[i].real() >= 0 && burst[i].imag() < 0) {
            out[i] = 3;
        } else if (burst[i].real() < 0 && burst[i].imag() < 0) {
            out[i] = 2;
        } else {
            out[i] = 1;
        }

        int phase = (gr::fast_atan2f(burst[i]) + M_PI) * 180 / M_PI;
        int offset = 45 - (phase % 90);
        if (abs(offset) <= 22) {
            n_ok++;
        }
    }

    *level = sum / n_symbols;
    *confidence = (int)(100. * n_ok / n_symbols);
}

// The 32 symbol UW with at most 5 symbol errors (scaled from 2 of 12 for legacy bursts)
bool tch_c_demod_impl::uw_ok(const int* demodulated_burst, const int* uw)
{
    int diffs = 0;
    for (int i = 0; i < TCH_C_UW_LENGTH; i++) {
        int diff = abs(demodulated_burst[i] - uw[i]);
        if (diff == 3) {
            diff = 1;
        }
        diffs += diff;
    }
    // At 3 dB Es/N0 (where TCH_C2/C8 still decode) 5 of the 32 UW symbols are wrong on
    // average; the prober has checked the UW coherently already, and the CRC protects
    // against false frames
    return diffs <= 10;
}

void tch_c_demod_impl::decode_gray_qpsk(int* demodulated_burst, size_t n_symbols)
{
    // Coherent Gray coded QPSK: symbol k carries coded bits (c[2k], c[2k+1]),
    // c[2k] on I and c[2k+1] on Q, bit 0 -> positive, bit 1 -> negative axis.
    // demod_qpsk numbers the quadrants counter-clockwise starting at (+I,+Q):
    // 0: (+I,+Q) -> 00, 1: (-I,+Q) -> 10, 2: (-I,-Q) -> 11, 3: (+I,-Q) -> 01
    static const int gray[4] = { 0, 2, 3, 1 };

    for (size_t i = 0; i < n_symbols; i++) {
        demodulated_burst[i] = gray[demodulated_burst[i]];
    }
}

void tch_c_demod_impl::map_symbols_to_bits(const int* demodulated_burst,
                                           size_t n_symbols,
                                           std::vector<uint8_t>& bits)
{
    bits.clear();
    for (size_t i = 0; i < n_symbols; i++) {
        bits.push_back((demodulated_burst[i] >> 1) & 1);
        bits.push_back(demodulated_burst[i] & 1);
    }
}

/*
 * Positions and values of the known symbols of a TCH_C burst: the 32 symbol UW and the
 * 16 symbol pilots. Known symbols are +-(1+j)/sqrt(2); sign is -1 for (-1-j)/sqrt(2).
 * Returns the number of known symbols.
 */
static int
tch_c_known_symbols(const tch_c_format& f, const int* uw, int* pos, float* sign)
{
    int n = 0;
    for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
        pos[n] = k;
        sign[n++] = uw[k] ? -1 : 1;
    }
    for (int p = 0; p < f.n_pilots; p++) {
        for (int k = 0; k < TCH_C_PILOT_LENGTH; k++) {
            pos[n] = f.pilot_start[p] + k;
            sign[n++] = TCH_C_PILOT[k] ? -1 : 1;
        }
    }
    return n;
}

static const int TCH_C_MAX_KNOWN = 32 + TCH_C_MAX_PILOTS * TCH_C_PILOT_LENGTH;

/*
 * 16APSK bursts (three codewords) are only decoded up to this noise variance on the
 * known symbols: every 16APSK burst passing the CRC on test recordings had less than 0.2
 * (rate 2/3 needs about 0.11), except three which were not 16APSK bursts at all but
 * failed decodes passing the CRC-10 by chance. Most TCH_C8 bursts with a valid UW which
 * fail the CRC are above this.
 */
static const float APSK16_MAX_NOISE = 0.25f;

/*
 * Phase and frequency correction for TCH_C bursts using the known symbols (UW and
 * pilots). The phase is modelled as phi0 + omega * k over the burst.
 * Returns omega in radians per symbol.
 */
float tch_c_demod_impl::known_symbol_correction(const tch_c_format& f,
                                                const gr_complex* x,
                                                gr_complex* y)
{
    // The UW tells uplink from downlink: use the one which correlates better
    gr_complex acc_dl = 0, acc_ul = 0;
    for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
        acc_dl += x[k] * (TCH_C_UW_DL[k] ? -1.f : 1.f);
        acc_ul += x[k] * (TCH_C_UW_UL[k] ? -1.f : 1.f);
    }
    const int* uw = std::abs(acc_dl) >= std::abs(acc_ul) ? TCH_C_UW_DL : TCH_C_UW_UL;

    int pos[TCH_C_MAX_KNOWN];
    float sign[TCH_C_MAX_KNOWN];
    const int n_known = tch_c_known_symbols(f, uw, pos, sign);

    // Remove the modulation
    gr_complex z[TCH_C_MAX_KNOWN];
    const gr_complex s0_conj = gr_complex(M_SQRT1_2, -M_SQRT1_2);
    for (int i = 0; i < n_known; i++) {
        z[i] = x[pos[i]] * sign[i] * s0_conj;
    }

    // Coarse frequency: maximise |sum z_k exp(-j omega k)| on a grid.
    // +-3 kHz in 25 Hz steps. A finer grid is not needed, see below.
    const float step = 2 * M_PI * 25. / f.symbols_per_second;
    const int n_grid = 120;
    float best_omega = 0, best_mag = -1;
    for (int g = -n_grid; g <= n_grid; g++) {
        const float omega = g * step;
        // The known symbols come in runs of consecutive symbols: step a phasor
        const gr_complex rot = std::polar(1.f, -omega);
        gr_complex acc = 0;
        gr_complex r = 0;
        for (int i = 0; i < n_known; i++) {
            if (i == 0 || pos[i] != pos[i - 1] + 1) {
                r = std::polar(1.f, -omega * pos[i]);
            }
            acc += z[i] * r;
            r *= rot;
        }
        const float mag = std::norm(acc);
        if (mag > best_mag) {
            best_mag = mag;
            best_omega = omega;
        }
    }

    // Fine frequency: least squares slope through the phases of the known blocks (UW,
    // pilots) after the coarse correction. With two blocks (TCH_C1) this is the phase
    // difference between UW and pilot, unambiguous for residuals below ~160 Hz.
    const int n_blocks = 1 + f.n_pilots;
    gr_complex h_block[1 + TCH_C_MAX_PILOTS] = {};
    float center[1 + TCH_C_MAX_PILOTS], weight[1 + TCH_C_MAX_PILOTS];
    for (int b = 0; b < n_blocks; b++) {
        const int first = b == 0 ? 0 : TCH_C_UW_LENGTH + (b - 1) * TCH_C_PILOT_LENGTH;
        const int len = b == 0 ? TCH_C_UW_LENGTH : TCH_C_PILOT_LENGTH;
        for (int i = first; i < first + len; i++) {
            h_block[b] += z[i] * std::polar(1.f, -best_omega * pos[i]);
        }
        center[b] = pos[first] + (len - 1) / 2.f;
        weight[b] = len;
    }
    float phase[1 + TCH_C_MAX_PILOTS];
    phase[0] = std::arg(h_block[0]);
    for (int b = 1; b < n_blocks; b++) {
        phase[b] = phase[b - 1] + std::arg(h_block[b] * std::conj(h_block[b - 1]));
    }
    float sw = 0, sx = 0, sp = 0;
    for (int b = 0; b < n_blocks; b++) {
        sw += weight[b];
        sx += weight[b] * center[b];
        sp += weight[b] * phase[b];
    }
    float sxx = 0, sxp = 0;
    for (int b = 0; b < n_blocks; b++) {
        sxx += weight[b] * (center[b] - sx / sw) * (center[b] - sx / sw);
        sxp += weight[b] * (center[b] - sx / sw) * (phase[b] - sp / sw);
    }
    const float omega = best_omega + sxp / sxx;

    // Phase at symbol 0
    gr_complex h = 0;
    for (int i = 0; i < n_known; i++) {
        h += z[i] * std::polar(1.f, -omega * pos[i]);
    }
    const float phi0 = std::arg(h);

    for (int k = 0; k < f.n_symbols; k++) {
        y[k] = x[k] * std::polar(1.f, -(phi0 + omega * k));
    }

    return omega;
}

/*
 * Pick the n TCH_C symbols from a burst with a (possibly fractional) number of samples
 * per symbol (probes: one), interpolating linearly. Returns the number of symbols.
 */
int tch_c_demod_impl::decimate(
    const gr_complex* in, int size, float sps, int n, gr_complex* out)
{
    int j;
    for (j = 0; j < n; j++) {
        float t = j * sps;
        int i = (int)t;
        if (i + 1 >= size) {
            break;
        }
        float mu = t - i;
        out[j] = in[i] + mu * (in[i + 1] - in[i]);
    }
    return j;
}

/*
 * Soft bits of the data symbols of a phase corrected TCH_C burst. Gain and noise variance
 * (relative to the symbol energy, written to noise) come from the known (QPSK) symbols.
 * Returns the gain; writes f.bits_per_symbol LLRs per
 * data symbol, L = log P(0) / P(1). QPSK: c[2k] from I, c[2k+1] from Q. 16APSK: max-log
 * over the 16 points for the 4 bits of the point index (MSB first).
 */
float tch_c_demod_impl::soft_bits(
    const tch_c_format& f, const gr_complex* y, const int* uw, float* llr, float* noise)
{
    const gr_complex s0 = gr_complex(M_SQRT1_2, M_SQRT1_2);
    int pos[TCH_C_MAX_KNOWN];
    float sign[TCH_C_MAX_KNOWN];
    const int n_known = tch_c_known_symbols(f, uw, pos, sign);

    float gain = 0;
    for (int i = 0; i < n_known; i++) {
        gain += (y[pos[i]] * std::conj(sign[i] * s0)).real();
    }
    gain /= n_known;
    if (gain <= 0) {
        gain = 1e-6;
    }

    float sigma2 = 0;
    for (int i = 0; i < n_known; i++) {
        sigma2 += std::norm(y[pos[i]] / gain - sign[i] * s0);
    }
    sigma2 = std::max(sigma2 / n_known, 1e-3f);
    *noise = sigma2;

    int p = 0;
    if (f.bits_per_symbol == 2) {
        const float scale = 2 * M_SQRT2 / sigma2 / gain;
        for (int s = 0; s < f.n_segments; s++) {
            for (int k = f.segment[s][0]; k < f.segment[s][0] + f.segment[s][1]; k++) {
                llr[p++] = scale * y[k].real();
                llr[p++] = scale * y[k].imag();
            }
        }
    } else {
        for (int s = 0; s < f.n_segments; s++) {
            for (int k = f.segment[s][0]; k < f.segment[s][0] + f.segment[s][1]; k++) {
                const gr_complex z = y[k] / gain;
                float d[16];
                for (int i = 0; i < 16; i++) {
                    d[i] =
                        std::norm(z - gr_complex(TCH_C_APSK16_I[i] * TCH_C_APSK16_SCALE,
                                                 TCH_C_APSK16_Q[i] * TCH_C_APSK16_SCALE));
                }
                for (int m = 0; m < 4; m++) {
                    float d0 = 1e30f, d1 = 1e30f;
                    for (int i = 0; i < 16; i++) {
                        if ((i >> (3 - m)) & 1) {
                            d1 = std::min(d1, d[i]);
                        } else {
                            d0 = std::min(d0, d[i]);
                        }
                    }
                    llr[p++] = (d1 - d0) / sigma2;
                }
            }
        }
    }
    return gain;
}

/*
 * Decode the codewords of a TCH_C burst. The scrambler runs on across the codewords and
 * one CRC-10 at the end covers all of them. Returns true if the CRC matches.
 */
static bool decode_tch_c(const tch_c_format& f,
                         tch_c_decoder& decoder,
                         const float* llr,
                         uint8_t* msg)
{
    const int K = f.code.info_bits;
    const int L = f.code.coded_bits;
    for (int j = 0; j < f.n_codewords; j++) {
        decoder.decode_raw(llr + j * L, msg + j * K);
    }
    tch_c_decoder::descramble(msg, f.n_codewords * K, 0);
    return tch_c_decoder::check_crc(msg, f.n_codewords * K);
}

// Mean squared distance of the data symbols to the nearest constellation point
static float
tch_c_fit(const tch_c_format& f, const gr_complex* y, float gain, int bits_per_symbol)
{
    float sum = 0;
    int n = 0;
    for (int s = 0; s < f.n_segments; s++) {
        for (int k = f.segment[s][0]; k < f.segment[s][0] + f.segment[s][1]; k++) {
            const gr_complex z = y[k] / gain;
            float best;
            if (bits_per_symbol == 2) {
                best = std::norm(gr_complex(std::fabs(z.real()), std::fabs(z.imag())) -
                                 gr_complex(M_SQRT1_2, M_SQRT1_2));
            } else {
                best = 1e30f;
                for (int i = 0; i < 16; i++) {
                    best = std::min(
                        best,
                        std::norm(z -
                                  gr_complex(TCH_C_APSK16_I[i] * TCH_C_APSK16_SCALE,
                                             TCH_C_APSK16_Q[i] * TCH_C_APSK16_SCALE)));
                }
            }
            sum += best;
            n++;
        }
    }
    return sum / n;
}

tch_c_decoder& tch_c_demod_impl::decoder_for(const tch_c_format& f)
{
    if (f.type == 8) {
        return f.bits_per_symbol == 4 ? d_c8a_decoder : d_c8q_decoder;
    }
    return f.type == 2 ? d_c2_decoder : d_c1_decoder;
}

/*
 * TCH_C bursts: demodulate, decode and check the CRC-10. If the UW is found but the
 * CRC fails, other demodulation variants (and for TCH_C8 the other modulation) are tried
 * and the first one passing the CRC is used.
 */
void tch_c_demod_impl::handle_burst(int channel,
                                    const gr_complex* burst,
                                    size_t burst_size,
                                    float sample_rate,
                                    double center_frequency,
                                    uint64_t timestamp,
                                    uint64_t sub_id,
                                    float noise,
                                    float magnitude,
                                    int tch_type,
                                    int max_variants)
{
    // Timing, pilots and data segments; for TCH_C8 shared by QPSK and 16APSK
    const tch_c_format& layout = tch_c_format_by_type(tch_type);
    const int n_symbols = layout.n_symbols;
    const float sps = sample_rate / layout.symbols_per_second;

    // Carrier recovery variants: 0 = fit to the known symbols, 1 = PLL, 2 = fit followed
    // by the PLL. On a test recording the fit alone passes the CRC on 99.8% of the TCH_C1
    // bursts that pass at all. The decision directed QPSK PLL does not work for 16APSK,
    // so TCH_C8 uses the fit only. Most probes are no burst at all and every variant is
    // one more chance for a failed decode to pass the CRC-10, so only the fit is tried
    // unless the prober asks for more (`variants`, for those which clearly have a UW).
    const int n_variants = tch_type == 8 ? 1 : std::max(1, std::min(3, max_variants));

    d_n_handled_bursts++;

    if (decimate(burst, burst_size, sps, n_symbols, d_decimated_burst) < n_symbols) {
        return;
    }

    const tch_c_format* chosen_format = nullptr;
    int chosen = -1;
    bool crc_ok = false;
    bool uplink = false;
    int confidence = 0;
    float level = 0;
    float chosen_gain = 1;
    double frequency_offset = 0;

    for (int v = 0; v < n_variants; v++) {
        double df = 0;
        if (v == 0 || v == 2) {
            float omega =
                known_symbol_correction(layout, d_decimated_burst, d_burst_after_pll);
            df = omega / (2 * M_PI) * layout.symbols_per_second;
        } else {
            memcpy(d_burst_after_pll, d_decimated_burst, n_symbols * sizeof(gr_complex));
        }
        if (v == 1 || v == 2) {
            float total_phase = qpsk_first_order_pll(
                d_burst_after_pll, n_symbols, d_alpha, d_burst_after_pll);
            df += total_phase / (n_symbols / (double)layout.symbols_per_second) /
                  (2 * M_PI);
        }

        int v_confidence;
        float v_level;
        demod_qpsk(
            d_burst_after_pll, n_symbols, d_demodulated_burst, &v_level, &v_confidence);

        // The UW is QPSK in every format
        bool dl_uw_ok = uw_ok(d_demodulated_burst, TCH_C_UW_DL);
        bool ul_uw_ok = uw_ok(d_demodulated_burst, TCH_C_UW_UL);
        if (!dl_uw_ok && !ul_uw_ok) {
            if (v == 0) {
                // Not a TCH_C burst (most bursts reaching this are not), so do not spend
                // time on the other variants. They are only meant to fix the CRC.
                break;
            }
            continue;
        }
        bool v_uplink = ul_uw_ok && !dl_uw_ok;
        const int* uw = v_uplink ? TCH_C_UW_UL : TCH_C_UW_DL;

        // Candidate formats: TCH_C8 is QPSK or 16APSK, try the better fitting one first
        const tch_c_format* formats[2] = { &layout, nullptr };
        int n_formats = 1;
        if (tch_type == 8) {
            float noise;
            float gain = soft_bits(TCH_C8Q, d_burst_after_pll, uw, d_llr, &noise);
            bool apsk_first = tch_c_fit(TCH_C8A, d_burst_after_pll, gain, 4) <
                              tch_c_fit(TCH_C8Q, d_burst_after_pll, gain, 2);
            formats[0] = apsk_first ? &TCH_C8A : &TCH_C8Q;
            formats[1] = apsk_first ? &TCH_C8Q : &TCH_C8A;
            n_formats = 2;
        }

        for (int fi = 0; fi < n_formats; fi++) {
            const tch_c_format& f = *formats[fi];
            float noise;
            float gain = soft_bits(f, d_burst_after_pll, uw, d_llr, &noise);
            if (f.bits_per_symbol == 4 && noise > APSK16_MAX_NOISE) {
                continue;
            }
            bool ok = decode_tch_c(f, decoder_for(f), d_llr, d_msg);

            // Keep the first try passing the UW check, unless a later one passes the CRC
            if (chosen < 0 || ok) {
                chosen = v;
                chosen_format = &f;
                chosen_gain = gain;
                crc_ok = ok;
                uplink = v_uplink;
                confidence = v_confidence;
                level = v_level;
                frequency_offset = df;
                memcpy(d_best_symbols, d_demodulated_burst, n_symbols * sizeof(int));
                memcpy(
                    d_best_corrected, d_burst_after_pll, n_symbols * sizeof(gr_complex));
                memcpy(d_best_msg, d_msg, sizeof(d_best_msg));
            }
            if (ok) {
                break;
            }
        }
        if (crc_ok) {
            break;
        }
    }

    if (chosen < 0) {
        return;
    }
    const tch_c_format& f = *chosen_format;

    uint64_t id = sub_id / 10;
    if (d_channel_id[channel] != id) {
        d_channel_id[channel] = id;
        d_n_access_ok_bursts++;
    }
    d_n_access_ok_sub_bursts++;

    if (f.bits_per_symbol == 2) {
        decode_gray_qpsk(d_best_symbols, n_symbols);
        map_symbols_to_bits(d_best_symbols, n_symbols, d_bits);
    } else {
        // 16APSK: the 4 bits of the nearest point for every symbol
        d_bits.clear();
        for (int k = 0; k < n_symbols; k++) {
            const gr_complex z = d_best_corrected[k] / chosen_gain;
            int best = 0;
            float best_d = 1e30f;
            for (int i = 0; i < 16; i++) {
                float d =
                    std::norm(z - gr_complex(TCH_C_APSK16_I[i] * TCH_C_APSK16_SCALE,
                                             TCH_C_APSK16_Q[i] * TCH_C_APSK16_SCALE));
                if (d < best_d) {
                    best_d = d;
                    best = i;
                }
            }
            for (int m = 3; m >= 0; m--) {
                d_bits.push_back((best >> m) & 1);
            }
        }
    }

    const int payload_bits = f.n_codewords * f.code.info_bits - 10;

    pmt::pmt_t pdu_meta = pmt::make_dict();
    pmt::pmt_t pdu_vector = pmt::init_u8vector(d_bits.size(), d_bits);
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("timestamp"), pmt::mp(timestamp));
    pdu_meta = pmt::dict_add(pdu_meta,
                             pmt::mp("center_frequency"),
                             pmt::mp(center_frequency + frequency_offset));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("id"), pmt::mp(sub_id));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("confidence"), pmt::mp(confidence));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("level"), pmt::mp(level));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("noise"), pmt::mp(noise));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("magnitude"), pmt::mp(magnitude));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("n_symbols"), pmt::mp(n_symbols));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("tch_type"), pmt::mp(f.type));
    pdu_meta =
        pmt::dict_add(pdu_meta, pmt::mp("bits_per_symbol"), pmt::mp(f.bits_per_symbol));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("uplink"), pmt::from_bool(uplink));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("crc_ok"), pmt::from_bool(crc_ok));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("variant"), pmt::mp(chosen));
    pdu_meta = pmt::dict_add(
        pdu_meta, pmt::mp("payload"), pmt::init_u8vector(payload_bits, d_best_msg));

    pmt::pmt_t out_msg = pmt::cons(pdu_meta, pdu_vector);
    message_port_pub(pmt::mp("pdus"), out_msg);
}

/*
 * Report the probe as handled to the prober (probe_handled), so it sees the backlog of
 * the demodulator.
 */
void tch_c_demod_impl::frame_handled(int channel, const pmt::pmt_t& meta)
{
    pmt::pmt_t burst_id = pmt::dict_ref(meta, pmt::mp("burst_id"), pmt::PMT_NIL);
    if (!pmt::is_null(burst_id)) {
        message_port_pub(d_frame_handled_ports[channel], burst_id);
    }
}

void tch_c_demod_impl::enqueue(int channel, const pmt::pmt_t& msg)
{
    uint64_t timestamp =
        pmt::to_uint64(pmt::dict_ref(pmt::car(msg), pmt::mp("timestamp"), pmt::PMT_NIL));
    d_pending.emplace(timestamp, std::make_pair(channel, msg));
}

/*
 * The prober publishes the UL slots later than the DL slots and the tracks one after
 * the other: take what is waiting and always handle the oldest probe, so the frame
 * sorter gets the frames in time order.
 */
void tch_c_demod_impl::handler(int channel, pmt::pmt_t msg)
{
    enqueue(channel, msg);

    while (!d_pending.empty()) {
        for (size_t i = 0; i < d_input_ports.size(); i++) {
            while ((msg = delete_head_nowait(d_input_ports[i]))) {
                enqueue(i, msg);
            }
        }

        auto oldest = d_pending.begin();
        channel = oldest->second.first;
        msg = oldest->second.second;
        d_pending.erase(oldest);

        demodulate(channel, msg);
        frame_handled(channel, pmt::car(msg));
    }
}

void tch_c_demod_impl::demodulate(int channel, pmt::pmt_t msg)
{
    pmt::pmt_t samples = pmt::cdr(msg);
    size_t burst_size = pmt::length(samples);
    const gr_complex* burst =
        (const gr_complex*)pmt::c32vector_elements(samples, burst_size);

    pmt::pmt_t meta = pmt::car(msg);
    double center_frequency =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
    float sample_rate =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("sample_rate"), pmt::PMT_NIL));
    uint64_t sub_id = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("id"), pmt::PMT_NIL));
    double uw_start =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("uw_start"), pmt::PMT_NIL));
    float noise = pmt::to_float(pmt::dict_ref(meta, pmt::mp("noise"), pmt::PMT_NIL));
    float magnitude =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("magnitude"), pmt::PMT_NIL));
    uint64_t timestamp =
        pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    int tch_type = pmt::to_long(pmt::dict_ref(meta, pmt::mp("tch_type"), pmt::PMT_NIL));
    bool probe = pmt::to_bool(pmt::dict_ref(meta, pmt::mp("probe"), pmt::PMT_F));
    int max_variants = pmt::to_long(
        pmt::dict_ref(meta, pmt::mp("variants"), pmt::from_long(probe ? 1 : 3)));

    timestamp += uw_start * 1000000000ULL / (int)sample_rate;

    handle_burst(channel,
                 burst,
                 burst_size,
                 sample_rate,
                 center_frequency,
                 timestamp,
                 sub_id,
                 noise,
                 magnitude,
                 tch_type,
                 max_variants);
}

} /* namespace iridium */
} /* namespace gr */
