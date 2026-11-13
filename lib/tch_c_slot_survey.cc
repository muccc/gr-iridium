/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tch_c_slot_survey.h"
#include <volk/volk.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace gr {
namespace iridium {

namespace {

/*
 * The energy window: 1 from half a symbol before the centre of UW symbol 0 to half a
 * symbol after the centre of UW symbol 31, cosine tapers of half a symbol on both
 * sides. The energy of the data symbols after the UW stays mostly outside.
 */
const double WIN_START = -0.5;
const double WIN_END = 31.5;
const double WIN_TAPER = 0.5;

// Root raised cosine pulse, t in symbols, peak 1 - alpha + 4 alpha / pi
double rrc(double t, double alpha)
{
    if (std::fabs(t) < 1e-9) {
        return 1 - alpha + 4 * alpha / M_PI;
    }
    if (std::fabs(std::fabs(4 * alpha * t) - 1) < 1e-9) {
        return alpha / M_SQRT2 *
               ((1 + 2 / M_PI) * std::sin(M_PI / (4 * alpha)) +
                (1 - 2 / M_PI) * std::cos(M_PI / (4 * alpha)));
    }
    return (std::sin(M_PI * t * (1 - alpha)) +
            4 * alpha * t * std::cos(M_PI * t * (1 + alpha))) /
           (M_PI * t * (1 - 16 * alpha * alpha * t * t));
}

// Raised cosine spectrum, f in units of the symbol rate, peak 1
double rc_spectrum(double f, double alpha)
{
    f = std::fabs(f);
    if (f <= (1 - alpha) / 2) {
        return 1;
    }
    if (f <= (1 + alpha) / 2) {
        return 0.5 * (1 + std::cos(M_PI / alpha * (f - (1 - alpha) / 2)));
    }
    return 0;
}

// Parabolic peak through (-1, a), (0, b), (1, c): offset in -1..1 and peak value
void parabola(float a, float b, float c, float* offset, float* peak)
{
    const float den = a - 2 * b + c;
    if (!(den < 0)) {
        // Not a maximum: take the largest sample
        if (a > b && a >= c) {
            *offset = -1;
            *peak = a;
        } else if (c > b) {
            *offset = 1;
            *peak = c;
        } else {
            *offset = 0;
            *peak = b;
        }
        return;
    }
    float d = 0.5f * (a - c) / den;
    d = std::max(-1.f, std::min(1.f, d));
    *offset = d;
    *peak = b - 0.25f * (a - c) * d;
}

int next_pow2(int n)
{
    int N = 1;
    while (N < n) {
        N *= 2;
    }
    return N;
}

const int* uw_of(bool uplink) { return uplink ? TCH_C_UW_UL : TCH_C_UW_DL; }

} // namespace

int tch_c_slot_survey::n_positions(int type)
{
    return type == 8 ? 30 : type == 2 ? 120 : 240;
}

double tch_c_slot_survey::nominal_frequency(int type, int position)
{
    if (type == 8) {
        return BAND_START + (8 * position + 4) * CHANNEL_WIDTH;
    }
    if (type == 2) {
        return BAND_START + (2 * position + 1) * CHANNEL_WIDTH;
    }
    return BAND_START + (position + 0.5) * CHANNEL_WIDTH;
}

tch_c_slot_survey::tch_c_slot_survey(double sample_rate, double center_frequency)
    : d_sample_rate(sample_rate), d_center_frequency(center_frequency), d_noise_db(0)
{
    init_type(d_t[0], TCH_C1);
    init_type(d_t[1], TCH_C2);
    init_type(d_t[2], TCH_C8Q);
}

tch_c_slot_survey::~tch_c_slot_survey() {}

void tch_c_slot_survey::init_type(type_tables& tt, const tch_c_format& f)
{
    const double alpha = RRC_ALPHA;
    const double Rs = f.symbols_per_second;
    const double sps = d_sample_rate / Rs;
    tt.f = &f;
    tt.sps = sps;
    tt.D = (int)std::ceil(TEMPLATE_SPAN * sps);
    tt.L = (int)std::floor(tt.D + (TCH_C_UW_LENGTH - 1 + TEMPLATE_SPAN) * sps) + 1;
    tt.N = next_pow2(tt.L);
    tt.K = (int)std::floor((1 + alpha) / 2 * Rs * tt.N / d_sample_rate);
    tt.h = std::max(1, (int)std::lround(sps / 4));
    const int N = tt.N, L = tt.L, K = tt.K;

    tt.fft_c.reset(new gr::fft::fft_complex_fwd(N, 1));
    tt.fft_p.reset(new gr::fft::fft_complex_fwd(N, 1));
    std::fill(tt.fft_c->get_inbuf(), tt.fft_c->get_inbuf() + N, gr_complex(0));
    std::fill(tt.fft_p->get_inbuf(), tt.fft_p->get_inbuf() + N, gr_complex(0));
    tt.p_shifted.resize(N);
    tt.u_rot.resize(L);

    // Energy window
    tt.w.resize(L);
    tt.sum_w2 = 0;
    for (int n = 0; n < L; n++) {
        const double s = (n - tt.D) / sps;
        double v = 0;
        if (s >= WIN_START && s <= WIN_END) {
            v = 1;
        } else if (s > WIN_START - WIN_TAPER && s < WIN_START) {
            v = 0.5 * (1 + std::cos(M_PI * (WIN_START - s) / WIN_TAPER));
        } else if (s > WIN_END && s < WIN_END + WIN_TAPER) {
            v = 0.5 * (1 + std::cos(M_PI * (s - WIN_END) / WIN_TAPER));
        }
        tt.w[n] = v;
        tt.sum_w2 += v * v;
    }

    // Band weights
    tt.H.resize(2 * K + 1);
    tt.sum_h = 0;
    for (int k = -K; k <= K; k++) {
        tt.H[k + K] = rc_spectrum(k * d_sample_rate / N / Rs, alpha);
        tt.sum_h += tt.H[k + K];
    }

    // UW templates: RRC pulses at the exact (fractional) symbol centres, unit energy
    for (int d = 0; d < 2; d++) {
        const int* uw = uw_of(d == 1);
        std::vector<gr_complex> u(L);
        double energy = 0;
        for (int n = 0; n < L; n++) {
            const double s = (n - tt.D) / sps;
            double v = 0;
            const int k0 = std::max(0, (int)std::ceil(s - TEMPLATE_SPAN));
            const int k1 =
                std::min(TCH_C_UW_LENGTH - 1, (int)std::floor(s + TEMPLATE_SPAN));
            for (int k = k0; k <= k1; k++) {
                v += (uw[k] ? -1 : 1) * rrc(s - k, alpha);
            }
            u[n] = gr_complex(v * M_SQRT1_2, v * M_SQRT1_2);
            energy += v * v;
        }
        const float scale = 1 / std::sqrt(energy);
        dir_tables& dt = tt.dir[d];
        dt.u_conj.resize(L);
        for (int n = 0; n < L; n++) {
            u[n] *= scale;
            dt.u_conj[n] = std::conj(u[n]);
        }

        // sum |FFT(u w)|^2 H: E(b) of a noise free UW of unit energy
        gr_complex* in = tt.fft_p->get_inbuf();
        for (int n = 0; n < L; n++) {
            in[n] = u[n] * tt.w[n];
        }
        tt.fft_p->execute();
        double u2h = 0;
        for (int k = -K; k <= K; k++) {
            u2h += std::norm(tt.fft_p->get_outbuf()[(k + N) % N]) * tt.H[k + K];
        }
        dt.u2h = u2h;

        /*
         * Response of the UW correlation to a noise free UW at d bins from a bin:
         * R(d) = |sum |u[n]|^2 e^-j2pi d n / N|, R(0) = 1. Tabulated for |d| <= 1.5 and
         * as q(d) = (R(1 - d) - R(1 + d)) / R(d), which the survey inverts to find the
         * peak between the bins.
         */
        dt.r_table.resize(R_RANGE * R_STEPS + 2);
        for (size_t i = 0; i < dt.r_table.size(); i++) {
            const double d = (double)i / R_STEPS;
            std::complex<double> acc = 0;
            for (int n = 0; n < L; n++) {
                acc += (double)std::norm(u[n]) * std::polar(1.0, -2 * M_PI * d * n / N);
            }
            dt.r_table[i] = std::abs(acc);
        }
        const int nq = 257;
        dt.q_table.resize(nq);
        for (int i = 0; i < nq; i++) {
            const double d = -0.5 + (double)i / (nq - 1);
            dt.q_table[i] = (dt.response(1 - d) - dt.response(1 + d)) / dt.response(d);
        }
        for (int i = 1; i < nq; i++) {
            if (!(dt.q_table[i] > dt.q_table[i - 1])) {
                throw std::runtime_error("tch_c_slot_survey: UW response not monotonic");
            }
        }
    }
    std::fill(tt.fft_p->get_inbuf(), tt.fft_p->get_inbuf() + N, gr_complex(0));

    // Polyphase matched filter for extract(): taps m = -M .. M + 1 around floor(t),
    // phase q at t = floor(t) + q / PHASES, normalised to a sum of 1
    tt.M = (int)std::ceil(EXTRACT_SPAN * sps);
    tt.n_taps = 2 * tt.M + 2;
    tt.poly.resize(PHASES);
    tt.rotated.resize(PHASES);
    tt.rotated_stamp.assign(PHASES, -1);
    tt.stamp = 0;
    tt.phasor.resize(tt.n_taps);
    for (int q = 0; q < PHASES; q++) {
        const double mu = (double)q / PHASES;
        tt.poly[q].resize(tt.n_taps);
        tt.rotated[q].resize(tt.n_taps);
        double sum = 0;
        for (int i = 0; i < tt.n_taps; i++) {
            const double s = (i - tt.M - mu) / sps;
            const double v = std::fabs(s) <= EXTRACT_SPAN ? rrc(s, alpha) : 0;
            tt.poly[q][i] = v;
            sum += v;
        }
        for (int i = 0; i < tt.n_taps; i++) {
            tt.poly[q][i] /= sum;
        }
    }
}

/*
 * C[b] = FFT_N(seg conj(u))[b] is the UW correlation at frequency b fs / N, the energy
 * E(b) = sum_k |FFT_N(seg w)[b + k]|^2 H[k] the band energy under the UW. For
 * seg = A u e^(j 2 pi b n / N): C = A (u has unit energy), E = A^2 U2H with
 * U2H = sum_k |FFT_N(u w)[k]|^2 H[k], so rho2 = |C|^2 U2H / E = 1.
 */
int tch_c_slot_survey::survey(const gr_complex* seg,
                              int type,
                              bool uplink,
                              double doppler_ring_hz,
                              double search_hz,
                              float threshold,
                              std::vector<candidate>& out,
                              double offset_hz)
{
    type_tables& tt = t(type);
    const dir_tables& dt = tt.dir[uplink ? 1 : 0];
    const int N = tt.N, L = tt.L, K = tt.K;
    const double fs = d_sample_rate;
    const double Rs = tt.f->symbols_per_second;

    // C = FFT(seg conj(u)), P = |FFT(seg w)|^2 (the input buffers stay zero beyond L)
    volk_32fc_x2_multiply_32fc(tt.fft_c->get_inbuf(), seg, &dt.u_conj[0], L);
    volk_32fc_32f_multiply_32fc(tt.fft_p->get_inbuf(), seg, &tt.w[0], L);
    tt.fft_c->execute();
    tt.fft_p->execute();
    const gr_complex* C = tt.fft_c->get_outbuf();
    const gr_complex* Pc = tt.fft_p->get_outbuf();
    volk_32fc_magnitude_squared_32f(&tt.p_shifted[0], Pc + N / 2, N / 2);
    volk_32fc_magnitude_squared_32f(&tt.p_shifted[N / 2], Pc, N / 2);

    // Noise floor: the 10% percentile of (a subsample of) the inner 80% of the band.
    // For exponentially distributed bin powers the mean is 1 / -ln(0.9) times that.
    const int step = std::max(1, N / 2048);
    d_noise_tmp.clear();
    for (int i = N / 10; i < N - N / 10; i += step) {
        d_noise_tmp.push_back(tt.p_shifted[i]);
    }
    const size_t nth = d_noise_tmp.size() / 10;
    std::nth_element(d_noise_tmp.begin(), d_noise_tmp.begin() + nth, d_noise_tmp.end());
    const float p_noise = std::max(d_noise_tmp[nth] / 0.105361f, 1e-30f);
    d_noise_db = 10 * std::log10(p_noise / tt.sum_w2 / fs);

    const int sb = std::max(1, (int)std::ceil(search_hz * N / fs));
    const double half_band = (1 + RRC_ALPHA) / 2 * Rs;
    const int n_pos = n_positions(type);
    int n_looked_at = 0;

    for (int pos = 0; pos < n_pos; pos++) {
        const double f_nom = nominal_frequency(type, pos);
        const double f = f_nom * (1 + doppler_ring_hz / F_RING) + offset_hz;
        const double rel = f - d_center_frequency;
        if (std::fabs(rel) + half_band + search_hz > fs / 2) {
            continue;
        }
        const int b0 = (int)std::lround(rel * N / fs);
        if (std::abs(b0) + sb + K + 2 > N / 2) {
            continue;
        }
        n_looked_at++;

        // Peak of |C| within +-sb bins, with one more bin on each side for the parabola
        float m[2 * 64 + 3];
        const int nb = std::min(2 * sb + 3, (int)(sizeof(m) / sizeof(m[0])));
        const int first = b0 - (nb - 1) / 2;
        int best = 1;
        for (int i = 0; i < nb; i++) {
            m[i] = std::abs(C[(first + i + N) % N]);
            if (i >= 1 && i < nb - 1 && m[i] > m[best]) {
                best = i;
            }
        }
        const int bp = first + best;
        float offset, peak;
        interpolate_peak(dt, m[best - 1], m[best], m[best + 1], &offset, &peak);

        float energy;
        volk_32f_x2_dot_prod_32f(
            &energy, &tt.p_shifted[N / 2 + bp - K], &tt.H[0], 2 * K + 1);
        const float rho2 = energy > 0 ? peak * peak * dt.u2h / energy : 0;
        if (!(rho2 >= threshold)) {
            continue;
        }
        const double f_est = (bp + offset) / (double)N;

        candidate c;
        c.position = pos;
        c.nominal_frequency = f_nom;
        c.frequency = f;
        c.estimated_frequency = d_center_frequency + f_est * fs;
        c.rho2 = rho2;
        c.energy = energy;
        c.magnitude = 10 * std::log10(std::max(energy / (p_noise * tt.sum_h), 1e-6f));
        c.bin = bp;
        out.push_back(c);
    }
    return n_looked_at;
}

void tch_c_slot_survey::interpolate_peak(
    const dir_tables& dt, float a, float b, float c, float* offset, float* peak) const
{
    // The bins hold R(-1 - d), R(-d), R(1 - d) for a peak at d (R is even): look up
    // d from q = (c - a) / b, increasing in d
    const float q = b > 0 ? (c - a) / b : 0;
    const std::vector<float>& t = dt.q_table;
    const int n = t.size();
    float d;
    if (!(q > t[0])) {
        d = -0.5f;
    } else if (!(q < t[n - 1])) {
        d = 0.5f;
    } else {
        const int i = std::upper_bound(t.begin(), t.end(), q) - t.begin() - 1;
        const float frac = (q - t[i]) / (t[i + 1] - t[i]);
        d = -0.5f + (i + frac) / (n - 1);
    }
    *offset = d;
    *peak = b / dt.response(d);
}

float tch_c_slot_survey::dir_tables::response(double d) const
{
    // Linear interpolation in r_table, |d| up to R_RANGE
    const double x = std::min(std::fabs(d), (double)R_RANGE) * R_STEPS;
    const int i = std::min((int)x, (int)r_table.size() - 2);
    const double frac = x - i;
    return r_table[i] + frac * (r_table[i + 1] - r_table[i]);
}

float tch_c_slot_survey::correlation_magnitude(const gr_complex* seg,
                                               const type_tables& tt,
                                               int shift)
{
    gr_complex r;
    volk_32fc_x2_dot_prod_32fc(&r, seg + shift, &tt.u_rot[0], tt.L);
    return std::abs(r);
}

double tch_c_slot_survey::fine_timing(
    const gr_complex* seg, int type, bool uplink, const candidate& c, float* rho2)
{
    type_tables& tt = t(type);
    const dir_tables& dt = tt.dir[uplink ? 1 : 0];

    // The template at the estimated frequency: sum seg[n + shift] conj(u[n]) e^-jwn
    const double w =
        2 * M_PI * (c.estimated_frequency - d_center_frequency) / d_sample_rate;
    lv_32fc_t inc = (lv_32fc_t)std::polar(1.0, -w);
    lv_32fc_t phase = 1;
    volk_32fc_s32fc_x2_rotator2_32fc(&tt.u_rot[0], &dt.u_conj[0], &inc, &phase, tt.L);

    const int h = tt.h;
    float offset, peak;
    parabola(correlation_magnitude(seg, tt, -h),
             correlation_magnitude(seg, tt, 0),
             correlation_magnitude(seg, tt, h),
             &offset,
             &peak);

    // Refine around the first estimate with a quarter of the step
    const int h2 = std::max(1, h / 4);
    const int center = (int)std::lround(offset * h);
    float offset2, peak2;
    parabola(correlation_magnitude(seg, tt, center - h2),
             correlation_magnitude(seg, tt, center),
             correlation_magnitude(seg, tt, center + h2),
             &offset2,
             &peak2);
    double tau = center + offset2 * h2;
    tau = std::max((double)-h, std::min((double)h, tau));

    if (rho2) {
        *rho2 = c.energy > 0 ? peak2 * peak2 * dt.u2h / c.energy : 0;
    }
    return tau;
}

void tch_c_slot_survey::extract_span(int type, double p_uw, long* first, long* end) const
{
    const type_tables& tt = t(type);
    *first = (long)std::floor(p_uw) - tt.M;
    *end = (long)std::floor(p_uw + tt.f->n_symbols * tt.sps) + 1 + tt.M + 2;
}

int tch_c_slot_survey::extract(const gr_complex* x,
                               double p_uw,
                               int type,
                               double freq_hz,
                               bool uplink,
                               gr_complex* out)
{
    type_tables& tt = t(type);
    const int n = tt.f->n_symbols + 1;
    const double w = 2 * M_PI * (freq_hz - d_center_frequency) / d_sample_rate;

    // e^-jw(i - M) for the taps, by a recursion renormalised in blocks
    {
        const std::complex<double> step = std::polar(1.0, -w);
        for (int i = 0; i < tt.n_taps; i += 64) {
            std::complex<double> z = std::polar(1.0, -w * (i - tt.M));
            const int end = std::min(tt.n_taps, i + 64);
            for (int k = i; k < end; k++) {
                tt.phasor[k] = gr_complex(z);
                z *= step;
            }
        }
    }
    tt.stamp++;

    // e^-jw n0 of the symbols: n0 advances by floor(sps) or floor(sps) + 1
    const long sps_floor = (long)std::floor(tt.sps);
    const std::complex<double> step0 = std::polar(1.0, -w * sps_floor);
    const std::complex<double> step1 = std::polar(1.0, -w * (sps_floor + 1));
    std::complex<double> z;
    long n_prev = 0;
    for (int j = 0; j < n; j++) {
        const double tj = p_uw + j * tt.sps;
        long n0 = (long)std::floor(tj);
        int q = (int)std::lround((tj - n0) * PHASES);
        if (q == PHASES) {
            q = 0;
            n0++;
        }
        if (tt.rotated_stamp[q] != tt.stamp) {
            volk_32fc_32f_multiply_32fc(
                &tt.rotated[q][0], &tt.phasor[0], &tt.poly[q][0], tt.n_taps);
            tt.rotated_stamp[q] = tt.stamp;
        }
        gr_complex r;
        volk_32fc_x2_dot_prod_32fc(&r, x + n0 - tt.M, &tt.rotated[q][0], tt.n_taps);
        const long d = n0 - n_prev;
        if (j % 64 == 0 || (d != sps_floor && d != sps_floor + 1)) {
            z = std::polar(1.0, -w * n0);
        } else {
            z *= d == sps_floor ? step0 : step1;
        }
        n_prev = n0;
        out[j] = r * gr_complex(z);
    }

    // Phase of the UW correlation to 0: UW symbols with bit 0 at +45 degrees
    const int* uw = uw_of(uplink);
    gr_complex acc = 0;
    for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
        acc += out[k] * (uw[k] ? -1.f : 1.f);
    }
    acc *= gr_complex(M_SQRT1_2, -M_SQRT1_2);
    if (std::abs(acc) > 0) {
        const gr_complex rot = std::conj(acc) / std::abs(acc);
        for (int j = 0; j < n; j++) {
            out[j] *= rot;
        }
    }
    return n;
}

} // namespace iridium
} // namespace gr
