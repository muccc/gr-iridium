/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_SLOT_SURVEY_H
#define INCLUDED_IRIDIUM_TCH_C_SLOT_SURVEY_H

#include "tch_c_format.h"
#include <gnuradio/fft/fft.h>
#include <gnuradio/gr_complex.h>
#include <volk/volk_alloc.hh>
#include <algorithm>
#include <memory>
#include <vector>

namespace gr {
namespace iridium {

/*
 * Survey of TCH_C slots on the wideband stream: at a predicted time, which raster
 * channels of a TCH_C type carry a burst? And extraction of the found bursts as symbol
 * rate samples for tch_c_demod, without burst_downmix.
 *
 * Sample index conventions (x is the wideband stream, sample i at time i / fs):
 *
 * - The centre of a symbol is the sample index (fractional) where the transmitted RRC
 *   pulse of that symbol peaks. p_uw is the centre of UW symbol 0, the time of the
 *   `timestamp` of a TCH_C frame.
 *
 * - survey() and fine_timing() work on a window seg = x + i0 of window_length(type)
 *   samples, i0 = lround(p_uw) - uw_offset(type): seg[uw_offset(type)] is the sample
 *   nearest to the centre of UW symbol 0 (the survey rounds a fractional prediction to
 *   the nearest sample, at most half a sample off: 1/80 of a symbol at 10 MS/s for
 *   TCH_C8, a correlation loss well below 1%). The window starts 6 symbols before the
 *   centre of UW symbol 0 and ends 6 symbols after the centre of UW symbol 31.
 *   fine_timing() also reads timing_margin(type) samples before and after the window:
 *   seg[-margin .. window_length + margin).
 *
 * - The UW position found is p_uw = i0 + uw_offset(type) + fine timing offset.
 *
 * - extract() takes the fractional p_uw directly and reads x[first, end) as returned by
 *   extract_span().
 *
 * What tch_c_demod expects in a cpdu (see tch_c_demod_impl::demodulate()/decimate()):
 * - vector: c32 samples, burst[j * sps] the centre of symbol j, j = 0 .. n_symbols - 1
 *   with sps = sample_rate / symbols_per_second. decimate() interpolates linearly
 *   between burst[i] and burst[i + 1], so it needs n_symbols + 1 samples even at one
 *   sample per symbol. burst[0] is the centre of UW symbol 0: burst_downmix cuts the
 *   burst at the (integer) correlation peak; the fractional part of the peak only goes
 *   into the `uw_start` meta and the timestamp, the samples are not shifted.
 * - meta: `sample_rate` (float, here the symbol rate: sps = 1), `center_frequency`
 *   (double, Hz, the frequency the burst was mixed down from; the demod adds the
 *   frequency its fit finds, +-3 kHz), `timestamp` (uint64 ns of burst[0]),
 *   `uw_start` (float, samples at sample_rate added to the timestamp: 0.0 here),
 *   `id` (uint64), `noise`, `magnitude` (float, printed), `tch_type` (long) and
 *   optionally `burst_id` (uint64, sent back on frame_handled).
 * - scale and phase: free. The demod estimates phase, frequency and gain from the UW
 *   and the pilots (the PLL variants start with the phase of the UW correlation).
 *   burst_downmix rotates the burst so the UW correlation has phase 0, i.e. UW and
 *   pilot symbols with bit 0 lie at +45 degrees; extract() does the same. Its RRC taps
 *   sum to 1 like those of gr::filter::firdes::root_raised_cosine(1.0, ...) in
 *   burst_downmix, so the `level` the demod prints is comparable.
 *
 * Not thread safe (scratch buffers); one instance per thread.
 */
class tch_c_slot_survey
{
public:
    // Ring alert channel and TCH_C channel raster
    static constexpr double F_RING = 1626270833.33;
    static constexpr double CHANNEL_WIDTH = 10e6 / 240;
    static constexpr double BAND_START = 1616e6;
    static constexpr float RRC_ALPHA = 0.2f;
    // Template: RRC pulses truncated to +-6 symbols
    static constexpr int TEMPLATE_SPAN = 6;
    // Matched filter of extract(): +-8 symbols like burst_downmix (+-2.5 symbols is too
    // short for the roll-off of 0.2 and lets the neighbouring channels in)
    static constexpr double EXTRACT_SPAN = 8;
    // Windows of the channel search: at least this long, and a long one
    static constexpr double CHANNEL_WINDOW_S = 2.5e-3;
    static constexpr double CHANNEL_WINDOW_LONG_S = 5e-3;

    struct candidate {
        int position;               // raster index: TCH_C1 c, TCH_C2 j, TCH_C8 m
        double nominal_frequency;   // Hz, raster frequency without Doppler
        double frequency;           // Hz, predicted carrier (with Doppler and offset)
        double estimated_frequency; // Hz, peak of the UW correlation (see survey())
        float rho2;                 // UW correlation power over in-band energy, 0..~1
        float energy;               // E(b), in-band energy under the UW (FFT units)
        float magnitude;            // dB, in-band PSD over noise(): (S + N) / N
        int bin;                    // FFT bin of the peak, -N/2 .. N/2 - 1
    };

    tch_c_slot_survey(double sample_rate, double center_frequency);
    ~tch_c_slot_survey();

    tch_c_slot_survey(const tch_c_slot_survey&) = delete;
    tch_c_slot_survey& operator=(const tch_c_slot_survey&) = delete;

    double sample_rate() const { return d_sample_rate; }
    double center_frequency() const { return d_center_frequency; }

    // Per TCH_C type (1, 2, 8)
    double samples_per_symbol(int type) const { return t(type).sps; }
    int window_length(int type) const { return t(type).L; }
    int uw_offset(int type) const { return t(type).D; }
    int fft_size(int type) const { return t(type).N; }
    int timing_margin(int type) const { return t(type).h + std::max(1, t(type).h / 4); }
    double bin_width(int type) const { return d_sample_rate / t(type).N; }

    // Raster positions of a type: 240, 120, 30
    static int n_positions(int type);
    static double nominal_frequency(int type, int position);

    /*
     * Survey one slot: seg[0 .. window_length(type)) (see above). Looks at every raster
     * position whose band, at f = f_nominal * (1 + doppler_ring_hz / F_RING) +
     * offset_hz, lies inside the captured band, in the FFT bins within +-search_hz of f
     * (at least +-1 bin, at most +-64 bins). Appends the positions with rho2 >=
     * threshold to out, ordered by position. Returns the number of positions looked
     * at.
     *
     * Per position: the largest |C| of the bins and its two neighbours; the offset of
     * the peak from that bin and its height come from the tabulated response of the UW
     * correlation to a noise free UW (exact for one, unlike a parabola).
     * rho2 = |C_peak|^2 U2H / E(b): 1 for a noise free UW, see the .cc. E(b) is the
     * band energy under an energy window from half a symbol before UW symbol 0 to half
     * a symbol after UW symbol 31, so the data after the UW hardly lowers rho2 (clean
     * bursts: 0.99). Noise only: mean ~0.05 (one frequency cell ~1/37, the search over
     * the bins about doubles it), ~2e-4 of the positions reach 0.25, none 0.35 in
     * 54,000. On a real recording (sat-147) baseline frames: median 0.95, minimum 0.68.
     *
     * Downlink: doppler_ring_hz = Doppler of the ring alerts, offset_hz = 0. Uplink
     * (pre-compensated by the terminal): doppler_ring_hz = -Doppler, offset_hz = the
     * mean residual, search_hz to cover the rest.
     *
     * TCH_C8 (4.17 us symbols): call it for several times around the prediction (e.g.
     * -2/0/+2 us) if the prediction is not good to ~1 us.
     */
    int survey(const gr_complex* seg,
               int type,
               bool uplink,
               double doppler_ring_hz,
               double search_hz,
               float threshold,
               std::vector<candidate>& out,
               double offset_hz = 0);

    // Noise floor of the last survey in dBFS/Hz (10 log10(sigma^2 / fs), the units of
    // fft_burst_tagger's `noise`): the 10% percentile of the inner 80% of the spectrum,
    // so it holds with many channels busy. The tagger's baseline includes the traffic
    // and reads several dB higher in a busy band.
    float noise() const { return d_noise_db; }

    /*
     * Fine timing of a candidate of survey() on the same seg (reads timing_margin(type)
     * samples beyond both ends): UW correlation at the estimated frequency at -h, 0, +h
     * samples (h = lround(sps / 4)), parabolic peak, refined once with steps of h / 4,
     * clamped to +-h. Returns the offset of the UW in samples (fractional):
     * p_uw = i0 + uw_offset(type) + offset. rho2 (optional): rho2 at the peak.
     * Accuracy: noise free ~0.005 symbol, at Es/N0 9 dB ~0.03 symbol.
     */
    double fine_timing(
        const gr_complex* seg, int type, bool uplink, const candidate& c, float* rho2);

    /*
     * Search of the UW of one type on single channels over a range of times, for bursts
     * whose time is known only to a fraction of a ms (uplink). set_window() takes the
     * FFT of x[0, len) once (len <= window_capacity(true), zero padded; the shorter
     * FFT if len <= window_capacity()), then
     * search_channel() per channel: matched filter in the frequency domain around f_hz,
     * IFFT to ~3 samples per symbol, and per output sample (lag) the UW correlation over
     * a frequency search of +-search_hz (zero padded DFT over the 32 UW symbols).
     * Returns false if no lag with p in [p_first, p_last] (sample indices of x, centre
     * of UW symbol 0) reaches threshold; else the best one in out. Its p_uw is good to
     * half an output sample (~0.16 symbol): refine it with fine_timing() on
     * seg = x + lround(p_uw) - uw_offset(type) (candidate: estimated_frequency = the
     * found frequency, energy = 0). The edges of the window (6 symbols of the matched
     * filter, and the UW length after p_last) are not searched.
     *
     * rho2 here: |sum_k y_k s_k*|^2 / (32 sum_k |y_k|^2) over the matched filter
     * outputs y_k at the UW symbols s_k, 1 for a clean UW, ~1/32 per cell for noise.
     * A search of +-0.35 ms and +-3.5 kHz has ~4,000 cells: noise reaches ~0.3.
     */
    struct channel_candidate {
        double p_uw;      // sample index of x (fractional)
        double frequency; // Hz, absolute
        float rho2;
        float magnitude; // dB, mean power at the UW symbols over the noise: (S + N) / N
    };
    int window_capacity(bool long_window = false) const { return d_nw[long_window]; }
    void set_window(const gr_complex* x, int len);
    // Noise floor of the last set_window() in dBFS/Hz, like noise()
    float window_noise() const { return d_window_noise_db; }
    bool search_channel(int type,
                        bool uplink,
                        double f_hz,
                        double search_hz,
                        double p_first,
                        double p_last,
                        float threshold,
                        channel_candidate* out);

    /*
     * The samples extract() reads: x[first, end) with p_uw relative to x.
     */
    void extract_span(int type, double p_uw, long* first, long* end) const;

    /*
     * Matched filter (RRC alpha 0.2, +-EXTRACT_SPAN symbols) at the centres p_uw + k *
     * sps of symbols k = 0 .. n_symbols, after mixing freq_hz (absolute) down to 0 Hz.
     * Writes n_symbols + 1 samples (one per symbol) to out, rotated so the UW
     * correlation has phase 0. Fractional sps and p_uw: polyphase taps at 1/128
     * sample, the taps of each phase rotated once per call (at 12 MS/s sps is 400,
     * 200, 50: one phase; at 10 MS/s 333.33, 166.67, 41.67: three phases). x must hold
     * x[first, end) of extract_span(). Returns the number of samples written.
     *
     * The PDU for tch_c_demod: these samples, `sample_rate` = symbols_per_second of
     * the type, `center_frequency` = freq_hz, `uw_start` = 0, `timestamp` = the time
     * of p_uw (rounded to ns), see the top of this file.
     */
    int extract(const gr_complex* x,
                double p_uw,
                int type,
                double freq_hz,
                bool uplink,
                gr_complex* out);

private:
    static const int PHASES = 128; // polyphase resolution of extract()

    // Tabulated UW correlation response over the frequency offset in bins
    static const int R_RANGE = 2;
    static const int R_STEPS = 128;

    struct dir_tables {
        volk::vector<gr_complex> u_conj; // conj of the unit energy UW template
        float u2h;                       // sum |FFT(u w)|^2 H
        std::vector<float> r_table;      // R(d), d = i / R_STEPS bins
        std::vector<float> q_table;      // q(d), d = -0.5 .. 0.5
        float response(double d) const;
    };

    struct type_tables {
        const tch_c_format* f;
        double sps;
        int L, D, N, K, h;
        dir_tables dir[2];
        volk::vector<float> w; // energy window under the UW
        volk::vector<float> H; // 2K + 1: raised cosine spectrum
        float sum_h;
        float sum_w2;
        std::unique_ptr<gr::fft::fft_complex_fwd> fft_c, fft_p;
        volk::vector<float> p_shifted;  // |FFT(seg w)|^2, bin -N/2 at index 0
        volk::vector<gr_complex> u_rot; // scratch for fine_timing
        // extract(): taps for m = -M .. M + 1 around floor(t), per phase
        int M;
        int n_taps;
        std::vector<volk::vector<float>> poly;
        std::vector<volk::vector<gr_complex>> rotated;
        std::vector<int> rotated_stamp;
        volk::vector<gr_complex> phasor;
        int stamp;
    };

    double d_sample_rate;
    double d_center_frequency;
    type_tables d_t[3];
    float d_noise_db;
    std::vector<float> d_noise_tmp;

    // Channel search: window FFT, IFFT per channel, DFT over the UW symbols
    static const int CS_DFT = 128;
    int d_nw[2];
    int d_wi; // window of the last set_window()
    int d_window_len;
    float d_window_noise; // mean noise power per bin of the window FFT
    float d_window_noise_db;
    std::unique_ptr<gr::fft::fft_complex_fwd> d_fft_w[2];
    std::unique_ptr<gr::fft::fft_complex_rev> d_ifft_ch[2][3];
    std::unique_ptr<gr::fft::fft_complex_fwd> d_fft_uw;
    std::vector<float> d_g; // matched filter spectrum, scratch

    type_tables& t(int type) { return d_t[type == 8 ? 2 : type == 2 ? 1 : 0]; }
    const type_tables& t(int type) const
    {
        return d_t[type == 8 ? 2 : type == 2 ? 1 : 0];
    }
    void init_type(type_tables& tt, const tch_c_format& f);
    void interpolate_peak(const dir_tables& dt,
                          float a,
                          float b,
                          float c,
                          float* offset,
                          float* peak) const;
    float correlation_magnitude(const gr_complex* seg, const type_tables& tt, int shift);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_SLOT_SURVEY_H */
