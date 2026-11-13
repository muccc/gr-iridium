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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "iridium.h"

#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>

#include "burst_downmix_impl.h"

#include <gnuradio/filter/firdes.h>
#include <inttypes.h>
#include <volk/volk.h>

namespace gr {
namespace iridium {


void write_data_c(const gr_complex* data, size_t len, char* name, int num)
{
    char filename[256];
    sprintf(filename, "/tmp/signals/%s-%d.cfile", name, num);
    FILE* fp = fopen(filename, "wb");
    fwrite(data, sizeof(gr_complex), len, fp);
    fclose(fp);
}

void write_data_f(const float* data, size_t len, char* name, int num)
{
    char filename[256];
    sprintf(filename, "/tmp/signals/%s-%d.f32", name, num);
    FILE* fp = fopen(filename, "wb");
    fwrite(data, sizeof(float), len, fp);
    fclose(fp);
}

float sinc(float x)
{
    if (x == 0) {
        return 1;
    }
    return sin(M_PI * x) / (M_PI * x);
}

float raised_cosine_h(float t, float Ts, float alpha)
{
    if (fabs(t) == Ts / (2 * alpha)) {
        return M_PI / (4 * Ts) * sinc(1 / (2 * alpha));
    }

    return 1 / Ts * sinc(t / Ts) * cos(M_PI * alpha * t / Ts) /
           (1 - powf((2 * alpha * t / Ts), 2));
}

burst_downmix::sptr burst_downmix::make(int sample_rate,
                                        int search_depth,
                                        size_t hard_max_queue_len,
                                        const std::vector<float>& input_taps,
                                        const std::vector<float>& start_finder_taps,
                                        bool handle_multiple_frames_per_burst,
                                        bool wait_for_demod)
{
    return gnuradio::get_initial_sptr(
        new burst_downmix_impl(sample_rate,
                               search_depth,
                               hard_max_queue_len,
                               input_taps,
                               start_finder_taps,
                               handle_multiple_frames_per_burst,
                               wait_for_demod));
}


/*
 * The private constructor
 */
burst_downmix_impl::burst_downmix_impl(int output_sample_rate,
                                       int search_depth,
                                       size_t hard_max_queue_len,
                                       const std::vector<float>& input_taps,
                                       const std::vector<float>& start_finder_taps,
                                       bool handle_multiple_frames_per_burst,
                                       bool wait_for_demod)
    : gr::block("burst_downmix",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(0, 0, 0)),
      d_symbols_per_second(::iridium::SYMBOLS_PER_SECOND),
      d_output_sample_rate(output_sample_rate),
      d_output_samples_per_symbol(d_output_sample_rate / d_symbols_per_second),
      d_samples_per_symbol((float)d_output_sample_rate / d_symbols_per_second),
      d_burst_sample_rate(d_output_sample_rate),
      d_max_burst_size(0),
      d_search_depth(search_depth),
      d_pre_start_samples(pre_start_samples()),
      d_n_dropped_bursts(0),
      d_debug_id(-1),

      // Take the FFT over the (short) preamble + 10 symbols from the unique word (UW)
      // (Frames with a 64 symbol preamble will use 26 symbols of the preamble)
      d_cfo_est_fft_size(pow(
          2,
          int(log(d_output_samples_per_symbol * (::iridium::PREAMBLE_LENGTH_SHORT + 10)) /
              log(2)))),

      d_fft_over_size_facor(16),
      d_sync_search_len(sync_search_len()),
      d_hard_max_queue_len(hard_max_queue_len),
      d_handle_multiple_frames_per_burst(handle_multiple_frames_per_burst),
      d_debug(false),

      d_input(NULL),
      d_frame(NULL),
      d_tmp_a(NULL),
      d_tmp_b(NULL),
      d_dl_preamble_reversed_conj_fft(NULL),
      d_ul_preamble_reversed_conj_fft(NULL),

      d_magnitude_f(NULL),
      d_magnitude_filtered_f(NULL),
      d_cfo_est_window_f(NULL),

      d_corr_fft(NULL),
      d_corr_dl_ifft(NULL),
      d_corr_ul_ifft(NULL),

      d_input_fir(input_taps),
      d_start_finder_fir(start_finder_taps),
      d_rrc_fir(gr::filter::firdes::root_raised_cosine(
          1.0, d_output_sample_rate, d_symbols_per_second, RRC_ALPHA, RRC_TAPS)),
      d_cfo_est_fft(d_cfo_est_fft_size),
      d_wait_for_demod(wait_for_demod),
      d_burst_id(0),
      d_n_frames(0)
{
    d_dl_preamble_reversed_conj = generate_sync_word(::iridium::direction::DOWNLINK);
    d_ul_preamble_reversed_conj = generate_sync_word(::iridium::direction::UPLINK);

    initialize_cfo_est_fft();

    initialize_correlation_filter();

    message_port_register_in(pmt::mp("cpdus"));
    message_port_register_out(pmt::mp("cpdus"));
    message_port_register_out(pmt::mp("burst_handled"));

    set_msg_handler(pmt::mp("cpdus"), [this](pmt::pmt_t msg) { this->handler(msg); });

    if (d_debug) {
        std::cout << "Start filter size:" << d_start_finder_fir.ntaps()
                  << " Search depth:" << d_search_depth << "\n";
    }
}

/*
 * Our virtual destructor.
 */
burst_downmix_impl::~burst_downmix_impl()
{
    if (d_input) {
        volk_free(d_input);
    }
    if (d_frame) {
        volk_free(d_frame);
    }
    if (d_tmp_a) {
        volk_free(d_tmp_a);
    }
    if (d_tmp_b) {
        volk_free(d_tmp_b);
    }
    if (d_magnitude_f) {
        volk_free(d_magnitude_f);
    }
    if (d_magnitude_filtered_f) {
        volk_free(d_magnitude_filtered_f);
    }
    if (d_cfo_est_window_f) {
        free(d_cfo_est_window_f);
    }
    free_correlation_filter();
}

void burst_downmix_impl::free_correlation_filter(void)
{
    volk_free(d_dl_preamble_reversed_conj_fft);
    d_dl_preamble_reversed_conj_fft = NULL;
    volk_free(d_ul_preamble_reversed_conj_fft);
    d_ul_preamble_reversed_conj_fft = NULL;
    delete d_corr_fft;
    d_corr_fft = NULL;
    delete d_corr_dl_ifft;
    d_corr_dl_ifft = NULL;
    delete d_corr_ul_ifft;
    d_corr_ul_ifft = NULL;
}

// Start the burst this many samples before the detected rise of the magnitude
int burst_downmix_impl::pre_start_samples() const
{
    return int(0.1e-3 * d_output_sample_rate);
}

// Number of samples searched for the sync word
int burst_downmix_impl::sync_search_len() const
{
    return (::iridium::PREAMBLE_LENGTH_LONG + ::iridium::UW_LENGTH + 8) *
           d_samples_per_symbol;
}

/*
 * Bursts can arrive at a rate which is not a multiple of the symbol rate. Run the filters
 * and the sync word correlation at the actual rate with a fractional number of samples
 * per symbol then. The demodulator interpolates the symbols accordingly. Only done when
 * the rate changes.
 */
void burst_downmix_impl::set_burst_sample_rate(float sample_rate)
{
    if (sample_rate == d_burst_sample_rate) {
        return;
    }

    d_burst_sample_rate = sample_rate;
    d_samples_per_symbol = sample_rate / d_symbols_per_second;
    d_pre_start_samples = pre_start_samples();
    d_sync_search_len = sync_search_len();

    d_rrc_fir.set_taps(gr::filter::firdes::root_raised_cosine(
        1.0, sample_rate, d_symbols_per_second, RRC_ALPHA, RRC_TAPS));

    d_dl_preamble_reversed_conj = generate_sync_word(::iridium::direction::DOWNLINK);
    d_ul_preamble_reversed_conj = generate_sync_word(::iridium::direction::UPLINK);
    free_correlation_filter();
    initialize_correlation_filter();
}

/*
 * The correlation template for the start of a burst: the 16 symbol preamble followed by
 * the 12 symbol UW.
 * Built from raised cosine pulses at the (possibly fractional) symbol positions,
 * truncated to RC_TAPS samples.
 */
volk::vector<gr_complex>
burst_downmix_impl::generate_sync_word(::iridium::direction direction)
{
    gr_complex s1 = gr_complex(-1, -1);
    gr_complex s0 = -s1;
    std::vector<gr_complex> sync_word;
    int i;

    if (direction == ::iridium::direction::DOWNLINK) {
        for (i = 0; i < ::iridium::PREAMBLE_LENGTH_SHORT; i++) {
            sync_word.push_back(s0);
        }
        for (i = 0; i < ::iridium::UW_LENGTH; i++) {
            sync_word.push_back(::iridium::UW_DL[i] ? s1 : s0);
        }
    } else if (direction == ::iridium::direction::UPLINK) {
        for (i = 0; i < ::iridium::PREAMBLE_LENGTH_SHORT; i += 2) {
            sync_word.push_back(s1);
            sync_word.push_back(s0);
        }
        for (i = 0; i < ::iridium::UW_LENGTH; i++) {
            sync_word.push_back(::iridium::UW_UL[i] ? s1 : s0);
        }
    }

    const int half_rc_size = (RC_TAPS - 1) / 2;
    const float Ts = 1. / d_symbols_per_second;
    const int len = (int)((sync_word.size() - 1) * d_samples_per_symbol) + 1;
    volk::vector<gr_complex> sync_word_padded(len);

    for (int n = 0; n < len; n++) {
        gr_complex v = 0;
        for (size_t k = 0; k < sync_word.size(); k++) {
            // The span of the 51 tap filter the legacy template was made with:
            // half_rc_size samples after the symbol and one less before it
            float x = n - k * d_samples_per_symbol;
            if (x > -half_rc_size && x <= half_rc_size) {
                v +=
                    sync_word[k] * raised_cosine_h(x / d_burst_sample_rate, Ts, 0.4) * Ts;
            }
        }
        sync_word_padded[n] = v;
    }

    if (d_debug) {
        std::cout << "Sync Word Unpadded: ";
        for (gr_complex s : sync_word) {
            std::cout << s << ", ";
        }
        std::cout << std::endl;

        std::cout << "Sync Word Padded: ";
        for (gr_complex s : sync_word_padded) {
            std::cout << s << ", ";
        }
        std::cout << std::endl;
    }

    std::reverse(sync_word_padded.begin(), sync_word_padded.end());
    volk_32fc_conjugate_32fc(
        &sync_word_padded[0], &sync_word_padded[0], sync_word_padded.size());
    return sync_word_padded;
}

void burst_downmix_impl::initialize_cfo_est_fft(void)
{
    // See find_cfo_peak()
    const int fine_size = d_cfo_est_fft_size * d_fft_over_size_facor;
    d_cfo_fine_phase.resize(fine_size);
    for (int i = 0; i < fine_size; i++) {
        d_cfo_fine_phase[i] = std::polar(1.0, -2 * M_PI * i / fine_size);
    }
    d_cfo_fine_step.resize(d_cfo_est_fft_size);
    d_cfo_fine_row.resize(d_cfo_est_fft_size);
    for (int i = 0; i < d_cfo_est_fft_size; i++) {
        d_cfo_fine_step[i] = d_cfo_fine_phase[i];
    }

    // Compute window and move it into aligned buffer
    std::vector<float> window =
        fft::window::build(fft::window::WIN_BLACKMAN, d_cfo_est_fft_size, 0);
    d_cfo_est_window_f =
        (float*)volk_malloc(sizeof(float) * d_cfo_est_fft_size, volk_get_alignment());
    memcpy(d_cfo_est_window_f, &window[0], sizeof(float) * d_cfo_est_fft_size);

    if (d_debug) {
        printf("fft_length=%d (%d)\n",
               d_cfo_est_fft_size,
               d_output_samples_per_symbol * (::iridium::PREAMBLE_LENGTH_SHORT + 10));
    }
}

/*
 * The peak of the spectrum of the windowed samples in the input buffer of d_cfo_est_fft,
 * at d_fft_over_size_facor times the resolution of that FFT: what a zero padded FFT of
 * that many times the size would give, but only the bins around the peak of the short
 * FFT are computed. Returns the index of the peak in [-N/2 .. (N/2)-1] notation and
 * leaves the magnitudes of the bin before, at and after the peak in d_cfo_peak.
 */
int burst_downmix_impl::find_cfo_peak(void)
{
    const int size = d_cfo_est_fft_size;
    const int over = d_fft_over_size_facor;
    const int fine_size = size * over;

    d_cfo_est_fft.execute();
    volk_32fc_magnitude_squared_32f(d_magnitude_f, d_cfo_est_fft.get_outbuf(), size);
    const float coarse_max = *std::max_element(d_magnitude_f, d_magnitude_f + size);

    // The peak is within one bin of a local maximum of the short FFT, which shows at
    // least 78% of it (the scalloping loss of the Blackman window is 1.1 dB)
    const int n = 2 * over + 3;
    float magnitude[2 * 64 + 3];
    const gr_complex* in = d_cfo_est_fft.get_inbuf();
    int max_index = 0;
    d_cfo_peak[1] = -1;
    for (int coarse = 0; coarse < size; coarse++) {
        if (d_magnitude_f[coarse] < 0.75f * coarse_max ||
            d_magnitude_f[coarse] < d_magnitude_f[(coarse + size - 1) % size] ||
            d_magnitude_f[coarse] < d_magnitude_f[(coarse + 1) % size]) {
            continue;
        }

        // Fine bins from one coarse bin below to one above, and one more on each side
        const int first = coarse * over - over - 1;
        for (int i = 0; i < size; i++) {
            d_cfo_fine_row[i] =
                d_cfo_fine_phase[((long)(first + fine_size) * i) % fine_size];
        }
        for (int i = 0; i < n; i++) {
            gr_complex bin;
            volk_32fc_x2_dot_prod_32fc(&bin, in, &d_cfo_fine_row[0], size);
            magnitude[i] = std::norm(bin);
            volk_32fc_x2_multiply_32fc(
                &d_cfo_fine_row[0], &d_cfo_fine_row[0], &d_cfo_fine_step[0], size);
        }

        const int peak = std::max_element(magnitude + 1, magnitude + n - 1) - magnitude;
        if (magnitude[peak] <= d_cfo_peak[1]) {
            continue;
        }
        max_index = fft_unshift_index((first + peak + fine_size) % fine_size, fine_size);

        // The spectrum does not wrap around for the interpolation
        d_cfo_peak[0] =
            max_index == -fine_size / 2 ? magnitude[peak] : magnitude[peak - 1];
        d_cfo_peak[1] = magnitude[peak];
        d_cfo_peak[2] =
            max_index == fine_size / 2 - 1 ? magnitude[peak] : magnitude[peak + 1];
    }
    return max_index;
}

void burst_downmix_impl::initialize_correlation_filter(void)
{
    // Based on code from synchronizer_v4_impl.cc in gr-burst

    // Make the FFT size a power of two
    int corr_fft_size_target = d_sync_search_len + d_dl_preamble_reversed_conj.size() - 1;
    d_corr_fft_size = pow(2, (int)(std::ceil(log(corr_fft_size_target) / log(2))));

    // TODO: We could increase the search size for free
    // d_sync_search_len = d_corr_fft_size - d_dl_preamble_reversed_conj.size() + 1;

    if (d_debug) {
        std::cout << "Conv FFT size:" << d_corr_fft_size << std::endl;
    }

    // Allocate space for the pre transformed filters
    d_dl_preamble_reversed_conj_fft = (gr_complex*)volk_malloc(
        d_corr_fft_size * sizeof(gr_complex), volk_get_alignment());
    d_ul_preamble_reversed_conj_fft = (gr_complex*)volk_malloc(
        d_corr_fft_size * sizeof(gr_complex), volk_get_alignment());

    // Temporary FFT to pre transform the filters
    fft::fft_complex_fwd fft_engine(d_corr_fft_size);
    memset(fft_engine.get_inbuf(), 0, sizeof(gr_complex) * d_corr_fft_size);

    int sync_word_len = d_dl_preamble_reversed_conj.size();

    // Transform the filters
    memcpy(fft_engine.get_inbuf(),
           &d_dl_preamble_reversed_conj[0],
           sizeof(gr_complex) * sync_word_len);
    fft_engine.execute();
    memcpy(d_dl_preamble_reversed_conj_fft,
           fft_engine.get_outbuf(),
           sizeof(gr_complex) * d_corr_fft_size);

    memcpy(fft_engine.get_inbuf(),
           &d_ul_preamble_reversed_conj[0],
           sizeof(gr_complex) * sync_word_len);
    fft_engine.execute();
    memcpy(d_ul_preamble_reversed_conj_fft,
           fft_engine.get_outbuf(),
           sizeof(gr_complex) * d_corr_fft_size);

    // Update the size of the work FFTs
    // TODO: This could be moved to the initialization list
    d_corr_fft = new fft::fft_complex_fwd(d_corr_fft_size);
    d_corr_dl_ifft = new fft::fft_complex_rev(d_corr_fft_size);
    d_corr_ul_ifft = new fft::fft_complex_rev(d_corr_fft_size);

    // The inputs need to zero, as we might not use it completely
    memset(d_corr_fft->get_inbuf(), 0, sizeof(gr_complex) * d_corr_fft_size);
}

void burst_downmix_impl::update_buffer_sizes(size_t burst_size)
{
    if (burst_size > d_max_burst_size) {
        d_max_burst_size = burst_size;
        if (d_input) {
            volk_free(d_input);
        }
        d_input = (gr_complex*)volk_malloc(d_max_burst_size * sizeof(gr_complex),
                                           volk_get_alignment());

        if (d_frame) {
            volk_free(d_frame);
        }
        d_frame = (gr_complex*)volk_malloc(d_max_burst_size * sizeof(gr_complex),
                                           volk_get_alignment());

        if (d_tmp_a) {
            volk_free(d_tmp_a);
        }
        // Room for the padding of the matched filter at both ends
        d_tmp_a = (gr_complex*)volk_malloc((d_max_burst_size + 4096) * sizeof(gr_complex),
                                           volk_get_alignment());

        if (d_tmp_b) {
            volk_free(d_tmp_b);
        }
        // Room for the padding of the matched filter at both ends
        d_tmp_b = (gr_complex*)volk_malloc((d_max_burst_size + 4096) * sizeof(gr_complex),
                                           volk_get_alignment());

        if (d_magnitude_f) {
            volk_free(d_magnitude_f);
        }
        d_magnitude_f =
            (float*)volk_malloc(d_max_burst_size * sizeof(float), volk_get_alignment());

        if (d_magnitude_filtered_f) {
            volk_free(d_magnitude_filtered_f);
        }
        d_magnitude_filtered_f =
            (float*)volk_malloc(d_max_burst_size * sizeof(float), volk_get_alignment());
    }
}

/*
 * The burst is shifted to the rough CFO estimate, low pass filtered and decimated on
 * demand: a burst can hold several frames, and only the parts looked at are filtered.
 */

// Shift the first n samples of the burst
void burst_downmix_impl::shift_input(size_t n)
{
    if (n > d_input_valid) {
        d_input_r.rotateN(
            d_input + d_input_valid, d_burst + d_input_valid, n - d_input_valid);
        d_input_valid = n;
    }
}

// Make the samples [from, to) of d_frame valid. Samples before from are not asked for
// again.
void burst_downmix_impl::filter_input(size_t from, size_t to)
{
    to = std::min(to, d_frame_size);
    d_frame_valid = std::max(d_frame_valid, from);
    if (to <= d_frame_valid) {
        return;
    }

    shift_input((to - 1) * d_decimation + d_input_fir.ntaps());
    d_input_fir.filterNdec(d_frame + d_frame_valid,
                           d_input + d_frame_valid * d_decimation,
                           to - d_frame_valid,
                           d_decimation);
    d_frame_valid = to;
}

size_t burst_downmix_impl::get_input_queue_size() { return nmsgs(pmt::mp("cpdus")); }

uint64_t burst_downmix_impl::get_n_dropped_bursts() { return d_n_dropped_bursts; }

void burst_downmix_impl::debug_id(uint64_t id) { d_debug_id = id; }

// Maps an index in [-N/2 .. (N/2)-1] notation to [0 .. N-1] notation
int burst_downmix_impl::fft_shift_index(int index, int fft_size)
{
    // Clamp the input to [-N/2 .. (N/2)-1]
    index = std::max(index, -fft_size / 2);
    index = std::min(index, fft_size / 2 - 1);

    if (index < 0) {
        index += fft_size;
    }
    return index;
}

// Maps an index in [0 .. N-1] notation to [-N/2 .. (N/2)-1] notation
int burst_downmix_impl::fft_unshift_index(int index, int fft_size)
{
    // Clamp the input to [0 .. N-1]
    index = std::max(index, 0);
    index = std::min(index, fft_size - 1);

    if (index >= fft_size / 2) {
        index -= fft_size;
    }
    return index;
}

float burst_downmix_impl::interpolate(float alpha, float beta, float gamma)
{
    const float correction = 0.5 * (alpha - gamma) / (alpha - 2 * beta + gamma);
    return correction;
}

int burst_downmix_impl::process_next_frame(float sample_rate,
                                           double center_frequency,
                                           uint64_t timestamp,
                                           uint64_t sub_id,
                                           size_t burst_size,
                                           int start,
                                           float noise,
                                           float magnitude)
{
    /*
     * Use the center frequency to make some assumptions about the burst.
     */
    int max_frame_length = 0;
    int min_frame_length = 0;

    if (center_frequency > ::iridium::SIMPLEX_FREQUENCY_MIN) {
        // Simplex transmissions and broadcast frames might have a 64 symbol preamble.
        // We ignore that and cut away the extra 48 symbols.
        // Frames above this frequency must be downlink and simplex frames.
        // XXX: If the SDR is not configured well, there might be aliasing from low
        // frequencies in this region.
        max_frame_length =
            ::iridium::MAX_FRAME_LENGTH_SIMPLEX * d_output_samples_per_symbol;
        min_frame_length =
            (::iridium::MIN_FRAME_LENGTH_SIMPLEX)*d_output_samples_per_symbol;
    } else {
        max_frame_length =
            ::iridium::MAX_FRAME_LENGTH_NORMAL * d_output_samples_per_symbol;
        min_frame_length =
            (::iridium::MIN_FRAME_LENGTH_NORMAL)*d_output_samples_per_symbol;
    }

    if (burst_size - start < min_frame_length) {
        return 0;
    }

    /*
     * Find the fine CFO estimate using an FFT over the preamble and the first symbols
     * of the unique word.
     * The signal gets squared to remove the BPSK modulation from the unique word.
     */

    if (burst_size - start < d_cfo_est_fft_size) {
        // There are not enough samples available to run the FFT.
        return 0;
    }

    /*
     * Only the samples searched for the sync word are RRC filtered until it is found,
     * then the rest of the frame: not the whole burst, which can hold several frames.
     */
    const int available = burst_size - start;
    const int half_rrc_size = (d_rrc_fir.ntaps() - 1) / 2;
    const int search_len = std::min(d_sync_search_len, available);
    // The UW follows the preamble
    const int uw_offset = ::iridium::PREAMBLE_LENGTH_SHORT * d_output_samples_per_symbol;
    int n_search = available;   // samples shifted to the fine CFO estimate
    int n_filtered = available; // samples RRC filtered
    int n_frame = available;    // samples needed for the longest frame
    if (!d_debug) {
        n_filtered = search_len;
        n_frame = std::min(available,
                           search_len + uw_offset + max_frame_length + half_rrc_size);
        n_search = n_frame;
    }
    filter_input(start, start + n_search);

    volk_32fc_x2_multiply_32fc(
        d_tmp_a, d_frame + start, d_frame + start, d_cfo_est_fft_size);
    volk_32fc_32f_multiply_32fc(
        d_cfo_est_fft.get_inbuf(), d_tmp_a, d_cfo_est_window_f, d_cfo_est_fft_size);
    const int max_index = find_cfo_peak();
    if (d_debug) {
        printf("max_index=%d\n", max_index);
    }

    // Interpolate the result of the FFT to get a finer resolution.
    // see
    // http://www.dsprelated.com/dspbooks/sasp/Quadratic_Interpolation_Spectral_Peaks.html
    // TODO: The window should be Gaussian and the output should be put on a log scale
    // https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html

    const float alpha = d_cfo_peak[0];
    const float beta = d_cfo_peak[1];
    const float gamma = d_cfo_peak[2];
    const float interpolated_index = max_index + interpolate(alpha, beta, gamma);

    // Normalize the result.
    // Divide by two to remove the effect of the squaring operation before.
    float center_offset =
        interpolated_index / (d_cfo_est_fft_size * d_fft_over_size_facor) / 2;

    if (d_debug) {
        printf("interpolated_index=%f center_offset=%f (%f)\n",
               interpolated_index,
               center_offset,
               center_offset * d_output_sample_rate);
    }


    /*
     * Shift the burst again using the result of the FFT.
     */
    float phase_inc = 2 * M_PI * -center_offset;
    d_r.set_phase_incr(exp(gr_complex(0, phase_inc)));
    d_r.set_phase(gr_complex(1, 0));
    // With some room at the start and the end, so the RRC can run
    gr_complex* shifted = d_tmp_b + half_rrc_size;
    memset(d_tmp_b, 0, half_rrc_size * sizeof(gr_complex));
    d_r.rotateN(shifted, d_frame + start, n_search);
    memset(shifted + n_search, 0, half_rrc_size * sizeof(gr_complex));
    center_frequency += center_offset * sample_rate;

    if (d_debug) {
        write_data_c(
            shifted, n_search, (char*)"signal-filtered-deci-cut-start-shift", sub_id);
    }

    /*
     * Apply the RRC filter.
     */
    d_rrc_fir.filterN(d_tmp_a, d_tmp_b, n_filtered);

    if (d_debug) {
        write_data_c(d_tmp_a,
                     n_filtered,
                     (char*)"signal-filtered-deci-cut-start-shift-rrc",
                     sub_id);
    }

    /*
     * Use a correlation to find the start of the sync word.
     * Uses an FFT to perform the correlation.
     */

    // The burst might be shorter than the search window
    memcpy(d_corr_fft->get_inbuf(), d_tmp_a, sizeof(gr_complex) * search_len);
    memset(d_corr_fft->get_inbuf() + search_len,
           0,
           sizeof(gr_complex) * (d_corr_fft_size - search_len));
    d_corr_fft->execute();

    // We use the initial FFT for both correlations (DL and UL)
    volk_32fc_x2_multiply_32fc(d_corr_dl_ifft->get_inbuf(),
                               d_corr_fft->get_outbuf(),
                               &d_dl_preamble_reversed_conj_fft[0],
                               d_corr_fft_size);
    volk_32fc_x2_multiply_32fc(d_corr_ul_ifft->get_inbuf(),
                               d_corr_fft->get_outbuf(),
                               &d_ul_preamble_reversed_conj_fft[0],
                               d_corr_fft_size);
    d_corr_dl_ifft->execute();
    d_corr_ul_ifft->execute();

    int corr_offset_dl;
    int corr_offset_ul;

    float correction_dl = 0;
    float correction_ul = 0;

    // Find the peaks of the correlations
    volk_32fc_magnitude_squared_32f(
        d_magnitude_f, d_corr_dl_ifft->get_outbuf(), d_corr_fft_size);
    float* max_dl_p = std::max_element(d_magnitude_f, d_magnitude_f + d_corr_fft_size);
    corr_offset_dl = max_dl_p - d_magnitude_f;
    if (corr_offset_dl > 0) {
        correction_dl = interpolate(d_magnitude_f[corr_offset_dl - 1],
                                    d_magnitude_f[corr_offset_dl],
                                    d_magnitude_f[corr_offset_dl + 1]);
    }
    float max_dl = *max_dl_p;

    volk_32fc_magnitude_squared_32f(
        d_magnitude_f, d_corr_ul_ifft->get_outbuf(), d_corr_fft_size);
    float* max_ul_p = std::max_element(d_magnitude_f, d_magnitude_f + d_corr_fft_size);
    corr_offset_ul = max_ul_p - d_magnitude_f;
    if (corr_offset_ul > 0) {
        correction_ul = interpolate(d_magnitude_f[corr_offset_ul - 1],
                                    d_magnitude_f[corr_offset_ul],
                                    d_magnitude_f[corr_offset_ul + 1]);
    }
    float max_ul = *max_ul_p;

    gr_complex corr_result;
    int corr_offset;
    float correction;
    ::iridium::direction direction;

    if (max_dl > max_ul) {
        direction = ::iridium::direction::DOWNLINK;
        corr_offset = corr_offset_dl;
        correction = correction_dl;
        corr_result = d_corr_dl_ifft->get_outbuf()[corr_offset];
    } else {
        direction = ::iridium::direction::UPLINK;
        corr_offset = corr_offset_ul;
        correction = correction_ul;
        corr_result = d_corr_ul_ifft->get_outbuf()[corr_offset];
    }


    if (d_debug) {
        printf("Conv max index = %d\n", corr_offset);
    }

    // Careful: The correlation might have found the start of the sync word
    // before the first sample => preamble_offset might be negative
    int preamble_offset = corr_offset - d_dl_preamble_reversed_conj.size() + 1;
    int uw_start = preamble_offset + uw_offset;

    // If the UW starts at an offset < 0, we will not be able to demodulate the signal
    if (uw_start < 0) {
        return 0;
    }

    size_t frame_size = std::min((int)burst_size - start, uw_start + max_frame_length);
    int consumed_samples = frame_size;

    // The rest of the frame
    if (n_filtered < (int)frame_size) {
        filter_input(start, start + n_frame);
        d_r.rotateN(shifted + n_search, d_frame + start + n_search, n_frame - n_search);
        memset(shifted + n_frame, 0, half_rrc_size * sizeof(gr_complex));
        d_rrc_fir.filterN(
            d_tmp_a + n_filtered, d_tmp_b + n_filtered, frame_size - n_filtered);
    }

    /*
     * Rotate the phase so the demodulation has a starting point.
     */
    d_r.set_phase_incr(exp(gr_complex(0, 0)));
    d_r.set_phase(std::conj(corr_result / abs(corr_result)));
    d_r.rotateN(d_tmp_b, d_tmp_a, frame_size);

    if (d_debug) {
        write_data_c(d_tmp_b,
                     frame_size,
                     (char*)"signal-filtered-deci-cut-start-shift-rrc-rotate",
                     sub_id);
    }

    /*
     * Align the burst so the first sample of the burst is the first symbol
     * of the UW after the RRC filter.
     *
     */

    // Update the amount of available samples after filtering
    frame_size = std::max((int)frame_size - uw_start, 0);
    start += uw_start;

    if (d_debug) {
        write_data_c(d_tmp_b + uw_start,
                     frame_size,
                     (char*)"signal-filtered-deci-cut-start-shift-rrc-rotate-cut",
                     sub_id);
    }

    timestamp += start * 1000000000ULL / (int)sample_rate;
    /*
     * Done :)
     */
    pmt::pmt_t pdu_meta = pmt::make_dict();
    pmt::pmt_t pdu_vector = pmt::init_c32vector(frame_size, d_tmp_b + uw_start);

    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("sample_rate"), pmt::mp(sample_rate));
    pdu_meta =
        pmt::dict_add(pdu_meta, pmt::mp("center_frequency"), pmt::mp(center_frequency));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("direction"), pmt::mp((int)direction));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("uw_start"), pmt::mp(correction));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("timestamp"), pmt::mp(timestamp));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("id"), pmt::mp(sub_id));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("noise"), pmt::mp(noise));
    pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("magnitude"), pmt::mp(magnitude));
    if (d_wait_for_demod) {
        // The demodulator reports this on its frame_handled port
        pdu_meta = pmt::dict_add(pdu_meta, pmt::mp("burst_id"), pmt::mp(d_burst_id));
    }

    if (d_debug) {
        printf("center_frequency=%f, uw_start=%u\n", center_frequency, uw_start);
    }

    pmt::pmt_t out_msg = pmt::cons(pdu_meta, pdu_vector);
    message_port_pub(pmt::mp("cpdus"), out_msg);
    d_n_frames++;

    return consumed_samples;
}

/*
 * Acknowledge a burst to the burst source. With wait_for_demod the acknowledgement
 * carries the number of frames published for the burst: (id . n_frames). The burst
 * source then also waits for the demodulator to report these frames as handled.
 */
void burst_downmix_impl::burst_done(uint64_t id)
{
    if (d_wait_for_demod) {
        message_port_pub(pmt::mp("burst_handled"),
                         pmt::cons(pmt::from_uint64(id), pmt::from_long(d_n_frames)));
    } else {
        message_port_pub(pmt::mp("burst_handled"), pmt::mp(id));
    }
}

/*
 * Search for the start of the burst by looking at the magnitude.
 * Look at most d_search_depth far.
 */
int burst_downmix_impl::find_start(size_t burst_size, uint64_t id)
{
    auto& fir = d_start_finder_fir;
    int fir_size = fir.ntaps();
    int half_fir_size = (fir_size - 1) / 2;

    // The burst might be shorter than d_search_depth.
    filter_input(0, burst_size);
    int N = std::min(d_search_depth, (int)burst_size - (fir_size - 1));

    volk_32fc_magnitude_squared_32f(d_magnitude_f, d_frame, N + fir_size - 1);

    if (d_debug) {
        write_data_f(d_magnitude_f, N + fir_size - 1, (char*)"signal-mag", id);
    }

    fir.filterN(d_magnitude_filtered_f, d_magnitude_f, N);

    if (d_debug) {
        write_data_f(d_magnitude_filtered_f, N, (char*)"signal-mag-filter", id);
    }

    float* max = std::max_element(d_magnitude_filtered_f, d_magnitude_filtered_f + N);
    float threshold = *max * 0.28;
    if (d_debug) {
        std::cout << "Threshold:" << threshold << " Max:" << *max << "("
                  << (max - d_magnitude_filtered_f) << ")\n";
    }

    int start;
    for (start = 0; start < N; start++) {
        if (d_magnitude_filtered_f[start] >= threshold) {
            break;
        }
    }

    if (start > 0) {
        start = std::max(start + half_fir_size - d_pre_start_samples, 0);
    }
    return start;
}

void burst_downmix_impl::handler(pmt::pmt_t msg)
{
    /*
     * Extract the burst and meta data from the cpdu
     */
    pmt::pmt_t samples = pmt::cdr(msg);
    size_t burst_size = pmt::length(samples);
    const gr_complex* burst =
        (const gr_complex*)pmt::c32vector_elements(samples, burst_size);

    pmt::pmt_t meta = pmt::car(msg);
    float relative_frequency =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("relative_frequency"), pmt::PMT_NIL));
    double center_frequency =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
    float sample_rate =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("sample_rate"), pmt::PMT_NIL));
    uint64_t id = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("id"), pmt::PMT_NIL));
    uint64_t timestamp =
        pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    float noise = pmt::to_float(pmt::dict_ref(meta, pmt::mp("noise"), pmt::PMT_NIL));
    float magnitude =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("magnitude"), pmt::PMT_NIL));

    d_burst_id = id;
    d_n_frames = 0;

    if (id == d_debug_id) {
        d_debug = true;
    }

    if (d_debug) {
        printf("---------------> id:%" PRIu64 " len:%zu\n", id, burst_size);
        float absolute_frequency = center_frequency + relative_frequency * sample_rate;
        printf("relative_frequency=%f, absolute_frequency=%f\n",
               relative_frequency,
               absolute_frequency);
        printf("sample_rate=%f\n", sample_rate);
    }

    if (d_hard_max_queue_len && get_input_queue_size() >= d_hard_max_queue_len) {
        std::cerr << "Warning: Dropping burst as hard queue length is reached!"
                  << std::endl;
        d_n_dropped_bursts++;
        burst_done(id);
        return;
    }

    // This burst might be larger than the one before.
    // Update he buffer sizes if needed.
    update_buffer_sizes(burst_size + 10000);

    if (d_debug) {
        write_data_c(burst, burst_size, (char*)"signal", id);
    }


    /*
     * Shift the center frequency of the burst to the provided rough CFO estimate.
     */
    float phase_inc = 2 * M_PI * -relative_frequency;
    d_input_r.set_phase_incr(exp(gr_complex(0, phase_inc)));
    d_input_r.set_phase(gr_complex(1, 0));
    center_frequency += relative_frequency * sample_rate;


    /*
     * Apply the initial low pass filter and decimate the burst.
     */
    int decimation = std::max(1, (int)(std::lround(sample_rate) / d_output_sample_rate));

#if 0
      // Option for additional padding. Probably not needed.
      int input_fir_pad_size = (d_input_fir.ntaps() - 1) / 2;
      memmove(d_tmp_a + input_fir_pad_size, d_tmp_a, sizeof(gr_complex) * burst_size);
      memset(d_tmp_a, 0, sizeof(gr_complex) * input_fir_pad_size);
      memset(d_tmp_a + input_fir_pad_size + burst_size, 0, sizeof(gr_complex) * input_fir_pad_size);

      burst_size = burst_size / decimation;
#else
    burst_size = (burst_size - d_input_fir.ntaps() + 1) / decimation;

    timestamp += d_input_fir.ntaps() / 2 * 1000000000ULL / (int)sample_rate;
#endif

    // See filter_input()
    d_burst = burst;
    d_decimation = decimation;
    d_frame_size = burst_size;
    d_frame_valid = 0;
    d_input_valid = 0;

    sample_rate /= decimation;

    // Only does something if the rate is not a multiple of the symbol rate
    set_burst_sample_rate(sample_rate);

    if (d_debug) {
        printf("---------------> id:%" PRIu64 " len:%lu\n",
               id,
               burst_size / d_output_sample_rate);
        filter_input(0, burst_size);
        write_data_c(d_frame, burst_size, (char*)"signal-filtered-deci", id);
    }

    int start = find_start(burst_size, id);

    if (d_debug) {
        std::cout << "Start:" << start << "\n";
        write_data_c(d_frame + start,
                     burst_size - start,
                     (char*)"signal-filtered-deci-cut-start",
                     id);
    }

    if (d_handle_multiple_frames_per_burst) {
        int handled_samples;
        int sub_id = id;
        do {
            const int n_frames = d_n_frames;
            handled_samples = process_next_frame(sample_rate,
                                                 center_frequency,
                                                 timestamp,
                                                 sub_id,
                                                 burst_size,
                                                 start,
                                                 noise,
                                                 magnitude);
            start += handled_samples;
            // This is OK as ids are incremented by 10 by the burst tagger. Only frames
            // sent on use one.
            if (d_n_frames != n_frames) {
                sub_id++;
            }
        } while (d_handle_multiple_frames_per_burst && handled_samples > 0);
    } else {
        process_next_frame(sample_rate,
                           center_frequency,
                           timestamp,
                           id,
                           burst_size,
                           start,
                           noise,
                           magnitude);
    }

    burst_done(id);

    if (d_debug_id >= 0) {
        d_debug = false;
    }
}

int burst_downmix_impl::work(int noutput_items,
                             gr_vector_const_void_star& input_items,
                             gr_vector_void_star& output_items)
{
    return 0;
}

} /* namespace iridium */
} /* namespace gr */
