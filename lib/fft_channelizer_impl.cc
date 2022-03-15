/* -*- c++ -*- */
/*
 * Copyright 2022 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Many of this is following directly what is presented
 * here: https://media.ccc.de/v/gpn18-15-channelizing-with-gnuradio
 *
 * No special windowing is performed. This leads to a lot of noise
 * in the output. First tests indicate no significant drop in
 * overall performance though.
 */

#include "fft_channelizer_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace iridium {

constexpr bool is_powerof2(int v) {
    return v && ((v & (v - 1)) == 0);
}

using input_type = gr_complex;
using output_type = gr_complex;
fft_channelizer::sptr fft_channelizer::make(int fft_size, int decimation)
{
    return gnuradio::make_block_sptr<fft_channelizer_impl>(fft_size, decimation);
}


/*
 * The private constructor
 */
fft_channelizer_impl::fft_channelizer_impl(int fft_size, int decimation)
    : gr::sync_decimator("fft_channelizer",
                         gr::io_signature::make(
                             1 /* min inputs */, 1 /* max inputs */, sizeof(input_type)),
                         gr::io_signature::make(decimation + 1 /* min outputs */,
                                                decimation + 1 /*max outputs */,
                                                sizeof(output_type)),
                         decimation),
    d_decimation(decimation),
    d_fft_size(fft_size),
    d_ifft_size(fft_size / decimation),
    d_fft(fft_size, 1),
    d_ifft(fft_size / decimation, 1),
    d_inverse_overlap(4),
    d_output_step((fft_size - fft_size / decimation) / decimation)
{
    // FFT size and decimation must be a power of two so the output sizes are
    // also a power of two.
    if(!is_powerof2(d_fft_size)) {
        throw std::runtime_error("FFT size must be a power of two");
    }

    if(!is_powerof2(d_decimation)) {
        throw std::runtime_error("Decimation must be a power of two");
    }

    // Make sure the output channels are aligned to integer bins of the input FFT.
    // E.g. fft size of 1024 and decimation of 16 work, but decimation of 32 does not.
    if((d_fft_size - d_ifft_size) % d_decimation != 0) {
        throw std::runtime_error("FFT size and decimation can not be matched.");
    }

    // Make sure we don't have phase shifts between two iterations.
    if(d_output_step % d_inverse_overlap != 0) {
        throw std::runtime_error("FFT size and decimation can not be matched.");
    }

    set_output_multiple(d_ifft_size - d_ifft_size / d_inverse_overlap);

    // We want to keep the last 1/d_inverse_overlap block in the history buffer
    set_history(d_fft_size / d_inverse_overlap + 1);

    d_out = (fftwf_complex *)volk_malloc(d_fft_size * sizeof(fftwf_complex), volk_get_alignment());
}

/*
 * Our virtual destructor.
 */
fft_channelizer_impl::~fft_channelizer_impl() {}


int fft_shift(int N, int f)
{
    if(f < N / 2) {
        return f + N / 2;
    } else {
        return f - N / 2;
    }
}

int positive_offset(int N, int step, int channel, int n)
{
    return fft_shift(N, channel * step + n / 2);
}

int negative_offset(int N, int step, int channel, int n)
{
    return fft_shift(N, channel * step);
}

int fft_channelizer_impl::work(int noutput_items,
                               gr_vector_const_void_star& input_items,
                               gr_vector_void_star& output_items)
{
    const input_type* in = reinterpret_cast<const input_type*>(input_items[0]);

    int ninput_items = noutput_items * d_decimation;

    // GNURadio supplies us with some history in the front of the buffer. We use
    // this for the overlap of the FFT.
    for(int i = 0; i < ninput_items; i += d_fft_size - d_fft_size / d_inverse_overlap) {
        // Initial forward FFT. No window is applied.
#if 0
        fftwf_complex * fft_in = (fftwf_complex*)&in[i];
        auto it = d_fft_plans.find(fft_in);

        if(it == d_fft_plans.end()) {
            fftwf_plan plan = fftwf_plan_dft_1d(d_fft_size, fft_in, d_out, FFTW_FORWARD, FFTW_MEASURE);
            d_fft_plans.insert({fft_in, plan});
            fftwf_execute(plan);
        } else {
            fftwf_execute(it->second);
        }
#else
        memcpy(d_fft.get_inbuf(), &in[i], d_fft_size * sizeof(gr_complex));
        d_fft.execute();
#endif

        // Cycle through each output and perform the smaller reverse FFT at the appropriate
        // location of the larger forward FFTs output.
        for(int j = 0; j < output_items.size(); j++) {

            // Construct the input to the reverse FFT. We need to copy the two halves
            // as the output of the initial FFT still is [DC, Positive Freq, Negative Freq].
            memcpy(&d_ifft.get_inbuf()[0], &d_fft.get_outbuf()[positive_offset(d_fft_size, d_output_step, j, d_ifft_size)], d_ifft_size / 2 * sizeof(gr_complex));
            memcpy(&d_ifft.get_inbuf()[d_ifft_size / 2], &d_fft.get_outbuf()[negative_offset(d_fft_size, d_output_step, j, d_ifft_size)], d_ifft_size / 2 * sizeof(gr_complex));

            d_ifft.execute();

            output_type* out = reinterpret_cast<output_type*>(output_items[j]);

            volk_32fc_s32fc_multiply_32fc(&out[i / d_decimation], &d_ifft.get_outbuf()[d_ifft_size / d_inverse_overlap], 1./d_fft_size, (d_ifft_size - d_ifft_size / d_inverse_overlap));
        }
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace iridium */
} /* namespace gr */
