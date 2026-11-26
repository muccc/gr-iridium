/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_DEMOD_IMPL_H
#define INCLUDED_IRIDIUM_TCH_C_DEMOD_IMPL_H

#include "tch_c_format.h"
#include "upstream_done_counter.h"
#include <iridium/tch_c_demod.h>
#include <map>

namespace gr {
namespace iridium {

class tch_c_demod_impl : public tch_c_demod
{
private:
    float d_alpha; // PLL loop gain

    uint64_t d_n_handled_bursts;
    uint64_t d_n_access_ok_bursts;
    uint64_t d_n_access_ok_sub_bursts;
    std::vector<uint64_t> d_channel_id;

    tch_c_decoder d_c1_decoder;
    tch_c_decoder d_c2_decoder;
    tch_c_decoder d_c8q_decoder;
    tch_c_decoder d_c8a_decoder;

    gr_complex d_decimated_burst[TCH_C_MAX_SYMBOLS];
    gr_complex d_burst_after_pll[TCH_C_MAX_SYMBOLS];
    int d_demodulated_burst[TCH_C_MAX_SYMBOLS];
    float d_llr[TCH_C_MAX_DATA_BITS];
    uint8_t d_msg[tch_c_decoder::MAX_SCRAMBLED];
    uint8_t d_best_msg[tch_c_decoder::MAX_SCRAMBLED];
    int d_best_symbols[TCH_C_MAX_SYMBOLS];
    gr_complex d_best_corrected[TCH_C_MAX_SYMBOLS];
    std::vector<uint8_t> d_bits;

    std::vector<pmt::pmt_t> d_input_ports;
    std::vector<pmt::pmt_t> d_frame_handled_ports;
    upstream_done_counter d_upstream_done;

    // Sharding: with shard_count parallel instances each decodes only the probes whose
    // arrival sequence number is this instance's, splitting the decode across threads.
    int d_shard_index;
    int d_shard_count;
    uint64_t d_probe_seq;

    // Frames taken from the input ports, by timestamp: (channel, frame)
    std::multimap<uint64_t, std::pair<int, pmt::pmt_t>> d_pending;
    void enqueue(int channel, const pmt::pmt_t& msg);
    void handler(int channel, pmt::pmt_t msg);
    void demodulate(int channel, pmt::pmt_t msg);
    void frame_handled(int channel, const pmt::pmt_t& meta);
    void handle_burst(int channel,
                      const gr_complex* burst,
                      size_t burst_size,
                      float sample_rate,
                      double center_frequency,
                      uint64_t timestamp,
                      uint64_t sub_id,
                      float noise,
                      float magnitude,
                      int tch_type,
                      int max_variants);
    int decimate(const gr_complex* in, int size, float sps, int n, gr_complex* out);
    float
    known_symbol_correction(const tch_c_format& f, const gr_complex* x, gr_complex* y);
    float qpsk_first_order_pll(const gr_complex* x, int size, float alpha, gr_complex* y);
    void demod_qpsk(const gr_complex* burst,
                    size_t n_symbols,
                    int* out,
                    float* level,
                    int* confidence);
    bool uw_ok(const int* demodulated_burst, const int* uw);
    float soft_bits(const tch_c_format& f,
                    const gr_complex* y,
                    const int* uw,
                    float* llr,
                    float* noise);
    tch_c_decoder& decoder_for(const tch_c_format& f);
    void decode_gray_qpsk(int* demodulated_burst, size_t n_symbols);
    void map_symbols_to_bits(const int* demodulated_burst,
                             size_t n_symbols,
                             std::vector<uint8_t>& bits);

public:
    tch_c_demod_impl(int shard_index = 0, int shard_count = 1);
    ~tch_c_demod_impl();

    uint64_t get_n_handled_bursts();
    uint64_t get_n_access_ok_bursts();
    uint64_t get_n_access_ok_sub_bursts();
    void set_upstream_count(int n) { d_upstream_done.set_upstream_count(n); }
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_DEMOD_IMPL_H */
