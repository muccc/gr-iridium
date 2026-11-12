/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
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
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_DECODER_H
#define INCLUDED_IRIDIUM_TCH_C_DECODER_H

#include <cstdint>

namespace gr {
namespace iridium {

// Channel coding parameters of a TCH_C burst format
struct tch_c_code {
    int info_bits;     // K: information bits incl. the CRC-10
    int coded_bits;    // L
    int p[4];          // CTC interleaver P0..P3
    int puncture_y[4]; // parity Y transmitted for couples with index % 4 == i
};

/*
 * Decoder for the TCH_C channel coding:
 * L coded bits -> channel de-interleaver -> depuncture ->
 * duo-binary circular turbo decoder (DVB-RCS / 802.16 CTC, K/2 couples) ->
 * DVB PRBS descrambler -> K - 10 payload bits + CRC-10.
 */
class tch_c_decoder
{
public:
    static const int MAX_COUPLES = 808;
    static const int MAX_CODED = 2424;
    static const int MAX_INFO = 2 * MAX_COUPLES;
    // Bursts with several codewords: the scrambler runs on across them
    static const int MAX_SCRAMBLED = 3 * MAX_INFO;

    tch_c_decoder(const tch_c_code& code, int iterations = 8);

    int info_bits() const { return d_n_info; }
    int coded_bits() const { return d_n_coded; }
    int payload_bits() const { return d_n_info - 10; }

    /*
     * llr: coded_bits() LLRs in transmission order, L = log P(bit=0) / P(bit=1).
     * msg: info_bits() descrambled bits out (payload followed by the CRC-10).
     * Returns true if the CRC-10 matches.
     */
    bool decode(const float* llr, uint8_t* msg);

    /*
     * Turbo decoding of one codeword only: info_bits() bits out, still scrambled.
     * For bursts with several codewords (TCH_C8), see descramble() and check_crc().
     */
    void decode_raw(const float* llr, uint8_t* bits);

    // XOR n bits with the DVB PRBS, starting offset bits into the sequence
    static void descramble(uint8_t* bits, int n, int offset);

    // CRC-10 over the first n - 10 bits compared with the last 10
    static bool check_crc(const uint8_t* msg, int n);

    static uint16_t crc10(const uint8_t* bits, int n);

private:
    int d_iterations;
    int d_n_couples;
    int d_n_coded;
    int d_n_info;
    int d_puncture_y[4];

    int d_next_state[4][8];
    int d_parity[4][8];
    int d_ilv_src[MAX_COUPLES];
    int d_chan_src[MAX_CODED];

    // Working memory: couple metrics for the 4 values of A and B
    float d_chan[MAX_COUPLES][4];
    float d_ext[MAX_COUPLES][4];
    float d_in[MAX_COUPLES][4];
    float d_out1[MAX_COUPLES][4];
    float d_out2[MAX_COUPLES][4];
    float d_gamma[MAX_COUPLES][8][4];
    float d_alpha[MAX_COUPLES][8];
    float d_beta[MAX_COUPLES][8];

    void siso(const float (*in)[4], const float* parity_llr, float (*out)[4]);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_DECODER_H */
