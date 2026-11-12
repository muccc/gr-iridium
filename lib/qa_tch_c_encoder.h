/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_QA_TCH_C_ENCODER_H
#define INCLUDED_IRIDIUM_QA_TCH_C_ENCODER_H

#include "tch_c_format.h"
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace gr {
namespace iridium {

/*
 * Transmit side of the TCH_C channel coding, for the unit tests: CRC-10, scrambling,
 * duo-binary circular turbo encoding (DVB-RCS / IEEE 802.16 CTC with the parameters
 * of tch_c_format.h), puncturing and channel interleaving. Written independently of
 * tch_c_decoder.
 */
class qa_tch_c_encoder
{
public:
    explicit qa_tch_c_encoder(const tch_c_format& f) : d_f(f) {}

    /*
     * payload: n_codewords * K - 10 bits. Returns the coded bits of all data symbols in
     * transmission order (n_codewords * L bits, then the padding bits).
     */
    std::vector<uint8_t> encode_burst(const std::vector<uint8_t>& payload) const
    {
        const int K = d_f.code.info_bits;
        const int n_info = d_f.n_codewords * K;
        if ((int)payload.size() != n_info - 10) {
            throw std::invalid_argument("wrong payload length");
        }

        // CRC-10 over the payload, then scrambling across all codewords
        std::vector<uint8_t> msg(payload);
        const uint16_t crc = crc10(payload);
        for (int i = 9; i >= 0; i--) {
            msg.push_back((crc >> i) & 1);
        }
        const std::vector<uint8_t> prbs = dvb_prbs(n_info);
        for (int i = 0; i < n_info; i++) {
            msg[i] ^= prbs[i];
        }

        std::vector<uint8_t> out;
        for (int j = 0; j < d_f.n_codewords; j++) {
            const std::vector<uint8_t> cw = encode_codeword(
                std::vector<uint8_t>(msg.begin() + j * K, msg.begin() + (j + 1) * K));
            out.insert(out.end(), cw.begin(), cw.end());
        }
        out.insert(out.end(), d_f.padding_bits, 0);
        return out;
    }

    // One codeword: K information bits in, L transmitted bits out
    std::vector<uint8_t> encode_codeword(const std::vector<uint8_t>& info) const
    {
        const int N = d_f.code.info_bits / 2;
        const int L = d_f.code.coded_bits;

        // Couples (A, B) as values 2A + B
        std::vector<int> u(N), v(N);
        for (int k = 0; k < N; k++) {
            u[k] = 2 * info[2 * k] + info[2 * k + 1];
        }

        // Interleaver: couple j of the second encoder is couple src(j), with A and B
        // swapped for odd source couples
        const int* p = d_f.code.p;
        const int P[4] = { 0, N / 2 + p[1], p[2], N / 2 + p[3] };
        for (int j = 0; j < N; j++) {
            const int m = (p[0] * j + P[j % 4] + 1) % N;
            v[j] = (m % 2) ? ((u[m] & 1) << 1 | u[m] >> 1) : u[m];
        }

        const std::vector<uint8_t> y1 = rsc_parity(u);
        const std::vector<uint8_t> y2 = rsc_parity(v);

        // Puncturing: A, B for every couple, Y1 and Y2 for the selected ones
        std::vector<uint8_t> coded;
        for (int k = 0; k < N; k++) {
            coded.push_back(u[k] >> 1);
            coded.push_back(u[k] & 1);
            if (d_f.code.puncture_y[k % 4]) {
                coded.push_back(y1[k]);
                coded.push_back(y2[k]);
            }
        }
        if ((int)coded.size() != L) {
            throw std::logic_error("puncturing does not give L bits");
        }

        // Channel interleaver: written column by column into 2 floor(sqrt(L) / 2) + 1
        // columns, read row by row
        const int ncol = 2 * (int)std::floor(std::sqrt((double)L) / 2) + 1;
        std::vector<int> start(ncol);
        int sum = 0;
        for (int c = 0; c < ncol; c++) {
            start[c] = sum;
            sum += (L - c + ncol - 1) / ncol;
        }
        std::vector<uint8_t> out(L);
        for (int t = 0; t < L; t++) {
            out[t] = coded[start[t % ncol] + t / ncol];
        }
        return out;
    }

    static uint16_t crc10(const std::vector<uint8_t>& bits)
    {
        // x^10 + x^9 + x^6 + x^3 + x^2 + x + 1, MSB first, initial value 0
        uint16_t r = 0;
        for (uint8_t b : bits) {
            const int fb = ((r >> 9) & 1) ^ b;
            r = (r << 1) & 0x3ff;
            if (fb) {
                r ^= 0x24f;
            }
        }
        return r;
    }

    // DVB energy dispersal sequence 1 + X^14 + X^15, seed 100101010000000
    static std::vector<uint8_t> dvb_prbs(int n)
    {
        uint16_t reg = 0x4a80; // 100101010000000, first stage in the MSB
        std::vector<uint8_t> out;
        for (int i = 0; i < n; i++) {
            const int b = ((reg >> 1) ^ reg) & 1;
            out.push_back(b);
            reg = (reg >> 1) | (b << 14);
        }
        return out;
    }

private:
    const tch_c_format& d_f;

    /*
     * Parity Y of the 8 state duo-binary RSC (feedback 1 + D + D^3, Y = 1 + D^2 + D^3,
     * B also fed into stages 2 and 3), circular: it starts in the state it ends in.
     */
    static std::vector<uint8_t> rsc_parity(const std::vector<int>& u)
    {
        auto run = [&](int s, std::vector<uint8_t>* y) {
            for (int x : u) {
                const int A = x >> 1, B = x & 1;
                const int s1 = s >> 2 & 1, s2 = s >> 1 & 1, s3 = s & 1;
                const int fb = A ^ B ^ s1 ^ s3;
                if (y) {
                    y->push_back(fb ^ s2 ^ s3);
                }
                s = fb << 2 | (s1 ^ B) << 1 | (s2 ^ B);
            }
            return s;
        };
        for (int s0 = 0; s0 < 8; s0++) {
            if (run(s0, nullptr) == s0) {
                std::vector<uint8_t> y;
                run(s0, &y);
                return y;
            }
        }
        throw std::logic_error("no circulation state");
    }
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_QA_TCH_C_ENCODER_H */
