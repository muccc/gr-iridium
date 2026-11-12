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

#include "tch_c_decoder.h"

#include <algorithm>
#include <cmath>

namespace gr {
namespace iridium {

namespace {

// The second constituent encoder sees the bits of odd couples swapped: 2A+B -> 2B+A
const int SWAP_AB[4] = { 0, 2, 1, 3 };

const float NEG_INF = -1e9f;

// log(1 + exp(-d)) for d = 0 .. 8 in steps of 1/8 (beyond that it is < 3.4e-4)
struct correction_table {
    float t[65];
    correction_table()
    {
        for (int i = 0; i <= 64; i++) {
            t[i] = std::log1p(std::exp(-i / 8.f));
        }
    }
};
const correction_table CORRECTION;

// Jacobian logarithm log(exp(a) + exp(b)). The table is ~5x faster than
// computing the correction and decodes just as well; max-log alone loses ~0.2 dB.
inline float max_star(float a, float b)
{
    float m = std::max(a, b);
    float d = std::fabs(a - b) * 8;
    return d < 64 ? m + CORRECTION.t[(int)(d + 0.5f)] : m;
}

} // namespace

tch_c_decoder::tch_c_decoder(const tch_c_code& code, int iterations)
    : d_iterations(iterations),
      d_n_couples(code.info_bits / 2),
      d_n_coded(code.coded_bits),
      d_n_info(code.info_bits)
{
    for (int i = 0; i < 4; i++) {
        d_puncture_y[i] = code.puncture_y[i];
    }

    // 8 state duo-binary RSC, state = 4*S1 + 2*S2 + S3,
    // feedback 1+D+D^3, parity Y 1+D^2+D^3
    for (int u = 0; u < 4; u++) {
        int A = u >> 1, B = u & 1;
        for (int s = 0; s < 8; s++) {
            int s1 = (s >> 2) & 1, s2 = (s >> 1) & 1, s3 = s & 1;
            int fb = A ^ B ^ s1 ^ s3;
            d_next_state[u][s] = (fb << 2) | ((s1 ^ B) << 1) | (s2 ^ B);
            d_parity[u][s] = fb ^ s2 ^ s3;
        }
    }

    // Interleaved couple j is drawn from natural couple d_ilv_src[j]
    const int N = d_n_couples;
    const int P[4] = { 0, N / 2 + code.p[1], code.p[2], N / 2 + code.p[3] };
    for (int j = 0; j < N; j++) {
        d_ilv_src[j] = (code.p[0] * j + P[j % 4] + 1) % N;
    }

    // Channel interleaver: written column-wise into 2 * floor(sqrt(L) / 2) + 1 columns
    // (21 for C1, 29 for C2), read row-wise.
    // Transmitted bit p is coded bit d_chan_src[p].
    const int ncol = 2 * (int)std::floor(std::sqrt((double)d_n_coded) / 2) + 1;
    int start[64];
    int sum = 0;
    for (int c = 0; c < ncol; c++) {
        start[c] = sum;
        sum += (d_n_coded - c + ncol - 1) / ncol;
    }
    for (int p = 0; p < d_n_coded; p++) {
        d_chan_src[p] = start[p % ncol] + p / ncol;
    }
}

namespace {

// DVB energy dispersal PRBS 1+X^14+X^15, seed 100101010000000, restarted every burst
struct prbs_table {
    uint8_t bits[tch_c_decoder::MAX_SCRAMBLED];
    prbs_table()
    {
        int reg[15] = { 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < tch_c_decoder::MAX_SCRAMBLED; i++) {
            int b = reg[13] ^ reg[14];
            bits[i] = b;
            for (int k = 14; k > 0; k--) {
                reg[k] = reg[k - 1];
            }
            reg[0] = b;
        }
    }
};
const prbs_table PRBS;

} // namespace

void tch_c_decoder::descramble(uint8_t* bits, int n, int offset)
{
    for (int i = 0; i < n; i++) {
        bits[i] ^= PRBS.bits[offset + i];
    }
}

bool tch_c_decoder::check_crc(const uint8_t* msg, int n)
{
    uint16_t crc = crc10(msg, n - 10);
    for (int i = 0; i < 10; i++) {
        if (msg[n - 10 + i] != ((crc >> (9 - i)) & 1)) {
            return false;
        }
    }
    return true;
}

uint16_t tch_c_decoder::crc10(const uint8_t* bits, int n)
{
    // x^10 + x^9 + x^6 + x^3 + x^2 + x + 1, MSB first, init 0, no final XOR
    uint16_t r = 0;
    for (int i = 0; i < n; i++) {
        int fb = ((r >> 9) & 1) ^ bits[i];
        r = (r << 1) & 0x3FF;
        if (fb) {
            r ^= 0x24F;
        }
    }
    return r;
}

/*
 * One log-MAP pass over the tail-biting trellis.
 * in[k][v]: log probability of couple value v at step k (channel systematic + a priori).
 * Writes the extrinsic log probabilities to out.
 */
void tch_c_decoder::siso(const float (*in)[4], const float* parity_llr, float (*out)[4])
{
    for (int k = 0; k < d_n_couples; k++) {
        for (int s = 0; s < 8; s++) {
            for (int v = 0; v < 4; v++) {
                d_gamma[k][s][v] = in[k][v] + (d_parity[v][s] ? -parity_llr[k] : 0);
            }
        }
    }

    std::fill(&d_alpha[0][0], &d_alpha[0][0] + d_n_couples * 8, 0.f);
    std::fill(&d_beta[0][0], &d_beta[0][0] + d_n_couples * 8, 0.f);

    // No known start state: run around the ring a few times to let the metrics settle
    for (int pass = 0; pass < 3; pass++) {
        float a[8];
        std::copy(d_alpha[0], d_alpha[0] + 8, a);
        for (int k = 0; k < d_n_couples; k++) {
            float na[8];
            std::fill(na, na + 8, NEG_INF);
            for (int s = 0; s < 8; s++) {
                for (int v = 0; v < 4; v++) {
                    int ns = d_next_state[v][s];
                    na[ns] = max_star(na[ns], a[s] + d_gamma[k][s][v]);
                }
            }
            float m = *std::max_element(na, na + 8);
            for (int s = 0; s < 8; s++) {
                na[s] -= m;
            }
            std::copy(na, na + 8, d_alpha[(k + 1) % d_n_couples]);
            std::copy(na, na + 8, a);
        }
    }

    for (int pass = 0; pass < 3; pass++) {
        float b[8];
        std::copy(d_beta[0], d_beta[0] + 8, b);
        for (int k = d_n_couples - 1; k >= 0; k--) {
            float nb[8];
            std::fill(nb, nb + 8, NEG_INF);
            for (int s = 0; s < 8; s++) {
                for (int v = 0; v < 4; v++) {
                    nb[s] = max_star(nb[s], b[d_next_state[v][s]] + d_gamma[k][s][v]);
                }
            }
            float m = *std::max_element(nb, nb + 8);
            for (int s = 0; s < 8; s++) {
                nb[s] -= m;
            }
            std::copy(nb, nb + 8, d_beta[k]);
            std::copy(nb, nb + 8, b);
        }
    }

    for (int k = 0; k < d_n_couples; k++) {
        const float* beta_next = d_beta[(k + 1) % d_n_couples];
        float m[4] = { NEG_INF, NEG_INF, NEG_INF, NEG_INF };
        for (int s = 0; s < 8; s++) {
            for (int v = 0; v < 4; v++) {
                m[v] = max_star(m[v],
                                d_alpha[k][s] + d_gamma[k][s][v] +
                                    beta_next[d_next_state[v][s]]);
            }
        }
        for (int v = 0; v < 4; v++) {
            out[k][v] = m[v] - in[k][v];
        }
        float mx = *std::max_element(out[k], out[k] + 4);
        for (int v = 0; v < 4; v++) {
            out[k][v] -= mx;
        }
    }
}

bool tch_c_decoder::decode(const float* llr, uint8_t* msg)
{
    decode_raw(llr, msg);
    descramble(msg, d_n_info, 0);
    return check_crc(msg, d_n_info);
}

void tch_c_decoder::decode_raw(const float* llr, uint8_t* bits)
{
    // Undo the channel interleaver
    float coded[MAX_CODED];
    for (int p = 0; p < d_n_coded; p++) {
        coded[d_chan_src[p]] = llr[p];
    }

    // Depuncture: A, B for every couple, parity Y1, Y2 for every 4th
    float la[MAX_COUPLES], lb[MAX_COUPLES], ly1[MAX_COUPLES], ly2[MAX_COUPLES];
    int p = 0;
    for (int k = 0; k < d_n_couples; k++) {
        la[k] = coded[p];
        lb[k] = coded[p + 1];
        p += 2;
        if (d_puncture_y[k % 4]) {
            ly1[k] = coded[p];
            ly2[k] = coded[p + 1];
            p += 2;
        } else {
            ly1[k] = 0;
            ly2[k] = 0;
        }
    }

    // Couple log probabilities from the two bit LLRs
    for (int k = 0; k < d_n_couples; k++) {
        for (int v = 0; v < 4; v++) {
            d_chan[k][v] = ((v >> 1) ? -la[k] : 0) + ((v & 1) ? -lb[k] : 0);
            d_ext[k][v] = 0;
        }
    }

    for (int it = 0; it < d_iterations; it++) {
        for (int k = 0; k < d_n_couples; k++) {
            for (int v = 0; v < 4; v++) {
                d_in[k][v] = d_chan[k][v] + d_ext[k][v];
            }
        }
        siso(d_in, ly1, d_out1);

        // To encoder 2 order, swapping A/B of odd source couples
        for (int j = 0; j < d_n_couples; j++) {
            int m = d_ilv_src[j];
            for (int v = 0; v < 4; v++) {
                int u = (m % 2) ? SWAP_AB[v] : v;
                d_in[j][v] = d_chan[m][u] + d_out1[m][u];
            }
        }
        siso(d_in, ly2, d_out2);

        // Back to natural order
        for (int j = 0; j < d_n_couples; j++) {
            int m = d_ilv_src[j];
            for (int v = 0; v < 4; v++) {
                int u = (m % 2) ? SWAP_AB[v] : v;
                d_ext[m][v] = d_out2[j][u];
            }
        }
    }

    // Decide each couple on its a posteriori metric
    for (int k = 0; k < d_n_couples; k++) {
        int best = 0;
        float best_metric = NEG_INF;
        for (int v = 0; v < 4; v++) {
            float metric = d_chan[k][v] + d_out1[k][v] + d_ext[k][v];
            if (metric > best_metric) {
                best_metric = metric;
                best = v;
            }
        }
        bits[2 * k] = best >> 1;
        bits[2 * k + 1] = best & 1;
    }
}

} // namespace iridium
} // namespace gr
