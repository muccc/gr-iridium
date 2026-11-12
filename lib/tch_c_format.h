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

#ifndef INCLUDED_IRIDIUM_TCH_C_FORMAT_H
#define INCLUDED_IRIDIUM_TCH_C_FORMAT_H

#include "tch_c_decoder.h"

namespace gr {
namespace iridium {

/*
 * TCH_C burst formats: TCH_C1 QPSK rate 4/5, TCH_C2 QPSK rate 2/3, TCH_C8 QPSK or 16APSK
 * rate 2/3. All of them last 8.27 ms: TCH_Cn runs at n * 30 ksym/s. Every burst starts
 * with the 32 symbol UW, followed by data segments interleaved with 16 symbol pilots
 * (0x5555). UW and pilots are QPSK.
 *
 * The data carry n_codewords codewords back to back, followed by padding_bits zero bits.
 * The scrambler runs on across the codewords and one CRC-10 at the end protects all of
 * them.
 */
struct tch_c_format {
    int type; // n of TCH_Cn
    int symbols_per_second;
    int n_symbols;
    int bits_per_symbol; // 2: QPSK, 4: 16APSK
    int n_codewords;
    int padding_bits;
    int n_pilots;
    int pilot_start[8];
    int n_segments;
    int segment[16][2]; // data symbols: start, length
    tch_c_code code;
};

/*
 * Unique words: 32 symbols, one bit per symbol on the diagonal (0 -> (1+j)/sqrt(2),
 * 1 -> -(1+j)/sqrt(2)), written as quadrants 0 / 2 like the legacy UWs.
 * DL 0x1B565EE9, UL 0x48147313
 */
const int TCH_C_UW_LENGTH = 32;
const int TCH_C_UW_DL[] = { 0, 0, 0, 2, 2, 0, 2, 2, 0, 2, 0, 2, 0, 2, 2, 0,
                            0, 2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 0, 2, 0, 0, 2 };
const int TCH_C_UW_UL[] = { 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0,
                            0, 2, 2, 2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 0, 2, 2 };
// Pilot: 16 symbols (0x5555) at the pilot_start positions of the format
const int TCH_C_PILOT[] = { 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2 };

const int TCH_C_MAX_SYMBOLS = 1984;
const int TCH_C_MAX_PILOTS = 8;
const int TCH_C_PILOT_LENGTH = 16;
const int TCH_C_MAX_DATA_BITS = 1824 * 4;

const tch_c_format TCH_C1 = { 1,
                              30000,
                              248,
                              2,
                              1,
                              0,
                              1,
                              { 102 },
                              2,
                              { { 32, 70 }, { 118, 130 } },
                              { 320, 400, { 13, 106, 108, 2 }, { 1, 0, 0, 0 } } };

const tch_c_format TCH_C2 = { 2,
                              60000,
                              496,
                              2,
                              1,
                              0,
                              2,
                              { 102, 358 },
                              4,
                              { { 32, 70 }, { 118, 130 }, { 248, 110 }, { 374, 122 } },
                              { 576, 864, { 7, 34, 32, 2 }, { 1, 0, 1, 0 } } };

// TCH_C8: QPSK (TCH_C8Q) and 16APSK (TCH_C8A) share timing, pilots and data segments
#define TCH_C8_LAYOUT                                                                  \
    8, { 102, 358, 606, 854, 1102, 1350, 1598, 1846 }, 16,                             \
    {                                                                                  \
        { 32, 70 }, { 118, 130 }, { 248, 110 }, { 374, 122 }, { 496, 110 },            \
            { 622, 122 }, { 744, 110 }, { 870, 122 }, { 992, 110 }, { 1118, 122 },     \
            { 1240, 110 }, { 1366, 122 }, { 1488, 110 }, { 1614, 122 }, { 1736, 110 }, \
            { 1862, 122 }                                                              \
    }

const tch_c_format TCH_C8Q = {
    8, 240000, 1984,          2,
    2, 0,      TCH_C8_LAYOUT, { 1216, 1824, { 23, 112, 4, 116 }, { 1, 0, 1, 0 } }
};

const tch_c_format TCH_C8A = {
    8, 240000, 1984,          4,
    3, 24,     TCH_C8_LAYOUT, { 1616, 2424, { 19, 2, 16, 6 }, { 1, 0, 1, 0 } }
};

#undef TCH_C8_LAYOUT

/*
 * 16APSK constellation: index = 8 c[4k] + 4 c[4k+1] + 2 c[4k+2] +
 * c[4k+3]. The two low bits select the quadrant like QPSK, the two high bits the point:
 * outer ring at 45, 15 or 75 degrees, or the inner ring. Scaled to the QPSK pilots
 * (16384): outer radius 18577 / 16384, inner radius 6193 / 16384.
 */
const float TCH_C_APSK16_I[16] = { 13136,  13136,  -13136, -13136, 17945, 17945,
                                   -17945, -17945, 4808,   4808,   -4808, -4808,
                                   4379,   4379,   -4379,  -4379 };
const float TCH_C_APSK16_Q[16] = {
    13136, -13136, 13136, -13136, 4808, -4808, 4808, -4808,
    17945, -17945, 17945, -17945, 4379, -4379, 4379, -4379
};
const float TCH_C_APSK16_SCALE = 1.f / 16384;

// The layout (rate, pilots, segments) of a TCH_C type; for TCH_C8 that of TCH_C8Q
inline const tch_c_format& tch_c_format_by_type(int type)
{
    return type == 8 ? TCH_C8Q : type == 2 ? TCH_C2 : TCH_C1;
}

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_FORMAT_H */
