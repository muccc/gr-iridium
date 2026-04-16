#pragma once
namespace iridium {

enum class direction {
    UNDEF = 0,
    DOWNLINK = 1,
    UPLINK = 2,
    DOWNLINK_NEXT = 3,
    UPLINK_NEXT = 4,
};

const int SYMBOLS_PER_SECOND = 25000;
const int SYMBOLS_PER_SECOND_NEXT = 30000;
const int UW_LENGTH = 12;

const int SIMPLEX_FREQUENCY_MIN = 1626000000;

const int PREAMBLE_LENGTH_SHORT = 16;
const int PREAMBLE_LENGTH_LONG = 64;

const int MIN_FRAME_LENGTH_NORMAL = 131; // IBC frame
const int MAX_FRAME_LENGTH_NORMAL = 191;

const int MIN_FRAME_LENGTH_SIMPLEX = 80; // Single page IRA
const int MAX_FRAME_LENGTH_SIMPLEX = 444;

const int UW_DL[] = { 0, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 2 };
const int UW_UL[] = { 2, 2, 0, 0, 0, 2, 0, 0, 2, 0, 2, 2 };

// NEXT/Certus-C1 unique words (symbol-level).
// UW_UL_NEXT is experimental and may change; it is the current best
// guess and has not yet been verified against real NEXT uplink
// traffic.
const int UW_DL_NEXT[] = { 2, 2, 0, 2, 2, 0, 2, 0, 2, 0, 2, 2 };
const int UW_UL_NEXT[] = { 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0 };
} // namespace iridium
