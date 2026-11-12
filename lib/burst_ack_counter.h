/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H
#define INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H

#include <pmt/pmt.h>
#include <unordered_map>
#include <algorithm>
#include <cstdint>

namespace gr {
namespace iridium {

/*
 * The burst sources can publish every burst to several downmixers (one per burst type).
 * Each of them acknowledges the burst on burst_handled. A burst only counts as handled
 * once all of them have acknowledged it, so the source sees the backlog of the slowest
 * one.
 *
 * Downmixers running with wait_for_demod acknowledge with (id . n_frames): the number of
 * frames they published for the burst. The demodulators report every handled frame
 * with its burst id on frame_handled, and the burst counts as handled only once all of
 * its frames are. This way the source also sees the backlog of the demodulators.
 *
 * A burst on the border of two channels can arrive twice with the same id; it then
 * needs twice the acknowledgements and completes twice.
 */
class burst_ack_counter
{
public:
    explicit burst_ack_counter(int acks_per_burst)
        : d_acks_per_burst(std::max(acks_per_burst, 1))
    {
    }

    // An acknowledgement from a downmixer. Returns the number of bursts completed.
    int burst_handled(const pmt::pmt_t& msg)
    {
        if (!pmt::is_pair(msg)) {
            if (d_acks_per_burst == 1) {
                return 1;
            }
            return update(pmt::to_uint64(msg), 1, 0);
        }
        return update(pmt::to_uint64(pmt::car(msg)), 1, pmt::to_long(pmt::cdr(msg)));
    }

    // A demodulator has handled a frame. Returns the number of bursts completed.
    int frame_handled(const pmt::pmt_t& msg)
    {
        return update(pmt::to_uint64(msg), 0, -1);
    }

private:
    struct pending {
        int acks = 0;   // acknowledgements from downmixers
        int frames = 0; // frames published but not yet handled by a demodulator
    };

    const int d_acks_per_burst;
    std::unordered_map<uint64_t, pending> d_pending;

    int update(uint64_t id, int acks, int frames)
    {
        pending& p = d_pending[id];
        p.acks += acks;
        // A frame can be handled before its downmixer acknowledged the burst
        p.frames += frames;
        int completed = 0;
        while (p.acks >= d_acks_per_burst && p.frames == 0) {
            p.acks -= d_acks_per_burst;
            completed++;
        }
        if (p.acks == 0 && p.frames == 0) {
            d_pending.erase(id);
        }
        return completed;
    }
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H */
