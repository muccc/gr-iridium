/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H
#define INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H

#include <pmt/pmt.h>
#include <algorithm>
#include <cstdint>
#include <map>

namespace gr {
namespace iridium {

/*
 * The downmixer acknowledges every burst of a source on burst_handled.
 *
 * Downmixers running with wait_for_demod acknowledge with (id . n_frames): the number of
 * frames they published for the burst. The demodulators report every handled frame
 * with its burst id on frame_handled, and the burst counts as handled only once all of
 * its frames are. This way the source also sees the backlog of the demodulators.
 *
 * A burst on the border of two channels can arrive twice with the same id; it then
 * needs twice the acknowledgements and completes twice.
 *
 * The source also registers every burst it publishes. backlog() is the time between the
 * oldest burst not handled yet and the newest one. The frame sorter only finds duplicates
 * and sorts frames which reach it less than 2 s (of the time of the bursts) apart, and a
 * number of unhandled bursts does not limit that: a source which can wait (reading from
 * a file) should also do so while the backlog is longer than MAX_BACKLOG.
 */
class burst_ack_counter
{
public:
    burst_ack_counter() {}

    static constexpr uint64_t MAX_BACKLOG = 1000000000; // ns

    // The source has published a burst
    void published(uint64_t id, uint64_t timestamp)
    {
        pending& p = d_pending[id];
        p.published++;
        p.timestamp = timestamp;
        d_newest = std::max(d_newest, timestamp);
    }

    // Burst ids are handed out in the order of the bursts
    uint64_t backlog() const
    {
        if (d_pending.empty()) {
            return 0;
        }
        const uint64_t oldest = d_pending.begin()->second.timestamp;
        return d_newest > oldest ? d_newest - oldest : 0;
    }

    // An acknowledgement from a downmixer. Returns the number of bursts completed.
    int burst_handled(const pmt::pmt_t& msg)
    {
        if (!pmt::is_pair(msg)) {
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
        int published = 0; // times the source published it and it is not completed yet
        int acks = 0;      // acknowledgements from the downmixer
        int frames = 0;    // frames published but not yet handled by a demodulator
        uint64_t timestamp = 0;
    };

    std::map<uint64_t, pending> d_pending;
    uint64_t d_newest = 0; // timestamp of the newest burst published

    int update(uint64_t id, int acks, int frames)
    {
        pending& p = d_pending[id];
        p.acks += acks;
        // A frame can be handled before its downmixer acknowledged the burst
        p.frames += frames;
        int completed = 0;
        while (p.acks >= 1 && p.frames == 0) {
            p.acks--;
            p.published--;
            completed++;
        }
        if (p.published <= 0 && p.acks == 0 && p.frames == 0) {
            d_pending.erase(id);
        }
        return completed;
    }
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_BURST_ACK_COUNTER_H */
