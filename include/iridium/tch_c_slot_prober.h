/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_H
#define INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_H

#include <gnuradio/sync_block.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief Looks for TCH_C bursts where the TDMA frame timing says they can be
 * \ingroup iridium
 *
 * Takes the wideband output stream of fft_burst_tagger (and its time_base tag), the
 * legacy frames of iridium_qpsk_demod on the message port frames (the ring alerts give
 * the frame timing and Doppler of the satellites) and the TCH_C frames of tch_c_demod on
 * hits. Publishes candidate bursts at one sample per symbol on cpdus, for the probes port
 * of a tch_c_demod, which acknowledges them on its probe_handled port (connect it to
 * frame_handled).
 *
 * Reading from a file (drop_overflow false) work() waits while outstanding_limit probes
 * (0: no limit) are not acknowledged or the oldest of them is more than 1 s older than
 * the newest; otherwise new probes are dropped then.
 *
 * Environment: IRIDIUM_PROBE_DEBUG=1 logs every accepted ring alert to stderr.
 */
class IRIDIUM_API tch_c_slot_prober : virtual public gr::sync_block
{
public:
    typedef std::shared_ptr<tch_c_slot_prober> sptr;

    /*!
     * \param center_frequency Center frequency of the input stream in Hz
     * \param sample_rate Sample rate of the input stream
     * \param types TCH_C types to look for (1, 2, 8)
     * \param history_s Seconds of samples kept to probe the past of a new satellite
     * \param threshold Minimum normalized UW correlation of a candidate
     * \param outstanding_limit Maximum number of probes not acknowledged yet
     * \param drop_overflow Drop probes instead of waiting (live input)
     * \param uplink Also probe uplink slots
     */
    static sptr make(double center_frequency,
                     double sample_rate,
                     const std::vector<int>& types,
                     float history_s = 0.75,
                     float threshold = 0.35,
                     int outstanding_limit = 500,
                     bool drop_overflow = false,
                     bool uplink = false);

    virtual uint64_t get_n_probes() = 0;
    virtual uint64_t get_n_dropped_probes() = 0;
    virtual uint64_t get_n_ring_alerts() = 0;
    virtual int get_n_outstanding() = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_SLOT_PROBER_H */
