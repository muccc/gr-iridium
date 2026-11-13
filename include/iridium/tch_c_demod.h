/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_TCH_C_DEMOD_H
#define INCLUDED_IRIDIUM_TCH_C_DEMOD_H

#include <gnuradio/block.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief Demodulates and decodes TCH_C1, TCH_C2 and TCH_C8 bursts.
 * \ingroup iridium
 *
 * Takes the frames of the burst_downmix blocks running with tch_c_type 1, 2 or 8.
 * Phase and frequency come from a fit to the known symbols (UW and pilots). The data
 * symbols are turned into soft bits (QPSK, or 16APSK for TCH_C8), turbo decoded,
 * descrambled and checked against the CRC-10. If the UW matches but the CRC fails,
 * other carrier recovery variants (and for TCH_C8 the other modulation) are tried.
 *
 * Output PDUs carry the hard decisions of all symbols (2 bits per symbol, 4 for
 * 16APSK) and in the meta data: tch_type, bits_per_symbol, uplink, crc_ok, variant and
 * the decoded payload, besides the fields of iridium_qpsk_demod.
 */
class IRIDIUM_API tch_c_demod : virtual public gr::block
{
public:
    typedef std::shared_ptr<tch_c_demod> sptr;

    /*!
     * \param n_channels Number of input ports (cpdus0 .. cpdusN-1), one per channel of
     * the channelizer. With 1 the input port is called cpdus.
     * \param probe_port Add the input port probes (bursts already at one sample per
     * symbol, e.g. from tch_c_slot_prober) and the output port probe_handled, where
     * the burst_id of every handled probe is reported. Probes are handled together with
     * the frames of the other ports, oldest first. A probe is decoded with the carrier
     * fit to the known symbols only.
     *
     * Meta data of a probe: timestamp (uint64, ns, of the first UW symbol), uw_start
     * (float, 0), center_frequency (double, Hz), sample_rate (float, the symbol rate of
     * the type), tch_type (long: 1, 2 or 8), id (uint64), burst_id (uint64), noise and
     * magnitude (float, as from fft_burst_tagger), probe (PMT_T). The samples (c32
     * vector, n_symbols + 1 of them) are the matched filter output at the symbol
     * centers, the first one at the first UW symbol. Carrier phase and gain do not
     * matter, the carrier frequency may be off by up to +-3 kHz.
     */
    static sptr make(int n_channels, bool probe_port = false);

    virtual uint64_t get_n_handled_bursts() = 0;
    virtual uint64_t get_n_access_ok_bursts() = 0;
    virtual uint64_t get_n_access_ok_sub_bursts() = 0;

    /*!
     * \brief Number of message connections into this block (default 1). The block
     * only finishes once all of them are done, instead of after the first one.
     */
    virtual void set_upstream_count(int n) = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TCH_C_DEMOD_H */
