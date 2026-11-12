/* -*- c++ -*- */
/*
 * Copyright 2020 Free Software Foundation, Inc.
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_IRIDIUM_TAGGED_BURST_TO_PDU_H
#define INCLUDED_IRIDIUM_TAGGED_BURST_TO_PDU_H

#include <gnuradio/sync_block.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief <+description of block+>
 * \ingroup iridium
 *
 */
class IRIDIUM_API tagged_burst_to_pdu : virtual public gr::sync_block
{
public:
    typedef std::shared_ptr<tagged_burst_to_pdu> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of iridium::tagged_burst_to_pdu.
     *
     * To avoid accidental use of raw pointers, iridium::tagged_burst_to_pdu's
     * constructor is in a private implementation
     * class. iridium::tagged_burst_to_pdu::make is the public interface for
     * creating new instances.
     *
     * \param acks_per_burst Number of blocks acknowledging every burst on the
     *        burst_handled port, one per connected burst_downmix block. A burst is
     *        outstanding until all of them have acknowledged it.
     *        Downmixers running with wait_for_demod also make it wait until the
     *        demodulators have reported all frames of the burst on frame_handled.
     */
    static sptr make(int max_burst_size,
                     float relative_center_frequency,
                     float relative_span,
                     float d_relative_sample_rate,
                     double sample_offset,
                     int outstanding_limit,
                     bool drop_overflow,
                     int acks_per_burst = 1);

    virtual uint64_t get_n_dropped_bursts() = 0;
    virtual int get_output_queue_size() = 0;
    virtual int get_output_max_queue_size() = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_TAGGED_BURST_TO_PDU_H */
