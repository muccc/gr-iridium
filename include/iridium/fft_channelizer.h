/* -*- c++ -*- */
/*
 * Copyright 2022 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_FFT_CHANNELIZER_H
#define INCLUDED_IRIDIUM_FFT_CHANNELIZER_H

#include <gnuradio/sync_decimator.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief <+description of block+>
 * \ingroup iridium
 *
 */
class IRIDIUM_API fft_channelizer : virtual public gr::sync_decimator
{
public:
    typedef std::shared_ptr<fft_channelizer> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of iridium::fft_channelizer.
     *
     * To avoid accidental use of raw pointers, iridium::fft_channelizer's
     * constructor is in a private implementation
     * class. iridium::fft_channelizer::make is the public interface for
     * creating new instances.
     *
     * \param acks_per_burst Number of blocks acknowledging every burst on the
     *        burst_handled port, one per burst_downmix block connected to a PDU port.
     *        A burst is outstanding until all of them have acknowledged it.
     *        Downmixers running with wait_for_demod also make it wait until the
     *        demodulators have reported all frames of the burst on frame_handled.
     */
    static sptr make(int fft_size,
                     int decimation,
                     bool activate_streams,
                     int pdu_ports,
                     int max_burst_size,
                     int outstanding_limit,
                     bool drop_overflow,
                     int acks_per_burst = 1);

    virtual uint64_t get_n_dropped_bursts() = 0;
    virtual int get_output_queue_size() = 0;
    virtual int get_output_max_queue_size() = 0;
};
} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FFT_CHANNELIZER_H */
