/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_FRAME_SORTER_CPP_H
#define INCLUDED_IRIDIUM_FRAME_SORTER_CPP_H

#include <gnuradio/block.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief <+description of block+>
 * \ingroup iridium
 *
 */
class IRIDIUM_API frame_sorter : virtual public gr::block
{
public:
    typedef std::shared_ptr<frame_sorter> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of iridium::frame_sorter.
     *
     * To avoid accidental use of raw pointers, iridium::frame_sorter's
     * constructor is in a private implementation
     * class. iridium::frame_sorter::make is the public interface for
     * creating new instances.
     */
    /*!
     * \param reorder_window_ns a frame is held until a frame this many nanoseconds
     * newer arrives, then emitted. Must exceed the largest latency difference between
     * the demodulators feeding the sorter: the TCH_C slot prober back-probes up to its
     * history, so with TCH_C the window must cover that history plus a margin, else its
     * late frames are emitted out of time order. Default 2 s (the legacy-only latency).
     */
    static sptr make(int64_t reorder_window_ns = 2000000000);

    /*!
     * \brief Number of message connections into this block (default 1). The block
     * only finishes once all of them are done, instead of after the first one.
     */
    virtual void set_upstream_count(int n) = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_FRAME_SORTER_CPP_H */
