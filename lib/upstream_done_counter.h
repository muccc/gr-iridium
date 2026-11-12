/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_UPSTREAM_DONE_COUNTER_H
#define INCLUDED_IRIDIUM_UPSTREAM_DONE_COUNTER_H

#include <gnuradio/block.h>
#include <atomic>

namespace gr {
namespace iridium {

/*
 * GNU Radio considers a block with only message ports finished as soon as the first of
 * its upstream message neighbours is done: it then stops processing messages, and the
 * messages that the other neighbours still send are lost. At the end of a recording
 * this drops the output of slow downmixers whenever a faster one finishes first.
 *
 * With this installed, the block only finishes once it has been notified by all of its
 * set_upstream_count() upstream connections.
 */
class upstream_done_counter
{
public:
    void install(gr::block* block)
    {
        block->set_msg_handler(pmt::mp("system"), [this, block](pmt::pmt_t msg) {
            if (pmt::is_pair(msg) && pmt::eqv(pmt::car(msg), pmt::mp("done")) &&
                ++d_n_done < d_n_upstream) {
                return;
            }
            block->system_handler(msg);
        });
    }

    void set_upstream_count(int n) { d_n_upstream = n; }

private:
    std::atomic<int> d_n_upstream{ 1 };
    int d_n_done = 0;
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_UPSTREAM_DONE_COUNTER_H */
