#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 Free Software Foundation, Inc.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#


import numpy
from gnuradio import gr
from collections import namedtuple

Message = namedtuple('Message', ['timestamp', 'center_frequency', 'confidence', 'pmt'])

PMT_TIMESTAMP = gr.pmt.intern('timestamp')
PMT_CENTER_FREQUENCY = gr.pmt.intern('center_frequency')
PMT_CONFIDENCE = gr.pmt.intern('confidence')

class frame_sorter(gr.sync_block):
    """
    docstring for block frame_sorter
    """
    def __init__(self):
        gr.sync_block.__init__(self,
            name="frame_sorter",
            in_sig=None,
            out_sig=None)

        self._messages = []
        self.message_port_register_in(gr.pmt.intern('pdus'))
        self.message_port_register_out(gr.pmt.intern('pdus'))
        self.set_msg_handler(gr.pmt.intern('pdus'), self.handle_msg)

    def handle_msg(self, msg_pmt):
        meta = gr.pmt.car(msg_pmt)

        timestamp = gr.pmt.to_uint64(gr.pmt.dict_ref(meta, PMT_TIMESTAMP, gr.pmt.PMT_NIL))
        center_frequency = gr.pmt.to_double(gr.pmt.dict_ref(meta, PMT_CENTER_FREQUENCY, gr.pmt.PMT_NIL))
        confidence = gr.pmt.to_long(gr.pmt.dict_ref(meta, PMT_CONFIDENCE, gr.pmt.PMT_NIL))

        remove_count = 0
        insert_index = -1
        remove_index = None
        for idx, message in enumerate(self._messages):
            ts_delta = timestamp - message.timestamp
            if ts_delta > 2e9:
                self.message_port_pub(gr.pmt.intern('pdus'), message.pmt)
                remove_count += 1
            elif abs(ts_delta) <= 1000:
                if abs(message.center_frequency - center_frequency) < 10000:
                    if message.confidence < confidence:
                        remove_index = idx
                    else:
                        insert_index=None
            elif ts_delta < 0:
                break
            if ts_delta > 0 and insert_index is not None:
                insert_index=idx

        if insert_index is not None:
            self._messages.insert(insert_index+1, Message(timestamp, center_frequency, confidence, msg_pmt))
            if remove_index is not None and remove_index > insert_index:
                remove_index+=1
        if remove_index is not None:
            del self._messages[remove_index]
        if remove_count > 0:
            self._messages = self._messages[remove_count:]

    def stop(self):
        # Flush remaining messages
        for message in self._messages:
            self.message_port_pub(gr.pmt.intern('pdus'), message.pmt)
        self._messages = []
        return True

    def work(self, input_items, output_items):
        return len(input_items[0])

