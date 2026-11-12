#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2022 gr-iridium authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import pmt

try:
    from iridium import frame_sorter
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from iridium import frame_sorter


class qa_frame_sorter(gr_unittest.TestCase):

    FLUSH_TIMEOUT = 2 * 10**9 + 1
    SECOND = 1 * 10**9

    def setUp(self):
        self.tb = gr.top_block()
        self.max_timestamp = 0
        self.sorter = frame_sorter()
        self.d1 = blocks.message_debug()
        self.tb.msg_connect((self.sorter, 'pdus'), (self.d1, 'store'))
        self.tb.start()

    def tearDown(self):
        self.tb.stop()
        self.tb.wait()
        self.tb = None

    def test_instance(self):
        instance = frame_sorter()

    def send_frame(self, timestamp, center_frequency, confidence, tch_type=None, crc_ok=None):
        msg_meta = pmt.dict_add(pmt.make_dict(), pmt.intern('timestamp'), pmt.from_uint64(timestamp))
        msg_meta = pmt.dict_add(msg_meta, pmt.intern('center_frequency'), pmt.from_float(center_frequency))
        msg_meta = pmt.dict_add(msg_meta, pmt.intern('confidence'), pmt.from_long(confidence))
        if tch_type is not None:
            # Like tch_c_demod: TCH_C frames carry their type, legacy frames do not
            msg_meta = pmt.dict_add(msg_meta, pmt.intern('tch_type'), pmt.from_long(tch_type))
        if crc_ok is not None:
            msg_meta = pmt.dict_add(msg_meta, pmt.intern('crc_ok'), pmt.from_bool(crc_ok))
        msg = pmt.cons(msg_meta, pmt.init_u8vector(2, range(2)))
        if timestamp > self.max_timestamp:
            self.max_timestamp = timestamp
        self.sorter.to_basic_block()._post(pmt.intern("pdus"), msg)

    def flush(self):
        self.send_frame(self.max_timestamp + self.FLUSH_TIMEOUT, 0, 0)

    def expect(self, num):
        # Publish all frames remaining
        self.flush()

        self.waitFor(lambda: self.d1.num_messages() == num, timeout=1.0, poll_interval=0.01)

    def expect_frames(self, frames):
        # Publish all frames remaining
        self.flush()

        self.waitFor(lambda: self.d1.num_messages() == len(frames), timeout=1.0, poll_interval=0.01)

        for i in range(self.d1.num_messages()):
            f = pmt.to_python(self.d1.get_message(i))[0]
            self.assertEqual(frames[i], f)

    def test_001_deduplicate_time(self):
        # Less than one ms difference. Second one should win (higher confidence)
        self.send_frame(4, 1000, 99)
        self.send_frame(3 + 1000 * 1000, 1000, 100)

        self.expect_frames([
            {'timestamp': 3 + 1000 * 1000, 'center_frequency': 1000, 'confidence': 100},
        ])

    def test_002_sort_time(self):
        # Exactly one ms difference. Both should be sorted and published
        self.send_frame(4 + 1000 * 1000, 4, 100)
        self.send_frame(4, 4, 99)

        self.expect_frames([
            {'timestamp': 4, 'center_frequency': 4, 'confidence': 99},
            {'timestamp': 4 + 1000 * 1000, 'center_frequency': 4, 'confidence': 100},
        ])

    def test_003_deduplicate_frequency(self):
        # Less than 10 kHz difference. First one should win (higher confidence)
        self.send_frame(1000, 4001, 100)
        self.send_frame(1000, 14000, 99)

        self.expect_frames([
            {'timestamp': 1000, 'center_frequency': 4001, 'confidence': 100},
        ])

    def test_004_sort_frequency(self):
        # Exactly 10 kHz difference. Both should be sorted and published
        self.send_frame(3, 14000, 99)
        self.send_frame(3, 4000, 100)

        self.expect_frames([
            {'timestamp': 3, 'center_frequency': 4000, 'confidence': 100},
            {'timestamp': 3, 'center_frequency': 14000, 'confidence': 99},
        ])

    def test_005_deduplicate_time_2(self):
        # Less than one ms difference with a valid packet in between on a different frequency.
        # Second one should go through (no duplicate)
        # Third one should win (higher confidence)
        self.send_frame(4 + 4, 4, 99)
        self.send_frame(4 + 100, 100000, 100)
        self.send_frame(4 + 3 + 1000 * 1000, 4, 100)

        self.expect_frames([
            {'timestamp': 4 + 100, 'center_frequency': 100000, 'confidence': 100},
            {'timestamp': 4 + 3 + 1000 * 1000, 'center_frequency': 4, 'confidence': 100},
        ])

    def test_006_deduplicate_time_3(self):
        # Less than one ms difference. First one should win (higher confidence)
        self.send_frame(26457929700, 1626074240, 100)
        self.send_frame(26457930500, 1626074368, 97)

        self.expect_frames([
            {'timestamp': 26457929700, 'center_frequency': 1626074240, 'confidence': 100},
        ])

    def test_007_valid_tch_c_frame_drops_legacy_frame(self):
        # A legacy frame on top of a TCH_C frame with a valid CRC is a misreading of the
        # TCH_C burst and is dropped, whatever the confidence and the order
        self.send_frame(1000, 4000, 100)
        self.send_frame(1000 + 100, 4000, 60, tch_type=1, crc_ok=True)
        self.send_frame(5 * 10**6, 4000, 60, tch_type=1, crc_ok=True)
        self.send_frame(5 * 10**6 + 100, 4000, 100)

        self.expect_frames([
            {'timestamp': 1100, 'center_frequency': 4000, 'confidence': 60, 'tch_type': 1, 'crc_ok': True},
            {'timestamp': 5 * 10**6, 'center_frequency': 4000, 'confidence': 60, 'tch_type': 1, 'crc_ok': True},
        ])

    def test_007b_invalid_tch_c_frame_keeps_legacy_frame(self):
        # Without a valid CRC, legacy and TCH_C frames are not duplicates of each other
        self.send_frame(1000, 4000, 100)
        self.send_frame(1000 + 100, 4000, 60, tch_type=1, crc_ok=False)

        self.expect_frames([
            {'timestamp': 1000, 'center_frequency': 4000, 'confidence': 100},
            {'timestamp': 1100, 'center_frequency': 4000, 'confidence': 60, 'tch_type': 1, 'crc_ok': False},
        ])

    def test_008_crc_beats_confidence(self):
        # Two copies of a TCH_C frame: the one with a valid CRC wins, even with a lower
        # confidence, in either order
        self.send_frame(1000, 4000, 90, tch_type=1, crc_ok=False)
        self.send_frame(1000 + 100, 4001, 80, tch_type=1, crc_ok=True)
        self.send_frame(5 * 10**6, 4000, 80, tch_type=1, crc_ok=True)
        self.send_frame(5 * 10**6 + 100, 4001, 90, tch_type=1, crc_ok=False)

        self.expect_frames([
            {'timestamp': 1100, 'center_frequency': 4001, 'confidence': 80, 'tch_type': 1, 'crc_ok': True},
            {'timestamp': 5 * 10**6, 'center_frequency': 4000, 'confidence': 80, 'tch_type': 1, 'crc_ok': True},
        ])

    def test_009_confidence_between_valid_frames(self):
        # Both copies pass the CRC: the higher confidence wins as before
        self.send_frame(1000, 4000, 80, tch_type=1, crc_ok=True)
        self.send_frame(1000 + 100, 4001, 90, tch_type=1, crc_ok=True)

        self.expect_frames([
            {'timestamp': 1100, 'center_frequency': 4001, 'confidence': 90, 'tch_type': 1, 'crc_ok': True},
        ])


if __name__ == '__main__':
    gr_unittest.run(qa_frame_sorter)
