#!/usr/bin/env python3

import zmq
import bits_pb2

context = zmq.Context()
socket = context.socket(zmq.SUB)

socket.connect("tcp://localhost:2323")
socket.setsockopt_string(zmq.SUBSCRIBE, "")

frame = bits_pb2.Frame()

while True:
    msg = socket.recv()
    frame.ParseFromString(msg)

    data_bin = ''.join(['{:08b}'.format(x) for x in frame.data])
    print("RAW: %s %012.4f %010d N:%05.2f%+06.2f I:%011d %3d%% %.5f %3d %s" % (frame.fileinfo, frame.timestamp/1000000, frame.frequency, frame.magnitude, frame.noise, frame.id, frame.confidence, frame.level, frame.len, data_bin))
