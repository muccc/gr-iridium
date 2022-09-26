#!/usr/bin/env python3

# python3 -m pip install git+https://github.com/miili/ringnes.git

import zmq
import ringnes
import numpy
import signal, os
import time

dwim = False
remaining = 0

def handler(signum, frame):
    global dwim
    dwim = True

signal.signal(signal.SIGUSR1, handler)

context = zmq.Context()

# ~2 seconds at 10 MSPS
cap = int(2 * 10e6 * 8 / 4096) * 4096

#print("buffering %d samples" % (cap / 8))
ring = ringnes.Ringbuffer(capacity=cap)

socket = context.socket(zmq.SUB)
socket.connect("tcp://localhost:5000")
socket.setsockopt(zmq.SUBSCRIBE, b'')


while True:
    data = socket.recv()
    ring.put(memoryview(data))

    if dwim:
        if remaining <= 0:
            f = open("%d.fc32" % time.time(), 'wb')
            samples = numpy.frombuffer(ring, dtype=numpy.complex64)
            f.write(ring)
            #print("write", len(samples))
        # 1 second at 10 MSPS
        remaining += 10000000
        dwim = False

    if remaining > 0:
        print(len(data))
        samples = numpy.frombuffer(data, dtype=numpy.complex64)
        f.write(data)
        remaining -= len(samples)
        #print("write", len(samples), "remain", remaining)

        if remaining <= 0:
            f.close()

