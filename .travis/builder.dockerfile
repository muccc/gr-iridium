FROM debian:buster-slim

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gnuradio gnuradio-dev \
    libvolk1.4 libvolk1-dev \
    cmake build-essential \
    python python-dev python-scipy swig \
    libboost-dev libcppunit-dev

RUN mkdir /build
WORKDIR /build