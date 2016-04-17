# GNURadio Iridium Out Of Tree Module

This module provides blocks to build an Iridium burst detector and demodulator.

It also provides a sample application which can be used to detect and demodulate data from the Iridium satellite network.

## Prerequisites
A working [GNURadio](https://gnuradio.com) installation with the following components is necessary:

 - VOLK
 - FFT
 - Filter
 - Python

No other OOT module is needed.

## Installation
```
mkdir build
cd build
cmake ..
make
sudo make install
sudo ldconfig
```
## Sample Usage
### Online (with an SDR)

`iridium-extractor examples/hackrf.conf | grep "A:OK" > output.bits`

This will capture the complete Iridium band using a connected HackRF and demodulate detected bursts into frames.

The final `grep "A:OK"` filters the output for frames which have a valid unique word.

### Offline (file based)
`iridium-extractor -c 1626000000 -r 2000000 -f float --offline name-f1.626000e+09-s2.000000e+06-t20160401000000.cfile | grep "A:OK" > output.bits`

This processes the file in offline mode and supplies the needed options via the command line.

## Extracting Iridium Frames From Raw Data

To capture and demodulate Iridium frames use `iridium-extractor`. You can either process a file offline or stream data into the tool.

The extractor can read a configuration file. It also accepts arguments from the command line.
The `examples/` directory contains example configuration files for common use cases.

Some options are only available via the configuration file, others are only available via the command line. If no configuration file is used, a file name can be provided to read samples from a file. If no (configuration) file is specified, samples are read from stdin.

### Configuration File
Configuration files need to have a `.conf` file extension.
The configuration file is grouped into sections. Each section starts with a `[section-name]` line.
### `osmosdr-source` Section

If this section is present an `osmosdr-source` is instantiated

| Option Name      | Description                                |
|------------------|--------------------------------------------|
| `sample_rate`    | Sample rate at which the source should run |
| `center_freq`    | Center frequency for the source in Hz      |
| `gain`           | (RF)-Gain in dB                            |
| `if_gain`        | IF-Gain in dB                              |
| `bb_gain`        | BB-Gain in dB                              |
| `bandwidth`      | Base band filter bandwidth in Hz           |


### Command Line Options
Command line options can be used instead of a configuration file. If a configuration file is also specified, command line options take precedence.

#### `-o`, `--offline`: Offline Processing
By default, the extractor will drop samples or burst if the computing power available is not enough to keep up.

If you have an already recorded file, use the `-o`,`--offline` option to not drop samples. In this case the extractor will pause reading the file (or input stream) until it can process more samples again.

#### `-c`: Center Frequency
The center frequency for the source or the file in Hz.

#### `-r`: Sample Rate
The sample rate of the source or the file.

#### `-f`: Input File Format
| File Format                                        | `iridium-extractor` format option |
|----------------------------------------------------|-----------------------------------|
| complex uint8 (RTLSDR)                             | `rtl`                             |
| complex int8 (hackrf, rad1o with hackrf-transfer)  | `hackrf`                          |
| complex int16 (USRP with specrec from gr-analysis) | `sc16`                            |
| complex float (GNURadio, `uhd_rx_cfile`)           | `float`                           |

#### `-q`: Queue Length
For each channel (by default there is one), a queue is filled with samples where the detector has detected activity. By default each queue is 500 frames long. You can tweak the length of the queue(s) with this option.

