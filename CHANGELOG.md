# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Indicate support for GNU Radio 3.10 in MANIFEST.md.
- Limit size of detected bursts to 90 ms for better frame sorting.
- Improve suppression of long lasting interference.
- Install example configurations.
- Better center frequency estimates
- `frame_sorter` prefers frames with a valid CRC.

### Added
- Support for ZeroMQ sources.
- Support for USRP devices via UHD.
- Backpressure from the demodulator: a burst only counts as handled once its frames are
  demodulated (`wait_for_demod` of `burst_downmix`, `frame_handled` ports).
- Wideband detection in `fft_burst_tagger` for bursts without a strong spectral line
  (`--wideband-db`, off by default).
- Reception of TCH_C1, TCH_C2 and TCH_C8 bursts with `--tch-c` (off by default): a
  `tch_c_type` parameter for `burst_downmix` and a new `tch_c_demod` block (known-symbol
  carrier recovery, QPSK/16APSK soft demapping, turbo decoding, CRC-10). Output as
  `NXD:` (downlink) / `NXU:` (uplink) lines carrying the decoded payload, see `docs/tch_c.md`.
- `tagged_burst_to_pdu` and `fft_channelizer` take the number of acknowledgements per burst,
  for sources feeding several downmixers.
- Tracking of TCH_C channels across TDMA frames in `fft_burst_tagger` (`--track-frames`, on
  with `--tch-c`).

### Fixed
- Use single "system" gain in BladeRF example.
- Fixed a bug which only allowed multiples of 128 Hz as detected frequency.
- Add a missing boost include which prevented building on some systems.
- Improved duplicate suppression.
- Frames lost at the end of a recording when several downmixers feed one demodulator.

## [1.0.0] - 2022-04-30
### Added
- A changelog

[Unreleased]: https://github.com/muccc/gr-iridium/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/muccc/gr-iridium/releases/tag/v1.0.0
