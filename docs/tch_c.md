# TCH_C reception

Besides the classic Iridium bursts (25 ksym/s DE-QPSK, printed as `RAW:` lines),
`iridium-extractor` can demodulate and decode TCH_C bursts. They are traffic channel bursts
with coherent modulation and turbo coding: TCH_C1, TCH_C2 and TCH_C8. This is off by
default. Enable it with `--tch-c`.

There are two ways to find the bursts (`--tch-c-detect`): probing the TDMA slots of the
satellites, timed by their ring alerts (`probe`, the default, see Slot probing), and the
older energy path, burst detection followed by TCH_C downmixers (`energy`, see Receiver).
`both` runs the two side by side.

## Burst formats

All TCH_C bursts last 8.27 ms, and a channel carries one burst per 90 ms TDMA frame. TCH_Cn
runs at n × 30 ksym/s. Every burst starts with a 32 symbol unique word (UW). The data
segments that follow are interleaved with 16 symbol pilots (`0x5555`). UW and pilots use
one bit per symbol on the diagonal (0 → (1+j)/√2, 1 → −(1+j)/√2).

| | TCH_C1 | TCH_C2 | TCH_C8 QPSK | TCH_C8 16APSK |
| --- | --- | --- | --- | --- |
| Symbol rate | 30 ksym/s | 60 ksym/s | 240 ksym/s | 240 ksym/s |
| Symbols | 248 | 496 | 1,984 | 1,984 |
| Approx. bandwidth | ±18 kHz | ±36 kHz | ±144 kHz | ±144 kHz |
| Pilots | 1 | 2 | 8 | 8 |
| Data symbols | 200 | 432 | 1,824 | 1,824 |
| Codewords × info / coded bits | 1 × 320 / 400 | 1 × 576 / 864 | 2 × 1,216 / 1,824 | 3 × 1,616 / 2,424 (+ 24 padding bits) |
| Code rate | 4/5 | 2/3 | 2/3 | 2/3 |
| Payload bits | 310 | 566 | 2,422 | 4,838 |

UW: downlink `0x1B565EE9`, uplink `0x48147313`. The exact layouts (pilot positions, data
segments, turbo interleaver parameters, puncturing) are in `lib/tch_c_format.h`.

Coding chain (transmit order): payload + CRC-10 (x¹⁰+x⁹+x⁶+x³+x²+x+1, MSB first, initial value 0)
→ scrambling with the DVB energy dispersal sequence (1+X¹⁴+X¹⁵, restarted at every burst
and running on across the codewords) → duo-binary circular turbo code per codeword
(DVB-RCS / IEEE 802.16 CTC) → puncturing → block interleaver → QPSK (c[2k] on I, c[2k+1] on Q,
0 on the positive axis) or 16APSK (4+12 points, ring ratio 3). One CRC-10 at the end
protects all codewords of a burst.

## Receiver

The energy path (`--tch-c-detect energy` or `both`) and the parts shared with the prober:

- **Detection** (`fft_burst_tagger`). TCH_C bursts have no spectral line, so the single-bin
  detector misses most of them. The wideband detection (`--wideband-db`, `--wideband-width`)
  averages the relative magnitude over 30 kHz of bins and compares it to a much lower
  threshold (6 dB). After a wideband burst 6–12 ms long, the tagger tracks its channel
  (`--track-frames`): it opens a window on that frequency in each of the next TDMA frames,
  whether anything is detected there or not.
- **Routing.** With TCH_C2 or TCH_C8 enabled, the tagger classifies the width of every
  burst from its spectrum averaged over 9 ms: occupied 22 to 34 kHz off center (class 1), or
  also 50 to 130 kHz off center (class 2). The TCH_C2 downmixer only processes bursts of
  class 1 or 2, the TCH_C8 downmixer only those of class 2, and windows on tracked channels
  take the class of the burst that started the tracking. This halves the CPU time without
  losing a TCH_C2 or TCH_C8 burst on the test recordings.
- **Downmixing** (`burst_downmix` with `tch_c_type` 1, 2 or 8). The legacy downmixer and
  one downmixer per TCH_C type get every burst of a channel. A TCH_C downmixer runs at
  10 samples per symbol (TCH_C8: as many as the channel allows), and the matched filter is RRC
  α=0.2. It searches for the 32 symbol UW and cuts out exactly one burst. If the channel rate is
  not a multiple of the symbol rate, it runs at the actual rate with a fractional number of
  samples per symbol. In tracked windows it searches the whole window for the UW at several
  frequency offsets. This search and the search for the start of a burst run on every
  2nd to 4th sample (about 2.6 samples per symbol), taken straight from the input filter.

  Most bursts a TCH_C downmixer gets are not TCH_C bursts. So it first filters only the
  samples it searches for the UW and normalizes the UW correlation to the energy of the
  samples it covers (1 for a noise free UW). Below 0.25 it skips ahead by one burst length
  and looks again (bursts of consecutive time slots can arrive as one detection); only
  otherwise it filters the rest of the burst and passes it on. Other bursts and noise
  mostly reach about 0.2, bursts passing the CRC had at least 0.27, nearly all of them more
  than 0.55.
- **Demodulation and decoding** (`tch_c_demod`). Phase and frequency come from a fit to the
  known symbols (UW and pilots). Gain and noise variance give soft bits (max-log for
  16APSK). These are then de-interleaved, turbo decoded (log-MAP, up to 8 iterations),
  descrambled and checked against the CRC-10. If the UW matches but the CRC fails, a
  decision-directed PLL and the fit followed by the PLL are tried as well.
  For TCH_C8, the modulation that fits the data symbols better is tried first. 16APSK is
  only decoded if the noise variance on the known symbols is below 0.25 of the symbol
  energy (it needs about 0.11): most TCH_C8 bursts failing the CRC are above that.
- **Sorting** (`frame_sorter`). A TCH_C frame with a valid CRC removes legacy frames at the
  same time and frequency: those are the legacy demodulator misreading the TCH_C burst.
  Among copies of a TCH_C frame, one with a valid CRC beats one without.
- **Backpressure.** The demodulators report every handled frame back to the burst source
  (`frame_handled`), and the downmixers tell it how many frames they published per burst.
  A burst only counts as handled, and frees its place in the queue (`-q`), once all its
  frames are demodulated. Offline the burst source also waits while the oldest unhandled
  burst is more than 1 s older than the newest one, and the TCH_C demodulator handles the
  frames waiting on its inputs oldest first. This keeps the channels and the legacy and
  TCH_C paths in step, so the sorter's 2 s window still catches all duplicates; live, the
  queue limit drops bursts when the demodulators fall behind.

## Slot probing

The energy detection misses weak bursts, and a call shows up as a run of frames with gaps.
`tch_c_slot_prober` looks for TCH_C bursts where the TDMA timing says they can be instead
(`--tch-c-detect probe`, the default, or `both` with the energy path). Its bursts go
straight to `tch_c_demod`; the tagger runs as without `--tch-c`:

- **Timing.** The ring alerts (legacy downlink frames on the simplex ring channel) give the
  90 ms frame phase and the Doppler of every satellite in view, one track per satellite
  (`lib/tch_c_frame_timing.h`). From a track the first UW symbol of every downlink slot and
  its carrier are known to about 1 µs and 10 Hz.
- **Survey.** In every downlink slot of every track, one FFT of the slot times the UW
  correlates all raster positions of a type at once (`lib/tch_c_slot_survey.h`). Positions
  whose normalized UW correlation (rho2, 1 for a clean UW) reaches the threshold (TCH_C1
  0.6, TCH_C2 and C8 0.55) are cut out at one sample per symbol and go to `tch_c_demod`.
- **History.** The prober keeps `--tch-c-history` seconds of samples and probes a slot only
  after that delay, so a satellite found by its first ring alert is probed back to the start
  of the buffer, and the slot times are interpolated between ring alerts.
- **Uplink** (TCH_C1). Seen from a receiver near the terminal, the UL slots are at
  (t + phase(t)) mod 90 ms = W0 + 8.407 ms × slot, phase being the downlink frame phase of
  the satellite. W0 is twice the frame epoch in the stream's time base, the same for all
  satellites, so it is learned from the UL frames: until the first one, a sweep searches all
  of the frame on the channels with downlink traffic. Then each UL slot is searched on the
  channels whose downlink slot of the same number carried a frame (a terminal's UL burst is
  on its DL channel), ±0.35 ms and −1…+6.5 kHz around the carrier with the Doppler
  removed, and once a link is found ±50 µs around its last frame. `--no-tch-c-uplink`
  turns this off.

On `sat-147` (12 MS/s, 279 s, `-D 8`), distinct frames with a valid CRC:

| | Energy path only (before the prober) | `probe` | `both` |
| --- | --- | --- | --- |
| TCH_C1 DL | 30,766 | 38,530 | 38,689 |
| TCH_C2 / C8 QPSK / C8 16APSK DL | 207 / 307 / 228 | 451 / 488 / 260 | 456 / 513 / 266 |
| TCH_C1 UL | 274 | 283 | 297 |
| CPU (user) | 1,192 s | 510 s | 1,272 s |

`probe` alone finds 99.5% of what `both` finds at 40% of the CPU time; it misses bursts of
satellites before their first ring alert and some bursts the energy path happens to
decode. Among the slots inside active TCH_C1 downlink links still without a frame, 76% have
no UW at all (rho2 < 0.1; idle links send only every 4th frame or so); 7% have a UW at the
CRC limit (rho2 0.6–0.8, ~4 dB Es/N0). The UL bursts on that recording are all strong
(rho2 0.8–1) and the UL slots between them empty. A control run with all probe times
shifted into the bursts (`IRIDIUM_PROBE_OFFSET_US=4190`) gives no frame.

## Output

With `--tch-c`, every TCH_C burst that passes the CRC gives one line. Its fields are
those of `RAW:` lines plus the burst type, and the bits are the decoded payload:

`NXD: i-1623786748-t1 0000835.1945 1623787583 N:22.41-99.16 I:00000006820 100% 0.01579 310 NC1 1000100011...`

(bits shortened)

| Column | Example | Explanation |
| --- | --- | --- |
| 1 | `NXD:` | `NXD:` for downlink, `NXU:` for uplink (from the UW) |
| 2–6 | | As for `RAW:` lines (file info, time in ms, frequency in Hz, magnitude and noise, ID). The time is that of the first UW symbol. |
| 7 | `100%` | Share of the symbols within ±22° of the diagonals after phase correction |
| 8 | `0.01579` | Mean magnitude of the symbols |
| 9 | `310` | Number of payload bits: 310 for TCH_C1, 566 for TCH_C2, 2422 for TCH_C8 QPSK, 4838 for TCH_C8 16APSK |
| 10 | `NC1` | Burst type: `NC1` TCH_C1, `NC2` TCH_C2, `NC8` TCH_C8 QPSK, `NCA` TCH_C8 16APSK |
| 11 | `1000100011...` | Decoded and descrambled payload without the CRC-10, in transmit order (not symbol-reversed like the bits of `RAW:` lines) |

Bursts failing the CRC give no line. A CRC-10 lets about 1 in 1,024 wrong decodes through,
so expect occasional lines with a wrong payload from weak bursts.

## Options

| Option | Default | |
| --- | --- | --- |
| `--tch-c` | off | Enable TCH_C reception |
| `--tch-c-types` | `1,2,8` | TCH_C types to look for |
| `--wideband-db` | 6 with `--tch-c-detect energy` or `both`, else 0 | Threshold of the wideband detection, 0 disables it |
| `--wideband-width` | 30000 | Bandwidth in Hz averaged by the wideband detection |
| `--track-frames` | 4 with `--tch-c-detect energy` or `both`, else 0 | TDMA frames to look at after a TCH_C like burst, 0 disables tracking |
| `--tch-c-detect` | `probe` | `energy` (burst detection and TCH_C downmixers), `probe` (slot probing) or `both` |
| `--tch-c-history` | 0.75 | Seconds of samples the prober keeps |
| `--tch-c-probe-threshold` | 0 (per type) | Minimum rho2 of a probed downlink burst |
| `--no-tch-c-uplink` | | Do not probe TCH_C1 uplink slots |

## Cost and limitations

- With the energy path, CPU time is several times that of legacy-only processing. 20 s of
  a 10 MS/s test recording with `-D 8`, `--tch-c-detect energy`:

  | Options | CPU time | TCH_C1 / C2 / C8 lines | `RAW:` lines |
  | --- | --- | --- | --- |
  | (none) | 22 s | | 12,078 |
  | `--tch-c --tch-c-types 1` | 44 s | 365 / – / – | 12,910 |
  | `--tch-c --tch-c-types 1,2` | 52 s | 365 / 40 / – | 12,910 |
  | `--tch-c` | 63 s | 365 / 40 / 10 | 12,910 |

  The `RAW:` lines grow by 7% through the wideband detection; how many of them are valid
  legacy frames rather than misread TCH_C bursts was not checked. TCH_C1 is the most expensive
  part, as its downmixer looks at every burst. Use `--tch-c-types 1` or `1,2` on slow
  machines.
- The TCH_C1 downmixer still looks at every burst: its bandwidth is about that of a
  legacy burst, and neither the bandwidth nor the spectral line tell them apart reliably.
  It only filters the part of a burst where the UW could be and drops the burst if the
  correlation with the UW is poor (see Downmixing), which costs about half of what the
  legacy downmixer needs for a burst.
- For the energy path the detector is the limit, not the demodulation: the slot prober
  finds 25% more TCH_C1 frames on sat-147, about 5 dB below the detection threshold.
- `probe` misses what comes before the first ring alert of a satellite (a satellite whose
  ring alerts are not received is not probed at all), uplink frames before the UL timing
  is learned, and uplink other than TCH_C1. Without the wideband detection, 5-10% fewer
  legacy uplink frames are found on sat-147 than with `--tch-c-detect energy`.
- Downlink TCH_C1, C2 and C8 (QPSK and 16APSK) and uplink TCH_C1 decode on real recordings.
  Uplink TCH_C2 and TCH_C8 are untested: there were no such bursts in the test recordings.
