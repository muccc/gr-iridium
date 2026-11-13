# TCH_C reception

Besides the classic Iridium bursts (25 ksym/s DE-QPSK, printed as `RAW:` lines),
`iridium-extractor` can demodulate and decode TCH_C bursts. They are traffic channel bursts
with coherent modulation and turbo coding: TCH_C1, TCH_C2 and TCH_C8. This is off by
default. Enable it with `--tch-c`.

The bursts are found by probing the TDMA slots of the satellites, timed by their ring
alerts, rather than by the burst detector.

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

TCH_C bursts have no spectral line, so the burst detector of `fft_burst_tagger` misses most
of them, and a call shows up as a run of frames with gaps. `tch_c_slot_prober` looks for
them where the TDMA timing says they can be instead. It reads the stream of the tagger (the
tagger itself runs as without `--tch-c`) and sends its candidates straight to
`tch_c_demod`.

- **Timing.** The ring alerts (legacy downlink frames on the simplex ring channel) give the
  90 ms frame phase and the Doppler of every satellite in view, one track per satellite
  (`lib/tch_c_frame_timing.h`). From a track the first UW symbol of every downlink slot and
  its carrier are known to about 1 µs and 10 Hz.
- **Survey.** In every downlink slot of every track, one FFT of the slot times the UW
  correlates all raster positions of a type at once (`lib/tch_c_slot_survey.h`). Positions
  whose normalized UW correlation (rho2, 1 for a clean UW) reaches the threshold (TCH_C1
  0.6, TCH_C2 and C8 0.55) get a fine timing (TCH_C8, with its 4.17 µs symbols, is
  surveyed at −4, −2, 0, 2 and 4 µs around the prediction: a satellite with few ring alerts
  can be a few µs off, and its timing is refined on all known symbols, UW and pilots,
  within ±0.16 symbol: 16APSK is sensitive to timing errors of a few percent of a
  symbol), are cut out with a matched filter (RRC α=0.2,
  ±8 symbols) at one sample per symbol and go to `tch_c_demod`.
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
- **Demodulation and decoding** (`tch_c_demod`). Phase and frequency come from a fit to the
  known symbols (UW and pilots). Gain and noise variance give soft bits (max-log for
  16APSK). These are then de-interleaved, turbo decoded (log-MAP, up to 8 iterations),
  descrambled and checked against the CRC-10. For probes with a clear UW (rho2 ≥ 0.7) whose
  CRC fails, a decision-directed PLL and the fit followed by the PLL are tried as well.
  For TCH_C8, the modulation that fits the data symbols better is tried first. 16APSK is
  only decoded if the noise variance on the known symbols is below 0.25 of the symbol
  energy (it needs about 0.11): most TCH_C8 bursts failing the CRC are above that.
- **Sorting** (`frame_sorter`). A TCH_C frame with a valid CRC removes legacy frames at the
  same time and frequency: those are the legacy demodulator misreading the TCH_C burst.
  Among copies of a TCH_C frame, one with a valid CRC beats one without.
- **Backpressure.** Reading from a file, the prober waits while 500 (`-q`) probes are not
  demodulated yet or the oldest of them is more than 1 s older than the newest; live it
  drops probes then. Likewise the legacy burst sources wait for the legacy demodulator
  (`frame_handled`). This keeps the paths in step, so the sorter's 2 s window still catches
  all duplicates.

On `sat-147` (12 MS/s, 279 s, `-D 8`), distinct frames with a valid CRC:

| | Slot prober | Energy detection (`schneider/iridium-next-harder`, before the prober) |
| --- | --- | --- |
| TCH_C1 DL | 38,964 | 34,957 |
| TCH_C2 / C8 QPSK / C8 16APSK DL | 498 / 579 / 293 | 383 / 567 / 292 |
| TCH_C1 UL | 283 | 294 |
| Legacy frames parsed | 64,355 | 64,051 |
| CPU (user) / wall time | 532 s / 2:28 | 1,534 s / 4:43 |

The energy detection found 406 frames the prober misses (mostly TCH_C1 bursts at the CRC
limit); the prober finds 4,530 the energy detection missed.

Among the slots inside active TCH_C1 downlink links still without a frame, 76% have no UW
at all (rho2 < 0.1; idle links send only every 4th frame or so); 7% have a UW at the CRC
limit (rho2 0.6–0.8, ~4 dB Es/N0). The UL bursts on that recording are all strong (rho2
0.8–1) and the UL slots between them empty. A control run with all probe times shifted into
the bursts (`IRIDIUM_PROBE_OFFSET_US=4190`) gives no frame.

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
| `--tch-c-history` | 0.75 | Seconds of samples the prober keeps |
| `--tch-c-probe-threshold` | 0 (per type) | Minimum rho2 of a probed downlink burst |
| `--no-tch-c-uplink` | | Do not probe TCH_C1 uplink slots |

## Cost and limitations

- CPU time is about 1.5 times that of legacy-only processing (20 s of sat-147, `-D 8`:
  31 s without `--tch-c`, 47 s with it). Of the prober's own time, about half goes into
  the extraction of the candidates and a quarter into the slot surveys.
- Only satellites whose ring alerts are received are probed, and only from their first
  ring alert on (minus the history). Uplink frames come only from terminals near the
  receiver, and only once the UL timing is learned; uplink TCH_C2 and TCH_C8 are not
  probed.
- The TCH_C output is not reproducible bit for bit from run to run, not even with `-D 1`:
  the TCH_C frames reach the prober (which refines its tracks with them) at times that
  depend on the threads. The frames found differ by a few in thousands.
- Downlink TCH_C1, C2 and C8 (QPSK and 16APSK) and uplink TCH_C1 decode on real recordings.
  Uplink TCH_C2 and TCH_C8 are untested: there were no such bursts in the test recordings.
