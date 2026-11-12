# TCH_C reception

Besides the classic Iridium bursts (25 ksym/s DE-QPSK, printed as `RAW:` lines),
`iridium-extractor` can demodulate and decode TCH_C bursts. They are traffic channel bursts
with coherent modulation and turbo coding: TCH_C1, TCH_C2 and TCH_C8. This is off by
default. Enable it with `--tch-c`.

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
  frequency offsets.
- **Demodulation and decoding** (`tch_c_demod`). Phase and frequency come from a fit to the
  known symbols (UW and pilots). Gain and noise variance give soft bits (max-log for
  16APSK). These are then de-interleaved, turbo decoded (log-MAP, 8 iterations),
  descrambled and checked against the CRC-10. If the UW matches but the CRC fails, a
  decision-directed PLL and the fit followed by the PLL are tried as well.
  For TCH_C8, the modulation that fits the data symbols better is tried first.
- **Sorting** (`frame_sorter`). A TCH_C frame with a valid CRC removes legacy frames at the
  same time and frequency: those are the legacy demodulator misreading the TCH_C burst.
  Among copies of a TCH_C frame, one with a valid CRC beats one without.
- **Backpressure.** The demodulators report every handled frame back to the burst source
  (`frame_handled`), and the downmixers tell it how many frames they published per burst.
  A burst only counts as handled, and frees its place in the queue (`-q`), once all its
  frames are demodulated. Offline this keeps the channels in step, so the sorter's 2 s
  window still catches all duplicates; live, the queue limit drops bursts when the
  demodulators fall behind.

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
| `--wideband-db` | 6 with `--tch-c`, else 0 | Threshold of the wideband detection, 0 disables it |
| `--wideband-width` | 30000 | Bandwidth in Hz averaged by the wideband detection |
| `--track-frames` | 4 with `--tch-c`, else 0 | TDMA frames to look at after a TCH_C like burst, 0 disables tracking |

## Cost and limitations

- CPU time is several times that of legacy-only processing. 20 s of a 10 MS/s test recording
  with `-D 8`:

  | Options | CPU time | TCH_C1 / C2 / C8 lines | `RAW:` lines |
  | --- | --- | --- | --- |
  | (none) | 36 s | | 12,053 |
  | `--tch-c --tch-c-types 1` | 98 s | 359 / – / – | 12,909 |
  | `--tch-c --tch-c-types 1,2` | 135 s | 364 / 40 / – | 12,910 |
  | `--tch-c` | 200 s | 364 / 40 / 10 | 12,910 |

  The `RAW:` lines grow by 7% through the wideband detection; how many of them are valid
  legacy frames rather than misread TCH_C bursts was not checked. TCH_C8 is the most expensive
  part: its downmixer runs at up to 2.4 MS/s for every burst, and a 16APSK burst needs three
  turbo decodes. Use `--tch-c-types 1` or `1,2` on slow machines.
- The TCH_C1 downmixer still processes every burst: its bandwidth is about that of a
  legacy burst, and neither the bandwidth nor the spectral line tell them apart reliably.
- The detector is the limit, not the demodulation: probing every TDMA slot of detected calls
  directly in the samples finds decodable bursts that are about 5 dB below the detection
  threshold.
- Downlink TCH_C1, C2 and C8 (QPSK and 16APSK) and uplink TCH_C1 decode on real recordings.
  Uplink TCH_C2 and TCH_C8 are untested: there were no such bursts in the test recordings.
