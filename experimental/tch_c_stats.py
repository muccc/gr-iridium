#!/usr/bin/env python3
#
# This file is a part of gr-iridium
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""
Metrics of iridium-extractor output (.bits files) for the TCH_C reception.

usage: tch_c_stats.py [--window T0 T1] [--min-link-frames N] new.bits [baseline.bits]

Prints
- line counts per kind (NXD/NXU/RAW) and burst type
- distinct NX frames: lines within 1 ms and 5 kHz of an earlier line of the same direction
  and type are duplicates
- out-of-order lines: time more than 2 s before the latest time seen so far in the file
- distinct NX frames by source: ID >= 50,000,000,000 from the slot prober, else energy path
- with a baseline: frames only in the baseline / only in the new file / common
- link fill of the TCH_C1 downlink: frames grouped into links by carrier (+-20 kHz, follows
  Doppler) and time mod 90 ms (+-0.5 ms, follows the drift); share of the 90 ms frames between
  the first and the last frame of a link with a hit, and the lengths of the gaps
"""

import argparse
import bisect
import collections
import sys

import numpy as np

DEDUP_MS = 1.0
DEDUP_HZ = 5000.0
OUT_OF_ORDER_MS = 2000.0
PROBE_ID = 50000000000
FRAME_MS = 90.0
LINK_HZ = 20000.0
LINK_PHASE_MS = 0.5
LINK_EXPIRY_MS = 20000.0


class Frames:
    """NX frames of one (direction, type): time [ms], frequency [Hz] and ID, sorted by time"""

    def __init__(self):
        self.t = []
        self.f = []
        self.id = []

    def finish(self):
        order = np.argsort(np.array(self.t), kind="stable")
        self.t = np.array(self.t, dtype=float)[order]
        self.f = np.array(self.f, dtype=float)[order]
        self.id = np.array(self.id, dtype=np.int64)[order]


def parse(path, window):
    counts = collections.Counter()
    frames = collections.defaultdict(Frames)
    max_t = -1e300
    ooo = collections.Counter()
    max_back = collections.defaultdict(float)
    with open(path, errors="replace") as fp:
        for line in fp:
            kind = line[:4]
            if kind not in ("NXD:", "NXU:", "RAW:"):
                continue
            fields = line.split(None, 10)
            try:
                t = float(fields[2])
                f = float(fields[3])
                ident = int(fields[5][2:])
            except (IndexError, ValueError):
                counts[(kind, "bad")] += 1
                continue
            if window and not (window[0] <= t < window[1]):
                continue
            typ = fields[9] if kind != "RAW:" else ""
            counts[(kind, typ)] += 1
            # out of order: relative to the latest time of any line so far
            if t < max_t - OUT_OF_ORDER_MS:
                ooo[kind] += 1
            max_back[kind] = max(max_back[kind], max_t - t)
            max_t = max(max_t, t)
            if kind != "RAW:":
                fr = frames[(kind, typ)]
                fr.t.append(t)
                fr.f.append(f)
                fr.id.append(ident)
    for fr in frames.values():
        fr.finish()
    return counts, frames, ooo, max_back


def dedup(fr):
    """
    Distinct frames: a frame within 1 ms / 5 kHz of a kept frame is a duplicate of it.
    Returns the indices of the kept frames and, per kept frame, whether any of its copies came
    from the slot prober and whether any came from the energy path.
    """
    keep = []
    probe = []
    energy = []
    recent = collections.deque()  # (time, frequency, index into keep) of the last DEDUP_MS
    for i in range(len(fr.t)):
        t = fr.t[i]
        while recent and recent[0][0] < t - DEDUP_MS:
            recent.popleft()
        k = next((k for _, f, k in recent if abs(fr.f[i] - f) <= DEDUP_HZ), None)
        if k is None:
            k = len(keep)
            keep.append(i)
            probe.append(False)
            energy.append(False)
            recent.append((t, fr.f[i], k))
        if fr.id[i] >= PROBE_ID:
            probe[k] = True
        else:
            energy[k] = True
    return np.array(keep, dtype=int), np.array(probe, dtype=bool), np.array(energy, dtype=bool)


def match(t_a, f_a, t_b, f_b):
    """Boolean arrays: which frames of a have a partner in b and vice versa (sorted times)"""
    has_a = np.zeros(len(t_a), dtype=bool)
    has_b = np.zeros(len(t_b), dtype=bool)
    lo = np.searchsorted(t_b, t_a - DEDUP_MS, side="left")
    hi = np.searchsorted(t_b, t_a + DEDUP_MS, side="right")
    for i in range(len(t_a)):
        for j in range(lo[i], hi[i]):
            if abs(f_a[i] - f_b[j]) <= DEDUP_HZ:
                has_a[i] = True
                has_b[j] = True
    return has_a, has_b


def links(t, f, min_frames):
    """Group TCH_C1 DL frames (sorted by time) into links; returns lists of frame numbers"""
    active = []  # [t_last, f_last, frame numbers]
    done = []
    for ti, fi in zip(t, f):
        best = None
        best_d = LINK_PHASE_MS
        still = []
        for link in active:
            gap = ti - link[0]
            if gap > LINK_EXPIRY_MS:
                done.append(link[2])
                continue
            still.append(link)
            if abs(fi - link[1]) > LINK_HZ:
                continue
            n = round(gap / FRAME_MS)
            d = abs(gap - n * FRAME_MS)
            if n >= 1 and d < best_d:
                best, best_d = link, d
        active = still
        if best is None:
            active.append([ti, fi, [0]])
        else:
            n = round((ti - best[0]) / FRAME_MS)
            best[0] = ti
            best[1] = fi
            best[2].append(best[2][-1] + n)
    done.extend(link[2] for link in active)
    return [x for x in done if len(x) >= min_frames]


def fmt_counts(counts):
    print("Lines:")
    for kind in ("NXD:", "NXU:", "RAW:"):
        for (k, typ), n in sorted(counts.items()):
            if k == kind:
                print("  %s %-4s %8d" % (kind, typ, n))


def report(path, window, min_link_frames, baseline):
    counts, frames, ooo, max_back = parse(path, window)
    print("File: %s" % path)
    if window:
        print("Window: %.3f .. %.3f s" % (window[0] / 1e3, window[1] / 1e3))
    fmt_counts(counts)

    print("Out of order (more than %.0f s back):" % (OUT_OF_ORDER_MS / 1e3))
    for kind in ("NXD:", "NXU:", "RAW:"):
        if kind in max_back:
            print("  %s %6d lines, max backward step %.1f ms" % (kind, ooo[kind], max_back[kind]))

    distinct = {}
    print("Distinct NX frames (1 ms / 5 kHz); source: probe = ID >= %d" % PROBE_ID)
    print("  %-9s %8s %9s %11s %11s %12s %6s"
          % ("", "lines", "distinct", "duplicates", "probe only", "energy only", "both"))
    tot = np.zeros(6, dtype=int)
    for key in sorted(frames):
        fr = frames[key]
        keep, probe, energy = dedup(fr)
        distinct[key] = keep
        row = np.array([len(fr.t), len(keep), len(fr.t) - len(keep),
                        np.sum(probe & ~energy), np.sum(energy & ~probe), np.sum(probe & energy)])
        tot += row
        print("  %s %-4s %8d %9d %11d %11d %12d %6d" % ((key[0], key[1]) + tuple(row)))
    print("  %-9s %8d %9d %11d %11d %12d %6d" % (("total",) + tuple(tot)))

    if baseline:
        _, bframes, _, _ = parse(baseline, window)
        print("Comparison with %s (distinct frames):" % baseline)
        print("  %-9s %10s %10s %10s" % ("", "only base", "only new", "common"))
        tot = np.zeros(3, dtype=int)
        for key in sorted(set(frames) | set(bframes)):
            if key in frames:
                fr = frames[key]
                k = distinct[key]
                t_n, f_n = fr.t[k], fr.f[k]
            else:
                t_n = f_n = np.zeros(0)
            if key in bframes:
                fr = bframes[key]
                k = dedup(fr)[0]
                t_b, f_b = fr.t[k], fr.f[k]
            else:
                t_b = f_b = np.zeros(0)
            has_n, has_b = match(t_n, f_n, t_b, f_b)
            row = np.array([np.sum(~has_b), np.sum(~has_n), np.sum(has_n)])
            tot += row
            print("  %s %-4s %10d %10d %10d" % ((key[0], key[1]) + tuple(row)))
        print("  %-9s %10d %10d %10d" % (("total",) + tuple(tot)))

    key = ("NXD:", "NC1")
    if key in frames:
        fr = frames[key]
        k = distinct[key]
        ls = links(fr.t[k], fr.f[k], min_link_frames)
        print("Link fill, TCH_C1 downlink (links with >= %d frames):" % min_link_frames)
        if not ls:
            print("  no links")
            return
        fill = np.array([len(set(x)) / (x[-1] - x[0] + 1) for x in ls])
        n_frames = sum(len(x) for x in ls)
        print("  %d links with %d frames (%d frames in shorter links)"
              % (len(ls), n_frames, len(k) - n_frames))
        pct = (10, 25, 50, 75, 90)
        print("  fill percentiles: " +
              "  ".join("p%d %.3f" % (p, v) for p, v in zip(pct, np.percentile(fill, pct))))
        print("  mean fill %.3f, frames weighted %.3f"
              % (np.mean(fill), n_frames / sum(x[-1] - x[0] + 1 for x in ls)))
        gaps = np.concatenate([np.diff(x) - 1 for x in ls])
        gaps = gaps[gaps > 0]
        edges = [1, 2, 3, 4, 5, 9, 17, 33, 65, 10**9]
        hist = np.histogram(gaps, bins=edges)[0]
        print("  gaps (missing frames): " + "  ".join(
            "%s:%d" % (str(a) if b - a == 1 else "%d-%s" % (a, b - 1 if b < 10**9 else ""), n)
            for a, b, n in zip(edges[:-1], edges[1:], hist)))
        print("  missing frames inside links: %d" % np.sum(gaps))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("new", help=".bits file")
    ap.add_argument("baseline", nargs="?", help="baseline .bits file to compare with")
    ap.add_argument("--window", nargs=2, type=float, metavar=("T0", "T1"),
                    help="only lines with T0 <= time < T1 (s)")
    ap.add_argument("--min-link-frames", type=int, default=10,
                    help="shortest link in the link fill report (default 10)")
    args = ap.parse_args()
    window = (args.window[0] * 1e3, args.window[1] * 1e3) if args.window else None
    report(args.new, window, args.min_link_frames, args.baseline)


if __name__ == "__main__":
    sys.exit(main())
