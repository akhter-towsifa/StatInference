#!/usr/bin/env python

# Iteratively coarsens a compute_equal_yield_binning.py JSON in place,
# targeting only the (channel, category) combos where a PROTECTED background
# (TT, DY) fails DatacardMaker's negative-bin resolution. Per policy: TT/DY
# never get tolerance flags -- a TT/DY failure means the binning is too fine
# for that (channel, category) and should be coarsened there, not papered
# over. Non-protected backgrounds are left to their existing tolerance flags
# in the datacard config; this script does not touch those failures.
#
# Reuses DatacardMaker (real create_datacards.py logic) to detect failures,
# and compute_equal_yield_binning.compute_edges to recompute coarser edges
# for just the failing (channel, category) entries, halving the requested
# bin count each round until the protected processes pass or --min-bins is
# reached.

import argparse
import json
import os
import sys

file_dir = os.path.dirname(os.path.abspath(__file__))
pkg_dir = os.path.dirname(file_dir)
base_dir = os.path.dirname(pkg_dir)
pkg_dir_name = os.path.split(pkg_dir)[1]
if base_dir not in sys.path:
    sys.path.append(base_dir)
__package__ = pkg_dir_name

import ROOT

from StatInference.dc_make.maker import DatacardMaker
from StatInference.dc_make.compute_equal_yield_binning import compute_edges

PROTECTED = {"TT", "DY"}

parser = argparse.ArgumentParser(description="Auto-coarsen a hist_bins JSON where TT/DY fail.")
parser.add_argument("--shapes-dir", required=True, type=str, help="directory containing <era>.root (create_datacards.py --input)")
parser.add_argument("--binning", required=True, type=str, help="hist_bins JSON to fix in place")
parser.add_argument("--config", required=True, type=str)
parser.add_argument("--signal-process", required=True, type=str, help="signal used to recompute coarser quantile edges")
parser.add_argument("--min-bins", default=1, type=int)
parser.add_argument("--max-rounds", default=10, type=int)
args = parser.parse_args()


def find_protected_failures():
    maker = DatacardMaker(args.config, args.shapes_dir, hist_bins=args.binning)
    for era, channel, category in maker.ECC():
        for name, p in maker.processes.items():
            if name not in maker.channel_processes[channel] or not p.is_signal:
                continue
            try:
                maker.addProcess(name, era, channel, category)
            except RuntimeError:
                pass  # only care about background failures here

    failures = set()
    for era, channel, category in maker.ECC():
        for name, p in maker.processes.items():
            if name not in maker.channel_processes[channel] or p.is_signal:
                continue
            try:
                maker.addProcess(name, era, channel, category)
            except RuntimeError:
                if name in PROTECTED:
                    failures.add((channel, category))
    return failures


def load():
    with open(args.binning) as f:
        return json.load(f)


def save(entries):
    with open(args.binning, "w") as f:
        json.dump(entries, f, indent=2)


def find_entry(entries, channel, category):
    for e in entries:
        if e["channels"] == [channel] and e["categories"] == [category]:
            return e
    raise RuntimeError(f"No binning entry for {channel}/{category} in {args.binning}")


f_shapes = ROOT.TFile(os.path.join(args.shapes_dir, "Run3_2022EE.root"), "READ")

for round_i in range(args.max_rounds):
    failures = find_protected_failures()
    if not failures:
        print(f"Round {round_i}: no TT/DY failures remain.")
        break

    print(f"Round {round_i}: TT/DY failing in {sorted(failures)}")
    entries = load()
    for channel, category in failures:
        entry = find_entry(entries, channel, category)
        n_current = len(entry["bins"]) - 1
        if n_current <= args.min_bins:
            raise RuntimeError(
                f"{channel}/{category} already at min-bins={args.min_bins} and TT/DY still fails -- "
                "cannot coarsen further automatically."
            )
        n_new = max(args.min_bins, n_current // 2)
        hist = f_shapes.Get(f"{channel}/{category}/{args.signal_process}")
        edges, degenerate = compute_edges(hist, n_new)
        entry["bins"] = edges
        print(f"  {channel}/{category}: {n_current} -> {len(edges) - 1} bins")
    save(entries)
else:
    raise RuntimeError(f"Did not converge within {args.max_rounds} rounds.")

f_shapes.Close()
print(f"Done. Final binning written to {args.binning}")
