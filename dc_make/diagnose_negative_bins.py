#!/usr/bin/env python

# Diagnostic: reports every (channel, category, process) combination that
# fails DatacardMaker's negative-bin resolution, in one pass, instead of
# discovering them one at a time via repeated create_datacards.py runs that
# each stop at the first failure. Reuses the real DatacardMaker/addProcess
# logic (not a reimplementation) so results match create_datacards.py exactly.
# Does not write any output; --output is unused by DatacardMaker.addProcess.

import argparse
import os
import sys

file_dir = os.path.dirname(os.path.abspath(__file__))
pkg_dir = os.path.dirname(file_dir)
base_dir = os.path.dirname(pkg_dir)
pkg_dir_name = os.path.split(pkg_dir)[1]
if base_dir not in sys.path:
    sys.path.append(base_dir)
__package__ = pkg_dir_name

from StatInference.dc_make.maker import DatacardMaker

parser = argparse.ArgumentParser(description="Find all negative-bin failures in one pass.")
parser.add_argument("--input", required=True, type=str)
parser.add_argument("--config", required=True, type=str)
parser.add_argument("--hist-bins", required=False, type=str, default=None)
args = parser.parse_args()

maker = DatacardMaker(args.config, args.input, hist_bins=args.hist_bins)

failures = []

# Signals first (populates signal_hists_by_key, needed by getRelevantBins for
# the background pass below), same order as DatacardMaker.createDatacards().
for era, channel, category in maker.ECC():
    for name, p in maker.processes.items():
        if name not in maker.channel_processes[channel] or not p.is_signal:
            continue
        try:
            maker.addProcess(name, era, channel, category)
        except RuntimeError as e:
            failures.append((channel, category, name, str(e).splitlines()[-1]))

for era, channel, category in maker.ECC():
    for name, p in maker.processes.items():
        if name not in maker.channel_processes[channel] or p.is_signal:
            continue
        try:
            maker.addProcess(name, era, channel, category)
        except RuntimeError as e:
            failures.append((channel, category, name, str(e).splitlines()[-1]))

print(f"\n{len(failures)} failures:")
by_process = {}
for channel, category, name, msg in failures:
    print(f"  {channel}/{category}/{name}: {msg}")
    by_process.setdefault(name, set()).add(f"{channel}/{category}")

print("\nBy process:")
for name, combos in sorted(by_process.items()):
    print(f"  {name}: {len(combos)} combo(s) -- {sorted(combos)}")
