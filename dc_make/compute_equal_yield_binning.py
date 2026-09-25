#!/usr/bin/env python

# Computes, per (channel, category), a set of bin edges over the DNN-score
# axis such that one chosen signal process's yield is split as evenly as
# possible across the requested number of bins (standard equal-population
# quantile binning), snapped to the source histogram's existing bin
# boundaries -- edges cannot be finer than what the source histogram has.
#
# Reads full-resolution shapes from StatInference/dc_make/select_region_shapes.py's
# output (channel/category/process TDirectories). Writes a JSON file in the
# format StatInference/dc_make/binner.py's Binner class expects: a list of
# {"bins": [...], "channels": [...], "categories": [...]} entries, directly
# usable as `create_datacards.py --hist-bins <this file>.json`.

import argparse
import json

import ROOT

DEFAULT_CHANNELS = ["eTau", "muTau", "tauTau"]
DEFAULT_CATEGORIES = ["res1b", "res2b", "boosted"]


def compute_edges(hist, n_bins):
    # Rounded to 4 decimals to match StatInference/common/tools.py:get_new_bin's
    # own `round(old_axis.GetBinLowEdge(...), 4)` comparison -- without this,
    # the raw float32-promoted edge values (e.g. 0.07999999821186066 instead
    # of 0.08) mismatch against that rounded comparison and rebinAndFill
    # raises "Incompatible bin edges" even though the edges are the same bin
    # boundaries.
    n_orig = hist.GetNbinsX()
    orig_edges = [round(hist.GetXaxis().GetBinLowEdge(i), 4) for i in range(1, n_orig + 2)]
    contents = [hist.GetBinContent(i) for i in range(1, n_orig + 1)]
    total = sum(contents)

    if total <= 0:
        return [orig_edges[0], orig_edges[-1]], True

    edges = [orig_edges[0]]
    running = 0.0
    target_idx = 1
    for i in range(n_orig):
        running += contents[i]
        while target_idx < n_bins and running >= target_idx * total / n_bins:
            edge = orig_edges[i + 1]
            if edge > edges[-1]:
                edges.append(edge)
            target_idx += 1
    if edges[-1] != orig_edges[-1]:
        edges.append(orig_edges[-1])
    return edges, False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute equal-signal-yield quantile bin edges per (channel, category)."
    )
    parser.add_argument("--input", required=True, type=str, help="input shape file (channel/category/process)")
    parser.add_argument("--output", required=True, type=str, help="output JSON file")
    parser.add_argument("--signal-process", required=True, type=str, help="signal histogram name to split on")
    parser.add_argument("--n-bins", default=30, type=int, help="requested number of bins (default: 30)")
    parser.add_argument("--channels", nargs="+", default=DEFAULT_CHANNELS, help="channels to process")
    parser.add_argument("--categories", nargs="+", default=DEFAULT_CATEGORIES, help="categories to process")
    args = parser.parse_args()

    f_in = ROOT.TFile(args.input, "READ")

    entries = []
    for channel in args.channels:
        for category in args.categories:
            hist_name = f"{channel}/{category}/{args.signal_process}"
            hist = f_in.Get(hist_name)
            if not hist:
                raise RuntimeError(f"Cannot find histogram {hist_name} in {args.input}")

            edges, degenerate = compute_edges(hist, args.n_bins)
            n_actual = len(edges) - 1
            flag = " (WARNING: zero signal yield -- single bin)" if degenerate else ""
            print(f"{channel}/{category}: requested {args.n_bins} bins, got {n_actual}{flag}")

            entries.append({
                "bins": edges,
                "channels": [channel],
                "categories": [category],
            })

    f_in.Close()

    with open(args.output, "w") as f_out:
        json.dump(entries, f_out, indent=2)
    print(f"Wrote {args.output}")
