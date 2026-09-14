#!/usr/bin/env python

# Selects one QCD region (default OS_Iso, the signal region) out of a merged
# HistMergerTask shape file and flattens it from
#   <channel>/<QCD region>/<category>/<process...> to
#   <channel>/<category>/<process...>
# i.e. the hist_name = f"{channel}/{category}/{process.hist_name}" convention
# StatInference/dc_make/maker.py expects.

import argparse

import ROOT

parser = argparse.ArgumentParser(description='Select a QCD region and flatten channel/region/category shapes.')
parser.add_argument('--input', required=True, type=str, help="input merged shape file")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--region', default='OS_Iso', type=str, help="QCD region to keep (default: OS_Iso)")
args = parser.parse_args()


def IsDir(key):
    return ROOT.gROOT.GetClass(key.GetClassName()).InheritsFrom("TDirectory")


def IsHist(key):
    return ROOT.gROOT.GetClass(key.GetClassName()).InheritsFrom("TH1")


f_in = ROOT.TFile(args.input, 'READ')
f_out = ROOT.TFile(args.output, 'RECREATE', '', 209)

for channel_key in sorted(f_in.GetListOfKeys(), key=lambda k: k.GetName()):
    if not IsDir(channel_key):
        continue
    channel_name = channel_key.GetName()
    channel_dir = f_in.Get(channel_name)

    region_dir = channel_dir.Get(args.region)
    if not region_dir:
        print(f"Warning: region '{args.region}' not found under channel '{channel_name}', skipping")
        continue

    for category_key in sorted(region_dir.GetListOfKeys(), key=lambda k: k.GetName()):
        if not IsDir(category_key):
            continue
        category_name = category_key.GetName()
        category_dir = region_dir.Get(category_name)

        out_dir = f_out.mkdir(f"{channel_name}/{category_name}")
        print(f"{channel_name}/{args.region}/{category_name} -> {channel_name}/{category_name}")
        for hist_key in sorted(category_dir.GetListOfKeys(), key=lambda k: k.GetName()):
            if not IsHist(hist_key):
                continue
            hist = category_dir.Get(hist_key.GetName()).Clone()
            out_dir.WriteTObject(hist, hist_key.GetName(), "Overwrite")

f_in.Close()
f_out.Close()
print(f"Wrote {args.output}")
