#!/usr/bin/env python3

"""
Given a suitable xml template file, this module will inject data.
"""
#TODO: should be implemented in a proper way, like with jinja2.
# Change it!

import sys
import datetime as dt
import argparse as ap
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Phylo
import pandas as pd
from pyimpasr.beast_injector import *
from pyimpasr.datebinfunc import datestr2float

if __name__ == "__main__":
            
    ap = ap.ArgumentParser()
    ap.add_argument("--sample", required=True, type=str)
    ap.add_argument("--fixed_tree", required=True, type=str)
    ap.add_argument("--template", required=True, type=str)
    ap.add_argument("--clock_rate", required=True, type=str)
    args = ap.parse_args()


    # Get records
    pfn = args.sample
    pref = pfn.split("/")[-1]
    fn = pfn + "/" + pref + ".hfix.m.fa"
    assert fn.endswith(".fa")
    with open(fn) as f:
        records = [tup for tup in SimpleFastaParser(f)]

    # Get taxa data
    taxa_data = []
    pref = pfn.split("/")[-1]
    csvfn = pfn + "/" + pref + ".csv"
    df = pd.read_csv(csvfn)
    datestrs = []
    for l in df.values:
        tax = l[1]
        datestr = l[2]
        loc = "WG"[l[0]=="UK-WLS"]
        date = datestr2float(datestr)
        datestrs.append(datestr)
        taxa_data.append((tax,date,loc,""))
    
    with open(args.fixed_tree) as f:
        fixed_tree = [l.rstrip("\n") for l in f][0]    

    with open(args.template) as f:
        for l in f:
            l = l.rstrip()
            if l.startswith("!@"):
                func = funcd[l.strip()]
                if l.startswith("!@TAXON_DATES_AND_LOCS@!"):
                    print(func(taxa_data))
                elif l.startswith("!@SINGLE_ALIGNMENT@!"):
                    print(func(records))
                elif l.startswith("!@FIXED_TREE@!"):
                    print(func(fixed_tree))
                elif l.startswith("!@CLOCK_RATE@!"):
                    print(func(args.clock_rate))
                else:
                    print(func(partition_names))
            else:
                print(l)
