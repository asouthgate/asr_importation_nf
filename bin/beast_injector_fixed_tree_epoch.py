#!/usr/bin/env python3

"""
Given a suitable xml template file, this module will inject data.
"""
#TODO: should be implemented in a proper way, like with jinja2.
# Change it!

import sys
import datetime as dt
import argparse as ap
import pathlib
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Phylo
import pandas as pd
from pyimpasr.beast_injector import *
from pyimpasr.datebinfunc import datestr2float

if __name__ == "__main__":
            
    ap = ap.ArgumentParser()
    ap.add_argument("--aln", required=True, type=str)
    ap.add_argument("--fixed_tree", required=True, type=str)
    ap.add_argument("--xml_template", required=True, type=str)
    ap.add_argument("--clock_rate", required=True, type=str)
    ap.add_argument("--dates", required=True, type=str)
    args = ap.parse_args()


    # Get records
    with open(args.aln) as f:
        records = [tup for tup in SimpleFastaParser(f)]

    # Get taxa data
    taxa_data = []
    df = pd.read_csv(args.dates, sep="\t")
    datestrs = []
    for l in df.values:
        tax = l[0]
        datestr = l[1]
        loc = "WG"[tax.split("/")[0]=="Wales"]
        date = datestr2float(datestr)
        datestrs.append(datestr)
        taxa_data.append((tax,date,loc,""))
    
    with open(args.fixed_tree) as f:
        fixed_tree = [l.rstrip("\n") for l in f][0]    

    with open(pathlib.Path(__file__).parent.absolute() / "pyimpasr" / args.xml_template) as f:
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
