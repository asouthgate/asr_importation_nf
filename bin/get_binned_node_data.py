#!/usr/bin/env python3

"""
Put node data into bins.
"""

import multiprocessing
import sys
import datetime as dt
import argparse as ap
import pandas as pd
from Bio import Phylo
from io import StringIO
import pyimpasr.datebinfunc as dbf

if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument("--root_date", required=True)
    parser.add_argument("--leaf_only", action="store_true", default=False, required=False)
    parser.add_argument("--data", required=True)
    args = parser.parse_args()


    print("tree,date,n_wales_imp,n_wales_nodes")
    all_df = pd.read_csv(args.data)

    # Check to make sure the max height for every tree is the same
    gb = all_df.groupby("tree")
    maxheights = [max(df['node_height']) for g, df in gb]
    assert len(set(maxheights)) == 1

    # Create bins
    first_date = dt.datetime.strptime(args.root_date, "%Y-%m-%d")
    last_date = dbf.height2date(maxheights[0], first_date)

    bins = dbf.cal_bins(14, first_date, last_date)


    for g, df in all_df.groupby("tree"):
        df = df.loc[df['node_label'] == "W"]
        transdf = df.loc[df['parent_label'] == "G"]

        if args.leaf_only:
            transdf = transdf.loc[df['node_is_terminal'] == True]
            
        trdates = [dbf.height2date(h, first_date) for h in transdf['node_height'].values]
        wdates = [dbf.height2date(h, first_date) for h in df['node_height'].values]

        # Defensive check to ensure dates are within range

        for date in trdates:
            assert date <= last_date, "%s cannot be greater than %s" % (str(date), str(last_date))
            assert date >= first_date, "%s cannot be less than %s" % (str(date), str(first_date))    

        for date in wdates:
            assert date <= last_date, "%s cannot be greater than %s" % (str(date), str(last_date))
            assert date >= first_date, "%s cannot be less than %s" % (str(date), str(first_date))    

        binned_trdates = dbf.histogram(trdates, bins)
        binned_wdates = dbf.histogram(wdates, bins)

        for bi in range(len(bins)):
            print(g, bins[bi], binned_trdates[bi], binned_wdates[bi], sep=',')




