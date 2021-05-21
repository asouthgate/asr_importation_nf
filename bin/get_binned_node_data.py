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
    parser.add_argument("--last_date", required=True)
    parser.add_argument("--leaf_only", action="store_true", default=False, required=False)
    parser.add_argument("--data", required=True)
    args = parser.parse_args()

    first_date = dt.datetime.strptime(args.root_date, "%Y-%m-%d")
    last_date = dt.datetime.strptime(args.last_date, "%Y-%m-%d")

    bins = dbf.cal_bins(14, first_date, last_date)

    print("tree,date,n_wales_imp,n_wales_nodes")
    all_df = pd.read_csv(args.data)
    for g, df in all_df.groupby("tree"):
        df = df.loc[df['node_label'] == "W"]
        transdf = df.loc[df['parent_label'] == "G"]

        if args.leaf_only:
            transdf = transdf.loc[df['node_is_terminal'] == True]
            
        trdates = [height2date(h, first_date) for h in transdf['node_height'].values]
        wdates = [height2date(h, first_date) for h in df['node_height'].values]

        binned_trdates = dbf.histogram(trdates, bins)
        binned_wdates = dbf.histogram(wdates, bins)

        for bi in range(len(bins)):
            print(g, bins[bi], binned_trdates[bi], binned_wdates[bi], sep=',')




