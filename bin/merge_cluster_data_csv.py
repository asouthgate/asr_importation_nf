#!/usr/bin/env python3

"""
This file takes two csv files, the first a (central_sample_id,cluster) cluster mapping, 
and the second (central_sample_id,...) metadata, and merges them into one csv.

Usage:
    python3 merge_cluster_data_csv.py clusters.csv metadata.csv
"""

import sys
import datetime as dt
import pandas as pd
import numpy as np

if __name__ == "__main__":

    # Read cluster designations
    df = pd.read_csv(sys.argv[1], dtype=str)[['central_sample_id','cluster']]

    # Read phylo metadata
    metadf = pd.read_csv(sys.argv[2], sep=',', dtype=str)
    metadf['adm1'] = [c.split("/")[0] for c in metadf['sequence_name']]
    metadf['central_sample_id'] = [c.split("/")[1] for c in metadf['sequence_name']]
    metadf['collection_date'] = metadf['sample_date']

    metadf = metadf[['adm1','central_sample_id', 'collection_date']]

    metadf['date'] = metadf['collection_date']
    metadf = metadf[pd.notnull(metadf['date'])]
    metadf = metadf[metadf['date'] != "None"]
    metadf['date'] = [dt.datetime.strptime(d, "%Y-%m-%d") for d in metadf['date']]

    metadf = metadf[['adm1', 'central_sample_id', 'date']]

    mergdf = df.merge(metadf, how='left', on='central_sample_id')

    mergdf.to_csv(sys.argv[3], index=False)




