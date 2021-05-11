#!/usr/bin/env python3

import sys
import pandas as pd
import pyimpasr.datebinfunc as dbf
import datetime
import datetime as dt
import numpy as np

# Sample stratified across time for each month, and also sample stratified by wales and not Wales
# (Optionally) downsample Wales to check sensitivity

N_SAMPLES = int(sys.argv[1])
#metadata_file = "res/cog_2021-01-31_all_metadata.wlineages.csv"
metadata_file = sys.argv[2]
df = pd.read_csv(metadata_file, dtype=str)
df = df[['adm1','central_sample_id','collection_date','uk_lineage']]
df = df[np.logical_not(df['central_sample_id'].str.contains("MILK"))]
df = df[np.logical_not(df['central_sample_id'].str.contains("QEUH"))]
df = df[np.logical_not(df['central_sample_id'].str.contains("CAMC"))]
df = df[np.logical_not(df['central_sample_id'].str.contains("ALDP"))]

# Drop any na
df = df.dropna()

# Get wales df
wales_df = df.loc[df['adm1'] == "UK-WLS"]
wales_lineages = set(wales_df['uk_lineage'])
# Get those sharing lineages with welsh sequences
lin_in_wal_df = df.loc[df['uk_lineage'].isin(wales_lineages)]

# Assign to date bins
all_dates = [dt.datetime.strptime(d, "%Y-%m-%d") for d in df['collection_date']]
wales_dates = [dt.datetime.strptime(d, "%Y-%m-%d") for d in wales_df['collection_date']]
bindaysize = 28
bins = dbf.cal_bins(bindaysize, min(wales_dates), max(wales_dates)+datetime.timedelta(days=bindaysize))
liw_dates = [dt.datetime.strptime(d, "%Y-%m-%d") for d in lin_in_wal_df['collection_date'].values]

lin_in_wal_df['date_bin'] = [bins[dbf.find_bin(d, bins)] for d in liw_dates]
wales_df['date_bin'] = [bins[dbf.find_bin(d, bins)] for d in wales_dates]
df['date_bin'] = [bins[dbf.find_bin(d, bins)] for d in all_dates]


bin_strat_props = np.array([len(subdf) for g, subdf in lin_in_wal_df.groupby("date_bin") ], dtype=float)
bin_strat_props /= sum(bin_strat_props)
# Now sample
sampled = set()
gi = 0
for g, subdf in lin_in_wal_df.groupby("date_bin"):
    p = bin_strat_props[gi]
    n = int(p * N_SAMPLES)
    gi += 1
    wls_subdf = subdf.loc[subdf['adm1'] == "UK-WLS"] 
    n_wls_samples = int(n * len(wls_subdf)/len(subdf))
    wls_samp = wls_subdf.sample(n_wls_samples, replace=False)
    nonwls_subdf = subdf.loc[subdf['adm1'] != "UK-WLS"]
    nonwls_samp = nonwls_subdf.sample(n-n_wls_samples,replace=False)
    for samp in nonwls_samp['central_sample_id']: sampled.add(samp)
    for samp in wls_samp['central_sample_id']: sampled.add(samp)
    print(g, n_wls_samples, n-n_wls_samples, n)


finaldf = df.loc[df['central_sample_id'].isin(sampled)]
finaldf.to_csv(sys.argv[3], index=False)
