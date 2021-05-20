#!/usr/bin/env python3

"""
This script takes a .fa.gz or .fa file and a list of headers.
Sequences with headers present in the list are outputted, others are not.

Usage:
    python3 extract_by_token_list_and_rename_unique.py test.fa.gz test_headers.txt
"""

import sys
import gzip
from Bio import SeqIO

if __name__ == "__main__":

    with open(sys.argv[2]) as f:
        tokens = set([l.strip() for l in f])

    retrieved = set()

    handle = sys.argv[1]
    if handle.endswith(".gz"):
        handle = gzip.open(handle, "rt")
    for r in SeqIO.parse(handle, "fasta"): 
        patt = r.description.rstrip()
        if patt in tokens:
            retrieved.add(patt)
            print(">"+patt)
            print(str(r.seq))
            tokens.remove(patt)
                

    for token in tokens:
        if token not in retrieved:
            sys.stderr.write("%s NOT RETRIEVED!\n" % token)

    if type(handle) != str:
        handle.close()


