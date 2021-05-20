#!/usr/bin/env python3

from Bio import SeqIO
import Bio.SeqIO.FastaIO as fio
import sys
import numpy as np
import random
import re

def convert_nonstandard_to_gaps(record_pairs):
    """ Force record strings to have characters only in ACGTN

    Args:
        pairs: (header, sequence) tuples for records
        
    Returns:
        pairs: (header, new_sequence) tuples for records

    """
    def fix(c):
        if c.upper() in "ACGT":
            return c
        else:
            return "-"
    return [(h, "".join([fix(c) for c in s])) for h,s in record_pairs]

def mask_spurious_indels(record_pairs,threshold=0.9):
    """ Remove any indel sites from an alignment if frequency < threshold

    Args:
        record_pairs: (h, seq) tuples

    Returns:
        maked_record_pairs: (h, new_seq) tuples, where new_seq is masked
        
    """
    
    def mask(string, inds):
        return "".join([c for i,c in enumerate(string) if i not in inds])

    inds = []
    for i in range(len(pairs[0][1])):
        chars = [s[i] for h,s in record_pairs]
        if chars.count("-")/len(chars) > threshold:
            inds.append(i)

    inds = set(inds)
    return [(h,mask(s,inds)) for h,s in record_pairs]

def trim_to_cds(pairs): 
    """
    Trims record sequences to CDS
    
    Args:
        pairs: (header, seq) fasta record tuples

    Returns:
        new_pairs: (header, new_seq) fasta record_tuples

    """
    # Find the indices of the start and end strings
    start = "ATGGAGAGCCTTG"
    end   = "TAATCTCACATAG"
    linds = []
    rinds = []
    # Randomly sample 10; find at least 10 with these strings present
    c = 0
    error_count = 0
    while c < 10:
        sys.stderr.write("%d\n" % c)
        pair = random.choice(pairs)
        randstr = pair[1].upper()
        startit = [f.start(0) for f in re.finditer(start, randstr)]
        endit = [f.end() for f in re.finditer(end, randstr)]
        if startit != [] and endit != []:
            linds.append(min(startit))
            rinds.append(max(endit))
            c += 1
        error_count += 1
        if error_count > 1000: raise Exception("Cannot find start and end.")
    start_position = min(linds)
    end_position = max(rinds)
    return [(h,s[start_position:end_position]) for h,s in pairs]
    
def gap_condition(s, p=0.1):
    """ Return False if the proportion of gaps is greater than p. """
    if (s.count("-")/len(s)) > p:
        return False
    return True

if __name__ == "__main__":
    import gzip

    fname = sys.argv[1]
    if fname.endswith(".gz"):
        handle = gzip.open(fname, "rt")
    else:
        handle = open(fname)

    pairs = [(h,s) for h,s in fio.SimpleFastaParser(handle)]

    handle.close()

    sys.stderr.write("Converting non-standard bases to gaps\n")
    pairs = convert_nonstandard_to_gaps(pairs)
    sys.stderr.write("Masking spurious indels\n")
    pairs = mask_spurious_indels(pairs)
    sys.stderr.write("Trimming to cds\n")
    pairs = trim_to_cds(pairs)
    sys.stderr.write("Masking spurious indels\n")
    pairs = mask_spurious_indels(pairs)

    for h,s in pairs:
        print(">"+h)
        print(s)
