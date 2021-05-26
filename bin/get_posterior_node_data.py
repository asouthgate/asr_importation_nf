#!/usr/bin/env python3

""" 
This module extracts and collects node data from a 
decorated tree with ancestral states.

Usage:

    python3 get_posterior_node_data.py --trees trees.newick --nthread x
    
"""

import multiprocessing
import sys
from functools import partial
from io import StringIO
import datetime as dt
import pyimpasr.datebinfunc as dbf
import argparse as ap
from Bio import Phylo

def get_node_data(tree):
    """Get data on each node: depth, parent_depth, state, parent_state, is_leaf

    Args:
        tree: Biopython Phylo tree object (decorated data)

    Returns:
        subtrees: list of subtrees
    """

    def get_locstate(node):
        return node.comment.split('"')[-2]

    results = []
    depths = tree.depths()
    assert depths[tree.root] == 0
    nodes = [n for n in tree.get_nonterminals()]
    sys.stderr.write("Got %d non-terminal nodes\n" % len(nodes))
    for clade in nodes:
        for cn in clade:
            results.append((depths[cn], depths[clade], get_locstate(cn), get_locstate(clade), int(cn.is_terminal())))
    return results

def work(pair, leaf_only=False):
    hi, handle = pair
    tree = Phylo.read(handle, "newick")
    res = get_node_data(tree)
    return (hi, *res)

if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument("--trees", required=True)
    parser.add_argument("--nproc", type=int, required=True)
    parser.add_argument("--leaf_only", action="store_true", default=False)
    args = parser.parse_args()

    print("tree,node_height,parent_height,node_label,parent_label,node_is_terminal")
    with open(args.trees) as f:
        handles = [(fi,StringIO(treel)) for fi, treel in enumerate(f)]

    with multiprocessing.Pool(processes=args.nproc) as pool:
        func = partial(work, leaf_only=args.leaf_only)
        results = pool.map(func, handles)
        for quint in results:
            hi = quint[0]
            for tup in quint[1:]:
                print(hi, *tup, sep=',')



