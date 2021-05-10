#!/usr/bin/env python3

"""
Extract a mapping of tip label \t tip state from the tree labels themselves.
"""

if __name__ == "__main__":

    import argparse as ap
    from pyimpasr.tree_state_extraction import tip_state_mapping_from_tree

    parser = ap.ArgumentParser()
    parser.add_argument("newick_file", help="Newick file to retrieve tip label, state csv from")
    args = parser.parse_args()
    
    print("tip\tstate")
    for l, state in tip_state_mapping_from_tree(args.newick_file):
        print(l.rstrip(), state, sep='\t')

