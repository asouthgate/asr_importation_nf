import sys
from pyimpasr.treetools import get_treefile_tip_labels

def get_taxon_states(newick):
    """
    Args:
        newick: .newick tree file
    
    Yields:
        mapping: tuples (label, state)

    """
    tip_labels = get_treefile_tip_labels(newick)
    if len(tip_labels) > len(set(tip_labels)):
        raise Exception("Invalid tree: non-unique tip labels")
    mapping = []
    UK_locations = ["Wales", "England", "Northern_Ireland", "Scotland"]
    for l in tip_labels:
        loc = l.split("/")[0]
        if "'" in loc: 
            # Remove quotes from tip labels if they are there
            loc = loc.replace("'", "")
        if loc in UK_locations:
            # If the location is a UK location, use it
            loc = loc.lower() 
        else: 
            # Otherwise, give it a world location state
            loc = "world"
        yield (l, state)


if __name__ == "__main__":
    import argparse as ap
    ap = ArgumentParser()
    ap.add_argument("newick_file", help="Newick file to retrieve tip label, state csv from")
    args = ap.parse_args()
    
    print("tip\tstate")
    for l, state in get_taxon_states(args.newick_file)
        print(l.rstrip(), loc, sep='\t')

