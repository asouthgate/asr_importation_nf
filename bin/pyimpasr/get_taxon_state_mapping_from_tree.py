import get_treefile_tip_labels
from Bio import Phylo

def get_treefile_tip_labels(newick_fname):
    """
    Extract tip labels from a newick file.

    Args:
        newick_fname: newick file name

    Returns:
        labs: tip label names
    """
    treegen = Phylo.parse(newick_fname, "newick")
    tree = next(treegen)
    labs = [n.name for n in tree.get_terminals()]
    return labs

def get_treefile_tip_state_mapping(newick):
    """
    Extract tip label: state mapping from newick file.

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


