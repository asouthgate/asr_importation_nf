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

def cal_state(label):
    """ Calculate a state from a tip label. """
    UK_locations = ["Wales", "England", "Northern_Ireland", "Scotland"]
    state = label.split("/")[0]
    if "'" in state: 
        # Remove quotes from tip labels if they are there
        state = state.replace("'", "")
    if state in UK_locations:
        # If the location is a UK location, use it
        state = state.lower() 
    elif "PHWC" in label:
        state = "wales"
    else: 
        # Otherwise, give it a world location state
        state = "world"
    return state

def tip_state_mapping_from_tree(newick_fname):
    """
    Extract tip label: state mapping from newick file.

    Args:
        newick_fname: .newick tree file name
    
    Yields:
        mapping: tuples (label, state)

    """
    tip_labels = get_treefile_tip_labels(newick_fname)
    if len(tip_labels) > len(set(tip_labels)):
        raise Exception("Invalid tree: non-unique tip labels")
    mapping = []
    for l in tip_labels:
        state = cal_state(l)
        yield (l, state)


