from Bio import Phylo
import sys

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
