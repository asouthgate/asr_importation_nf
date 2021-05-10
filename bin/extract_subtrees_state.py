import sys
from Bio import Phylo
import sys
from pyimpasr.tree_state_extraction import extract_transition_clades

"""
This script can be used to extract taxon sets which are putative imported mono clades.
"""


if __name__ == "__main__":

    # Parse tree
    treegen = Phylo.parse(sys.argv[1], "newick")
    sys.stderr.write("parsing tree file\n")
    tree = next(treegen)
    sys.stderr.write("parsed!\n")

    # Get result
    res = [(l,[t for t in l.get_terminals()]) for l in extract_transition_clades(tree)]

    # Select those groups that are sufficiently large
    bigenough = sorted([(l, terminals) for l, terminals in res if len(terminals) >= int(sys.argv[2])], key=lambda x:len(x[1]))[::-1]


    # Output those that are big enough
    odir = sys.argv[3] + "/"

    if len(bigenough) > 10000: assert False, "Too big. Something gone wrong %d" % len(bigenough)

    for li, pair in enumerate(bigenough):
        l, terminals = pair
        with open((odir + "extracted_subtree_%s") % str(li), "w") as of:
            for k in terminals: of.write(k.name+"\n")
        print(li, l.parental_loc, sep=",")

