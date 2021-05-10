import sys
from Bio import Phylo
import sys

def extract_transition_clades(tree):
    """Extract subtrees with wales states and non wales parents

    Args:
        tree: Biopython Phylo tree object (decorated data)

    Returns:
        subtrees: list of subtrees
    """
    extracted_clades = set()

    nodes = [n for n in tree.get_nonterminals()]
    sys.stderr.write("Got %d non-terminal nodes\n" % len(nodes))
    for clade in nodes:
        for cn in clade:
            if clade.comment != "wales":
                if cn.comment == "wales":
                    cn.parental_loc = clade.comment
                    extracted_clades.add(cn)

    n_extracted = len(extracted_clades)
    sys.stderr.write("Got %d extracted clades\n" % n_extracted)
    if not n_extracted: 
        raise Exception("No clades extracted!")
    return [st for st in extracted_clades]

treegen = Phylo.parse(sys.argv[1], "newick")
sys.stderr.write("parsing tree file\n")
tree = next(treegen)
sys.stderr.write("parsed!\n")

res = [(l,[t for t in l.get_terminals()]) for l in extract_transition_clades(tree)]

bigenough = sorted([(l, terminals) for l, terminals in res if len(terminals) >= int(sys.argv[2])], key=lambda x:len(x[1]))[::-1]

odir = sys.argv[3] + "/"

if len(bigenough) > 10000: assert False, "Too big. Something gone wrong %d" % len(bigenough)

for li, pair in enumerate(bigenough):
    l, terminals = pair
    with open((odir + "extracted_subtree_%s") % str(li), "w") as of:
        for k in terminals: of.write(k.name+"\n")
    print(li, l.parental_loc, sep=",")

