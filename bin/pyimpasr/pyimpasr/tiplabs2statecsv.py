import sys
from get_tree_tax_labels_biopython import get_tree_tip_labels

tip_labels = get_tree_tip_labels(sys.argv[1])

print("tip\tstate")
for l in tip_labels:
    loc = l.split("/")[0].replace("'", "")
    if loc in ["Wales", "England", "Northern_Ireland", "Scotland"]: loc = loc.lower()
    else: loc = "world"
    # Additionally, if it is only a heron, use that
    if "PHWC" in l:
        loc = "wales"
    print(l.rstrip(), loc, sep='\t')
#    print("'" + l.rstrip() + "'", loc, sep='\t')
