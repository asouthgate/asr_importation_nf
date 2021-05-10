import argparse as ap

if __name__ == "__main__":
    args = ap.argumentParser()
    subparsers = args.add_subparsers()
    g1 = subparsers.add_parser("tiplabs2statecsv")
    g1.add_argument("-i", required=True, help="input treefile")
    g2 = subparsers.add_parser("extract_transition_subtrees")
    g2.add_argument("-i", required=True, help="input newick treefile with labelled nodes")
    
