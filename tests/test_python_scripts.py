import unittest
import os
import sys
import pathlib
cpath = pathlib.Path(__file__).parent.absolute() 
sys.path.append( str(cpath / "../bin/") )
sys.path.append( str(cpath) )
import pyimpasr.tree_state_extraction as tse

class CsvFunctionsTester(unittest.TestCase):
    def test_tip_state_mapping_from_tree(self):
        newick_fname = "tests/test_data/test.newick"
        val_tip_state_mapping = sorted([
                            ("Wales/A5314Z/2021", "wales"),
                            ("Wales/A5A4Z/2021", "wales"),
                            ("PHWC1231G", "wales"),
                            ("England/GHAE8/2020", "england"),
                            ("ageaga2021/2021", "world")
                            ])
        testfunc = tse.tip_state_mapping_from_tree
        res = sorted([tup for tup in testfunc(newick_fname)])
        for j, tup in enumerate(res):
            self.assertEqual(val_tip_state_mapping[j], tup)
        
