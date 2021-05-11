import unittest
import os
import sys
import pathlib
cpath = pathlib.Path(__file__).parent.absolute() 
sys.path.append( str(cpath / "../bin/") )
sys.path.append( str(cpath) )
import datetime as dt
import pyimpasr.tree_state_extraction as tse
import pyimpasr.datebinfunc as dbf

class pyimpasrTester(unittest.TestCase):
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

    def test_datebinfunc_cal_bins_hist(self):
        get_date = lambda x: dt.datetime.strptime(x, "%Y-%m-%d")
        start_date = get_date("2020-01-01")
        end_date = get_date("2020-01-15")
        bins = dbf.cal_bins(5, start_date, end_date)
        val_bins = [get_date(x) for x in ["2020-01-01", "2020-01-06", "2020-01-11", "2020-01-16"]]
        for j, vb in enumerate(val_bins):
            self.assertEqual(vb, bins[j])

        datestrs_to_bin = [
                        "2019-12-15", "2020-01-01", "2020-01-01", "2020-01-11", "2020-01-12",
                        "2020-01-05", "2020-01-09", "2020-01-02", "2020-01-03"
                            ]
        dates_to_bin = [get_date(x) for x in datestrs_to_bin] 
        hist = dbf.histogram(dates_to_bin, bins)
        val_hist = [1+1+1, 1+1+1, 1+1, 1]
        for vi, v in enumerate(val_hist):
            self.assertEqual(v, hist[vi])

        
    
        
        
        
