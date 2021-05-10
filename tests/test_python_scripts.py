import unittest
import os
import sys
import pathlib
cpath = pathlib.Path(__file__).parent.absolute() 
sys.path.append( str(cpath / "../bin/") )
import pyimpasr

class CsvFunctionsTester(unittest.TestCase):
    def test_get_taxon_states
