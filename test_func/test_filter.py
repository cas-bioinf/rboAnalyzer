import os
import unittest

from Bio.Blast import NCBIXML
from rna_blast_analyze.BR_core.BA_support import blast_hsps2list
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits


class TestFilterBlastHits(unittest.TestCase):
    def setUp(self):
        bixml = os.path.abspath(os.path.dirname(__file__) + '/test_data/web_multi_hit.xml')
        with open(bixml, 'r') as b:
            self.blast = blast_hsps2list([i for i in NCBIXML.parse(b)][0])

    def test_filter_by_eval(self):
        filtered = filter_by_eval(self.blast, '<', 5)
        self.assertEqual(len(filtered), 12)

    def test_filter_by_bits(self):
        filtered = filter_by_bits(self.blast, '>', 30)
        self.assertEqual(len(filtered), 12)
