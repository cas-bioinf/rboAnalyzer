import os
import json
import unittest

from Bio.Blast import NCBIXML
from rna_blast_analyze.BR_core.BA_support import blast_hit_getter_from_hits, blast_hit_getter_from_subseq
from rna_blast_analyze.BR_core.BA_support import blast_hsps2list
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict


class TestFilterBlastHits(unittest.TestCase):
    def setUp(self):
        bixml = os.path.abspath(os.path.dirname(__file__) + '/test_data/web_multi_hit.xml')
        with open(bixml, 'r') as b:
            self.blast = blast_hsps2list([i for i in NCBIXML.parse(b)][0])

        jfile = os.path.abspath(os.path.dirname(__file__) + '/test_data/RF00001_output.json')
        with open(jfile, 'r') as j:
            self.res = blastsearchrecomputefromdict(json.load(j))

    def test_filter_by_eval_hits(self):
        filtered = filter_by_eval(self.blast, blast_hit_getter_from_hits, '<', 5)
        self.assertEqual(len(filtered), 12)

    def test_filter_by_bits_hits(self):
        filtered = filter_by_bits(self.blast, blast_hit_getter_from_hits, '>', 30)
        self.assertEqual(len(filtered), 12)

    def test_filter_by_eval_subseq(self):
        filtered = filter_by_eval(self.res.hits, blast_hit_getter_from_subseq, '<', 10**-30)
        self.assertEqual(len(filtered), 4)

    def test_filter_by_bits_subseq(self):
        filtered = filter_by_bits(self.res.hits, blast_hit_getter_from_subseq, '>', 50)
        self.assertEqual(len(filtered), 4)


