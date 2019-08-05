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
        filtered = filter_by_eval(self.blast, blast_hit_getter_from_hits, [('<', 5)])
        for f in filtered:
            self.assertTrue(f[1].expect < 5)

    def test_filter_by_bits_hits(self):
        filtered = filter_by_bits(self.blast, blast_hit_getter_from_hits, [('>', 30)])
        for f in filtered:
            self.assertTrue(f[1].bits > 30)

    def test_filter_by_eval_subseq(self):
        filtered = filter_by_eval(self.res.hits, blast_hit_getter_from_subseq, [('<', 10**-30)])
        for f in filtered:
            self.assertTrue(_get_eval(f) < 10**-30)

    def test_filter_by_bits_subseq(self):
        filtered = filter_by_bits(self.res.hits, blast_hit_getter_from_subseq, [('>', 50)])
        for f in filtered:
            self.assertTrue(_get_bits(f) > 50)

    def test_filter_by_2_bits(self):
        filtered = filter_by_bits(self.res.hits, blast_hit_getter_from_subseq, [('>', 30), ('<', 50)])
        for f in filtered:
            self.assertTrue(30 < _get_bits(f) < 50)

    def test_filter_by_2_eval(self):
        filtered = filter_by_eval(self.res.hits, blast_hit_getter_from_subseq, [('>', 10**-30), ('<', 5)])
        for f in filtered:
            self.assertTrue(10**-30 < _get_eval(f) < 5)


def _get_bits(f):
    return f.source.annotations['blast'][1].bits


def _get_eval(f):
    return f.source.annotations['blast'][1].expect