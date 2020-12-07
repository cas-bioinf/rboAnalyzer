import unittest
import tempfile
import os

from rna_blast_analyze.BR_core.viennaRNA import run_rnaplot

def loadfile_len(f):
    with open(f) as a:
        return len(a.read())

class TestRNAplotCall(unittest.TestCase):
    def setUp(self):
        self.seq = "GCCUCAUAGCUCAGAGGUUUAGAGCACUGGUCUUGUAAACCAGGGGUCGUGAGUUCGAGUCUCACUGGGGCCU"
        self.str = "(((((...((((.........)))).(((((.......))))).....(((((.......))))).))))).."

        fd, fname = tempfile.mkstemp()
        os.close(fd)
        self.file = fname

    def tearDown(self):
        os.remove(self.file)

    def test_basic(self):
        run_rnaplot(self.seq, self.str, outfile=self.file)
        self.assertTrue(loadfile_len(self.file))

    def test_plot_ps(self):
        run_rnaplot(self.seq, self.str, outfile=self.file, format="ps")
        self.assertTrue(loadfile_len(self.file))

    def test_plot_svg(self):
        run_rnaplot(self.seq, self.str, outfile=self.file, format="svg")
        self.assertTrue(loadfile_len(self.file))

    def test_plot_gml(self):
        run_rnaplot(self.seq, self.str, outfile=self.file, format="gml")
        self.assertTrue(loadfile_len(self.file))

    def test_dummy(self):
        self.assertFalse(loadfile_len(self.file))

