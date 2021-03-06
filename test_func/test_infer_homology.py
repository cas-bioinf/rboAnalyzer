import unittest
import os
import json

from rna_blast_analyze.BR_core import infer_homology
from rna_blast_analyze.BR_core import convert_classes

fwd = os.path.dirname(__file__)
test_dir = 'test_data'


class Test1WithCMfile(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(fwd, test_dir, 'RF00001_output.json'), 'r') as ff:
            ll = json.load(ff)
            self._bsdata = convert_classes.blastsearchrecomputefromdict(ll)
            self.data = self._bsdata.copy()

    def clean_object(self):
        self.data.query.annotations = dict()
        for h in self.data.hits:
            h.extension.annotations = dict()
        self.data.args.blast_query = os.path.join(fwd, test_dir, 'RF00001.fasta')
        self.data.cm_file = os.path.join(fwd, test_dir, 'RF00001.cm')

    def test_simple(self):
        self.clean_object()
        cm_file, _ = infer_homology.find_and_extract_cm_model(self.data.args, self.data)
        pred, sel, _ = infer_homology.infer_homology(self.data, self.data.args, cm_file)

    def test_rfam(self):
        self.clean_object()
        self.data.args.use_rfam = True
        cm_file, _ = infer_homology.find_and_extract_cm_model(self.data.args, self.data)
        pred, sel, cmfile = infer_homology.infer_homology(self.data, self.data.args, cm_file)
        os.remove(cmfile)

    def test_with_cm_file(self):
        self.clean_object()
        self.data.args.use_rfam = False
        cm_file, _ = infer_homology.find_and_extract_cm_model(self.data.args, self.data)
        pred, sel, _ = infer_homology.infer_homology(self.data, self.data.args, cm_file)


if __name__ == '__main__':
    unittest.main()