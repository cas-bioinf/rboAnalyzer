import unittest
import os
import json

import add_path
from rna_blast_analyze.BR_core import infer_homology
from rna_blast_analyze.BR_core import convert_classes

cwd = os.getcwd()
test_dir = 'test_data'


class Test1WithCMfile(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(cwd, test_dir, 'RF00001_output.json'), 'r') as ff:
            ll = json.load(ff)
            self._bsdata = convert_classes.blastsearchrecomputefromdict(ll)
            self.data = self._bsdata.copy()

    def clean_object(self):
        self.data.query.annotations = dict()
        for h in self.data.hits:
            h.subs[h.ret_keys[0]].annotations = dict()

    def test_simple(self):
        self.clean_object()
        pred, sel = infer_homology.infer_homology(self.data, self.data.args)

    def test_rfam(self):
        self.clean_object()
        self.data.args.use_rfam = True
        self.data.args.blast_query = os.path.join(cwd, test_dir, 'RF00001.fasta')
        pred, sel = infer_homology.infer_homology(self.data, self.data.args)

    def test_with_cm_file(self):
        self.clean_object()
        self.data.args.use_rfam = False
        self.data.cm_file = os.path.join(cwd, test_dir, 'RF00001.cm')
        pred, sel = infer_homology.infer_homology(self.data, self.data.args)


if __name__ == '__main__':
    unittest.main()