import json
import os
import unittest
from argparse import Namespace

from Bio.Alphabet import RNAAlphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

from rna_blast_analyze.BR_core import convert_classes
from rna_blast_analyze.BR_core.BA_methods import Subsequences, HitList, BlastSearchRecompute
from rna_blast_analyze.BR_core.parser_to_bio_blast import blast_parse_txt


class TestRecrusive(unittest.TestCase):
    def recrusive_compare(self, a, b):
        # assume that list and tuple are ok to interchange
        # because after json, everything is list
        try:
            assert(type(a), type(b), 'rec compare: type not same {} {}'.format(a, b))
        except AssertionError:
            if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
                pass
            else:
                self.assertEqual(type(a), type(b), 'rec compare: type not same {} {}'.format(a, b))

        if self.comparable(a):
            self.assertEqual(a, b)
        else:
            if isinstance(a, dict) and isinstance(b, dict):
                for key in a.keys():
                    self.recrusive_compare(a[key], b[key])
            elif hasattr(a, '__dict__'):
                self.assertEqual(set(vars(a).keys()), set(vars(b).keys()), 'rec compare: params not same')

                for key in vars(a).keys():
                    self.recrusive_compare(getattr(a, key), getattr(b, key))
            else:
                if isinstance(a, (tuple, list)):
                    for i, j in zip(a, b):
                        self.recrusive_compare(i, j)
                else:
                    raise TypeError('cant compare {}'.format(type(a)))

    @staticmethod
    def comparable(o):
        return isinstance(o, (str, int, float, type(None), set))


tc = TestRecrusive('__init__')


class TestBioSeq(unittest.TestCase):
    def test_convert_simple(self):
        s = Seq('CAGTACGTGAC')
        encoded = convert_classes.seq2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.seqfromdict(encoded)
        self.assertEqual(s, decoded, 'sequences dont match')

    def test_convert_alphabet(self):
        s = Seq('ACUACUCAGGAC', RNAAlphabet())
        encoded = convert_classes.seq2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.seqfromdict(encoded)
        self.assertEqual(s, decoded)


class TestBioSeqRecord(unittest.TestCase):
    def test_convert_simple(self):
        s = SeqRecord(
            Seq('ACGUACGUAC'),
            id='aaa'
        )
        encoded = convert_classes.seqrecord2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.seqrecordfromdict(encoded)
        tc.recrusive_compare(s, decoded)

    def test_compare_all(self):
        s = SeqRecord(
            Seq('CAGUGCAGU'),
            letter_annotations={
                'as': '123456789',
                'bb': '.((..))..'
            },
            id='asdoq',
            name='some name',
            description='longer description goes here',
            dbxrefs=['ref1', 'ref2', 'ref3'],
            features=[
                SeqFeature(
                    location=FeatureLocation(1, 3),
                    id='f1',
                ),
                SeqFeature(
                    id='f2',
                    type='bb',
                )
            ],
            annotations={
                'aaa': 123,
                'bbb': ['dumlist', 2],
                'ccc': None
            }
        )
        encoded = convert_classes.seqrecord2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.seqrecordfromdict(encoded)
        tc.recrusive_compare(s, decoded)


class TestBlastTxt(unittest.TestCase):
    def test_convert_simple(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'blast_parse_hits_txt_standalone.txt'), 'r') as f:
            for r in blast_parse_txt(f):
                blast_outputs.append(r)
        encoded_original = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded_original)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)

    def test_convert_multi_queries(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'blast_parse_multi_query_web.txt'), 'r') as f:
            for r in blast_parse_txt(f):
                blast_outputs.append(r)
        encoded = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)

    def test_convert_multi_hsps(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'blast_parse_web_multi_hit.txt'), 'r') as f:
            for r in blast_parse_txt(f):
                blast_outputs.append(r)
        encoded = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)


class TestBlastXML(unittest.TestCase):
    def test_convert_simple(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'web_multi_hit.xml'), 'r') as f:
            for r in NCBIXML.parse(f):
                blast_outputs.append(r)
        encoded_original = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded_original)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)

    def test_convert_multi_queries(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'web_multi_hit.xml'), 'r') as f:
            for r in NCBIXML.parse(f):
                blast_outputs.append(r)
        encoded = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)

    def test_convert_multi_hsps(self):
        blast_outputs = []
        with open(os.path.join('test_data', 'web_multi_hit.xml'), 'r') as f:
            for r in NCBIXML.parse(f):
                blast_outputs.append(r)
        encoded = [convert_classes.blasttodict(i) for i in blast_outputs]
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = [convert_classes.blastfromdict(i) for i in encoded]
        for orig_r, rec_r in zip(blast_outputs, decoded):
            tc.recrusive_compare(orig_r, rec_r)


class TestSubsequences(unittest.TestCase):
    def test_basic(self):
        s = Subsequences(SeqRecord(Seq('ACGUTGU'), id='aa'))
        encoded = convert_classes.subsequences2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.subsequencesfromdict(encoded)
        tc.recrusive_compare(s, decoded)

    def test_full(self):
        s = Subsequences(
            SeqRecord(
                Seq('ACGUACGUGAC'),
                id='qq'
            )
        )
        s.subs = {'asd': SeqRecord(Seq('ACGUR'), id='asd')}
        s.ret_keys = ['1', 'asd', '4t']
        s.query_name = 'myquery'
        s.best_start = 12
        s.best_end = 324
        s.templates = {'t0': 'ACGU'}
        encoded = convert_classes.subsequences2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.subsequencesfromdict(encoded)
        tc.recrusive_compare(s, decoded)


class TestHitList(unittest.TestCase):
    def test_hitlist(self):
        s = HitList()
        a = Subsequences(SeqRecord(Seq('ACGUTGU'), id='aa'))
        b = Subsequences(SeqRecord(Seq('ACGAUCGUGAC'), id='bb'))
        s.append(a)
        s.append(b)
        encoded = convert_classes.hitlist2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.hitlistfromdict(encoded)
        tc.recrusive_compare(s, decoded)


class TestBlastRecompute(unittest.TestCase):
    def test_blastrecompute(self):
        s = BlastSearchRecompute()
        h = HitList()
        a = Subsequences(SeqRecord(Seq('ACGUTGU'), id='aa'))
        b = Subsequences(SeqRecord(Seq('ACGAUCGUGAC'), id='bb'))
        h.append(a)
        h.append(b)
        s.hits = h
        s.query = SeqRecord(Seq('ACGUGUGCA'), id='query')
        s.args = Namespace(**{'aa': 'asdq', 'bb': 'acoi'})
        encoded = convert_classes.blastsearchrecompute2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.blastsearchrecomputefromdict(encoded)
        tc.recrusive_compare(s, decoded)

    def test_blastrecompute_with_blast_data(self):
        # load blast data
        blast_outputs = []
        with open(os.path.join('test_data', 'blast_parse_hits_txt_standalone.txt'), 'r') as f:
            for r in blast_parse_txt(f):
                blast_outputs.append(r)

        s = BlastSearchRecompute()
        h = HitList()
        a = Subsequences(
            SeqRecord(
                Seq('ACGUTGU'),
                id='aa',
                annotations={
                    'blast': (0, blast_outputs[0].alignments[0].hsps[0])
                }
            )
        )
        b = Subsequences(
            SeqRecord(
                Seq('ACGAUCGUGAC'),
                id='bb',
                annotations={
                    'blast': (1, blast_outputs[0].alignments[0].hsps[1])
                }
            )
        )
        h.append(a)
        h.append(b)
        s.hits = h
        s.query = SeqRecord(Seq('ACGUGUGCA'), id='query')
        s.args = Namespace(**{'aa': 'asdq', 'bb': 'acoi'})
        encoded = convert_classes.blastsearchrecompute2dict(s)
        encoded_json = json.dumps(encoded)
        encoded = json.loads(encoded_json)
        decoded = convert_classes.blastsearchrecomputefromdict(encoded)
        tc.recrusive_compare(s, decoded)

if __name__ == '__main__':
    unittest.main()