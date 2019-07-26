import unittest
import os
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.expand_by_LOCARNA import locarna_worker
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict

class TestLocarna(unittest.TestCase):
    def setUp(self):
        self.json_file = os.path.abspath(os.path.dirname(__file__) + '/test_data/RF00001_output.json')
        f = open(self.json_file, 'r')
        mydata = json.load(f)
        f.close()
        self.data = blastsearchrecomputefromdict(mydata)

    def test_locarna_worker(self):
        a = locarna_worker((
            self.data.hits[0].source,
            str(self.data.query.seq),
            self.data.args.locarna_params,
            self.data.args.locarna_anchor_length
        ))
        self.assertIsNotNone(a.extension)

    def test_locarna_worker_2(self):

        a = locarna_worker((
            self.data.hits[0].source,
            str(self.data.query.seq) + 'CGATGCTAGT',
            self.data.args.locarna_params,
            self.data.args.locarna_anchor_length
        ))
        self.assertIsNotNone(a.extension)

    def test_locarna_worker_3(self):

        a = locarna_worker((
            self.data.hits[0].source,
            'CGATGCTAGT' + str(self.data.query.seq),
            self.data.args.locarna_params,
            self.data.args.locarna_anchor_length
        ))
        self.assertIsNotNone(a.extension)

    def test_locarna_worker_4(self):
        aq = self.data.query
        nq = SeqRecord(
            Seq('ASCIAPsoeqwofidasf ijajsdofij'),
            id=aq.id
        )

        a = locarna_worker((
            self.data.hits[0].source,
            nq,
            self.data.args.locarna_params,
            self.data.args.locarna_anchor_length
        ))
        self.assertIsNone(a.extension)
        self.assertIn(
            'Could not find match between provided query sequence and QUERY sequence from BLAST output.',
            a.source.annotations['msgs']
        )

    def test_locarna_worker_5(self):
        s = self.data.hits[0].source
        ns = SeqRecord(
            Seq('CAGUCGUAGCUGAUCGUAGCUGUCGAUCGUAGCUGUACGUG'),
            id=s.id,
            annotations=s.annotations,
            letter_annotations=s.letter_annotations
        )

        a = locarna_worker((
            ns,
            str(self.data.query.seq),
            self.data.args.locarna_params,
            self.data.args.locarna_anchor_length
        ))
        self.assertIsNone(a.extension)
        self.assertIn(
            'Could not find match between subject sequence from the database and subject sequence from BLAST output.'
            ' This can be DATABASE problem. Try to recreate the database and/or check if the accessions in the BLAST'
            ' output correspond with accessions in the provided database (i.e. the sequences are the same ones).',
            a.source.annotations['msgs']
        )
