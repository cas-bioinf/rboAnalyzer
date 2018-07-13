import unittest
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rna_blast_analyze.BR_core import BA_support


class TestNR(unittest.TestCase):
    def setUp(self):
        self.test_cases = []
        for i in range(100):
            seqs = []
            for j in range(1000):
                seqs.append(
                    ''.join([random.choice('AGTC') for x in range(6)])
                )
            l = set(seqs)
            seqs_records = []
            for j, s in enumerate(seqs):
                seqs_records.append(
                    SeqRecord(
                        Seq(s),
                        id=str(j)
                    )
                )
            self.test_cases.append((seqs_records, l))

    def test_nr(self):
        for seqs_list, check_n in self.test_cases:
            nr = BA_support.non_redundant_seqs(seqs_list)
            self.assertEqual(len(nr), len(check_n), 'length do not match')
            self.assertEqual({str(seq.seq) for seq in nr}, check_n)
