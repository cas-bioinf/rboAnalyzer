import os
import unittest
import re
import logging

from test_func.test_execution import fwd, test_data_dir
from test_func.pseudoargs_class import Pseudoargs
from rna_blast_analyze.BR_core.expand_by_BLAST import blast_wrapper_inner
from rna_blast_analyze.BR_core.expand_by_LOCARNA import locarna_anchored_wrapper_inner
from rna_blast_analyze.BR_core.expand_by_joined_pred_with_rsearch import joined_wrapper_inner
from Bio import SeqIO

logger = logging.getLogger()

blast_in = os.path.join(fwd, test_data_dir, 'edge_cases_out.txt')
blast_query = os.path.join(fwd, test_data_dir, 'edge_cases_q.fa')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'db1')
source_fasta = os.path.join(fwd, test_data_dir, 'blastdb','db1.fa')


class TestExecution(unittest.TestCase):
    def setUp(self):
        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp='(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True
        )

        with open(source_fasta, 'r') as fh:
            self.seqs = [s for s in SeqIO.parse(fh, format='fasta')]

    def test_logs_and_correct_seqs_simple(self):
        logger.setLevel('WARNING')
        with self.assertLogs('rna_blast_analyze.BR_core.extend_hits', level='WARNING') as l:
            _, out = blast_wrapper_inner(self.args, [])

            for o in l.output:
                self.assertIn('Sequence cannot be extended sufficiently', o.split(':')[-1])

            recs = out[0].res_2_record_list()

            with open(source_fasta, 'r') as fh:
                for seq in SeqIO.parse(fh, format='fasta'):
                    ss = re.sub('[a-z]', '', str(seq.seq))
                    if 'rc' == seq.description.split()[-1]:
                        strand = -1
                    else:
                        strand = 1

                    for rec in recs:
                        if seq.id in rec.id:
                            if strand == 1:
                                rec_s = str(rec.seq)
                            else:
                                rec_s = str(rec.seq.reverse_complement())
                            self.assertEqual(ss, rec_s)

    def test_logs_and_correct_seqs_locarna(self):
        logger.setLevel('WARNING')
        with self.assertLogs('rna_blast_analyze.BR_core.extend_hits', level='WARNING') as l:
            _, out = locarna_anchored_wrapper_inner(self.args, [])

            for o in l.output:
                self.assertIn('Sequence cannot be extended sufficiently', o.split(':')[-1])

            recs = out[0].res_2_record_list()

            with open(source_fasta, 'r') as fh:
                for seq in SeqIO.parse(fh, format='fasta'):
                    ss = re.sub('[a-z]', '', str(seq.seq))
                    if 'rc' == seq.description.split()[-1]:
                        strand = -1
                    else:
                        strand = 1

                    for rec in recs:
                        if seq.id in rec.id:
                            if strand == 1:
                                rec_s = str(rec.seq)
                            else:
                                rec_s = str(rec.seq.reverse_complement())
                            self.assertEqual(ss, rec_s)

    def test_logs_and_correct_seqs_joined(self):
        logger.setLevel('WARNING')
        with self.assertLogs('rna_blast_analyze.BR_core.extend_hits', level='WARNING') as l:
            _, out = joined_wrapper_inner(self.args, [])

            for o in l.output:
                self.assertIn('Sequence cannot be extended sufficiently', o.split(':')[-1])

            recs = out[0].res_2_record_list()

            with open(source_fasta, 'r') as fh:
                for seq in SeqIO.parse(fh, format='fasta'):
                    ss = re.sub('[a-z]', '', str(seq.seq))
                    if 'rc' == seq.description.split()[-1]:
                        strand = -1
                    else:
                        strand = 1

                    for rec in recs:
                        if seq.id in rec.id:
                            if strand == 1:
                                rec_s = str(rec.seq)
                            else:
                                rec_s = str(rec.seq.reverse_complement())
                            self.assertEqual(ss, rec_s)
