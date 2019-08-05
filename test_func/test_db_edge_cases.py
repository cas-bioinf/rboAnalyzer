import os
import unittest
import re

from test_func.test_execution import fwd, test_data_dir
from test_func.pseudoargs_class import Pseudoargs
from rna_blast_analyze.BA import lunch_with_args
from Bio import SeqIO
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
import glob

blast_in = os.path.join(fwd, test_data_dir, 'edge_cases_out.txt')
blast_query = os.path.join(fwd, test_data_dir, 'edge_cases_q.fa')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'db1')
source_fasta = os.path.join(fwd, test_data_dir, 'blastdb','db1.fa')
test_output_file = os.path.join(fwd, test_data_dir, 'rboAnalyzer_test.html' )


class TestExecution(unittest.TestCase):
    def setUp(self):
        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp=r'(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            html=test_output_file,
        )

        with open(source_fasta, 'r') as fh:
            self.seqs = [s for s in SeqIO.parse(fh, format='fasta')]

    def tearDown(self):
        files = glob.glob(blast_in + '.r-*')
        remove_files_with_try(
            [
                test_output_file,
            ] + files
        )

    def test_logs_and_correct_seqs_simple(self):
        with self.assertLogs('rboAnalyzer', level='WARNING') as l:
            args = self.args
            args.mode = 'simple'
            out = lunch_with_args(args)

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
        with self.assertLogs('rboAnalyzer', level='WARNING') as l:
            args = self.args
            args.mode = 'locarna'
            out = lunch_with_args(args)

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
        with self.assertLogs('rboAnalyzer', level='WARNING') as l:
            args = self.args
            args.mode = 'meta'
            out = lunch_with_args(args)

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
