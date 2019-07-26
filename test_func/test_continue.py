import json
import os
import tempfile
import unittest
from copy import copy

from unittest.mock import patch

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict, blastsearchrecompute2dict, seqrecord2dict
from test_func.pseudoargs_class import Pseudoargs
from test_func.test_execution import tab_output_equal, tab_output_equal_structures
from rna_blast_analyze.BA import lunch_with_args

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'RF00001_short.blastout')
blast_query = os.path.join(fwd, test_data_dir, 'RF00001.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb')
blast_db_fasta = os.path.join(fwd, test_data_dir, 'blastdb')
blast_output = os.path.join(fwd, test_data_dir, 'RF00001_output.json')
backup_file = blast_in + '.tmp_rboAnalyzer'
test_output_file = os.path.join(fwd, test_data_dir, 'test_output.html')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']


def get_input(text):
    return input(text)


class TestContinuation(unittest.TestCase):
    def setUp(self):
        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold',],
            blast_regexp='(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            mode='simple',
            html=test_output_file,
        )

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, fasta_structures = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        os.close(ff)
        self.fasta_structures = fasta_structures

    @patch('builtins.input', return_value='yes')
    def test_continuation(self, mock_input):

        with open(blast_output, 'r') as f, open(backup_file, 'w') as ff:
            data = blastsearchrecomputefromdict(json.load(f))
            data.args.blast_in = blast_in
            data.args.json = None
            data.args.html = test_output_file
            for hit in data.hits:
                del hit.extension.letter_annotations['rnafold']

            json.dump([blastsearchrecompute2dict(data), {}, {}], ff, indent=2)

        out = lunch_with_args(self.args)
        self.assertEqual(1, 1)

        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            j_obj = json.dumps(blastsearchrecompute2dict(out[0]), indent=2)
            with open(self.json, 'w') as ff:
                ff.write(j_obj)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures,
                ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
            )
            self.assertEqual(t, True)

            ts = tab_output_equal_structures(
                csvfile=self.csv,
                jsonfile=self.json,
                fastastructures=self.fasta_structures
            )
            self.assertTrue(ts)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.fasta_structures,
                    self.fasta,
                    test_output_file,
                ],
                'failed to delete: '
            )

    @patch('builtins.input', return_value='yes')
    def test_continuation2(self, mock_input):

        with open(blast_output, 'r') as f, open(backup_file, 'w') as ff:
            data = blastsearchrecomputefromdict(json.load(f))
            data.args.blast_in = blast_in
            data.args.json = None
            data.args.html = test_output_file
            data.args.prediction_method += ['centroid']

            new_structures = {'rnafold': []}
            for h in data.hits:
                n = copy(h.extension)
                n.letter_annotations['ss0'] = n.letter_annotations['rnafold']
                del n.letter_annotations['rnafold']

                new_structures['rnafold'].append(n)

            json.dump([blastsearchrecompute2dict(data), {k: [seqrecord2dict(s) for s in v] for k, v in new_structures.items()}, {}], ff, indent=2)

        out = lunch_with_args(self.args)
        self.assertEqual(1, 1)

        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            j_obj = json.dumps(blastsearchrecompute2dict(out[0]), indent=2)
            with open(self.json, 'w') as ff:
                ff.write(j_obj)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures,
                ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
            )
            self.assertEqual(t, True)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.fasta_structures,
                    self.fasta,
                    test_output_file
                ],
                'failed to delete: '
            )