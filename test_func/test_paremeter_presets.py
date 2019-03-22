import os
import tempfile
import unittest
import json
from subprocess import call

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from test_func.test_execution import fwd, test_data_dir, tab_output_equal, base_script, root


class TestDirectExecution_with_prediction(unittest.TestCase):
    # check if prediction works with empty parameters
    def setUp(self):
        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t13')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t14')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t15')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t16')
        os.close(ff)
        self.json = json_file

        self.cmd = base_script + [
            '--blast_in', os.path.join(fwd, test_data_dir, 'RF00001_short.blastout'),
            '--blast_query', os.path.join(fwd, test_data_dir, 'RF00001.fasta'),
            '--blast_db', os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--mode', 'simple',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--threads', '2',
            '--pm_param_file', os.path.join(fwd, test_data_dir, 'empty_params.json'),
            '--centroid_fast_preset',
            '--turbofold_fast_preset',
            '--enable_overwrite'
        ]

        def run(mm):
            bb = call(self.cmd + ['--prediction_method', mm], cwd=root)
            self.assertEqual(bb, 0)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
            )
            self.assertTrue(t)

            with open(self.json, 'r') as j:
                data = json.load(j)
                # check if prediction params are empty
                self.assertIn(data['args']['pred_params'][mm]['max_seqs_in_prediction'], [1, 2])

            remove_files_with_try(
                [
                    self.html,
                    self.csv,
                    self.json,
                    self.pandas_dump,
                ],
                ''
            )
        self.run = run

    def test_TurboFold_fast(self):
        self.run('TurboFold_fast')

    def test_centroid_homfold_fast(self):
        self.run('centroid_homfold_fast')


if __name__ == '__main__':
    unittest.main()
