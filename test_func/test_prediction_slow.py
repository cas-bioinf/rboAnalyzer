import os
import tempfile
import unittest
from subprocess import call

from test_func.pseudoargs_class import Pseudoargs
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from rna_blast_analyze.BR_core.expand_by_BLAST import blast_wrapper_inner
from test_func.test_execution import cwd, test_dir, tab_output_equal_structures, tab_output_equal

pm2test = [
    # 'TurboFold',
    'alifold_refold',
    'alifold_refold_rnafold_c',
    'alifold_unpaired_conserved_refold',
    # 'clustalo_alifold_rapidshapes',
    'dh_clustal_alifold_conserved_ss_rnafoldc',
    'dh_clustal_alifold_refold',
    'dh_clustal_alifold_refold_rnafoldc',
    'dh_tcoffee_alifold_conserved_ss_rnafoldc',
    'dh_tcoffee_alifold_refold',
    'dh_tcoffee_alifold_refold_rnafoldc',
    # 'muscle_alifold_rapidshapes',
    'muscle_alifold_refold',
    'muscle_alifold_refold_rnafold_c',
    'muscle_alifold_unpaired_conserved_refold',
    'pairwise_centroid_homfold',
    # 'rcoffee_alifold_rapidshapes',
    # 'rfam_rapidshapes',
    'rfam_rnafoldc',
    'rfam_subopt',
    'rnafold',
    'subopt_fold_clustal_alifold',
    'subopt_fold_muscle_alifold',
    'subopt_fold_query',
    'tcoffee_rcoffee_alifold_conserved_ss_rnafoldc',
    'tcoffee_rcoffee_alifold_refold',
    'tcoffee_rcoffee_alifold_refold_rnafoldc',
]


class TestPredictionMethods(unittest.TestCase):
    def setUp(self):
        self.args = Pseudoargs(
            os.path.join(cwd, test_dir, 'RF00001.fasta'),
            os.path.join(cwd, test_dir, 'RF00001.blastout'),
            os.path.join(cwd, test_dir, 'blastdb', 'RF00001-art.blastdb'),
            b_type='plain',
            prediction_method=pm2test,
            blast_regexp='(?<=\|)[A-Z0-9]*\.?\d*$'
        )

        ff, csv = tempfile.mkstemp()
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp()
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp()
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp()
        os.close(ff)
        self.json = json_file

        ff, fasta_structures = tempfile.mkstemp()
        os.close(ff)
        self.fasta_structures = fasta_structures

    def test_prediction(self):
        _, out = blast_wrapper_inner(self.args, [])
        self.assertEqual(1, 1)
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            out[0].to_json(self.json)
            out[0].to_html(self.html)
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal_structures(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastastructures=self.fasta_structures
            )
            self.assertTrue(t)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.pandas_dump,
                    self.fasta_structures,
                ],
                ''
            )


class TestDirectExecution_with_prediction(unittest.TestCase):
    def setUp(self):
        ff, csv = tempfile.mkstemp()
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp()
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp()
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp()
        os.close(ff)
        self.json = json_file

    def test_BA(self):
        a = [
            'python3',
            '../rna_blast_analyze/BA.py',
            '-blast_in', os.path.join(cwd, test_dir, 'RF00001.blastout'),
            '-blast_query', os.path.join(cwd, test_dir, 'RF00001.fasta'),
            '-blast_db', os.path.join(cwd, test_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--prediction_method'] + pm2test + [
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
        ]
        bb = call(a)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
        )
        self.assertTrue(t)

        remove_files_with_try(
            [
                self.csv,
                self.json,
                self.pandas_dump,
            ],
            ''
        )

    def test_BA_shell(self):
        a = [
            'python3',
            '../rna_blast_analyze/BA.py',
            '-blast_in', os.path.join(cwd, test_dir, 'RF00001.blastout'),
            '-blast_query', os.path.join(cwd, test_dir, 'RF00001.fasta'),
            '-blast_db', os.path.join(cwd, test_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '"(?<=\|)[A-Z0-9]*\.?\d*$"',
            '--b_type', 'plain',
            '--prediction_method'
            ] + pm2test + [
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
        ]
        bb = call(' '.join(a), shell=True)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
        )
        self.assertEqual(t, True)

        remove_files_with_try(
            [
                self.csv,
                self.json,
                self.pandas_dump,
            ],
            ''
        )

if __name__ == '__main__':
    unittest.main()
