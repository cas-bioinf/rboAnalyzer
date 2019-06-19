import os
import tempfile
import unittest
from subprocess import call

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from test_func.test_execution import fwd, test_data_dir, tab_output_equal, base_script, root


class TestDirectExecution_with_prediction(unittest.TestCase):
    def setUp(self):
        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t17')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t18')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t19')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t20')
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

    def test_TurboFold(self):
        self.run('TurboFold')

    def test_TurboFold_fast(self):
        self.run('TurboFold_fast')

    def test_clustalo_alifold_refold_rnafoldc(self):
        self.run('clustalo_alifold_refold_rnafoldc')

    def test_clustalo_alifold_unpaired_conserved_refold_rnafoldc(self):
        self.run('clustalo_alifold_unpaired_conserved_refold_rnafoldc')

    def test_muscle_alifold_refold_rnafoldc(self):
        self.run('muscle_alifold_refold_rnafoldc')

    def test_muscle_alifold_unpaired_conserved_refold_rnafoldc(self):
        self.run('muscle_alifold_unpaired_conserved_refold_rnafoldc')

    def test_centroid_homfold(self):
        self.run('centroid_homfold')

    def test_centroid_homfold_fast(self):
        self.run('centroid_homfold_fast')

    def test_rfam_centroid_homfold(self):
        self.run('rfam_centroid_homfold')

    def test_rfam_rnafoldc(self):
        self.run('rfam_rnafoldc')

    # ----only if mfold is installed----
    # As most of the testing is done internally, we expect that this is True
    def test_rfam_subopt(self):
        self.run('rfam_subopt')

    def test_subopt_fold_query(self):
        self.run('subopt_fold_query')

    def test_subopt_fold_clustal_alifold(self):
        self.run('subopt_fold_clustal_alifold')

    def test_subopt_fold_muscle_alifold(self):
        self.run('subopt_fold_muscle_alifold')
    # ----------------------------------

    def test_rnafold(self):
        self.run('rnafold')


if __name__ == '__main__':
    unittest.main()
