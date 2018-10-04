import os
import tempfile
import unittest
from subprocess import call

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from test_func.test_execution import fwd, test_data_dir, tab_output_equal, base_script, root


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

        self.cmd = base_script + [
            '-blast_in', os.path.join(fwd, test_data_dir, 'RF00001_short.blastout'),
            '-blast_query', os.path.join(fwd, test_data_dir, 'RF00001.fasta'),
            '-blast_db', os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--mode', 'simple',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
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

    def test_alifold_refold(self):
        self.run('alifold_refold')

    def test_alifold_refold_rnafold_c(self):
        self.run('alifold_refold_rnafold_c')

    def test_alifold_unpaired_conserved_refold(self):
        self.run('alifold_unpaired_conserved_refold')

    def test_dh_clustal_alifold_conserved_ss_rnafoldc(self):
        self.run('dh_clustal_alifold_conserved_ss_rnafoldc')

    def test_dh_clustal_alifold_refold(self):
        self.run('dh_clustal_alifold_refold')

    def test_dh_clustal_alifold_refold_rnafoldc(self):
        self.run('dh_clustal_alifold_refold_rnafoldc')

    def test_dh_tcoffee_alifold_conserved_ss_rnafoldc(self):
        self.run('dh_tcoffee_alifold_conserved_ss_rnafoldc')

    def test_dh_tcoffee_alifold_refold(self):
        self.run('dh_tcoffee_alifold_refold')

    def test_dh_tcoffee_alifold_refold_rnafoldc(self):
        self.run('dh_tcoffee_alifold_refold_rnafoldc')

    def test_muscle_alifold_refold(self):
        self.run('muscle_alifold_refold')

    def test_muscle_alifold_refold_rnafold_c(self):
        self.run('muscle_alifold_refold_rnafold_c')

    def test_muscle_alifold_unpaired_conserved_refold(self):
        self.run('muscle_alifold_unpaired_conserved_refold')

    def test_pairwise_centroid_homfold(self):
        self.run('pairwise_centroid_homfold')

    def test_rfam_rnafoldc(self):
        self.run('rfam_rnafoldc')

    # ----only if mfold is installed----
    # def test_rfam_subopt(self):
    #     self.run('rfam_subopt')
    #
    # def test_subopt_fold_query(self):
    #     self.run('subopt_fold_query')
    #
    # def test_subopt_fold_clustal_alifold(self):
    #     self.run('subopt_fold_clustal_alifold')
    #
    # def test_subopt_fold_muscle_alifold(self):
    #     self.run('subopt_fold_muscle_alifold')

    def test_rnafold(self):
        self.run('rnafold')

    def test_tcoffee_rcoffee_alifold_conserved_ss_rnafoldc(self):
        self.run('tcoffee_rcoffee_alifold_conserved_ss_rnafoldc')

    def test_tcoffee_rcoffee_alifold_refold(self):
        self.run('tcoffee_rcoffee_alifold_refold')

    def test_tcoffee_rcoffee_alifold_refold_rnafoldc(self):
        self.run('tcoffee_rcoffee_alifold_refold_rnafoldc')


if __name__ == '__main__':
    unittest.main()
