import os
import tempfile
import unittest
from subprocess import call

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from test_func.test_execution import fwd, test_data_dir, tab_output_equal, base_script, root
import json
import glob


class TestDirectExecution_with_prediction(unittest.TestCase):
    def setUp(self):
        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t25')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t26')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t27')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t28')
        os.close(ff)
        self.json = json_file

        ppfile = os.path.join(root, 'rna_blast_analyze', 'BR_core', 'prediction_parameters.json')
        np_fd, self.np_file = tempfile.mkstemp(prefix='rba_', suffix='_t29')
        with open(ppfile, 'r') as f, os.fdopen(np_fd, 'w') as np:
            prediction_parameters = json.load(f)

            for met in prediction_parameters:
                if 'cmscore_percent' in prediction_parameters[met]:
                    del prediction_parameters[met]['cmscore_percent']

                    # effectively set all NoHomologousSeq exceptions where applicable
                    prediction_parameters[met]['cmscore_tr'] = 1000

            json.dump(prediction_parameters, np)

        self.cmd = base_script + [
            '--blast_in', os.path.join(fwd, test_data_dir, 'RF00001_short.blastout'),
            '--blast_query', os.path.join(fwd, test_data_dir, 'RF00001.fasta'),
            '--blast_db', os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', r'(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--mode', 'simple',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--threads', '2',
            '--pm_param_file', self.np_file,
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
        self.run = run

    def tearDown(self):
        files = glob.glob(os.path.join(fwd, test_data_dir, 'RF00001_short.blastout') + '.r-*')
        remove_files_with_try(
            [
                self.html,
                self.csv,
                self.json,
                self.pandas_dump,
                self.np_file
            ] + files
        )

    def test_TurboFold(self):
        self.run('TurboFold')

    def test_TurboFold_fast(self):
        self.run('Turbo-fast')

    def test_clustalo_alifold_refold_rnafoldc(self):
        self.run('C-A-r-Rc')

    def test_clustalo_alifold_unpaired_conserved_refold_rnafoldc(self):
        self.run('C-A-U-r-Rc')

    def test_muscle_alifold_refold_rnafoldc(self):
        self.run('M-A-r-Rc')

    def test_muscle_alifold_unpaired_conserved_refold_rnafoldc(self):
        self.run('M-A-U-r-Rc')

    def test_centroid_homfold(self):
        self.run('centroid')

    def test_centroid_homfold_fast(self):
        self.run('centroid-fast')

    def test_rfam_centroid_homfold(self):
        self.run('rfam-centroid')

    def test_rfam_rnafoldc(self):
        self.run('rfam-Rc')

    # ----only if mfold is installed----
    # As most of the testing is done internally, we expect that this is True
    def test_rfam_subopt(self):
        self.run('rfam-sub')

    def test_subopt_fold_query(self):
        self.run('fq-sub')

    def test_subopt_fold_clustal_alifold(self):
        self.run('C-A-sub')

    def test_subopt_fold_muscle_alifold(self):
        self.run('M-A-sub')
    # -----------------------------------

    def test_rnafold(self):
        self.run('rnafold')


if __name__ == '__main__':
    unittest.main()
