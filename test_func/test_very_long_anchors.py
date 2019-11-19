import os
import tempfile
import unittest
from subprocess import call

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
from test_func.test_execution import tab_output_equal
import glob

# test SE
# test locarna
# test LoSe

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'Alig_long.xml')
blast_query = os.path.join(fwd, test_data_dir, 'Alig_long_query.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'Alig_long.blastdb')
# blast_db_fasta = os.path.join(fwd, test_data_dir, 'blastdb')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']
download_base = ['python3', '-m', 'rna_blast_analyze.download_blast_genomes']


class TestExecution(unittest.TestCase):
    def setUp(self):
        # generage blast_db
        r = call(download_base + ['-e', 'marek.schwarz@biomed.cas.cz', '-in', blast_in, '-o', blast_db], cwd=root)
        if r:
            raise ChildProcessError('Failed to download database.')

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t2')
        os.close(ff)
        self.html = html

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, fasta_structures = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        os.close(ff)
        self.fasta_structures = fasta_structures

    def tearDown(self):
        files = glob.glob(blast_in + '.r-*')
        db_files = glob.glob(blast_db + '*')
        remove_files_with_try(
            [
                self.csv,
                self.json,
                self.fasta_structures,
                self.fasta,
                self.html,
            ] + files + db_files
        )

    def test_BA_one_pred_method_locarna(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'locarna',
            '--b_type', 'xml',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
        )
        self.assertTrue(t)
