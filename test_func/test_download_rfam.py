import os
import tempfile
import unittest
import shutil

from rna_blast_analyze.BR_core.BA_support import remove_files_with_try
import glob
import subprocess

# test SE
# test locarna
# test LoSe

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'RF00001_short.blastout')
blast_query = os.path.join(fwd, test_data_dir, 'RF00001.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb')
blast_db_fasta = os.path.join(fwd, test_data_dir, 'blastdb')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']


class TestDownloadRfam(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix='rba_')

        fd, tmp_config_file = tempfile.mkstemp(prefix='rba_')
        with os.fdopen(fd, 'w') as f:
            f.write(
                "[TOOL_PATHS]\n"
                "[DATA]\n"
                "rfam_dir = " + self.tmpdir
            )
        self.tmp_config_file = tmp_config_file

    def tearDown(self):
        files = glob.glob(blast_in + '.r-*')
        remove_files_with_try(
            [self.tmp_config_file
             ] + files
        )
        shutil.rmtree(self.tmpdir)

    def test_download_simple(self):
        a = base_script + [
            '--download_rfam',
            '--config_file', self.tmp_config_file
        ]
        subprocess.call(a)
        files = os.listdir(self.tmpdir)
        self.assertIn('Rfam.cm.gz', files)
        self.assertIn('Rfam.cm', files)
        self.assertIn('Rfam.cm.i1f', files)
        self.assertIn('Rfam.cm.i1i', files)
        self.assertIn('Rfam.cm.i1m', files)
        self.assertIn('Rfam.cm.i1p', files)

    def test_downloaded(self):
        a = base_script + [
            '--download_rfam',
            '--config_file', self.tmp_config_file
        ]
        subprocess.call(a)

        files = os.listdir(self.tmpdir)
        self.assertIn('Rfam.cm.gz', files)
        self.assertIn('Rfam.cm', files)
        self.assertIn('Rfam.cm.i1f', files)
        self.assertIn('Rfam.cm.i1i', files)
        self.assertIn('Rfam.cm.i1m', files)
        self.assertIn('Rfam.cm.i1p', files)

        with tempfile.TemporaryFile(mode='w+', encoding='utf-8') as tmp:
            subprocess.call(a, stdout=tmp, stderr=tmp)

            tmp.seek(0)
            cmdout = tmp.read()
            self.assertIn('No new data. Nothing to do.', cmdout)

            files = os.listdir(self.tmpdir)
            self.assertIn('Rfam.cm.gz', files)
            self.assertIn('Rfam.cm', files)
            self.assertIn('Rfam.cm.i1f', files)
            self.assertIn('Rfam.cm.i1i', files)
            self.assertIn('Rfam.cm.i1m', files)
            self.assertIn('Rfam.cm.i1p', files)
