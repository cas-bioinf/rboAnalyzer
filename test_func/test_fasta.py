import os
import tempfile
import unittest
import subprocess

from test_func.test_execution import remove_files_with_try

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'RF00001_short.blastout')
blast_query = os.path.join(fwd, test_data_dir, 'RF00001.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb')
blast_db_fasta = os.path.join(fwd, test_data_dir, 'blastdb')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']
test_html_file = os.path.join(fwd, 'test_html.html')


class TestFasta(unittest.TestCase):
    def test_bad_fasta(self):
        run_with_fasta(
            '>acb\nACGTCAGTGCAT123CAGCAGT\n',
            self.assertEqual, 1
        )

    def test_short_query(self):
        run_with_fasta(
            '>acb\nCGATCGTAGCTGATGCTGAGCTGATGTGCATGCA\n',
            self.assertEqual, 1
        )

    def test_lowercase(self):
        run_with_fasta(
            '>X15126.1_3-120 X15126.1/3-120\n'
            'UACGGCGGCCcagucgauguagcgACCGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCC'
            'CUCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAGAGUAGGACACCGCCGAAC\n',
            self.assertEqual, 0
        )

    def test_bad_header(self):
        run_with_fasta(
            'acb\nCGATCGTAGCTGATGCTGAGCTGATGTGCATGCA\n',
            self.assertEqual, 1
        )

    def test_good_header2(self):
        run_with_fasta(
            '>A-Za-z0-9_|[]-. asd\n'
            'UACGGCGGCCcagucgauguagcgACCGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCC'
            'CUCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAGAGUAGGACACCGCCGAAC\n',
            self.assertEqual, 0
        )

    def test_bad_header2(self):
        run_with_fasta(
            '>A-Za-z0-9_|[]@$% asd\n'
            'UACGGCGGCCcagucgauguagcgACCGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCC'
            'CUCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAGAGUAGGACACCGCCGAAC\n',
            self.assertEqual, 1
        )

    def test_protein(self):
        run_with_fasta(
            '>acb\nEMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL\n',
            self.assertEqual, 1
        )

    def test_bad_description(self):
        run_with_fasta(
            '>A-Za-z0-9_|[] @!)(#&!@$^)!(*@&#_!@()*#&#*(!@#$^*\n'
            'UACGGCGGCCcagucgauguagcgACCGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCC'
            'CUCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAGAGUAGGACACCGCCGAAC\n',
            self.assertEqual, 0
        )

    def test_bare_sequence(self):
        run_with_fasta(
            'UACGGCGGCCcagucgauguagcgACCGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCC'
            'CUCCAGCGCCGAUGGUACUGCACUCGACAGGGUGUGGGAGAGUAGGACACCGCCGAAC\n',
            self.assertEqual, 1
        )


def run_with_fasta(fasta_example, func, expect):
    fd, fdpath = tempfile.mkstemp()
    with os.fdopen(fd, 'w') as f:
        f.write(fasta_example)

    a = base_script + [
        '--blast_in', blast_in,
        '--blast_query', fdpath,
        '--blast_db', blast_db,
        '--mode', 'simple',
        '--blast_regexp', r'"(?<=\|)[A-Z0-9]*\.?\d*$"',
        '--b_type', 'plain',
        '--prediction_method', 'rnafold',
        '--enable_overwrite',
        '--html', test_html_file
    ]
    bb = subprocess.call(' '.join(a), cwd=root, shell=True)
    func(bb, expect)

    remove_files_with_try(
        [
            fdpath,
            test_html_file,
        ],
        ''
    )
    print(test_html_file)