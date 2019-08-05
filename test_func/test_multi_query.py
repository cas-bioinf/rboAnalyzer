import os
import tempfile
import unittest
from subprocess import call
from copy import deepcopy
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try, iter2file_name
from test_func.test_execution import tab_output_equal

from Bio import SearchIO

# test SE
# test locarna
# test LoSe

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'RF00001_blastout.xml')
blast_query = os.path.join(fwd, test_data_dir, 'RF00001.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb')
blast_db_fasta = os.path.join(fwd, test_data_dir, 'blastdb')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']


class TestExecMultiQuery(unittest.TestCase):
    def setUp(self):
        ff, blast_xml = tempfile.mkstemp(prefix='rba_test_xml', suffix='_t21')
        blast = SearchIO.read(blast_in, 'blast-xml')
        blast1 = deepcopy(blast)
        blast2 = deepcopy(blast)
        dbblast = [blast1.hsp_filter(lambda x: x.bitscore_raw >= 95), blast2.hsp_filter(lambda x: 80 < x.bitscore_raw < 95)]
        with os.fdopen(ff, 'w') as b:
            SearchIO.write(dbblast, b, 'blast-xml')
        self.blast_xml = blast_xml

        ff, double_fasta = tempfile.mkstemp(prefix='rba_test_double_fasta', suffix='_t22')
        with os.fdopen(ff, 'w') as t, open(blast_query, 'r') as s:
            q = s.read()
            t.write('{}\n{}\n'.format(q, q))
        self.query_double = double_fasta

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t7.csv')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t8.html')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t9.pandas_dump')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t10.json')
        os.close(ff)
        self.json = json_file

    def tearDown(self):
        remove_files_with_try(
            [
                self.blast_xml,
                self.query_double,
            ],
            ''
        )

    def test_BA_one_pred_method_simple(self):
        a = base_script + [
            '--blast_in', self.blast_xml,
            '--blast_query', self.query_double,
            '--blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', r'(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'xml',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
            '--threads', '2',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        # we have two input query sequences so output files will be numbered

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 0),
            jsonfile=iter2file_name(self.json, True, 0),
            pdfile=iter2file_name(self.pandas_dump, True, 0),
        )
        self.assertTrue(t)

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 1),
            jsonfile=iter2file_name(self.json, True, 1),
            pdfile=iter2file_name(self.pandas_dump, True, 1),
        )
        self.assertTrue(t)

        remove_files_with_try(
            [
                iter2file_name(self.csv, True, 0),
                iter2file_name(self.json, True, 0),
                iter2file_name(self.pandas_dump, True, 0),
                iter2file_name(self.html, True, 0),
                iter2file_name(self.csv, True, 1),
                iter2file_name(self.json, True, 1),
                iter2file_name(self.pandas_dump, True, 1),
                iter2file_name(self.html, True, 1)
            ],
            ''
        )

    def test_BA_one_pred_method_locarna(self):
        a = base_script + [
            '--blast_in', self.blast_xml,
            '--blast_query', self.query_double,
            '--blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', r'(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'xml',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
            '--threads', '2',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 0),
            jsonfile=iter2file_name(self.json, True, 0),
            pdfile=iter2file_name(self.pandas_dump, True, 0),
        )
        self.assertTrue(t)

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 1),
            jsonfile=iter2file_name(self.json, True, 1),
            pdfile=iter2file_name(self.pandas_dump, True, 1),
        )
        self.assertTrue(t)

        remove_files_with_try(
            [
                iter2file_name(self.csv, True, 0),
                iter2file_name(self.json, True, 0),
                iter2file_name(self.pandas_dump, True, 0),
                iter2file_name(self.html, True, 0),
                iter2file_name(self.csv, True, 1),
                iter2file_name(self.json, True, 1),
                iter2file_name(self.pandas_dump, True, 1),
                iter2file_name(self.html, True, 1)
            ],
            ''
        )

    def test_BA_one_pred_method_joined(self):
        a = base_script + [
            '--blast_in', self.blast_xml,
            '--blast_query', self.query_double,
            '--blast_db', blast_db,
            '--mode', 'meta',
            '--blast_regexp', r'(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'xml',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
            '--threads', '2',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 0),
            jsonfile=iter2file_name(self.json, True, 0),
            pdfile=iter2file_name(self.pandas_dump, True, 0),
        )
        self.assertTrue(t)

        t = tab_output_equal(
            csvfile=iter2file_name(self.csv, True, 1),
            jsonfile=iter2file_name(self.json, True, 1),
            pdfile=iter2file_name(self.pandas_dump, True, 1),
        )
        self.assertTrue(t)

        remove_files_with_try(
            [
                iter2file_name(self.csv, True, 0),
                iter2file_name(self.json, True, 0),
                iter2file_name(self.pandas_dump, True, 0),
                iter2file_name(self.html, True, 0),
                iter2file_name(self.csv, True, 1),
                iter2file_name(self.json, True, 1),
                iter2file_name(self.pandas_dump, True, 1),
                iter2file_name(self.html, True, 1)
            ],
            ''
        )