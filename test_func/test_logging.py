import os
import tempfile
import unittest
from subprocess import call
from test_func.test_execution import blast_in, blast_query, blast_db, fwd, base_script, root
test_dir = 'test_data'
test_html_file = os.path.join(fwd, test_dir, 'test_rboAnalyzer.html')


class TestDirectExecution(unittest.TestCase):
    def setUp(self):
        ff, log = tempfile.mkstemp(prefix='rba_', suffix='_t12')
        os.close(ff)
        self.log = log
        self.cmd = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold',
            '--logfile', self.log,
            '--enable_overwrite',
            '--html', test_html_file
        ]

    def test_logger_warn(self):
        bb = call(self.cmd, cwd=root)
        self.assertEqual(bb, 0)

        self.assertTrue(os.path.isfile(self.log))

        with open(self.log, 'r') as l:
            a = l.read()
            self.assertNotRegex(a, '.+')

        os.remove(self.log)
        os.remove(test_html_file)

    def test_logger_info(self):
        bb = call(self.cmd + ['-v'], cwd=root)
        self.assertEqual(bb, 0)

        self.assertTrue(os.path.isfile(self.log))

        with open(self.log, 'r') as l:
            a = l.read()
            self.assertRegex(a, 'INFO')

        os.remove(self.log)
        os.remove(test_html_file)

    def test_logger_debug(self):
        bb = call(self.cmd + ['-vv'], cwd=root)
        self.assertEqual(bb, 0)

        self.assertTrue(os.path.isfile(self.log))

        with open(self.log, 'r') as l:
            a = l.read()
            self.assertRegex(a, 'DEBUG')

        os.remove(self.log)
        os.remove(test_html_file)