import json
import os
import tempfile
import unittest
from subprocess import call

import pandas as pd
from Bio import SeqIO

from rna_blast_analyze.BR_core.tools_versions import pred_method_required_tools
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try, parse_named_structure_file
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict
from rna_blast_analyze.BR_core.expand_by_BLAST import blast_wrapper_inner
from rna_blast_analyze.BR_core.expand_by_LOCARNA import locarna_anchored_wrapper_inner
from rna_blast_analyze.BR_core.expand_by_joined_pred_with_rsearch import joined_wrapper_inner
from test_func.pseudoargs_class import Pseudoargs

# test SE
# test locarna
# test LoSe

fwd = os.path.dirname(__file__)
test_dir = 'test_func'
test_data_dir = 'test_data'
blast_in = os.path.join(fwd, test_data_dir, 'RF00001_short.blastout')
blast_query = os.path.join(fwd, test_data_dir, 'RF00001.fasta')
blast_db = os.path.join(fwd, test_data_dir, 'blastdb', 'RF00001-art.blastdb')
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
base_script = ['python3', '-m', 'rna_blast_analyze.BA']


class TestExecution(unittest.TestCase):
    def setUp(self):
        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
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

        ff, fasta = tempfile.mkstemp()
        os.close(ff)
        self.fasta = fasta

        ff, fasta_structures = tempfile.mkstemp()
        os.close(ff)
        self.fasta_structures = fasta_structures

    def test_simple_extension(self):
        _, out = blast_wrapper_inner(self.args, [])
        self.assertEqual(1, 1)

        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            out[0].to_json(self.json)
            out[0].to_html(self.html)
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures
            )
            self.assertEqual(t, True)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.pandas_dump,
                    self.fasta_structures,
                    self.fasta
                ],
                ''
            )

    def test_locarna_extension(self):
        _, out = locarna_anchored_wrapper_inner(self.args, [])
        self.assertEqual(1, 1)
        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            out[0].to_json(self.json)
            out[0].to_html(self.html)
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures
            )
            self.assertEqual(t, True)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.pandas_dump,
                    self.fasta_structures,
                    self.fasta
                ],
                ''
            )

    def test_lose_extenstion(self):
        _, out = joined_wrapper_inner(self.args, [])
        self.assertEqual(1, 1)
        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            out[0].to_json(self.json)
            out[0].to_html(self.html)
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures
            )
            self.assertEqual(t, True)

            remove_files_with_try(
                [
                    self.csv,
                    self.json,
                    self.pandas_dump,
                    self.fasta_structures,
                    self.fasta
                ],
                ''
            )


class TestDirectExecution(unittest.TestCase):
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

    def test_BA_one_pred_method_simple(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_one_pred_method_locarna(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_one_pred_method_joined(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'joined',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold', 'subopt_fold_query',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_shell(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--blast_regexp', '"(?<=\|)[A-Z0-9]*\.?\d*$"',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
        ]
        bb = call(' '.join(a), shell=True, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_shell_relative_path(self):
        a = base_script + [
            '-blast_in', os.path.join(test_dir, test_data_dir, 'RF00001.blastout'),
            '-blast_query', os.path.join(test_dir, test_data_dir, 'RF00001.fasta'),
            '-blast_db', os.path.join(test_dir, test_data_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '"(?<=\|)[A-Z0-9]*\.?\d*$"',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
        ]
        bb = call(' '.join(a), shell=True, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_rfam_simple(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_rfam_locarna(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_rfam_joined(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'joined',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_cm_file_simple(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_cm_file_locarna(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )

    def test_BA_cm_file_joined(self):
        a = base_script + [
            '-blast_in', blast_in,
            '-blast_query', blast_query,
            '-blast_db', blast_db,
            '--mode', 'joined',
            '--blast_regexp', '(?<=\|)[A-Z0-9]*\.?\d*$',
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv'
        ]
        bb = call(a, cwd=root)
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
                self.html
            ],
            ''
        )


def tab_output_equal(csvfile=None, jsonfile=None, pdfile=None, fastafile=None, fastastructures=None):
    c = None
    j = None
    p = None
    ff = None
    fs = None
    if csvfile is not None:
        c = pd.read_csv(csvfile)['best_sequence'].tolist()
    if jsonfile is not None:
        with open(jsonfile, 'r') as f:
            j = blastsearchrecomputefromdict(json.load(f)).hits
            j = [i.subs[i.ret_keys[0]] for i in j]
            j = [str(i.seq) for i in j]
    if pdfile is not None:
        p = pd.read_pickle(pdfile)['best_sequence'].tolist()
    if fastafile is not None:
        with open(fastafile, 'r') as f:
            ff = [str(i.seq) for i in SeqIO.parse(f, format='fasta')]
    if fastastructures is not None:
        fs = [i for i in parse_named_structure_file(fastastructures)]
        fs = [str(i.seq) for i in fs]

    outputs = [c, j, p, ff, fs]
    outputs = [i for i in outputs if i is not None]

    # check length
    if not all(len(i) == len(outputs[0]) for i in outputs):
        raise AssertionError(
            'All output files have not same length'
        )

    # check sequences (all outputs should have same output sequences)
    for ll in zip(*outputs):
        if not all(i == ll[0] for i in ll):
            raise AssertionError(
                'All output sequences should be same.'
            )
    return True


def tab_output_equal_structures(csvfile=None, jsonfile=None, pdfile=None, fastastructures=None):

    names = pred_method_required_tools.keys()

    cc = None
    jj = None
    pp = None
    fsfs = None
    if csvfile is not None:
        c = pd.read_csv(csvfile)
        cc = dict()
        for n in names:
            if n in c.columns:
                cc[n] = c[n].tolist()
    if jsonfile is not None:
        with open(jsonfile, 'r') as f:
            j = blastsearchrecomputefromdict(json.load(f)).hits
            j = [i.subs[i.ret_keys[0]] for i in j]
            jj = dict()
            for s in j:
                for key in s.letter_annotations.keys():
                    if key not in jj:
                        jj[key] = []
                    jj[key].append(s.letter_annotations[key])
    if pdfile is not None:
        p = pd.read_pickle(pdfile)
        pp = dict()
        for n in names:
            if n in p.columns:
                pp[n] = p[n].tolist()
    if fastastructures is not None:
        fs = [i for i in parse_named_structure_file(fastastructures)]
        assert all([len(fs[0].letter_annotations) == k for k in [len(i.letter_annotations) for i in fs]])
        fsfs = dict()
        for s in fs:
            for key in s.letter_annotations.keys():
                if key not in fsfs:
                    fsfs[key] = []
                fsfs[key].append(s.letter_annotations[key])

    outputs = [cc, jj, pp, fsfs]
    outputs = [i for i in outputs if i is not None]

    # check length
    if not all(len(i) == len(outputs[0]) for i in outputs):
        raise AssertionError(
            'All output files have not same length'
        )

    # check structures (all outputs should have same output structures)
    keys = outputs[0].keys()
    for k in keys:
        for ss in zip(*[ll[k] for ll in outputs]):
            if not all(i == ss[0] for i in ss):
                raise AssertionError(
                    'All output structures should be same.'
                )
    return True


if __name__ == '__main__':
    unittest.main()
