import json
import os
import tempfile
import unittest
from subprocess import call

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.tools_versions import method_required_tools
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try, parse_one_rec_in_multiline_structure
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict, blastsearchrecompute2dict
from test_func.pseudoargs_class import Pseudoargs
from rna_blast_analyze.BR_core.output.htmloutput import write_html_output
from rna_blast_analyze.BA import lunch_with_args
import glob

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


class TestExecution(unittest.TestCase):
    def setUp(self):
        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t2')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t3')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, fasta_structures = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        os.close(ff)
        self.fasta_structures = fasta_structures

        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp='(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            html=self.html,
        )

    def tearDown(self):
        files = glob.glob(blast_in + '.r-*')
        remove_files_with_try(
            [
                self.csv,
                self.json,
                self.pandas_dump,
                self.fasta_structures,
                self.fasta,
                self.html,

            ] + files
        )

    def test_simple_extension(self):
        self.args.mode = 'simple'
        out = lunch_with_args(self.args)
        self.assertEqual(1, 1)

        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            j_obj = json.dumps(blastsearchrecompute2dict(out[0]), indent=2)
            with open(self.json, 'w') as ff:
                ff.write(j_obj)
            with open(self.html, 'wb') as h:
                h.write(write_html_output(out[0]))
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures,
                ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
            )
            self.assertEqual(t, True)

    def test_locarna_extension(self):
        self.args.mode = 'locarna'
        out = lunch_with_args(self.args)
        self.assertEqual(1, 1)
        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            j_obj = json.dumps(blastsearchrecompute2dict(out[0]), indent=2)
            with open(self.json, 'w') as ff:
                ff.write(j_obj)
            with open(self.html, 'wb') as h:
                h.write(write_html_output(out[0]))
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures,
                ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
            )
            self.assertEqual(t, True)

    def test_lose_extenstion(self):
        self.args.mode = 'meta'
        out = lunch_with_args(self.args)
        self.assertEqual(1, 1)
        # test_output
        for i in range(len(out)):
            out[0].to_csv(self.csv)
            j_obj = json.dumps(blastsearchrecompute2dict(out[0]), indent=2)
            with open(self.json, 'w') as ff:
                ff.write(j_obj)
            with open(self.html, 'wb') as h:
                h.write(write_html_output(out[0]))
            out[0].to_pandas_dump(self.pandas_dump)
            out[0].write_results_fasta(self.fasta)
            out[0].write_results_structures(self.fasta_structures)

            t = tab_output_equal(
                csvfile=self.csv,
                jsonfile=self.json,
                pdfile=self.pandas_dump,
                fastafile=self.fasta,
                fastastructures=self.fasta_structures,
                ref_seqs_file=os.path.join(fwd, test_data_dir, 'joined.json')
            )
            self.assertEqual(t, True)


class TestDirectExecution(unittest.TestCase):
    def setUp(self):
        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t7')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t8')
        os.close(ff)
        self.html = html

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t9')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t10')
        os.close(ff)
        self.json = json_file
        self.re = '(?<=\|)[A-Z0-9]*\.?\d*$'

    def tearDown(self):

        files = glob.glob(blast_in + '.r-*')
        remove_files_with_try(
            [
                self.csv,
                self.json,
                self.pandas_dump,
                self.html

            ] + files
        )

    def test_BA_one_pred_method_simple(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
        )
        self.assertTrue(t)

    def test_BA_one_pred_method_locarna(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertTrue(t)

    def test_BA_one_pred_method_joined(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'meta',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'joined.json')
        )
        self.assertTrue(t)

    def test_BA(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--prediction_method', 'rnafold', 'centroid',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertTrue(t)

    def test_BA_shell(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--blast_regexp', '"' + self.re + '"',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--enable_overwrite',
        ]
        bb = call(' '.join(a), shell=True, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertEqual(t, True)

    def test_BA_shell_relative_path(self):
        a = base_script + [
            '--blast_in', os.path.join(test_dir, test_data_dir, 'RF00001_short.blastout'),
            '--blast_query', os.path.join(test_dir, test_data_dir, 'RF00001.fasta'),
            '--blast_db', os.path.join(test_dir, test_data_dir, 'blastdb', 'RF00001-art.blastdb'),
            '--blast_regexp', '"' + self.re + '"',
            '--b_type', 'plain',
            '--prediction_method', 'rnafold',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--enable_overwrite',
        ]
        bb = call(' '.join(a), shell=True, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertEqual(t, True)

    def test_BA_rfam_simple(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
        )
        self.assertTrue(t)

    def test_BA_rfam_locarna(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertTrue(t)

    def test_BA_rfam_joined(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'meta',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--use_rfam',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'joined.json')
        )
        self.assertTrue(t)

    def test_BA_cm_file_simple(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'simple',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
        )
        self.assertTrue(t)

    def test_BA_cm_file_locarna(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'locarna',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'locarna.json')
        )
        self.assertTrue(t)

    def test_BA_cm_file_joined(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db,
            '--mode', 'meta',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--cm_file', os.path.join(fwd, test_data_dir, 'RF00001.cm'),
            '-vv',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'joined.json')
        )
        self.assertTrue(t)

    def test_BA_fasta_db(self):
        a = base_script + [
            '--blast_in', blast_in,
            '--blast_query', blast_query,
            '--blast_db', blast_db_fasta,
            '--db_type', 'fasta',
            '--mode', 'simple',
            '--blast_regexp', self.re,
            '--b_type', 'plain',
            '--html', self.html,
            '--json', self.json,
            '--csv', self.csv,
            '--pandas_dump', self.pandas_dump,
            '--prediction_method', 'rnafold',
            '--enable_overwrite',
        ]
        bb = call(a, cwd=root)
        self.assertEqual(bb, 0)

        t = tab_output_equal(
            csvfile=self.csv,
            jsonfile=self.json,
            pdfile=self.pandas_dump,
            ref_seqs_file=os.path.join(fwd, test_data_dir, 'simple.json')
        )
        self.assertTrue(t)


def tab_output_equal(csvfile=None, jsonfile=None, pdfile=None, fastafile=None, fastastructures=None, ref_seqs_file=None):
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
            j = [i.extension for i in j]
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

    if ref_seqs_file is not None:
        json_file = ref_seqs_file

        f = open(json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = blastsearchrecomputefromdict(mydata).hits
        bb = [i.extension for i in bb]
        bb = [str(i.seq) for i in bb]

        outputs += [bb]

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

    names = method_required_tools.keys()

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
            j = [i.extension for i in j]
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


def parse_named_structure_file(file):
    with open(file, 'r') as f:
        for sr in parse_one_rec_in_multiline_structure(f):
            cf = sr.strip().splitlines()
            cfr = SeqRecord(Seq(cf[1]), id=cf[0])
            cfr.annotations['sss'] = []
            for i, ll in enumerate(cf[2:]):
                structure, str_name = ll.split()
                cfr.letter_annotations[str_name] = structure
                cfr.annotations['sss'].append(str_name)

            yield cfr