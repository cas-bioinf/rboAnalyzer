import os
import tempfile
import unittest
try:
    import mock
except ImportError:
    from unittest import mock
import json
from Bio import SeqIO
import re
import glob

from test_func.pseudoargs_class import Pseudoargs
from rna_blast_analyze.BA import lunch_with_args
from rna_blast_analyze.BR_core import exceptions
from rna_blast_analyze.BR_core import convert_classes

from test_func.test_execution import remove_files_with_try
from rna_blast_analyze.BR_core.repredict_structures import repredict_structures_for_homol_seqs

from rna_blast_analyze.BR_core.infer_homology import infer_homology, find_and_extract_cm_model

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
ref_json_file = os.path.abspath(os.path.dirname(__file__) + '/test_data/RF00001_output.json')
base_script = ['python3', '-m', 'rna_blast_analyze.BA']


class TestPredictExceptions(unittest.TestCase):
    def setUp(self):
        f = open(ref_json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        self.data = bb

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t2')
        os.close(ff)
        self.html = blast_query + 'test_html.html'

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t3')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, allHits_fasta = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        with os.fdopen(ff, 'w') as f:
            SeqIO.write([i.extension for i in self.data.hits], f, 'fasta')

        self.fasta_structures = allHits_fasta

        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp=r'(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            html=self.html,
        )

        self.func_args = {
            'query': self.data.query,
            'seqs2predict_fasta': allHits_fasta,
            'pred_method_params': {},
            'all_hits_list': [i.extension for i in self.data.hits],
            'seqs2predict_list': [i.extension for i in self.data.hits],
            'use_cm_file': 'abc',
        }

    def tearDown(self):
        remove_files_with_try(
            [
                self.csv,
                self.html,
                self.json,
                self.fasta,
                self.fasta_structures
            ],
            ''
        )

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.rnafold", side_effect=exceptions.RNAfoldException('expected_error_rnafold', 'b'))
    def test_rnafold(self, callMock):
        self.func_args['prediction_method'] = 'rnafold'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ['expected_error_rnafold'])

    @mock.patch("rna_blast_analyze.BR_core.repredict_structures.run_cmemit", side_effect=exceptions.CmemitException('expected_error_cmemit', 'b'))
    def test_cmemit(self, callMock):
        self.func_args['prediction_method'] = 'rfam-centroid'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_cmemit"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.run_cmalign_on_fasta", side_effect=exceptions.CmemitException('expected_error_cmalign', 'b'))
    def test_cmalign(self, callMock):
        self.func_args['prediction_method'] = 'rfam-Rc'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_cmalign"])

    @mock.patch("rna_blast_analyze.BR_core.centroid_homfold.run_centroid_homfold", side_effect=exceptions.CentroidHomfoldException('expected_error_centroid', 'b'))
    def test_centroid_homfold(self, callMock):
        self.func_args['prediction_method'] = 'centroid-fast'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_centroid"])

    @mock.patch("rna_blast_analyze.BR_core.turbofold._turbofold_worker", side_effect=exceptions.TurboFoldException('expected_error_TurboFold', 'b'))
    def test_TurboFold(self, callMock):
        self.func_args['prediction_method'] = 'Turbo-fast'
        self.func_args['threads'] = 1
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        # check that sequences are returned does not have Turbofold structure and have error message that turbo failed
        for s in a:
            self.assertNotIn('ss0', s.letter_annotations)
            self.assertEqual(['expected_error_TurboFold'], s.annotations['msgs'])

    @mock.patch("rna_blast_analyze.BR_core.BA_support.run_muscle", side_effect=exceptions.MuscleException('expected_error_muscle', 'b'))
    def test_muscle(self, callMock):
        self.func_args['prediction_method'] = 'M-A-U-r-Rc'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_muscle"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.compute_clustalo_clasic", side_effect=exceptions.ClustaloException('expected_error_clustalo', 'b'))
    def test_clustalo(self, callMock):
        self.func_args['prediction_method'] = 'C-A-sub'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_clustalo"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.compute_alifold", side_effect=exceptions.RNAalifoldException('expected_error_alifold', 'b'))
    def test_clustalo(self, callMock):
        self.func_args['prediction_method'] = 'C-A-sub'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_alifold"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.rnafold", side_effect=exceptions.NoHomologousSequenceException())
    def test_nh_ext(self, callMock):
        self.func_args['prediction_method'] = 'rnafold'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertTrue(c[0].startswith('No sequence was inferred homologous'))

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.rnafold", side_effect=exceptions.AmbiguousQuerySequenceException())
    def test_nh_ext(self, callMock):
        self.func_args['prediction_method'] = 'rnafold'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertTrue(c[0].startswith('Query sequence contains ambiguous characters'))

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.compute_refold", side_effect=exceptions.RefoldException('expected_error_refold', 'b'))
    def test_refold_pl(self, callMock):
        self.func_args['prediction_method'] = 'C-A-r-Rc'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_refold"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.run_clustal_profile2seqs_align", side_effect=exceptions.RefoldException('expected_error_clustalo', 'b'))
    def test_clustalo_profile(self, callMock):
        self.func_args['prediction_method'] = 'C-A-r-Rc'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_clustalo"])

    @mock.patch("rna_blast_analyze.BR_core.predict_structures.run_hybrid_ss_min", side_effect=exceptions.HybridssminException('expected_error_hybrid_ss_min', 'b'))
    def test_hybrid_ss_min(self, callMock):
        self.func_args['prediction_method'] = 'fq-sub'
        a, b, c = repredict_structures_for_homol_seqs(**self.func_args)
        self.assertIsNone(a)
        self.assertIsNone(b)
        self.assertEqual(c, ["expected_error_hybrid_ss_min"])

    def test_wrong_pm(self):
        self.func_args['prediction_method'] = 'blabla'
        with self.assertRaises(AssertionError):
            a, b, c = repredict_structures_for_homol_seqs(**self.func_args)


class TestFailExceptions(unittest.TestCase):
    def setUp(self):
        f = open(ref_json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        self.data = bb

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t2')
        os.close(ff)
        self.html = blast_query + 'test_html.html'

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t3')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, allHits_fasta = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        with os.fdopen(ff, 'w') as f:
            SeqIO.write([i.extension for i in self.data.hits], f, 'fasta')

        self.fasta_structures = allHits_fasta

        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp=r'(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            html=self.html,
        )

        self.func_args = {
            'query': self.data.query,
            'seqs2predict_fasta': allHits_fasta,
            'pred_method_params': {},
            'all_hits_list': [i.extension for i in self.data.hits],
            'seqs2predict_list': [i.extension for i in self.data.hits],
            'use_cm_file': 'abc',
        }
        remove_files_with_try([blast_in + '.tmp_rboAnalyzer'])

    def tearDown(self):
        files = glob.glob(blast_in + '.r-*')
        remove_files_with_try(
            [
                self.csv,
                self.html,
                self.json,
                self.fasta,
                self.fasta_structures
            ] + files,
            ''
        )

    def test_exit_on_missing_rfam(self):
        op = os.path.join(test_data_dir, 'config_rfam_here.txt')
        with open(op, 'w') as f:
            f.write(
                "[TOOL_PATHS]\n[DATA]\nrfam_dir = test_func/test_data\n"
            )
        self.args.config_file = op
        self.args.download_rfam = False

        with self.assertRaises(SystemExit):
            lunch_with_args(lunch_with_args(self.args))

        remove_files_with_try([op])

    def test_exit_on_unqual_blast_simple(self):
        op = os.path.join(test_data_dir, 'two_fasta.fa')
        with open(op, 'w') as f:
            f.write(
                ">a\nCAGCATGCTAGCTGATGCTA\n>b\nAGCTGATCGTAGCTGCTAGTCGTA\n"
            )
        self.args.blast_query = op
        self.args.mode = 'simple'
        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op])

    def test_exit_on_unqual_blast_locarna(self):
        op = os.path.join(test_data_dir, 'two_fasta.fa')
        with open(op, 'w') as f:
            f.write(
                ">a\nCAGCATGCTAGCTGATGCTA\n>b\nAGCTGATCGTAGCTGCTAGTCGTA\n"
            )
        self.args.blast_query = op
        self.args.mode = 'locarna'
        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op])

    def test_exit_on_unqual_blast_meta(self):
        op = os.path.join(test_data_dir, 'two_fasta.fa')
        with open(op, 'w') as f:
            f.write(
                ">a\nCAGCATGCTAGCTGATGCTA\n>b\nAGCTGATCGTAGCTGCTAGTCGTA\n"
            )
        self.args.blast_query = op
        self.args.mode = 'meta'
        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op])

    def test_wrong_blast_simple(self):
        op = os.path.join(test_data_dir, 'wrong_blast.txt')
        with open(self.args.blast_in, 'r') as f, open(op, 'w') as o:
            blast_raw = f.read()
            blast_altered = blast_raw.replace(
                'CGGAGGGGGAACA-CCCGGTCCCATTCCGAACCC',
                'CGGAGGGGGAACA-CCCGGTKYLATTCCGAACCC'
            )
            if blast_raw == blast_altered:
                assert False
            o.write(blast_altered)

        self.args.blast_in = op
        self.args.mode = 'simple'

        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op] + glob.glob(op + '.r-*'))

    def test_wrong_blast_locarna(self):
        op = os.path.join(test_data_dir, 'wrong_blast.txt')
        with open(self.args.blast_in, 'r') as f, open(op, 'w') as o:
            blast_raw = f.read()
            blast_altered = blast_raw.replace(
                'CGGAGGGGGAACA-CCCGGTCCCATTCCGAACCC',
                'CGGAGGGGGAACA-CCCGGTKYLATTCCGAACCC'
            )
            if blast_raw == blast_altered:
                assert False
            o.write(blast_altered)

        self.args.blast_in = op
        self.args.mode = 'locarna'

        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op] + glob.glob(op + '.r-*'))

    def test_wrong_blast_meta(self):
        op = os.path.join(test_data_dir, 'wrong_blast.txt')
        with open(self.args.blast_in, 'r') as f, open(op, 'w') as o:
            blast_raw = f.read()
            blast_altered = blast_raw.replace(
                'CGGAGGGGGAACA-CCCGGTCCCATTCCGAACCC',
                'CGGAGGGGGAACA-CCCGGTKYLATTCCGAACCC'
            )
            if blast_raw == blast_altered:
                assert False
            o.write(blast_altered)

        self.args.blast_in = op
        self.args.mode = 'meta'

        with self.assertRaises(SystemExit):
            lunch_with_args(self.args)
        remove_files_with_try([op] + glob.glob(op + '.r-*'))


class TestExtensionError(unittest.TestCase):
    def setUp(self):
        f = open(ref_json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        self.data = bb

        ff, csv = tempfile.mkstemp(prefix='rba_', suffix='_t1')
        os.close(ff)
        self.csv = csv

        ff, html = tempfile.mkstemp(prefix='rba_', suffix='_t2')
        os.close(ff)
        self.html = blast_query + 'test_html.html'

        ff, pandas_dump = tempfile.mkstemp(prefix='rba_', suffix='_t3')
        os.close(ff)
        self.pandas_dump = pandas_dump

        ff, json_file = tempfile.mkstemp(prefix='rba_', suffix='_t4')
        os.close(ff)
        self.json = json_file

        ff, fasta = tempfile.mkstemp(prefix='rba_', suffix='_t5')
        os.close(ff)
        self.fasta = fasta

        ff, allHits_fasta = tempfile.mkstemp(prefix='rba_', suffix='_t6')
        with os.fdopen(ff, 'w') as f:
            SeqIO.write([i.extension for i in self.data.hits], f, 'fasta')

        self.fasta_structures = allHits_fasta

        self.args = Pseudoargs(
            blast_query,
            blast_in,
            blast_db,
            b_type='plain',
            prediction_method=['rnafold'],
            blast_regexp=r'(?<=\|)[A-Z0-9]*\.?\d*$',
            enable_overwrite=True,
            html=self.html,
        )

        self.func_args = {
            'query': self.data.query,
            'seqs2predict_fasta': allHits_fasta,
            'pred_method_params': {},
            'all_hits_list': [i.extension for i in self.data.hits],
            'seqs2predict_list': [i.extension for i in self.data.hits],
            'use_cm_file': 'abc',
        }

    def tearDown(self):
        files = glob.glob(os.path.join(fwd, test_data_dir, 'RF00001_short.blastout') + '.r-*')
        remove_files_with_try(
            [
                self.csv,
                self.html,
                self.json,
                self.fasta,
                self.fasta_structures
            ] + files
        )

    @mock.patch("rna_blast_analyze.BR_core.expand_by_LOCARNA.run_locarna", side_effect=exceptions.LocarnaException('expected_error_locarna', 'b'))
    def test_locarna_complete_fail(self, callMock):
        self.args.mode = 'locarna'
        with self.assertRaises(SystemExit):
            out = lunch_with_args(self.args)

    @mock.patch("rna_blast_analyze.BR_core.expand_by_LOCARNA.run_locarna", side_effect=exceptions.LocarnaException('expected_error_locarna', 'b'))
    def test_meta_complete_fail(self, callMock):
        self.args.mode = 'locarna'
        with self.assertRaises(SystemExit):
            out = lunch_with_args(self.args)

    def test_locarna_partial_fail(self):
        bb = self.data
        hit = bb.hits.pop(1)
        hit.extension = None
        bb.hits_failed.append(hit)

        fda, all_hits_fasta = tempfile.mkstemp(prefix='rba_', suffix='_22')
        os.close(fda)
        bb.write_results_fasta(all_hits_fasta)

        ih_model, analyzed_hits = find_and_extract_cm_model(bb.args, bb)

        # this part predicts homology - it is not truly part of repredict
        homology_prediction, homol_seqs, cm_file_rfam_user = infer_homology(
            analyzed_hits=bb, args=bb.args, cm_model_file=ih_model, multi_query=False, iteration=0
        )

        self.assertIsNone(cm_file_rfam_user)
        self.assertEqual(len(analyzed_hits.hits), len(homology_prediction))

        remove_files_with_try([
            ih_model,
            all_hits_fasta
        ])

    def test_locarna_with_one_seq_fail(self):
        op = os.path.join(test_data_dir, 'wrong_blast.txt')
        with open(self.args.blast_in, 'r') as f, open(op, 'w') as o:
            blast_raw = f.read()
            blast_altered = re.sub(
                'AGCGGAGGGGAAACCGCCCGGTCCCATTCCGAACCCGGAAGC',
                'AGCGGAGGGGAAACCGCAGTCGATGTTTCCGAACCCCAGTGC',
                blast_raw,
                1
            )
            if blast_raw == blast_altered:
                assert False
            o.write(blast_altered)

        self.args.blast_in = op
        self.args.mode = 'locarna'
        self.args.threads = 4

        out = lunch_with_args(self.args)

        self.assertEqual(len(out[0].hits_failed), 1)

        remove_files_with_try([op] + glob.glob(op + '.r-*'))

    def test_meta_with_one_seq_fail(self):
        op = os.path.join(test_data_dir, 'wrong_blast.txt')
        with open(self.args.blast_in, 'r') as f, open(op, 'w') as o:
            blast_raw = f.read()
            blast_altered = re.sub(
                'AGCGGAGGGGAAACCGCCCGGTCCCATTCCGAACCCGGAAGC',
                'AGCGGAGGGGAAACCGCAGTCGATGTTTCCGAACCCACGTCG',
                blast_raw,
                1
            )
            if blast_raw == blast_altered:
                assert False
            o.write(blast_altered)

        self.args.blast_in = op
        self.args.mode = 'meta'
        self.args.threads = 4

        out = lunch_with_args(self.args)

        self.assertEqual(len(out[0].hits_failed), 0)

        remove_files_with_try([op] + glob.glob(op + '.r-*'))
