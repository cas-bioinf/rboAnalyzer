import os
import pickle
import re
from copy import deepcopy
from multiprocessing import Pool
from random import shuffle
from subprocess import call
from tempfile import mkstemp
import itertools
import logging
import shlex

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.BA_support as BA_support
import rna_blast_analyze.BR_core.extend_hits
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.alifold4all import compute_refold
from rna_blast_analyze.BR_core.config import CONFIG, tools_paths
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core.locarna_clustal_like_2stockholm import parse_locarna_alignment
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.stockholm_parser import StockholmFeatureStock, trim_alignment_by_sequence
from rna_blast_analyze.BR_core import BA_verify
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.cmalign import get_cm_model, run_cmfetch, RfamInfo
from rna_blast_analyze.BR_core.validate_args import validate_args

ml = logging.getLogger(__name__)


def write_clustal_like_file_with_anchors(fid, seq_name, cseq, my_anchors):
    fid.write('CLUSTAL W --- Clustal like format for locarna\n')
    fid.write('{} {}\n'.format(seq_name, cseq))
    for anch in my_anchors:
        fid.write('{} {}\n'.format(anch[0], anch[1]))


def locarna_worker(pack):
    ml.debug(fname())
    one_expanded_hit, query_seq, stockholm_features, locarna_params, anchor_length = pack
    # read the aligned segment and use it as anchors for locarna
    # run locarna in local mode and put the query sequence with the extended sequence with the blast aligned
    # segment as anchor
    blast_entry = one_expanded_hit.annotations['blast'][1]

    anchors = LocarnaAnchor(
        to_rna(blast_entry.query),
        blast_entry.match,
        to_rna(blast_entry.sbjct),
        anchor_length=anchor_length
    )
    # extracted temp is my query

    # access the locarna aligner directly
    # CARNA for some reason wount accept its own designed format, but, it will eat a pseudo-clustal
    fd1, locarna_file1 = mkstemp(prefix='rba_', suffix='_20', dir=CONFIG.tmpdir)
    with os.fdopen(fd1, 'w') as fp_locarna_file_1:
        ql1, ql2 = anchors.anchor_whole_seq(str(query_seq), 'query')
        write_clustal_like_file_with_anchors(fp_locarna_file_1,
                                             'query',
                                             str(query_seq),
                                             (
                                                 ('#A1', ql1.split()[0]),
                                                 ('#A2', ql2.split()[0])
                                             ))

    fd2, locarna_file2 = mkstemp(prefix='rba_', suffix='_21', dir=CONFIG.tmpdir)
    with os.fdopen(fd2, 'w') as fp_locarna_file_2:
        sl1, sl2 = anchors.anchor_whole_seq(str(one_expanded_hit.seq), 'subject')
        write_clustal_like_file_with_anchors(fp_locarna_file_2,
                                             one_expanded_hit.id,
                                             str(one_expanded_hit.seq),
                                             (
                                                 ('#A1', sl1.split()[0]),
                                                 ('#A2', sl2.split()[0])
                                             ))

    loc_out_file = run_locarna(
        locarna_file1,
        locarna_file2,
        locarna_params
    )

    # debug only
    # call('cat {} >> /tmp/all_locarna_files.txt'.format(loc_out_file), shell=True)

    # read locarna alignment
    with open(loc_out_file, 'r') as f:
        locarna_alig = parse_locarna_alignment(f)

    if len(locarna_alig) != 2:
        raise Exception('locarna alignment must be length 2')

    loc_rep = create_report_object_from_locarna(one_expanded_hit, locarna_alig)

    os.remove(locarna_file1)
    os.remove(locarna_file2)
    os.remove(loc_out_file)
    return loc_rep


def locarna_anchored_wrapper(args_inner, shared_list=None):
    ml.debug(fname())
    ret_line, hits = locarna_anchored_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def locarna_anchored_wrapper_inner(args_inner, shared_list=None):
    ml.debug(fname())
    if not shared_list:
        shared_list = []

    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    if not validate_args(args_inner):
        print("There was an error with provided arguments. Please see the error message.")
        exit(1)

    stockholm_features = StockholmFeatureStock()
    stockholm_features.add_custom_parser_tags('GC', {'cA1': 'anchor letter tag',
                                                     'cA2': 'anchor number tag'})

    p_blast = BA_support.blast_in(args_inner.blast_in, b=args_inner.b_type)
    # this is done for each query

    ml_out_line = []
    all_analyzed = []
    for iteration, (bhp, query) in enumerate(
            itertools.zip_longest(p_blast, SeqIO.parse(args_inner.blast_query, 'fasta'))
    ):
        print('processing query: {}'.format(query.id))
        # check query and blast
        if bhp is None:
            raise ValueError('There is more query sequences in provided fastafile then there are BLAST outputs.')
        if query is None:
            raise ValueError('There is more BLAST outputs then query sequences provided.')

        BA_verify.verify_query_blast(blast=bhp, query=query)

        # select all
        all_blast_hits = BA_support.blast_hsps2list(bhp)

        # filter if needed
        if args_inner.filter_by_eval is not None:
            tmp = filter_by_eval(all_blast_hits, *args_inner.filter_by_eval)
            if len(tmp) == 0 and len(all_blast_hits) != 0:
                ml.error('The requested filter removed all BLAST hits. Nothing to do.')
                exit(0)
        elif args_inner.filter_by_bitscore is not None:
            tmp = filter_by_bits(all_blast_hits, *args_inner.filter_by_bitscore)
            if len(tmp) == 0 and len(all_blast_hits) != 0:
                ml.error('The requested filter removed all BLAST hits. Nothing to do.')
                exit(0)

        all_short = all_blast_hits

        # expand hits according to query + 10 nucleotides +-
        if args_inner.db_type == "blastdb":
            shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits(
                all_short,
                args_inner.blast_db,
                bhp.query_length,
                extra=args_inner.subseq_window_simple_ext,
                blast_regexp=args_inner.blast_regexp,
                skip_missing=args_inner.skip_missing,
                msgs=args_inner.logmsgs,
            )
        elif args_inner.db_type in ["fasta", "gb"]:
            shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits_from_fasta(
                all_short,
                args_inner.blast_db,
                bhp.query_length,
                extra=args_inner.subseq_window_simple_ext,
                blast_regexp=args_inner.blast_regexp,
                skip_missing=args_inner.skip_missing,
                msgs=args_inner.logmsgs,
                format=args_inner.db_type,
            )
        else:
            raise ValueError

        # check, if blast hits are non - overlapping, if so, add the overlapping hit info to the longer hit
        # reflect this in user output
        # shorts_expanded = merge_blast_hits(shorts_expanded)

        shorts_expanded = BA_support.rc_hits_2_rna(shorts_expanded)

        analyzed_hits = BlastSearchRecompute()
        analyzed_hits.args = args_inner
        all_analyzed.append(analyzed_hits)

        query_seq = query.seq.transcribe()

        # compute alignment here
        analyzed_hits.query = query

        if args_inner.threads == 1:
            result = []
            for oeh in shorts_expanded:
                result.append(
                    locarna_worker(
                        (
                            oeh,
                            query_seq,
                            stockholm_features,
                            args_inner.locarna_params,
                            args_inner.locarna_anchor_length
                        )
                    )
                )
        else:
            pack = []
            for oeh in shorts_expanded:
                pack.append(
                    (
                        oeh,
                        query_seq,
                        stockholm_features,
                        args_inner.locarna_params,
                        args_inner.locarna_anchor_length
                    )
                )
            pool = Pool(processes=args_inner.threads)
            result = pool.map(locarna_worker, pack)
            pool.close()

        for res in result:
            analyzed_hits.hits.append(res)

        # infer homology
        # consider adding query sequence to alignment and scoring it as a proof against squed alignments and datasets
        # write all hits to fasta
        fda, all_hits_fasta = mkstemp(prefix='rba_', suffix='_22', dir=CONFIG.tmpdir)
        os.close(fda)
        analyzed_hits.write_results_fasta(all_hits_fasta)

        # this part predicts homology - it is not truly part of repredict
        homology_prediction, homol_seqs, cm_file = infer_homology(analyzed_hits=analyzed_hits, args=args_inner)

        # add homology prediction to the data
        for hit, pred in zip(analyzed_hits.hits, homology_prediction):
            hit.hpred = pred

        out_line = []
        # multiple prediction params
        if args_inner.dev_pred:
            dp_list = []
            # acomodate more dev pred outputs
            dpfile = None
            if getattr(args_inner, 'dump', False):
                dpfile = args_inner.dump.strip('dump')
            if getattr(args_inner, 'pandas_dump', False):
                dpfile = args_inner.pandas_dump.strip('pandas_dump')
            if getattr(args_inner, 'json', False):
                dpfile = args_inner.json.strip('json')

            # optimization so the rfam cm file is used only once
            if cm_file is None and 'rfam' in ''.join(args_inner.prediction_method):
                best_model = get_cm_model(args_inner.blast_query, threads=args_inner.threads)
                rfam = RfamInfo()
                cm_file = run_cmfetch(rfam.file_path, best_model)

            for method in args_inner.prediction_method:
                # cycle the prediction method settings
                # get set of params for each preditcion
                selected_pred_params = [kk for kk in args_inner.pred_params if method in kk]
                shuffle(selected_pred_params)
                # for method_params in args_inner.pred_params:
                for i, method_params in enumerate(selected_pred_params):
                    ah = deepcopy(analyzed_hits)

                    random_flag = BA_support.generate_random_name(8, shared_list)
                    shared_list.append(random_flag)

                    pname = re.sub(' ', '', str(method))
                    flag = '|pred_params|' + random_flag

                    # rebuild the args only with actualy used prediction settings
                    ah.args.prediction_method = method
                    ah.args.pred_params = method_params

                    if getattr(args_inner, 'dump', False):
                        spa = args_inner.dump.split('.')
                        ah.args.dump = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                    if getattr(args_inner, 'pandas_dump', False):
                        spa = args_inner.pandas_dump.split('.')
                        ah.args.pandas_dump = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                    if getattr(args_inner, 'dill', False):
                        spa = args_inner.dill.split('.')
                        ah.args.dill = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                    if getattr(args_inner, 'pdf_out', False):
                        spa = args_inner.pdf_out.split('.')
                        ah.args.pdf_out = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                    if getattr(args_inner, 'json', False):
                        spa = args_inner.json.split('.')
                        ah.args.json = '.'.join(spa[:-1]) + flag + '.' + spa[-1]

                    wrapped_ending_with_prediction(
                        args_inner=ah.args,
                        analyzed_hits=ah,
                        all_hits_fasta=all_hits_fasta,
                        query=query,
                        pred_method=method,
                        method_params=method_params,
                        used_cm_file=cm_file,
                    )
                    success = True
                    out_line.append(to_tab_delim_line_simple(ah.args))

                    dp_list.append((i, method_params, success, flag, pname, random_flag, args_inner.pred_params))

            if dpfile is not None:
                with open(dpfile + 'devPredRep', 'wb') as devf:
                    pickle.dump(dp_list, devf)
        else:
            wrapped_ending_with_prediction(
                args_inner=args_inner,
                analyzed_hits=analyzed_hits,
                all_hits_fasta=all_hits_fasta,
                query=query,
                used_cm_file=cm_file,
            )
            out_line.append(to_tab_delim_line_simple(args_inner))

        ml_out_line.append('\n'.join(out_line))

        if cm_file is not None and os.path.exists(cm_file):
            os.remove(cm_file)

        os.remove(all_hits_fasta)

    return '\n'.join(ml_out_line), all_analyzed


def create_report_object_from_locarna(exp_hit, locarna_alig):
    """
    create object which will be appended to BlastSearchRecompute class
    This needs to be Subsequences object

    :param exp_hit:
    :param locarna_alig:
    :return:
    """
    ml.debug(fname())
    # chop alignment by seq
    query_ind = [i for i, j in enumerate(locarna_alig) if j.id == 'query']
    if len(query_ind) != 1:
        raise Exception('multiple hits with id "query" in the alignment')
    trimmed_locarna_alig = trim_alignment_by_sequence(
        locarna_alig,
        str(locarna_alig[query_ind[0]].seq),
        structure_annotation='SS_cons'
    )

    aligned_subsequence = BA_support.select_analyzed_aligned_hit(trimmed_locarna_alig, exp_hit.id)

    # add annotations from exp hit
    aligned_subsequence.annotations = exp_hit.annotations
    aligned_subsequence.name = exp_hit.name

    # also add annotations from locarna, mainly score
    aligned_subsequence.annotations.update(locarna_alig.annotations)

    # get the structure
    # by refold
    refold_structures = refold_stockholm(trimmed_locarna_alig, trimmed_locarna_alig.column_annotations['SS_cons'])

    # select refold structure for my seq
    def _select_refold_structure(refold_structures_, exp_hit_id):
        for seq_refold in refold_structures_:
            # seq_refold can be different due too pipe through CLUSTAL
            if seq_refold.id in exp_hit_id:
                return seq_refold

    seq_refold_structure = _select_refold_structure(refold_structures, exp_hit.id)

    aligned_subsequence.letter_annotations['ss0'] = seq_refold_structure.letter_annotations['ss0']
    aligned_subsequence.annotations['sss'] = ['ss0']

    # prepare seq_record for subsequences
    aligned_subsequence.description = ''
    hit = BA_support.Subsequences(exp_hit)

    hit.extension = aligned_subsequence

    # find the matching sequence
    pos_match = re.search(str(aligned_subsequence.seq), str(exp_hit.seq), flags=re.IGNORECASE)
    if not pos_match:
        raise Exception('Subsequnce not found in supersequence. Terminating.')

    hit.best_start, hit.best_end = compute_true_location_locarna(hit, pos_match)

    return hit


def compute_true_location_locarna(hit, match):
    start, end = match.span()
    if hit.source.annotations['blast'][1].strand == 1:
        if hit.source.annotations['trimmed_ss']:
            s = start + 1
        else:
            s = hit.source.annotations['super_start'] + start

        e = end - start + s + 1
    else:
        if hit.source.annotations['trimmed_ss']:
            s = start + 1
        else:
            s = hit.source.annotations['super_start'] + start

        e = end - start + s - 1

    return s, e


def to_rna(seq):
    rna_seq = re.sub('t', 'u', seq, flags=re.IGNORECASE)
    return rna_seq.upper()


def run_locarna(query_file, subject_file, locarna_params):
    """
    possible settings:
    struct_local : seq_local
    1 : 1   - breaks sequences and joins them !!!! DO NOT USE
    1 : 0   - tries to align globally, end to end - lot of gaps
    0 : 1   - wanted alignment - unaligned ends of input sequences are trimmed in output alignment
    0 : 0   - wanted alignment - returns whole alignment including unaligned ends

    :param query_file:
    :param subject_file:
    :param locarna_params:
    :return:
    """
    ml.debug('running locarna')
    if not os.path.isfile(query_file):
        raise FileNotFoundError('Provided file {} was not found'.format(query_file))
    if not os.path.isfile(subject_file):
        raise FileNotFoundError('Provided file {} was not found'.format(subject_file))

    cmd = '{} {} {} {} > {}'.format(
        shlex.quote('{}locarna'.format(CONFIG.locarna_path)),
        ' '.join([shlex.quote(i) for i in shlex.split(locarna_params)]),
        shlex.quote(query_file),
        shlex.quote(subject_file),
        shlex.quote(subject_file + '.loc_out'),
    )
    ml.debug(cmd)
    r = call(cmd, shell=True)
    if r:
        msgfail = 'call to locarna failed for files in1:{} in2:{} out:{}'.format(
            query_file, subject_file, subject_file + '.loc_out'
        )
        ml.error(msgfail)
        raise ChildProcessError(msgfail)

    return subject_file + '.loc_out'


def write_locarna_anchors_with_min_length(match_line, min_anchor_length=1):
    ml.debug(fname())
    h1 = []
    h2 = []
    pa = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    part_desig = 0

    for match in re.finditer('\|+', match_line, flags=re.IGNORECASE):
        if len(match.group()) < min_anchor_length:
            # skip the iterations below minimum length
            continue

        for i in range(match.span()[0] - len(h1)):
            h1.append('.')
            h2.append('.')

        c = 0
        part_desig += 1
        for _ in match.group():
            if c == 9:
                c = 0
                part_desig += 1
            c += 1

            h1.append(pa[part_desig])
            h2.append(str(c))

    for i in range(len(match_line) - len(h1)):
        h1.append('.')
        h2.append('.')

    anchor_l1 = ''.join(h1)
    anchor_l2 = ''.join(h2)
    return anchor_l1, anchor_l2


def squeeze_locarna_anchors_to_aligned_seq(aligned_seq, anchor_line1, anchor_line2):
    ml.debug(fname())
    out_seq = []
    out_al1 = []
    out_al2 = []
    for seq_pos, al1, al2 in zip(aligned_seq, anchor_line1, anchor_line2):
        if seq_pos == '-':
            continue
        out_seq.append(seq_pos)
        out_al1.append(al1)
        out_al2.append(al2)

    squeezed_seq = ''.join(out_seq)
    squeezed_al1 = ''.join(out_al1)
    squeezed_al2 = ''.join(out_al2)

    return squeezed_seq, squeezed_al1, squeezed_al2


class LocarnaAnchor(object):
    """
    while initiating LocarnaAnchor object U can specify minimal anchor length to be used
    If default (-1) is kept, then minimal anchor length for succesfull usage for locarna is infered
    and the number is returned in anchor_length parameter
    """
    def __init__(self, query, match, subject, anchor_length=-1):
        self.match = match
        self.query = query
        self.subject = subject
        # self.anchor_l1, self.anchor_l2 = write_locarna_anchors(self.match)
        # compute anchor length

        self.anchor_length = anchor_length
        if anchor_length < 0:
            while True:
                self.anchor_l1, self.anchor_l2 = write_locarna_anchors_with_min_length(self.match, self.anchor_length)
                if '[' in self.anchor_l1:
                    self.anchor_length += 1
                else:
                    break
        else:
            self.anchor_l1, self.anchor_l2 = write_locarna_anchors_with_min_length(self.match, self.anchor_length)

        assert len(self.anchor_l1) == len(self.anchor_l2) == len(self.query) == len(self.subject)

        if anchor_length < 0:
            print('infered anchor length {}'.format(self.anchor_length))

        self.squeezed_query, self.q_al1, self.q_al2 = squeeze_locarna_anchors_to_aligned_seq(self.query,
                                                                                             self.anchor_l1,
                                                                                             self.anchor_l2)
        self.squeezed_subject, self.s_al1, self.s_al2 = squeeze_locarna_anchors_to_aligned_seq(self.subject,
                                                                                               self.anchor_l1,
                                                                                               self.anchor_l2)

    def anchor_whole_seq(self, seq, seq_line):
        """
        write anchors to whole sequence
        seq input must be appropriate sequence for the seq_line chosen
        :param seq: str or SeqRecord object
        :param seq_line: str - 'query' | 'subject'
        :return:
        """
        if not isinstance(seq, str) | isinstance(seq, SeqRecord):
            raise TypeError("types 'str' and 'SeqRecord' accepted, but {} given.".format(type(seq)))

        if not (seq_line == 'query') | (seq_line == 'subject'):
            raise AttributeError("Expected 'query' or 'subject', but {} given.".format(seq_line))

        if isinstance(seq, SeqRecord):
            sequence = str(seq.seq)
        else:
            sequence = seq

        # find matching sequence in seq
        if seq_line == 'query':
            match_o = re.search(self.squeezed_query, sequence)
            if not match_o:
                raise Exception('no match for squeeze query with provided sequence')

            start = ''.join(['.' for i in range(match_o.span()[0])])
            end = ''.join(['.' for i in range(match_o.span()[1], len(sequence))])
            l1 = start + self.q_al1 + end + ' #1'
            l2 = start + self.q_al2 + end + ' #2'
        else:
            match_o = re.search(self.squeezed_subject, sequence)
            if not match_o:
                raise Exception('no match for squeeze query with provided sequence')

            start = ''.join(['.' for i in range(match_o.span()[0])])
            end = ''.join(['.' for i in range(match_o.span()[1], len(sequence))])
            l1 = start + self.s_al1 + end + ' #1'
            l2 = start + self.s_al2 + end + ' #2'

        return l1, l2


def refold_stockholm(stockholm_alig, consensus_structure):
    """
    compute refold.pl from Vienna RNA package
    :param stockholm_alig:
    :param consensus_structure:
    :return:
    """
    ml.debug(fname())
    # convert to clustal alignment
    fd, clust_tempfile = mkstemp(prefix='rba_', suffix='_23', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as f:
        stockholm_alig.write_clustal(f)

    # write fake alifold output with given consensus structure
    fd, alif_fake_file = mkstemp(prefix='rba_', suffix='_24', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as f:
        # the consensus sequence in alifold file is really not used for anything
        f.write('A'*len(consensus_structure) + '\n')
        f.write(consensus_structure + '\n')

    # compute refold
    # refold_path = locate_refold()
    refold_constrained_file = compute_refold(clust_tempfile, alif_fake_file)

    parsed_seqs = []
    with open(refold_constrained_file, 'r') as f:
        # read the file
        for seq in BA_support.parse_seq_str(f):
            parsed_seqs.append(seq)

    # cleanup
    os.remove(clust_tempfile)
    os.remove(alif_fake_file)
    os.remove(refold_constrained_file)

    return parsed_seqs
