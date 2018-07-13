# try different approach
# use structure of template as guide and blast alignment as seed
# then gradually add nucleotides and reward addition which is in possible structure
# rather then structure, use probability dot plot for that
# ? how to allow insertions and deletions?
# is it just structure rna alignment ?
# maybe do it over all extended sequences at once?

# its like locarna's anchored mode - try that
# just take blast, mark matched parts and run locarna for each of them
# locarna in anchored mode with exp window 10 provided best results

# todo structure prediction
#  - repredict structures for aligned sequences
#    - alifold + refold of aligned sequences - or other consensus prediction algorithms
#    - custom alignment ad identification of conserved unaligned regions?
#    - predict by shapes cast () - unable to handle large records
#
#  - how to infer homology of record from alignment?
#  - repredict after homology inference
#  - run cmalign after homology inference
#  - rscape analysis
# prediction programs: MCFold, ProbKnot?, RSpredict, RNAwolf, TurboFold

import os
import pickle
import re
from copy import deepcopy
from multiprocessing import Pool
from random import shuffle
from subprocess import call
from tempfile import mkstemp
import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.add_usr_local_bin
import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.alifold4all import compute_refold
from rna_blast_analyze.BR_core.blast_hits_merging import merge_blast_hits
from rna_blast_analyze.BR_core.config import CONFIG, tools_paths
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core.locarna_clustal_like_2stockholm import parse_locarna_alignment
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.stockholm_parser import StockholmFeatureStock, trim_alignment_by_sequence
from rna_blast_analyze.BR_core import BA_verify


def locarna_worker(pack):
    one_expanded_hit, query_seq, stockholm_features, locarna_params, anchor_length = pack
    # read the aligned segment and use it as anchors for locarna
    # run locarna in local mode and put the query sequence with the extended sequence with the blast aligned
    # segment as anchor
    blast_entry = one_expanded_hit.annotations['blast'][1]

    anchors = LocarnaAnchor(to_rna(blast_entry.query),
                            blast_entry.match,
                            to_rna(blast_entry.sbjct),
                            anchor_length=anchor_length)
    # extracted temp is my query

    # access the locarna aligner directly
    # CARNA for some reason wount accept its own designed format, but, it will eat a pseudo-clustal
    def write_clustal_like_file_with_anchors(fid, seq_name, cseq, my_anchors):
        fid.write('CLUSTAL W --- Clustal like format for locarna\n')
        fid.write('{} {}\n'.format(seq_name, cseq))
        for anch in my_anchors:
            fid.write('{} {}\n'.format(anch[0], anch[1]))

    fd1, locarna_file1 = mkstemp()
    with os.fdopen(fd1, 'w') as fp_locarna_file_1:
        ql1, ql2 = anchors.anchor_whole_seq(str(query_seq), 'query')
        write_clustal_like_file_with_anchors(fp_locarna_file_1,
                                             'query',
                                             str(query_seq),
                                             (
                                                 ('#A1', ql1.split()[0]),
                                                 ('#A2', ql2.split()[0])
                                             ))

    fd2, locarna_file2 = mkstemp()
    with os.fdopen(fd2, 'w') as fp_locarna_file_2:
        sl1, sl2 = anchors.anchor_whole_seq(str(one_expanded_hit.seq), 'subject')
        write_clustal_like_file_with_anchors(fp_locarna_file_2,
                                             one_expanded_hit.id,
                                             str(one_expanded_hit.seq),
                                             (
                                                 ('#A1', sl1.split()[0]),
                                                 ('#A2', sl2.split()[0])
                                             ))

    loc_out_file = run_locarna(locarna_file1,
                               locarna_file2,
                               locarna_params)

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
    ret_line, hits = locarna_anchored_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def locarna_anchored_wrapper_inner(args_inner, shared_list=None):
    if not shared_list:
        shared_list = []

    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    stockholm_features = StockholmFeatureStock()
    stockholm_features.add_custom_parser_tags('GC', {'cA1': 'anchor letter tag',
                                                     'cA2': 'anchor number tag'})

    if args_inner.logfile:
        fid = open(args_inner.logfile, 'w')
        fid.write('Program BA runned at: {}\n'.format(BA_support.print_time()))
        BA_support.print_parameters(args_inner, fid)

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
        all_short = BA_support.blast_hsps2list(bhp)
        # expand hits according to query + 10 nucleotides +-
        shorts_expanded, _ = BA_support.expand_hits(
            all_short,
            args_inner.blast_db,
            bhp.query_length,
            extra=args_inner.subseq_window_locarna,
            blast_regexp=args_inner.blast_regexp
        )

        # check, if blast hits are non - overlapping, if so, add the overlapping hit info to the longer hit
        # reflect this in user output
        shorts_expanded = merge_blast_hits(shorts_expanded)

        shorts_expanded = BA_support.rc_hits_2_rna(shorts_expanded)

        analyzed_hits = BlastSearchRecompute()
        analyzed_hits.args = args_inner
        all_analyzed.append(analyzed_hits)

        query_seq = query.seq.transcribe()

        # compute alignment here
        analyzed_hits.query = query
        pack = []
        for oeh in shorts_expanded:
            pack.append((oeh,
                         query_seq,
                         stockholm_features,
                         args_inner.locarna_params,
                         args_inner.locarna_anchor_length))
        pool = Pool(processes=args_inner.threads)
        result = pool.map(locarna_worker, pack)
        pool.close()

        for res in result:
            analyzed_hits.hits.append(res)

        # infer homology
        # consider adding query sequence to alignment and scoring it as a proof against squed alignments and datasets
        # write all hits to fasta
        fda, all_hits_fasta = mkstemp()
        analyzed_hits.write_results_fasta(all_hits_fasta)

        # this part predicts homology - it is not truly part of repredict
        print('infering homology')
        homology_prediction, homol_seqs = infer_homology(analyzed_hits=analyzed_hits, args=args_inner)

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
                    if getattr(args_inner, 'o_tbl', False):
                        spa = args_inner.o_tbl.split('.')
                        ah.args.o_tbl = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
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
                query=query
            )
            out_line.append(to_tab_delim_line_simple(args_inner))

        ml_out_line.append('\n'.join(out_line))

        os.remove(all_hits_fasta)

    if args_inner.logfile:
        fid.write('ended at: {}'.format(BA_support.print_time()))
        fid.close()

    return '\n'.join(ml_out_line), all_analyzed


def create_report_object_from_locarna(exp_hit, locarna_alig):
    """
    create object which will be appended to BlastSearchRecompute class
    This needs to be Subsequences object

    :param exp_hit:
    :param locarna_alig:
    :return:
    """

    # chop alignment by seq
    query_ind = [i for i, j in enumerate(locarna_alig) if j.id == 'query']
    if len(query_ind) != 1:
        raise Exception('multiple hits with id "query" in the alignment')
    trimmed_locarna_alig = trim_alignment_by_sequence(
        locarna_alig,
        str(locarna_alig[query_ind[0]].seq),
        structure_annotation='SS_cons'
    )

    # select analyzed hit, drop query
    # def _select_analyzed_hit(trimmed_locarna_alig):
    #     for seq_in_alig in trimmed_locarna_alig.get_unalined_seqs(keep_letter_ann=True):
    #         if seq_in_alig.id == exp_hit.id:
    #             # aligned_subsequence = exp_hit
    #             aligned_subsequence = seq_in_alig
    #             return aligned_subsequence

    aligned_subsequence = BA_support.select_analyzed_aligned_hit(trimmed_locarna_alig, exp_hit.id)

    # add annotations from exp hit
    aligned_subsequence.annotations = exp_hit.annotations
    aligned_subsequence.name = exp_hit.name

    # also add annotations from locarna, mainly score
    aligned_subsequence.annotations.update(locarna_alig.annotations)

    # get the structure
    # by refold
    # todo check where secondary structure is
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
    hit = BA_support.Subsequences(exp_hit)
    sub_id = aligned_subsequence.id[:-2].split(':')[-1]

    hit.subs[sub_id] = aligned_subsequence

    # find the matching sequence
    pos_match = re.search(str(aligned_subsequence.seq), str(exp_hit.seq), flags=re.IGNORECASE)
    if not pos_match:
        raise Exception('Subsequnce not found in supersequence. Terminating.')

    hit.ret_keys = [sub_id,
                    'ss0',
                    None]
    hit.best_start = hit.source.annotations['super_start'] + pos_match.span()[0]
    hit.best_end = hit.source.annotations['super_start'] + pos_match.span()[1] - 1
    return hit


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
    if not os.path.isfile(query_file):
        raise FileNotFoundError('Provided file {} was not found'.format(query_file))
    if not os.path.isfile(subject_file):
        raise FileNotFoundError('Provided file {} was not found'.format(subject_file))

    r = call(
        '{}locarna {} {} {} > {}'.format(
            CONFIG.locarna_path,
            locarna_params,
            query_file,
            subject_file,
            subject_file + '.loc_out',
        ),
        shell=True
    )
    if r:
        raise ChildProcessError('call to locarna failed for files '
                                'in1:{} in2:{} out:{}'.format(query_file,
                                                              subject_file,
                                                              subject_file + '.loc_out'))

    return subject_file + '.loc_out'


def write_locarna_anchors_with_min_length(match_line, min_anchor_length=1):
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
    # convert to clustal alignment
    fd, clust_tempfile = mkstemp()
    with os.fdopen(fd, 'w') as f:
        stockholm_alig.write_clustal(f)

    # write fake alifold output with given consensus structure
    fd, alif_fake_file = mkstemp()
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
