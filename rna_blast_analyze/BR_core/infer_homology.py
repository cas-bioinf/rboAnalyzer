import os
import sys
from tempfile import mkstemp
from copy import deepcopy
import logging

import rna_blast_analyze.BR_core.BA_support as BA_support
import rna_blast_analyze.BR_core.viennaRNA
from rna_blast_analyze.BR_core.cmalign import run_cmalign_on_fasta, read_cmalign_sfile, run_cmbuild, run_cmfetch, \
    RfamInfo, get_cm_model_table, select_best_matching_model_from_cmscan
from rna_blast_analyze.BR_core.stockholm_alig import StockholmAlig
from rna_blast_analyze.BR_core.stockholm_parser import read_st
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
import shutil
from rna_blast_analyze.BR_core import exceptions

ml = logging.getLogger('rboAnalyzer')


def find_and_extract_cm_model(args, analyzed_hits, rfam=None):
    if rfam is None:
        rfam = RfamInfo()

    ml.info('Infer homology - searching RFAM for best matching model.')
    cmscan_results = get_cm_model_table(args.blast_query, threads=args.threads)
    best_matching_cm_model = select_best_matching_model_from_cmscan(cmscan_results)
    analyzed_hits.best_matching_model = best_matching_cm_model

    if best_matching_cm_model is not None:
        ml.info('Infer homology - best matching RFAM model: {}'.format(analyzed_hits.best_matching_model['target_name']))

    try:
        if args.cm_file:
            # use provided cm file
            ml.info('Infer homology - using provided CM file: {}'.format(args.cm_file))
            fd, cm_model_file = mkstemp(prefix='rba_', suffix='_27', dir=CONFIG.tmpdir)
            os.close(fd)
            ml.debug('Making a copy of provided cm model file to: {}'.format(cm_model_file))
            shutil.copy(args.cm_file, cm_model_file)

        elif args.use_rfam:
            ml.info('Infer homology - using RFAM CM file as reference model.')
            if analyzed_hits.best_matching_model is None:
                ml.error('No RFAM model was matched with score > 0. Nothing to build homology to.')
                sys.exit(1)

            cm_model_file = run_cmfetch(rfam.file_path, analyzed_hits.best_matching_model['target_name'])
        else:
            ml.info('Infer homology - using RSEARCH to build model')
            # default to using RSEARCH
            cm_model_file = build_cm_model_rsearch(analyzed_hits.query, CONFIG.rsearch_ribosum)

        return cm_model_file, analyzed_hits
    except exceptions.SubprocessException as e:
        ml.error("Can't obtain covariance model.")
        ml.error(str(e))
        sys.exit(1)


def infer_homology(analyzed_hits, args, cm_model_file, multi_query=False, iteration=0):
    """
    This is wrapper for infer homology methods. It deals with different options for generating CM.
    :return:
    """
    ml.info('Infering homology...')
    ml.debug(fname())
    bits, eval, loc_score, alig_length = hit_cons_characteristic(analyzed_hits.hits)

    # always run cmscan on rfam for informative reasons
    #  but use inferred CM only if --use_rfam was given
    #  if CM provided, also run inference but use provided file
    # print explanation alongside this information

    # find and extract cm model
    # This code is moved to each extension method to allow fail-fast if model is found in RFAM
    # cm_model_file, analyzed_hits = find_and_extract_cm_model(args, analyzed_hits)

    # include query seq in fasta file to get relevant bit score
    fd_f, fd_fasta = mkstemp(prefix='rba_', suffix='_28', dir=CONFIG.tmpdir)
    with os.fdopen(fd_f, 'w') as f:
        for seq in [analyzed_hits.query] + analyzed_hits.res_2_record_list():
            f.write('>{}\n{}\n'.format(
                seq.id,
                str(seq.seq))
            )

    cm_msa, cm_align_scores = run_cmalign_with_scores(fd_fasta, cm_model_file, threads=args.threads)

    _add_rsearch_align_scores2anal_hits(analyzed_hits, cm_align_scores)

    # remove first 1 (query) from the prediction scores
    prediction = infer_hits_cm(cm_align_scores[1:].bit_sc)

    # write scores to a table, compute it for all data and run some correlation statistics
    if args.repredict_file:
        # note that the first score is for the query and act as a benchmark here
        cm_msa_conservation = alignment_sequence_conservation(cm_msa, gap_chars='-.')

        repredict_file = BA_support.iter2file_name(args.repredict_file, multi_query, iteration)
        with open(repredict_file, 'w') as f:
            _print_table_for_corelation(
                f,
                cm_align_scores.seq_name[1:],
                bits,
                eval,
                loc_score,
                alig_length,
                cm_msa_conservation[1:],
                cm_align_scores.bit_sc[1:],
                cm_msa_conservation[0],
                cm_align_scores.bit_sc[0]
            )

    BA_support.remove_one_file_with_try(fd_fasta)

    selected_hits = [hit.extension for b, hit in zip(prediction, analyzed_hits.hits) if b]

    if args.cm_file or args.use_rfam:
        r_cm_file = cm_model_file
    else:
        r_cm_file = None
        BA_support.remove_one_file_with_try(cm_model_file)

    return prediction, selected_hits, r_cm_file


def run_cmalign_with_scores(fasta_file, cm_file, threads=None):
    fd_sfile, cm_sfile_path = mkstemp(prefix='rba_', suffix='_29', dir=CONFIG.tmpdir)
    os.close(fd_sfile)
    if threads:
        cm_params = '--notrunc --cpu {} --sfile {}'.format(threads, cm_sfile_path)
    else:
        cm_params = '--notrunc --sfile {}'.format(cm_sfile_path)
    cm_msa_file = run_cmalign_on_fasta(fasta_file, cm_file, cmalign_params=cm_params)

    cm_msa = read_st(cm_msa_file)

    # combine the eval and cm_msa_conservation_score
    # the cmalign scores somehow, look into the scoring if those scores are accessible, maybe they are far better then
    # my made up msa_conservation
    # there is - by option --sfile
    cm_align_scores = read_cmalign_sfile(cm_sfile_path)
    # the bit score can be probably directly comparable with blast bit score
    # i can also leverage the fact, that the badly aligned sequences with cmalign have negative bitscore
    # so my score can be
    cm_align_scores.index = range(len(cm_align_scores.index))

    BA_support.remove_files_with_try([cm_sfile_path, cm_msa_file])

    return cm_msa, cm_align_scores


def build_cm_model_rsearch(query_seq, path2selected_sim_array):
    ml.debug(fname())
    query_structure = rna_blast_analyze.BR_core.viennaRNA.RNAfold(str(query_seq.seq))[0]

    # remove any annotations from query:
    qs_clean = deepcopy(query_seq)
    qs_clean.annotations = dict()
    qs_clean.letter_annotations = dict()

    # query_structure = RNA.fold(str(analyzed_hits.query.seq))[0]
    # build stockholm like file for use in cm mohdel build
    st_like = StockholmAlig()
    st_like.append(qs_clean)
    st_like.column_annotations['SS_cons'] = query_structure

    fds, stock_file = mkstemp(prefix='rba_', suffix='_30', dir=CONFIG.tmpdir)
    with os.fdopen(fds, 'w') as f:
        st_like.write_stockholm(f)

    # run actual cmbuild
    cm_model_file = run_cmbuild(stock_file, cmbuild_params='--rsearch {}'.format(
        path2selected_sim_array
    ))

    # cleanup
    BA_support.remove_one_file_with_try(stock_file)
    return cm_model_file


def _add_rsearch_align_scores2anal_hits(ahits, s_table):
    """
    Adds scores from homology inference to analyzed hits
    :param ahits:
    :param s_table:
    :return:
    """
    ml.debug(fname())
    # add result for query to the query
    ahits.query.annotations['cmstat'] = s_table.iloc[0]
    query_len = len(ahits.query)

    for hit, (i, row) in zip(ahits.hits, s_table[1:].iterrows()):
        assert hit.extension.id == row.seq_name
        hit.extension.annotations['cmstat'] = row
        hit.extension.annotations['homology_estimate'] = compute_homology(query_len, row.bit_sc)
    return


def compute_homology(lx, h_bit_sc):
    if h_bit_sc < 0:
        h_estimate = 'Not homologous'
    elif h_bit_sc / lx >= 0.5 and h_bit_sc >= 20:
        h_estimate = 'Homologous'
    else:
        h_estimate = 'Uncertain'
    return h_estimate


def infer_hits_cm(bit_sc, tr=0):
    ml.debug(fname())
    pred = []
    for i in bit_sc:
        if i > tr:
            pred.append(True)
        else:
            pred.append(False)
    return pred


def _print_table_for_corelation(file, names, bits, eval, loc_score, alig_length,
                                cm_msa_conservation, cm_bit_score, cm_msa_default, cm_bit_default):
    ml.debug(fname())
    # write header
    file.write('table for correlation analysis cm_msa_default={} cm_bit_default={}\n'.format(
        cm_msa_default, cm_bit_default
    ))
    # file.write('name bits eval loc_score alig_length msa_conservation cm_msa_conservation cm_bit_score\n')
    file.write('name bits eval loc_score alig_length cm_msa_conservation cm_bit_score\n')
    # for line in zip(names, bits, eval, loc_score, alig_length, msa_conservation, cm_msa_conservation, cm_bit_score):
    for line in zip(names, bits, eval, loc_score, alig_length, cm_msa_conservation, cm_bit_score):
        file.write(' '.join([str(i) for i in line]).strip() + '\n')


def alignment_sequence_conservation(msa, gap_chars='-'):
    """
    compute per column conservation, then compute sum of the conservation per sequence, but count only those columns
    which in aligned sequence are

    consider adding query sequence to alignment and scoring it as a proof against squed alignments and datasets
    :param msa:
    :return:
    """
    ml.debug(fname())
    column_cons = alignment_column_conservation(msa, gap_chars=gap_chars)

    conservation_score_v = []
    for aligned_seq in msa:
        msa_score = 0
        for i, aligned_pose in enumerate(aligned_seq.seq):
            if aligned_pose not in gap_chars:
                msa_score += column_cons[i]

        conservation_score_v.append(msa_score)
    return conservation_score_v


def alignment_column_conservation(msa, gap_chars='-'):
    """
    computes column non-weighted column conservation
    :param msa: multiple alignment iterable over columns and with get_alignment_length() method
    :param gap_chars: characters denoting the gaps in alignment
    :return: list of ints
    """
    ml.debug(fname())
    column_cons = []
    for i in range(msa.get_alignment_length()):
        al_col = msa[:, i]
        total_gaps = 0
        for ch in gap_chars:
            total_gaps += al_col.count(ch)

        column_cons.append(len(al_col) - total_gaps)
    return column_cons


def hit_cons_characteristic(sequence_hits):
    # use blast e-value score and locarna alignment score to select trusted alignments
    # can also use blast hit bits info
    ml.debug(fname())
    bits = []
    eval = []
    align_score = []
    alig_l = []
    for hit in sequence_hits:
        bits.append(hit.source.annotations['blast'][1].bits)
        eval.append(hit.source.annotations['blast'][1].expect)

        # locarna score is not defined for SE and others
        if 'score' in hit.extension.annotations:
            align_score.append(hit.extension.annotations['score'])
        else:
            align_score.append(None)
        alig_l.append(hit.source.annotations['blast'][1].align_length)

    return bits, eval, align_score, alig_l
