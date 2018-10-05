import os
from tempfile import mkstemp
import logging

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.cmalign import run_cmalign_on_fasta, read_cmalign_sfile, run_cmbuild, get_cm_model, run_cmfetch, RfamInfo
from rna_blast_analyze.BR_core.stockholm_alig import StockholmAlig
from rna_blast_analyze.BR_core.stockholm_parser import read_st
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
import shutil

ml = logging.getLogger(__name__)


def infer_homology(analyzed_hits, args):
    """
    This is wrapper for infer homology methods. It deals with different options for generating CM.
    :return:
    """
    ml.info('Infering homology...')
    ml.debug(fname())
    bits, eval, loc_score, alig_length = hit_cons_characteristic(analyzed_hits.hits)

    if args.cm_file:
        # use provided cm file
        ml.info('Infer homology - using provided CM file: {}'.format(args.cm_file))
        fd, cm_model_file = mkstemp()
        os.close(fd)
        ml.debug('Making a copy of provided cm model file to: {}'.format(cm_model_file))
        shutil.copy(args.cm_file, cm_model_file)

    elif args.use_rfam:
        ml.info('Infer homology - using RFAM')
        # find and extract cm model
        rfam = RfamInfo()
        best_matching_model_name = get_cm_model(args.blast_query, threads=args.threads)
        ml.info('Infer homology - best matching model: {}'.format(best_matching_model_name))
        cm_model_file = run_cmfetch(rfam.file_path, best_matching_model_name)
    else:
        ml.info('Infer homology - using RSEARCH to build model')
        # default to using RSEARCH
        # cm_model_file = build_cm_model_rsearch(analyzed_hits.query, args.ribosum)
        cm_model_file = build_cm_model_rsearch(analyzed_hits.query, CONFIG.rsearch_ribosum)

    # do i need to include query in fasta file?
    # yes - to get idea of cm_conservation score value - it does not affect bit score for others
    fd_f, fd_fasta = mkstemp()
    with os.fdopen(fd_f, 'w') as f:
        for seq in [analyzed_hits.query] + analyzed_hits.res_2_record_list():
            f.write('>{}\n{}\n'.format(seq.id,
                                       str(seq.seq)))

    fd_sfile, cm_sfile_path = mkstemp()
    os.close(fd_sfile)
    if args.threads:
        cm_params = '--notrunc --cpu {} --sfile {}'.format(args.threads, cm_sfile_path)
    else:
        cm_params = '--notrunc --sfile {}'.format(cm_sfile_path)
    cm_msa_file = run_cmalign_on_fasta(fd_fasta, cm_model_file, cmalign_params=cm_params)

    cm_msa = read_st(cm_msa_file)
    # note that the first score is for the query and act as a benchmark here
    cm_msa_conservation = alignment_sequence_conservation(cm_msa, gap_chars='-.')

    # combine the eval and cm_msa_conservation_score
    # the cmalign scores somehow, look into the scoring if those scores are accessible, maybe they are far better then
    # my made up msa_conservation
    # there is - by option --sfile
    cm_align_scores = read_cmalign_sfile(cm_sfile_path)
    # the bit score can be probably directly comparable with blast bit score
    # i can also leverage the fact, that the badly aligned sequences with cmalign have negative bitscore
    # so my score can be
    cm_align_scores.index = range(len(cm_align_scores.index))

    _add_rsearch_align_scores2anal_hits(analyzed_hits, cm_align_scores)

    # remove first 1 (query) from the prediction scores
    prediction = _infer_hits_cm(cm_align_scores[1:].bit_sc)

    # write scores to a table, compute it for all data and run some correlation statistics

    if args.repredict_file:
        with open(args.repredict_file, 'w') as f:
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

    os.remove(fd_fasta)
    os.remove(cm_sfile_path)
    os.remove(cm_model_file)
    os.remove(cm_msa_file)

    selected_hits = [hit.subs[hit.ret_keys[0]] for b, hit in zip(prediction, analyzed_hits.hits) if b]

    return prediction, selected_hits


def build_cm_model_rsearch(query_seq, path2selected_sim_array):
    ml.debug(fname())
    query_structure = BA_support.RNAfold(str(query_seq.seq))[1]
    # query_structure = RNA.fold(str(analyzed_hits.query.seq))[0]
    # build stockholm like file for use in cm mohdel build
    st_like = StockholmAlig()
    st_like.append(query_seq)
    st_like.column_annotations['SS_cons'] = query_structure

    fds, stock_file = mkstemp()
    with os.fdopen(fds, 'w') as f:
        st_like.write_stockholm(f)

    # run actual cmbuild
    cm_model_file = run_cmbuild(stock_file, cmbuild_params='--rsearch {}'.format(
        path2selected_sim_array
    ))

    # cleanup
    os.remove(stock_file)
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

    for hit, (i, row) in zip(ahits.hits, s_table[1:].iterrows()):
        assert hit.subs[hit.ret_keys[0]].id == row.seq_name
        hit.subs[hit.ret_keys[0]].annotations['cmstat'] = row
    return


def _infer_hits_cm(bit_sc, tr=-2.03):
    """ current best guess mcc 0.9
    best found threshold = -2.03
    """
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
    column_cons = _alignment_column_conservation(msa, gap_chars=gap_chars)

    # todo continue here
    conservation_score_v = []
    for aligned_seq in msa:
        msa_score = 0
        for i, aligned_pose in enumerate(aligned_seq.seq):
            if aligned_pose not in gap_chars:
                msa_score += column_cons[i]

        conservation_score_v.append(msa_score)
    return conservation_score_v


def _alignment_column_conservation(msa, gap_chars='-'):
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
        if 'score' in hit.subs[hit.ret_keys[0]].annotations:
            align_score.append(hit.subs[hit.ret_keys[0]].annotations['score'])
        else:
            align_score.append(None)
        alig_l.append(hit.source.annotations['blast'][1].align_length)

    return bits, eval, align_score, alig_l
