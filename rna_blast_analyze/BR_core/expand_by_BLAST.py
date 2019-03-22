import os
import pickle
import re
from copy import deepcopy
from random import shuffle
from tempfile import mkstemp
import itertools
import logging

from Bio import SeqIO

import rna_blast_analyze.BR_core.BA_support as BA_support
import rna_blast_analyze.BR_core.extend_hits
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.stockholm_alig import StockholmFeatureStock
from rna_blast_analyze.BR_core import BA_verify
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.cmalign import run_cmfetch, get_cm_model, RfamInfo
from rna_blast_analyze.BR_core.validate_args import validate_args

ml = logging.getLogger(__name__)


def blast_wrapper(args_inner, shared_list=None):
    ml.debug(fname())
    ret_line, _ = blast_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def trim_before(seqs):
    nseqs = []
    for seq in seqs:
        ann = seq.annotations
        tss = ann['trimmed_ss']
        tse = ann['trimmed_se']
        tes = ann['trimmed_es']
        tee = ann['trimmed_ee']

        if not tss and not tse and not tes and not tee:
            # |-----|----------|---|-------------|-----|
            s = ann['extended_start'] - ann['super_start']
            e = ann['extended_end'] - ann['super_end']
        elif tss and not tse and not tes and not tee:
            # |---:--|----------|---|-------------|-----|
            s = ann['extended_start'] - 1
            e = ann['extended_end']
        elif not tss and tse and not tes and not tee:
            # |-----|----------|---|-------------|---:--|
            s = ann['extended_start'] - ann['super_start']
            e = ann['extended_end'] - ann['super_start'] + 1
        elif tss and tse and not tes and not tee:
            # |--:---|----------|---|-------------|---:--|
            if ann['strand'] == -1:
                s = len(seq.seq) - ann['extended_end']
                e = len(seq.seq) - ann['extended_start'] + 1
            else:
                s = ann['extended_start'] - 1
                e = ann['extended_end']
        elif tss and not tse and tes and not tee:
            # |-----|----:------|---|-------------|-----|
            s = 0
            e = ann['extended_end']
        elif not tss and tse and not tes and tee:
            # |-----|----------|---|-------:------|-----|
            s = ann['extended_start'] - ann['super_start']
            e = ann['extended_end'] - ann['super_start'] + 1
        elif tss and tse and tes and tee:
            # |-----|----:------|---|------:-------|-----|
            s = 0
            e = len(seq.seq)
        elif tss and tse and tee and not tes:
            # |--:---|----------|---|------:-------|-----|
            s = ann['extended_start'] - 1
            e = len(seq.seq[s:])
        elif tss and tes and tse and not tee:
            # |-----|-----:-----|---|-------------|--:---|
            s = 0
            e = ann['extended_end']
        else:
            raise NotImplementedError('Unexpected combination when trimming extended sequence.')

        nseq = seq[s:e]
        nseq.annotations = ann
        nseqs.append(nseq)
    return nseqs


def compute_true_location_se(hit, ql):
    hh = hit.extension
    ann = hh.annotations
    if ann['strand'] == 1:
        if ann['trimmed_es']:
            s = 1
        else:
            s = ann['blast'][1].sbjct_start - ann['blast'][1].query_start + 1
        e = s + len(hh.seq) - 1
    else:
        if ann['trimmed_es']:
            s = 1
        else:
            s = ann['blast'][1].sbjct_end - (ql - ann['blast'][1].query_end)
        e = s + len(hh.seq) - 1
    return s, e


def blast_wrapper_inner(args_inner, shared_list=None):
    ml.debug(fname())
    if not shared_list:
        shared_list = []

    if not validate_args(args_inner):
        print("There was an error with provided arguments. Please see the error message.")
        exit(1)

    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

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
        ml.info('processing query: {}'.format(query.id))
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
        elif args_inner.db_type in ["fasta", "gb", "server"]:
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

        shorts_expanded = trim_before(shorts_expanded)

        shorts_expanded = BA_support.rc_hits_2_rna(shorts_expanded)

        analyzed_hits = BlastSearchRecompute()
        analyzed_hits.args = args_inner
        all_analyzed.append(analyzed_hits)

        query_seq = query.seq.transcribe()

        analyzed_hits.query = query

        def create_blast_only_report_object(exp_hit, query_len):
            # init new Subsequences object
            #  here the object source shadows the final hit
            hit = BA_support.Subsequences(exp_hit)

            # init new SeqRecord object
            ns = deepcopy(exp_hit)
            ann = ns.annotations
            tss = ann['trimmed_ss']
            tse = ann['trimmed_se']
            tes = ann['trimmed_es']
            tee = ann['trimmed_ee']

            ns.letter_annotations['ss0'] = '.' * len(ns.seq)
            ns.annotations['sss'] = ['ss0']
            ns.description = ''

            hit.extension = ns

            pos_match = re.search(str(ns.seq), str(ns.seq), flags=re.IGNORECASE)
            if not pos_match:
                raise Exception('Subsequnce not found in supersequence.')

            bl = ns.annotations['blast'][1]

            if bl.sbjct_start < bl.sbjct_end:
                bls = bl.sbjct_start - bl.query_start + 1
                ble = bl.sbjct_end + (query_len - bl.query_end)
            elif bl.sbjct_end < bl.sbjct_start:
                bls = bl.sbjct_end - (query_len - bl.query_end)
                ble = bl.sbjct_start + bl.query_start - 1
            else:
                raise Exception('Unknown strand option.')

            # if whole subject sequence too short, this assertion will fail
            if tss or tse or tes or tee:
                ml.warning('Skipping check ({}) - subject sequence too short.'.format(ns.id))
            else:
                assert len(ns.seq) == abs(bls - ble) + 1
                assert bls == ns.annotations['extended_start']
                assert ble == ns.annotations['extended_end']

            hit.best_start, hit.best_end = compute_true_location_se(hit, query_len)

            return hit

        # blast only extension
        for exp_hit in shorts_expanded:
            analyzed_hits.hits.append(create_blast_only_report_object(exp_hit, len(query_seq)))

        # assign Locarna score to None as it is not directly accessible from mlocarna
        for hit in analyzed_hits.hits:
            hit.extension.annotations['score'] = None

        # infer homology
        # consider adding query sequence to alignment and scoring it as a proof against squed alignments and datasets
        # write all hits to fasta
        fda, all_hits_fasta = mkstemp(prefix='rba_', suffix='_17', dir=CONFIG.tmpdir)
        with os.fdopen(fda, 'w') as fah:
            # analyzed_hits.write_results_fasta(fah)
            for hit in analyzed_hits.hits:
                if len(hit.extension.seq) == 0:
                    continue
                fah.write('>{}\n{}\n'.format(
                    hit.extension.id,
                    str(hit.extension.seq))
                )

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
                method = method
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
                used_cm_file=cm_file
            )
            out_line.append(to_tab_delim_line_simple(args_inner))

        ml_out_line.append('\n'.join(out_line))

        if cm_file is not None and os.path.exists(cm_file):
            os.remove(cm_file)

        os.remove(all_hits_fasta)

    return '\n'.join(ml_out_line), all_analyzed
