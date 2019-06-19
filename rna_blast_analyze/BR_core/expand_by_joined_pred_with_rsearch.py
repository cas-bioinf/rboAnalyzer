import os
import pickle
import re
from copy import deepcopy
from tempfile import mkstemp
import logging

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
from rna_blast_analyze.BR_core.expand_by_BLAST import blast_wrapper_inner
from rna_blast_analyze.BR_core.expand_by_LOCARNA import locarna_anchored_wrapper_inner
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.stockholm_alig import StockholmFeatureStock
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.cmalign import get_cm_model, run_cmfetch, RfamInfo
from rna_blast_analyze.BR_core.validate_args import validate_args

ml = logging.getLogger(__name__)


def joined_wrapper(args_inner, shared_list=None):
    ml.debug(fname())
    ret_line, _ = joined_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def joined_wrapper_inner(args_inner, shared_list=None):
    ml.debug(fname())
    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    if not validate_args(args_inner):
        print("There was an error with provided arguments. Please see the error message.")
        exit(1)

    blast_args = deepcopy(args_inner)
    locarna_args = deepcopy(args_inner)

    if args_inner.repredict_file is None:
        fd, repred_file = mkstemp(prefix='rba_', suffix='_18', dir=CONFIG.tmpdir)
        os.close(fd)
    else:
        repred_file = args_inner.repredict_file

    for i, args in enumerate([blast_args, locarna_args]):
        args.prediction_method = []
        args.pred_params = dict()
        args.dump = None
        args.dill = None
        args.pdf_out = None
        args.pandas_dump = None
        args.repredict_file = repred_file + str(i)
        args.dev_pred = False
        args.logfile = None
        args.json = None
        args.html = None

    _, simple_seq = blast_wrapper_inner(blast_args)
    _, locarna_seq = locarna_anchored_wrapper_inner(locarna_args)

    if not shared_list:
        shared_list = []

    stockholm_features = StockholmFeatureStock()
    stockholm_features.add_custom_parser_tags('GC', {'cA1': 'anchor letter tag',
                                                     'cA2': 'anchor number tag'})

    p_blast = BA_support.blast_in(args_inner.blast_in, b=args_inner.b_type)
    # this is done for each query

    assert len(p_blast) == len(simple_seq) == len(locarna_seq)

    if len(p_blast) > 1:
        multi_query = True
    else:
        multi_query = False

    ml_out_line = []
    all_analyzed = []
    for (iteration, bhp), b_hits, l_hits in zip(enumerate(p_blast), simple_seq, locarna_seq):
        # now, the blast, muscle and locarna are precomputed
        # only thing needed is to make a CMmodel based on rsearch and evaluate which expansion is the best one
        # maybe even not needed, as the homology haw been already computed when running hmology inference with rsearch
        # it suffice to compare values passed from cmsearch
        # but the model is uncalibrated
        # update -> use only blast and locarna
        #  also combine repredict file here

        analyzed_hits = BlastSearchRecompute()
        analyzed_hits.args = args_inner

        # add query from simple extension to analyzed hits of joined pred (need cm bit score for html out)
        analyzed_hits.query = query = b_hits.query
        analyzed_hits.best_matching_model = b_hits.best_matching_model

        all_analyzed.append(analyzed_hits)
        order_out = []
        for bh, lh in zip(b_hits.hits, l_hits.hits):
            hits = [bh, lh]

            bit_scores = [i.extension.annotations['cmstat']['bit_sc'] for i in hits]

            mb = max(bit_scores)
            bit_index = [i for i, j in enumerate(bit_scores) if j == mb][0]
            order_out.append(bit_index)

            analyzed_hits.hits.append(hits[bit_index])

        # build the repredict file here if needed
        if args_inner.repredict_file:
            b_repredict = BA_support.iter2file_name(blast_args.repredict_file, multi_query, iteration)
            l_repredict = BA_support.iter2file_name(blast_args.repredict_file, multi_query, iteration)
            o_repredict = BA_support.iter2file_name(args_inner.repredict_file, multi_query, iteration)
            with open(b_repredict, 'r') as barf, open(l_repredict, 'r') as larf, open(o_repredict, 'w') as reprf:
                """
                please note that order of files to merge must be same as the order of methods in previous for cycle
                ie same as the one in which order_out var is set
                """
                bb = (barf, larf)

                fl = bb[0].readline()
                reprf.write(fl)
                fl = bb[0].readline()
                reprf.write(fl)
                # dump first line of the other documents
                [[i.readline() for _ in range(1)] for i in bb[1:]]

                for o in order_out:
                    lll = [i.readline() for i in bb]
                    reprf.write(lll[o])

        # recreate needed data from selected hits
        homology_prediction = []
        homol_seqs = []
        for hit in analyzed_hits.hits:
            homology_prediction.append(hit.hpred)
            if hit.hpred:
                homol_seqs.append(
                    hit.extension
                )

            # add default prediction if it is not present
            if 'ss0' not in hit.extension.letter_annotations:
                if 'sss' not in hit.extension.annotations:
                    hit.extension.anotations['sss'] = []
                hit.extension.annotations['sss'] += ['ss0']
                hit.extension.letter_annotations['ss0'] = '.'*len(hit.extension.seq)

        fda, all_hits_fasta = mkstemp(prefix='rba_', suffix='_19', dir=CONFIG.tmpdir)
        with os.fdopen(fda, 'w') as fah:
            # analyzed_hits.write_results_fasta(fah)
            for hit in analyzed_hits.hits:
                if len(hit.extension.seq) == 0:
                    continue
                fah.write('>{}\n{}\n'.format(
                    hit.extension.id,
                    str(hit.extension.seq))
                )

        # remove description from hits and sources
        for hit in analyzed_hits.hits:
            hit.extension.description = ''

        if args_inner.cm_file:
            cm_file = args_inner.cm_file
        else:
            cm_file = None

        out_line = []
        # multiple prediction params
        if args_inner.dev_pred:
            dp_list = []
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
                        pred_method=method,
                        method_params=method_params,
                        used_cm_file=cm_file,
                        multi_query=multi_query,
                        iteration=iteration,
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
                used_cm_file=cm_file,
                multi_query=multi_query,
                iteration=iteration,
            )
            out_line.append(to_tab_delim_line_simple(args_inner))

        ml_out_line.append('\n'.join(out_line))

        if cm_file is not None and os.path.exists(cm_file) and args_inner.cm_file is None:
            os.remove(cm_file)

        # remove repredict files for each method
        for i in range(2):
            rpi = BA_support.iter2file_name(repred_file + str(i), multi_query, iteration)
            os.remove(rpi)

        os.remove(all_hits_fasta)

    if args_inner.repredict_file is None:
        os.remove(repred_file)

    return '\n'.join(ml_out_line), all_analyzed
