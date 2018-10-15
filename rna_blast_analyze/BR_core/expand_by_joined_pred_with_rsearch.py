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

ml = logging.getLogger(__name__)


def joined_wrapper(args_inner, shared_list=None):
    ml.debug(fname())
    ret_line, _ = joined_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def joined_wrapper_inner(args_inner, shared_list=None):
    ml.debug(fname())
    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    blast_args = deepcopy(args_inner)
    locarna_args = deepcopy(args_inner)

    if args_inner.repredict_file is None:
        fd, repred_file = mkstemp()
        os.close(fd)
    else:
        repred_file = args_inner.repredict_file

    for i, args in enumerate([blast_args, locarna_args]):
        args.prediction_method = []
        args.pred_params = dict()
        args.dump = None
        args.dill = None
        args.o_tbl = None
        args.pdf_out = None
        args.pandas_dump = None
        args.repredict_file = repred_file + str(i)
        args.dev_pred = False
        args.logfile = None

    _, blast_seq = blast_wrapper_inner(blast_args)
    _, locarna_seq = locarna_anchored_wrapper_inner(locarna_args)

    if not shared_list:
        shared_list = []

    stockholm_features = StockholmFeatureStock()
    stockholm_features.add_custom_parser_tags('GC', {'cA1': 'anchor letter tag',
                                                     'cA2': 'anchor number tag'})

    p_blast = BA_support.blast_in(args_inner.blast_in, b=args_inner.b_type)
    # this is done for each query

    ml_out_line = []
    all_analyzed = []
    for (iteration, bhp), b_hits, l_hits in zip(enumerate(p_blast), blast_seq, locarna_seq):
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

        all_analyzed.append(analyzed_hits)
        order_out = []
        for bh, lh in zip(b_hits.hits, l_hits.hits):
            hits = [bh, lh]

            bit_scores = [i.subs[i.ret_keys[0]].annotations['cmstat']['bit_sc'] for i in hits]

            mb = max(bit_scores)
            bit_index = [i for i, j in enumerate(bit_scores) if j == mb][0]
            order_out.append(bit_index)

            analyzed_hits.hits.append(hits[bit_index])

        # build the repredict file here if needed
        if args_inner.repredict_file:
            with open(blast_args.repredict_file, 'r') as barf, open(locarna_args.repredict_file, 'r') as larf, \
                    open(args_inner.repredict_file, 'w') as reprf:
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
                [[i.readline() for j in range(1)] for i in bb[1:]]

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
                    hit.subs[hit.ret_keys[0]]
                )

            # add default prediction if it is not present
            if 'ss0' not in hit.subs[hit.ret_keys[0]].letter_annotations:
                if 'sss' not in hit.subs[hit.ret_keys[0]].annotations:
                    hit.subs[hit.ret_keys[0]].anotations['sss'] = []
                hit.subs[hit.ret_keys[0]].annotations['sss'] += ['ss0']
                hit.subs[hit.ret_keys[0]].letter_annotations['ss0'] = '.'*len(hit.subs[hit.ret_keys[0]].seq)

        fda, all_hits_fasta = mkstemp()
        with os.fdopen(fda, 'w') as fah:
            # analyzed_hits.write_results_fasta(fah)
            for hit in analyzed_hits.hits:
                if len(hit.subs[hit.ret_keys[0]].seq) == 0:
                    continue
                fah.write('>{}\n{}\n'.format(hit.subs[hit.ret_keys[0]].id,
                                           str(hit.subs[hit.ret_keys[0]].seq)))

        # remove description from hits and sources
        for hit in analyzed_hits.hits:
            hit.subs[hit.ret_keys[0]].description = ''

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

    return '\n'.join(ml_out_line), all_analyzed
