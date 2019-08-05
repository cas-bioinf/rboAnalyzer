import os
from copy import deepcopy
from tempfile import mkstemp
import logging

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
from rna_blast_analyze.BR_core.expand_by_BLAST import extend_simple_core
from rna_blast_analyze.BR_core.expand_by_LOCARNA import extend_locarna_core
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger('rboAnalyzer')


def extend_meta_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model):
    ml.debug(fname())
    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    blast_args = deepcopy(args_inner)
    locarna_args = deepcopy(args_inner)
    b_all_short = deepcopy(all_short)
    l_all_short = deepcopy(all_short)

    if args_inner.repredict_file is None:
        fd, repred_file = mkstemp(prefix='rba_', suffix='_18', dir=CONFIG.tmpdir)
        os.close(fd)
    else:
        repred_file = args_inner.repredict_file

    for i, args in enumerate([blast_args, locarna_args]):
        args.prediction_method = []
        args.pred_params = dict()
        args.dump = None
        args.pdf_out = None
        args.pandas_dump = None
        args.repredict_file = repred_file + str(i)
        args.dev_pred = False
        args.logfile = None
        args.json = None
        args.html = None
        args.cm_file = ih_model

    analyzed_hits_simple = deepcopy(analyzed_hits)
    analyzed_hits_locarna = deepcopy(analyzed_hits)

    analyzed_hits_simple, _, _, _ = extend_simple_core(analyzed_hits_simple, query, blast_args, b_all_short, multi_query, iteration, ih_model)
    analyzed_hits_locarna, _, _, _ = extend_locarna_core(analyzed_hits_locarna, query, locarna_args, l_all_short, multi_query, iteration, ih_model)

    # add cmstat to query
    analyzed_hits.query = analyzed_hits_simple.query

    order_out = []

    b_dict = {BA_support.get_hit_n(h): h for h in analyzed_hits_simple.hits}
    l_dict = {BA_support.get_hit_n(h): h for h in analyzed_hits_locarna.hits}
    ok_keys = sorted(set(b_dict.keys()) | set(l_dict.keys()))
    for inum in ok_keys:
        bh = b_dict.get(inum, None)
        lh = l_dict.get(inum, None)

        hits = [bh, lh]
        # fallback to simple if locarna returned empty hit
        # deal with the situation when both ways returned empty hits
        filtered_hits = [h for h in hits if h is not None]
        if len(filtered_hits) == 1:
            msg = 'Only one extension method completed successfully for {}. ' \
                  'Choosing the successfully extended sequence to the output.'.format(filtered_hits[0].extension.id)
            ml.info(msg)
            if ml.getEffectiveLevel() < 20:
                print(msg)
            analyzed_hits.hits.append(filtered_hits[0])
            continue
        elif len(filtered_hits) == 0:
            # append empty extension
            analyzed_hits.hits_failed.append(lh)
            continue

        bit_scores = [i.extension.annotations['cmstat']['bit_sc'] for i in hits]

        mb = max(bit_scores)
        bit_index = [i for i, j in enumerate(bit_scores) if j == mb][0]
        order_out.append(bit_index)

        analyzed_hits.hits.append(hits[bit_index])

    # build failed hits
    b_dict_failed = {BA_support.get_hit_n(h): h for h in analyzed_hits_simple.hits_failed}
    l_dict_failed = {BA_support.get_hit_n(h): h for h in analyzed_hits_locarna.hits_failed}
    for inum in sorted(set(b_dict_failed) | set(l_dict_failed)):
        if inum not in ok_keys:
            if inum in b_dict_failed:
                analyzed_hits.hits_failed.append(b_dict_failed[inum])
            elif inum in l_dict_failed:
                analyzed_hits.hits_failed.append(l_dict_failed[inum])
            else:
                raise KeyError("Failed to find inum key in failed extensions. This should not happen.")

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
            hit.extension.letter_annotations['ss0'] = '.' * len(hit.extension.seq)

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
            hit.extension.letter_annotations['ss0'] = '.' * len(hit.extension.seq)

    # remove description from hits and sources
    for hit in analyzed_hits.hits:
        hit.extension.description = ''

    if args_inner.cm_file or args_inner.use_rfam:
        cm_file_rfam_user = ih_model
    else:
        cm_file_rfam_user = None
        BA_support.remove_one_file_with_try(ih_model)
    return analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user
