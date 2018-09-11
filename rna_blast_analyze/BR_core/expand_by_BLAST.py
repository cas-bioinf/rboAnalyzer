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
from random import shuffle
from tempfile import mkstemp
import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.add_usr_local_bin
import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.blast_hits_merging import merge_blast_hits
from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.stockholm_alig import StockholmFeatureStock
from rna_blast_analyze.BR_core import BA_verify


def blast_wrapper(args_inner, shared_list=None):
    ret_line, _ = blast_wrapper_inner(args_inner, shared_list=shared_list)
    return ret_line


def blast_wrapper_inner(args_inner, shared_list=None):
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
            extra=args_inner.subseq_window_simple_ext,
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

        analyzed_hits.query = query

        def create_blast_only_report_object(exp_hit, query_len):
            # get extension by blast indices
            # create new record by those indices
            s = exp_hit.annotations['extended_start'] - exp_hit.annotations['super_start']
            e = exp_hit.annotations['extended_end'] - exp_hit.annotations['super_end']

            nseq = exp_hit.seq[s:e]
            nr = SeqRecord(nseq, id=exp_hit.id)

            nr.letter_annotations['ss0'] = '.' * len(nr.seq)
            nr.annotations['sss'] = ['ss0']

            hit = BA_support.Subsequences(exp_hit)
            sub_id = nr.id[:-2].split(':')[1]
            hit.subs[sub_id] = nr

            pos_match = re.search(str(nr.seq), str(exp_hit.seq), flags=re.IGNORECASE)
            if not pos_match:
                raise Exception('Subsequnce not found in supersequence. Terminating.')

            bl = exp_hit.annotations['blast'][1]

            if bl.sbjct_start < bl.sbjct_end:
                bls = bl.sbjct_start - bl.query_start
                ble = bl.sbjct_end + (query_len - bl.query_end)
            elif bl.sbjct_end < bl.sbjct_start:
                bls = bl.sbjct_end - (query_len - bl.query_end)
                ble = bl.sbjct_start + bl.query_start
            else:
                raise Exception('unknown strand option')

            # if whole subject sequence too short, this assertion will fail
            if exp_hit.annotations['trimmed_start'] or  exp_hit.annotations['trimmed_end']:
                print('skipping check')
            else:
                assert len(nr.seq) == abs(bls - ble)
                assert bls + 1 == exp_hit.annotations['extended_start']
                assert ble == exp_hit.annotations['extended_end']

            hit.ret_keys = [sub_id,
                            'ss0',
                            None]
            hit.best_start = hit.source.annotations['super_start'] + pos_match.span()[0]
            hit.best_end = hit.source.annotations['super_start'] + pos_match.span()[1] - 1
            return hit

        # blast only extension
        for exp_hit in shorts_expanded:
            analyzed_hits.hits.append(create_blast_only_report_object(exp_hit, len(query_seq)))

        # assign Locarna score to None as it is not directly accessible from mlocarna
        for hit in analyzed_hits.hits:
            hit.subs[hit.ret_keys[0]].annotations['score'] = None

        # infer homology
        # consider adding query sequence to alignment and scoring it as a proof against squed alignments and datasets
        # write all hits to fasta
        fda, all_hits_fasta = mkstemp()
        with os.fdopen(fda, 'w') as fah:
            # analyzed_hits.write_results_fasta(fah)
            for hit in analyzed_hits.hits:
                if len(hit.subs[hit.ret_keys[0]].seq) == 0:
                    continue
                fah.write('>{}\n{}\n'.format(hit.subs[hit.ret_keys[0]].id,
                                           str(hit.subs[hit.ret_keys[0]].seq)))

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
