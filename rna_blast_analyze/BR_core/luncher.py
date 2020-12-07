import os
import sys
import pickle
import re
from copy import deepcopy
from random import shuffle
from tempfile import mkstemp
import logging
import json
from Bio import SeqIO

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core import validate_args
from rna_blast_analyze.BR_core.BA_methods import BlastSearchRecompute, to_tab_delim_line_simple
from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
from rna_blast_analyze.BR_core.infer_homology import find_and_extract_cm_model
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.cmalign import run_cmfetch, get_cm_model, RfamInfo
from rna_blast_analyze.BR_core.expand_by_BLAST import extend_simple_core
from rna_blast_analyze.BR_core.expand_by_LOCARNA import extend_locarna_core
from rna_blast_analyze.BR_core.expand_by_joined_pred_with_rsearch import extend_meta_core
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecompute2dict, blastsearchrecomputefromdict
from rna_blast_analyze.BR_core import exceptions

ml = logging.getLogger('rboAnalyzer')


def lunch_computation(args_inner, shared_list=None):
    ml.debug(fname())
    if not shared_list:
        shared_list = []

    # update params if different config is requested
    CONFIG.override(tools_paths(args_inner.config_file))

    p_blast = BA_support.blast_in(args_inner.blast_in, b=args_inner.b_type)
    query_seqs = [i for i in SeqIO.parse(args_inner.blast_query, 'fasta')]

    if len(p_blast) != len(query_seqs):
        ml.error('Number of query sequences in provided BLAST output file ({}) does not match number of query sequences'
                 ' in query FASTA file ({}).'.format(len(p_blast), len(query_seqs)))
        sys.exit(1)

    # check if BLAST does not contain unexpected sequence characters
    validate_args.check_blast(p_blast)

    # create list of correct length if needed
    all_saved_data = [None] * len(query_seqs)
    saved_file = '{}.r-{}'.format(args_inner.blast_in, args_inner.sha1[:10])
    with open(saved_file, 'r+') as f:
        _saved = json.load(f)
        if _saved is None:
            f.seek(0)
            f.truncate()
            json.dump(all_saved_data, f)
        else:
            msg = "Loading backup data."
            print('STATUS: ' + msg)
            ml.info(msg + ' file: ' + saved_file)
            all_saved_data = _saved

            for saved_data in all_saved_data:
                # we can have partially computed data
                if saved_data is None:
                    continue
                if saved_data['args']['sha1'] != args_inner.sha1:
                    msg = "Input argument hash does not match the saved argument hash. "
                    if saved_data['args']['sha1'][:10] == args_inner.sha1[:10]:
                        msg += "This is because of truncating hashes to first 10 characters. "
                        msg += "Please remove the '{}' file.".format(saved_file)
                        ml.error(msg)
                        sys.exit(1)
                    else:
                        msg += "Please remove the '{}' file.".format(saved_file)
                        sys.exit(1)

    if len(p_blast) > 1:
        multi_query = True
    else:
        multi_query = False

    # this is done for each query
    ml_out_line = []
    all_analyzed = []
    for iteration, (bhp, query, saved_data) in enumerate(zip(p_blast, query_seqs, all_saved_data)):
        if saved_data is None:
            print('STATUS: processing query: {}'.format(query.id))
            validate_args.verify_query_blast(blast=bhp, query=query)

            analyzed_hits = BlastSearchRecompute(args_inner, query, iteration)
            analyzed_hits.multi_query = multi_query

            # run cm model build
            # allows to fail fast if rfam was selected and we dont find the model
            try:
                ih_model, analyzed_hits = find_and_extract_cm_model(args_inner, analyzed_hits)
            except (exceptions.MissingCMexception, exceptions.SubprocessException):
                sys.exit(1)

            # select all
            all_blast_hits = BA_support.blast_hsps2list(bhp)

            if len(all_blast_hits) == 0:
                ml.error('No hits found in {} - {}. Nothing to do.'.format(args_inner.blast_in, bhp.query))
                continue

            # filter if needed
            if args_inner.filter_by_eval is not None:
                tmp = filter_by_eval(all_blast_hits, BA_support.blast_hit_getter_from_hits, args_inner.filter_by_eval)
                if len(tmp) == 0 and len(all_blast_hits) != 0:
                    ml.error('The requested filter removed all BLAST hits {} - {}. Nothing to do.'.format(args_inner.blast_in,
                                                                                                          bhp.query))
                    continue
            elif args_inner.filter_by_bitscore is not None:
                tmp = filter_by_bits(all_blast_hits, BA_support.blast_hit_getter_from_hits, args_inner.filter_by_bitscore)
                if len(tmp) == 0 and len(all_blast_hits) != 0:
                    ml.error('The requested filter removed all BLAST hits {} - {}. Nothing to do.'.format(args_inner.blast_in,
                                                                                                          bhp.query))
                    continue

            all_short = all_blast_hits

            # now this is different for each mode
            if args_inner.mode == 'simple':
                analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user = extend_simple_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model)
            elif args_inner.mode == 'locarna':
                analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user = extend_locarna_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model)
            elif args_inner.mode == 'meta':
                analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user = extend_meta_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model)
            else:
                raise ValueError('Unknown option - should be cached by argparse.')

            if len(analyzed_hits.hits) == 0:
                ml.error(
                    "Extension failed for all sequences. Please see the error message. You can also try '--mode simple'."
                )
                sys.exit(1)

            analyzed_hits.copy_hits()

            with open(args_inner.blast_in + '.r-' + args_inner.sha1[:10], 'r+') as f:
                all_saved_data = json.load(f)
                all_saved_data[iteration] = blastsearchrecompute2dict(analyzed_hits)
                f.seek(0)
                f.truncate()
                json.dump(all_saved_data, f, indent=2)

        else:
            print('STATUS: extended sequences loaded from backup file for query {}'.format(query.id))
            analyzed_hits = blastsearchrecomputefromdict(saved_data)

            # overwrite the saved args with current
            # this will update used prediction methods and other non essential stuff
            analyzed_hits.args = args_inner

            if analyzed_hits.args.cm_file:
                cm_file_rfam_user = analyzed_hits.args.cm_file
            else:
                cm_file_rfam_user = None

        all_analyzed.append(analyzed_hits)

        # write all hits to fasta
        fda, all_hits_fasta = mkstemp(prefix='rba_', suffix='_22', dir=CONFIG.tmpdir)
        os.close(fda)
        analyzed_hits.write_results_fasta(all_hits_fasta)

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
            if cm_file_rfam_user is None and 'rfam' in ''.join(args_inner.prediction_method):
                best_model = get_cm_model(args_inner.blast_query, threads=args_inner.threads)
                rfam = RfamInfo()
                cm_file_rfam_user = run_cmfetch(rfam.file_path, best_model)

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
                        used_cm_file=cm_file_rfam_user,
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
                used_cm_file=cm_file_rfam_user,
                multi_query=multi_query,
                iteration=iteration,
            )
            out_line.append(to_tab_delim_line_simple(args_inner))

        ml_out_line.append('\n'.join(out_line))

        if cm_file_rfam_user is not None and os.path.exists(cm_file_rfam_user):
            BA_support.remove_one_file_with_try(cm_file_rfam_user)

        BA_support.remove_one_file_with_try(all_hits_fasta)
    return '\n'.join(ml_out_line), all_analyzed