import logging
import os
import pickle
import gzip
import numpy as np
import pandas
import json
from tempfile import mkstemp
from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.BA_methods import HitList
from rna_blast_analyze.BR_core.BA_support import iter2file_name, add_loc_to_description
from rna_blast_analyze.BR_core.centroid_homfold import me_centroid_homfold, centroid_homfold_fast
from rna_blast_analyze.BR_core.cmalign import RfamInfo, run_cmfetch, extract_ref_from_cm, run_cmemit
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.infer_homology import infer_hits_cm
from rna_blast_analyze.BR_core.predict_structures import alifold_refold_prediction, \
    rnafold_wrap_for_predict, subopt_fold_query,\
    subopt_fold_alifold, cmmodel_rnafold_c, rfam_subopt_pred
from rna_blast_analyze.BR_core.predict_structures import find_nc_and_remove, check_lonely_bp
from rna_blast_analyze.BR_core.turbofold import turbofold_fast, turbofold_with_homologous
from rna_blast_analyze.BR_core.output.htmloutput import write_html_output
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecompute2dict
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core import exceptions
from hashlib import sha1

ml = logging.getLogger('rboAnalyzer')

safe_prediction_method = [
    'rnafold',
    'centroid-fast',
    'Turbo-fast',
    'rfam-Rc',
    'sub-fq'
]


def wrapped_ending_with_prediction(
    args_inner, analyzed_hits, pred_method=None, method_params=None, used_cm_file=None, multi_query=False, iteration=0,
):
    """
    wrapper for prediction of secondary structures
    :param args_inner: Namespace of input arguments
    :param analyzed_hits: BlastSearchRecompute object
    :param all_hits_fasta: fasta file with all extended sequences
    :param query: query sequence
    :param pred_method:
    :param method_params:
    :param used_cm_file: cmfile if cmfile is known (user given or computed)
    :return:
    """
    ml.debug(fname())
    exec_time = {}
    msg = 'Entering structure prediction..'
    if ml.level < 21:
        ml.info(msg)
    else:
        print(msg)
        ml.info(msg)

    if pred_method is None:
        pred_method = args_inner.prediction_method

    if isinstance(pred_method, str):
        pred_method = (pred_method,)

    if method_params is None:
        method_params = args_inner.pred_params

    # ======= filter if needed =======
    # do the filtering based on e-val or bitscore
    # homologous hits still gets used for prediction

    # annotate ambiguous bases
    query = BA_support.annotate_ambiguos_base(analyzed_hits.query)

    # copy the list before filtering
    all_hits_list = [i.extension for i in analyzed_hits.get_all_hits()]

    if args_inner.filter_by_eval is not None:
        hits2predict = filter_by_eval(
            analyzed_hits.get_all_hits(), BA_support.blast_hit_getter_from_subseq, args_inner.filter_by_eval
        )
        _hits = HitList()
        for h in hits2predict:
            _hits.append(h)
        analyzed_hits.hits = _hits
    elif args_inner.filter_by_bitscore is not None:
        hits2predict = filter_by_bits(
            analyzed_hits.get_all_hits(), BA_support.blast_hit_getter_from_subseq, args_inner.filter_by_bitscore
        )
        _hits = HitList()
        for h in hits2predict:
            _hits.append(h)
        analyzed_hits.hits = _hits
    else:
        analyzed_hits.hits = analyzed_hits.get_all_hits()

    # if used_cm_file is provided do not override it with CM from RFAM
    # if use_rfam flag was given, then used_cm_file is already the best_matching model
    # if analyzed_hits.best_matching_model is None - then we could not find the best matching model in RFAM
    #  and the rfam based methods should fail (i.e. not predict anything)
    delete_cm = False
    if used_cm_file is None and analyzed_hits.best_matching_model is not None:
        rfam = RfamInfo()
        used_cm_file = run_cmfetch(rfam.file_path, analyzed_hits.best_matching_model['target_name'])
        delete_cm = True

    fd, seqs2predict_fasta = mkstemp(prefix='rba_', suffix='_83', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as fah:
        for hit in analyzed_hits.hits:
            if len(hit.extension.seq) == 0:
                continue
            fah.write('>{}\n{}\n'.format(
                hit.extension.id,
                str(hit.extension.seq))
            )

    if not isinstance(method_params, dict):
        raise Exception('prediction method parameters must be python dict')

    # prediction methods present in analyzed_hits
    #  which might be loaded from intermediate file

    # check if structures of a method are predicted for all required hit
    # check if the prediction parameters of such method were same

    prediction_msgs = []
    # compute prediction methods which were not computed
    for pkey in set(pred_method):
        # add sha1 hashes
        nh = sha1()
        nh.update(str(sorted(method_params.get(pkey, {}).items())).encode())
        current_hash = nh.hexdigest()

        if all(pkey in h.extension.letter_annotations for h in analyzed_hits.hits) and \
                len(
                    {
                        h.extension.annotations.get('sha1', {}).get(pkey, None) for h in analyzed_hits.hits
                    } | {current_hash, }
                ) == 1:
            msg_skip = 'All structures already computed for {}. Skipping...'.format(pkey)
            ml.info(msg_skip)
            if ml.level > 20:
                print(msg_skip, flush=True)
            continue

        msg_run = 'Running: {}...'.format(pkey)
        ml.info(msg_run)

        if ml.level > 20:
            print(msg_run, flush=True)

        structures, etime, msgs = repredict_structures_for_homol_seqs(
            query,
            seqs2predict_fasta,
            args_inner.threads,
            prediction_method=pkey,
            pred_method_params=method_params,
            all_hits_list=all_hits_list,
            seqs2predict_list=[i.extension for i in analyzed_hits.hits],
            use_cm_file=used_cm_file,
            )

        exec_time[pkey] = etime

        if structures is None:
            msg = 'Structures not predicted with {} method'.format(pkey)
            ml.info(msg)
            if ml.level > 20:
                print('STATUS: ' + msg)

        else:
            for i, hit in enumerate(analyzed_hits.hits):
                assert str(hit.extension.seq) == str(structures[i].seq)
                hit.extension.annotations['sss'] += [pkey]

                hit.extension.annotations['msgs'] += structures[i].annotations.get('msgs', [])

                # expects "predicted" in annotations - for now, if not given, default is True, as not all prediction
                #  methods implement "predicted" in their output
                if structures[i].annotations.get('predicted', True):
                    hit.extension.letter_annotations[pkey] = structures[i].letter_annotations['ss0']

                if 'sha1' not in hit.extension.annotations:
                    hit.extension.annotations['sha1'] = dict()
                hit.extension.annotations['sha1'][pkey] = current_hash

                try:
                    del hit.extension.letter_annotations['ss0']
                except KeyError:
                    pass
                try:
                    hit.extension.annotations['sss'].remove('ss0')
                except ValueError:
                    pass

            analyzed_hits.update_hit_stuctures()

        # check if msgs are not empty
        if msgs:
            prediction_msgs.append('{}: {}'.format(pkey, '\n'.join(msgs)))

        analyzed_hits.msgs = prediction_msgs

        with open(args_inner.blast_in + '.r-' + args_inner.sha1[:10], 'r+') as f:
            all_saved_data = json.load(f)
            all_saved_data[iteration] = blastsearchrecompute2dict(analyzed_hits)
            f.seek(0)
            f.truncate()
            json.dump(all_saved_data, f, indent=2)

    # remove structures predicted by different methods (which might be saved from previous computation)
    for hit in analyzed_hits.hits:
        for pkey in set(hit.extension.letter_annotations.keys()):
            if pkey not in pred_method:
                del hit.extension.letter_annotations[pkey]
                try:
                    hit.extension.annotations['sss'].remove(pkey)
                except ValueError:
                    pass

    BA_support.remove_one_file_with_try(seqs2predict_fasta)

    if delete_cm:
        BA_support.remove_one_file_with_try(used_cm_file)

    add_loc_to_description(analyzed_hits)

    # write html if requested
    if args_inner.html:
        html_file = iter2file_name(args_inner.html, multi_query, iteration)
        ml.info('Writing html to {}.'.format(html_file))
        with open(html_file, 'w') as h:
            h.write(write_html_output(analyzed_hits))

    # write csv file if requested
    if args_inner.csv:
        csv_file = iter2file_name(args_inner.csv, multi_query, iteration)
        ml.info('Writing csv to {}.'.format(csv_file))
        analyzed_hits.to_csv(csv_file)

    # replace with json
    if args_inner.json:
        json_file = iter2file_name(args_inner.json, multi_query, iteration)
        ml.info('Writing json to {}.'.format(json_file))
        j_obj = json.dumps(blastsearchrecompute2dict(analyzed_hits), indent=2)
        if getattr(args_inner, 'zip_json', False):
            with open(json_file + '.gz', 'wb') as ff:
                ff.write(
                    gzip.compress(j_obj.encode())
                )
        else:
            with open(json_file, 'w') as ff:
                ff.write(j_obj)

    if args_inner.pandas_dump:
        pickle_file = iter2file_name(args_inner.pandas_dump, multi_query, iteration)
        ml.info('Writing pandas pickle to {}.'.format(pickle_file))
        pandas.to_pickle(analyzed_hits.pandas, pickle_file)

    if args_inner.dump:
        dump_file = iter2file_name(args_inner.dump, multi_query, iteration)
        ml.info('Writing dump files base: {}.'.format(dump_file))
        with open(dump_file, 'wb') as pp:
            pickle.dump(analyzed_hits, pp, pickle.HIGHEST_PROTOCOL)

        with open(dump_file + '.time_dump', 'wb') as pp:
            pickle.dump(exec_time, pp, pickle.HIGHEST_PROTOCOL)

    return analyzed_hits


def create_nr_trusted_hits_file_MSA_safe(
        sim_threshold_percent=None,
        all_hits=None,
        query=None,
        cmscore_tr=-2.03,
        cm_threshold_percent=None,
        check_unambiguous=False,
        len_diff=0.1,
):
    """
    create non redundant trusted hits file

    multiple at minimum (2) sequences are needed for profile alignment for some aligners
    so this function always return two or more sequences or raises exception

    :param sim_threshold_percent:   seq similarity threshold for homology exclusion
    :param all_hits:                list of hits
    :param query:                   blast query
    :param cmscore_tr:              threshold for homology inclusion in bits
    :param cm_threshold_percent:    threshold for homology inclusion in % of query bits
    :param check_unambiguous:       bool whether to check unambiguous seqs
    :param len_diff:                threshold for exclusion of hits lq=len(query) lq - diff*lq < len(seq) < lq + diff*lq
    :return:
    """
    ml.debug(fname())
    # i need to leave query, even if with ambiguous basepairs in
    # because it is used as an reference during distance computation and subsequence selection,
    # however i don't need to have all homologous seqs there
    if check_unambiguous:
        all_hits = BA_support.filter_ambiguous_seqs_from_list(all_hits)

    dist_table, homologous_seqs, msgs = _trusted_hits_selection_wrapper(
        all_hits,
        query,
        cmscore_tr,
        cm_threshold_percent,
        len_diff_=len_diff
    )

    if dist_table.size == 0:
        raise exceptions.NoHomologousSequenceException

    to_include = BA_support.select_sequences_from_similarity_rec(
        dist_table,
        sim_threshold_percent=sim_threshold_percent
    )
    nr_homolog_hits = [homologous_seqs[i] for i in to_include]

    # final checking of nr homologs
    # if sequence is filtered here, it is ambiguous basepair in query
    # removing it is fine if multiple homologous sequences are present
    # the problem will arise when only 1 homologous sequence will remain
    # if we added sequence in previous step, raise exception, else behave like in prev step
    # what if trusted hit is only one?

    msg = (
        'STATUS: Only one sequence remained under defined "pred_sim_threshold" parameter.\n'
        ' Mitigation: Adding the most dissimilar homologous sequence to the non redundant sequences list.'
    )
    if len(nr_homolog_hits) < 2 and not check_unambiguous:
        msgs.append(msg)
        ml.info(msg)
        if ml.level > 20:
            print(msg)

        dis_hom_index = dist_table[:, 0].argmin()
        nr_homolog_hits.append(SeqRecord(homologous_seqs[dis_hom_index].seq, id='dummy_seq_01'))
        del dis_hom_index

    elif len(nr_homolog_hits) < 2 and check_unambiguous:
        if len(BA_support.filter_ambiguous_seqs_from_list(nr_homolog_hits)) == 0:
            # this mean query contain ambiguous bases
            raise exceptions.NoHomologousSequenceException
        else:
            msgs.append(msg)
            ml.info(msg)
            if ml.level > 20:
                print(msg)

            dis_hom_index = dist_table[:, 0].argmin()
            nr_homolog_hits.append(SeqRecord(homologous_seqs[dis_hom_index].seq, id='dummy_seq_01'))
            del dis_hom_index
        homologous_seqs = BA_support.filter_ambiguous_seqs_from_list(homologous_seqs)

    elif len(nr_homolog_hits) >= 2 and not check_unambiguous:
        pass

    elif len(nr_homolog_hits) > 2 and check_unambiguous:
        nr_homolog_hits = BA_support.filter_ambiguous_seqs_from_list(nr_homolog_hits)
        homologous_seqs = BA_support.filter_ambiguous_seqs_from_list(homologous_seqs)

    elif len(nr_homolog_hits) == 2 and check_unambiguous:
        homologous_seqs = BA_support.filter_ambiguous_seqs_from_list(homologous_seqs)
        if len(BA_support.filter_ambiguous_seqs_from_list(nr_homolog_hits)) == 1:
            # this mean that query contains ambiguous base
            raise exceptions.NoHomologousSequenceException

    else:
        raise Exception()

    fd_h, nr_homo_hits_file = mkstemp(prefix='rba_', suffix='_58', dir=CONFIG.tmpdir)
    with os.fdopen(fd_h, 'w') as f:
        SeqIO.write(nr_homolog_hits, f, 'fasta')

    return nr_homo_hits_file, homologous_seqs, msgs


def create_nr_homolog_hits_file_MSA_unsafe(sim_threshold_percent=None, all_hits=None, query=None, cmscore_tr=0.0,
                                           cm_threshold_percent=None, len_diff=0.1):
    """
    create non redundant homologous hits file
    """
    ml.debug(fname())
    dist_table, homologous_seqs, msgs = _trusted_hits_selection_wrapper(
        all_hits,
        query,
        cmscore_tr,
        cm_threshold_percent,
        len_diff_=len_diff
    )
    if dist_table.size == 0:
        nr_homolog_hits = [query]
    else:
        # normal execution
        to_include = BA_support.select_sequences_from_similarity_rec(
            dist_table,
            sim_threshold_percent=sim_threshold_percent
        )
        nr_homolog_hits = [homologous_seqs[i] for i in to_include]

    fd_h, nr_homo_hits_file = mkstemp(prefix='rba_', suffix='_59', dir=CONFIG.tmpdir)
    with os.fdopen(fd_h, 'w') as f:
        SeqIO.write(nr_homolog_hits, f, 'fasta')

    return nr_homo_hits_file, homologous_seqs, msgs


def _extract_cmscore_from_hom_seqs(hom_seqs):
    """
    will return list of cm_scores as obtained from homology inference
     usefull for selecting sequences more relevant
    :param hom_seqs:
    :return:
    """
    return [i.annotations['cmstat']['bit_sc'] for i in hom_seqs]


def _trusted_hits_selection_wrapper(all_hits_, query_, cmscore_tr_, cm_threshold_percent_, len_diff_=0.1):
    """
    runs basic non_redundant sequences calculation (ie exact sequence match)
    selects homologous sequences from all hits list by cmscore threshold or by query sequence

    behaviour:
        will return distance array with similarities in % including query sequence and list of homologous sequences
        including query sequence

        if no sequence is homologous
        it will return empty array for distance matrix and list with query sequence
    """
    ml.debug(fname())
    msgs = []
    # trusted sequence selection
    # ========================================================
    assert (cmscore_tr_ == 0) or cm_threshold_percent_ is None

    score = _extract_cmscore_from_hom_seqs(all_hits_)

    if cm_threshold_percent_ is not None:
        selection_threshold = cm_threshold_percent_ * query_.annotations['cmstat'].bit_sc / 100
    else:
        selection_threshold = cmscore_tr_

    pred = infer_hits_cm(score, tr=selection_threshold)
    trusted_seqs_ = [i for i, j in zip(all_hits_, pred) if j]

    if len(trusted_seqs_) == 0:
        msg = 'STATUS: No estimated full-length sequences from BLAST output ' \
              'selected as reference for structure prediction.\n' \
              ' Using query sequence as reference.'
        msgs.append(msg)
        ml.info(msg)
        if ml.level > 20:
            print(msg)
        return np.empty(0), [query_], msgs

    # add query to trusted sequences
    trusted_seqs_query = [query_] + trusted_seqs_

    # make nr list of sequences -> faster alignment
    # better selection
    nr_trusted_seqs_query = BA_support.non_redundant_seqs(trusted_seqs_query)

    # check if the homologous sequence is not exact match as query
    #  (ie taking non redundant set would be only one sequence)
    if len(nr_trusted_seqs_query) == 1:
        msg = 'STATUS: All sequences selected as reference are exactly same as query sequence.'
        msgs.append(msg)
        ml.info(msg)
        if ml.level > 20:
            print(msg)
        return np.empty(0), [query_], msgs

    # select only sequences in some predifined length range to query
    # this is needed for longish ncRNAs
    #   tolerate 10 % length difference?
    ref_len = len(query_)
    nr_len_selected_trusted = [
        seq for seq in nr_trusted_seqs_query if ref_len * (1 - len_diff_) < len(seq) < ref_len * (1 + len_diff_)
    ]

    # this is to control if only one sequence remained after filtering for length difference
    if len(nr_len_selected_trusted) == 1:
        msg = \
            'No sequence satisfy the length difference condition ({}: {}-{})'.format(
                len_diff_,
                ref_len * (1 - len_diff_),
                ref_len * (1 + len_diff_)
            )
        msgs.append(msg)
        ml.info(msg)
        if ml.level > 20:
            print(msg)
        return np.empty(0), [query_], msgs

    # sanitize seq names (muscle has issues with too long names)
    san_hom_seqs, san_dict = BA_support.sanitize_fasta_names_in_seqrec_list(nr_len_selected_trusted)

    c_fd, trusted_sequence_file_ = mkstemp(prefix='rba_', suffix='_60', dir=CONFIG.tmpdir)
    with os.fdopen(c_fd, 'w') as f:
        SeqIO.write(san_hom_seqs, f, 'fasta')

    align_file = BA_support.run_muscle(trusted_sequence_file_, reorder=True)
    alig = AlignIO.read(align_file, format='clustal')
    distance_calc = DistanceCalculator(model='identity')
    dist_mat = distance_calc.get_distance(alig)
    # rebuild index from sanitized
    orig_index = [san_dict[i] for i in dist_mat.names]
    dist_mat_pd = pandas.DataFrame.from_records(dist_mat.matrix, index=orig_index)
    dist_table_ = (1 - dist_mat_pd.values) * 100

    BA_support.remove_files_with_try([align_file, trusted_sequence_file_])
    return dist_table_, trusted_seqs_query, msgs


def nonhomseqwarn(method_name):
    msg = \
        'No sequence was inferred homologous, need at least 1 for {} type of prediction.\n' \
        'Try changing the parameters for selection of reference sequences' \
        ' (e.g. the "cmscore_percent" and "pred_sim_threshold")\n' \
        'For more information see the docs/prediction_methods.md\n' \
        'Or try one of following prediction methods: {}.'.format(
            method_name,
            ', '.join(safe_prediction_method)
        )
    ml.warning(msg)
    if ml.level > 30:
        print(msg)
    return msg


def repredict_structures_for_homol_seqs(
        query, seqs2predict_fasta,
        threads=None,
        prediction_method=None,
        pred_method_params=None,
        all_hits_list=None,
        seqs2predict_list=None,
        use_cm_file=None,
):
    """Run RNA structure prediction based on chosen method and parameters.
    """

    default_sim_tr_perc = 90
    default_score_tr = 0.0
    query_max_len_diff = 0.1

    try:
        if 'default' == prediction_method:
            # do nothing
            return None, None, []

        elif 'rfam-Rc' == prediction_method:
            if use_cm_file is None:
                msg = "No CM model. Can't use {}.".format(prediction_method)
                ml.warning(msg)
                return None, None, [msg]
            else:
                structures, exec_time = cmmodel_rnafold_c(
                    seqs2predict_fasta,
                    use_cm_file,
                    threads=threads,
                    params=pred_method_params.get(prediction_method, {})
                )
                return structures, exec_time, []

        elif 'rfam-centroid' == prediction_method:
            # run cmscan if needed
            # run cmfetch
            # run cmemit -> homologous seqs
            # run centroid_homfold

            method_parameters = pred_method_params.get(prediction_method, {})
            if use_cm_file is None:
                msg = "No CM model. Can't use {}.".format(prediction_method)
                ml.warning(msg)
                return None, None, [msg]
            else:
                cep = method_parameters.get('cmemit', '')
                if '-u' not in cep:
                    cep += ' -u'
                if '-N' not in cep:
                    cep += ' -N {}'.format(method_parameters.get('n_seqs', 10))

                hf_file = run_cmemit(use_cm_file, params=cep)

                structures, exec_time = me_centroid_homfold(seqs2predict_fasta, hf_file, params=method_parameters)

                BA_support.remove_one_file_with_try(hf_file)
                return structures, exec_time, []

        elif 'rfam-sub' == prediction_method:
            if use_cm_file is None:
                msg = "No CM model. Can't use {}.".format(prediction_method)
                ml.warning(msg)
                return None, None, [msg]
            else:
                ref_structure = extract_ref_from_cm(use_cm_file)

                structures, exec_time = rfam_subopt_pred(
                    seqs2predict_fasta,
                    ref_structure,
                    params=pred_method_params.get(prediction_method, None),
                    threads=threads,
                )
                return structures, exec_time, []

        elif 'rnafold' == prediction_method:
            structures, exec_time = rnafold_wrap_for_predict(
                seqs2predict_fasta,
                params=pred_method_params.get(prediction_method, {}).get('RNAfold', '')
            )
            return structures, exec_time, []

        elif 'fq-sub' == prediction_method:
            a, qf = mkstemp(prefix='rba_', suffix='_55', dir=CONFIG.tmpdir)
            with os.fdopen(a, 'w') as fd:
                fd.write('>query\n{}\n'.format(str(query.seq)))

            structures, exec_time = subopt_fold_query(
                seqs2predict_fasta,
                qf,
                params=pred_method_params.get(prediction_method, None),
                threads=threads
            )
            BA_support.remove_one_file_with_try(qf)
            return structures, exec_time, []

        elif 'C-A-sub' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, homologous_seqs, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file

            f, homologous_sequence_file = mkstemp(prefix='rba_', suffix='_64', dir=CONFIG.tmpdir)
            with os.fdopen(f, 'w') as fh:
                SeqIO.write(homologous_seqs, fh, 'fasta')

            structures, exec_time = subopt_fold_alifold(
                seqs2predict_fasta,
                homologous_sequence_file,
                aligner='clustalo',
                params=method_parameters,
                threads=threads
            )
            BA_support.remove_one_file_with_try(homologous_sequence_file)
            del homologous_sequence_file
            del homologous_seqs
            return structures, exec_time, msgs

        elif 'M-A-sub' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, homologous_seqs, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file

            f, homologous_sequence_file = mkstemp(prefix='rba_', suffix='_65', dir=CONFIG.tmpdir)
            with os.fdopen(f, 'w') as fh:
                SeqIO.write(homologous_seqs, fh, 'fasta')

            structures, exec_time = subopt_fold_alifold(
                seqs2predict_fasta,
                homologous_sequence_file,
                aligner='muscle',
                params=method_parameters,
                threads=threads,
            )

            BA_support.remove_one_file_with_try(homologous_sequence_file)
            del homologous_sequence_file
            del homologous_seqs
            return structures, exec_time, msgs

        elif 'C-A-r-Rc' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            structures, exec_time = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='clustalo'
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file
            return structures, exec_time, msgs

        elif 'M-A-r-Rc' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures, exec_time = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='muscle'
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file
            return structures, exec_time, msgs

        elif 'C-A-U-r-Rc' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures, exec_time = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='clustalo'
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file
            return structures, exec_time, msgs

        elif 'M-A-U-r-Rc' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures, exec_time = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='muscle'
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file
            return structures, exec_time, msgs

        elif 'centroid' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_homolog_hits_file_MSA_unsafe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            raw_structures, exec_time = me_centroid_homfold(
                seqs2predict_fasta, nr_homo_hits_file,
                params=method_parameters
            )

            # check noncanonical
            if prediction_method in pred_method_params and pred_method_params[prediction_method]:
                allow_nc = pred_method_params[prediction_method].get('allow_noncanonical', False)
                allow_lp = pred_method_params[prediction_method].get('allow_lonely_pairs', False)
            else:
                allow_nc = False
                allow_lp = False
            if not allow_nc:
                for seq in raw_structures:
                    repstr = find_nc_and_remove(str(seq.seq), structure=seq.letter_annotations['ss0'])
                    seq.letter_annotations['ss0'] = repstr

            # check lonely basepairs
            if not allow_lp:
                for seq in raw_structures:
                    repstr = check_lonely_bp(seq.letter_annotations['ss0'])
                    seq.letter_annotations['ss0'] = repstr

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del nr_homo_hits_file
            return raw_structures, exec_time, msgs

        elif 'centroid-fast' == prediction_method:
            method_parameters = pred_method_params.get(prediction_method, {})
            if query.annotations['ambiguous']:
                raise exceptions.AmbiguousQuerySequenceException

            raw_structures, exec_time = centroid_homfold_fast(
                all_seqs=all_hits_list,
                query=query,
                all_seqs_fasta=seqs2predict_fasta,
                n=method_parameters.get('max_seqs_in_prediction', 10),
                centroid_homfold_params=method_parameters,
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff)
            )

            # check noncanonical
            if prediction_method in pred_method_params and pred_method_params[prediction_method]:
                allow_nc = pred_method_params[prediction_method].get('allow_noncanonical', False)
                allow_lp = pred_method_params[prediction_method].get('allow_lonely_pairs', False)
            else:
                allow_nc = False
                allow_lp = False
            if not allow_nc:
                for seq in raw_structures:
                    repstr = find_nc_and_remove(str(seq.seq), structure=seq.letter_annotations['ss0'])
                    seq.letter_annotations['ss0'] = repstr

            # check lonely basepairs
            if not allow_lp:
                for seq in raw_structures:
                    repstr = check_lonely_bp(seq.letter_annotations['ss0'])
                    seq.letter_annotations['ss0'] = repstr

            return raw_structures, exec_time, []

        elif 'TurboFold' == prediction_method:
            # set arbitrary sim_threshold_percent to 100, because we want to remove only identical sequences from prediction
            #  with TurboFold. The structure of redundant sequences will be set according to the one in prediction
            all_hits_filtered = BA_support.filter_ambiguous_seqs_from_list(all_hits_list)
            seqs2predict_filtered = BA_support.filter_ambiguous_seqs_from_list(seqs2predict_list)
            if len(seqs2predict_list) != len(seqs2predict_filtered):
                ml.warning('Some sequences contain ambiguous bases - they will not be predicted.')

            if query.annotations['ambiguous']:
                raise exceptions.AmbiguousQuerySequenceException()

            method_parameters = pred_method_params.get(prediction_method, {})

            nr_homo_hits_file, _, msgs = create_nr_homolog_hits_file_MSA_unsafe(
                all_hits=all_hits_filtered,
                query=query,
                sim_threshold_percent=100,
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
                )

            with open(nr_homo_hits_file, 'r') as nrf:
                nr_homo_hits = [seq for seq in SeqIO.parse(nrf, format='fasta')]

            nh = sha1()
            nh.update(str(sorted(method_parameters.items())).encode())
            nh_str = nh.hexdigest()

            structures_t, exec_time = turbofold_with_homologous(
                all_sequences=seqs2predict_filtered,
                nr_homologous=nr_homo_hits,
                params=method_parameters.get('TurboFold', {}),
                n=method_parameters.get('max_seqs_in_prediction', 3),
                cpu=threads,
                pkey=prediction_method,
                sha1val=nh_str,
            )

            structures = BA_support.rebuild_structures_output_from_pred(
                seqs2predict_list,
                structures_t
            )

            BA_support.remove_one_file_with_try(nr_homo_hits_file)
            del structures_t
            del nr_homo_hits
            del nr_homo_hits_file
            return structures, exec_time, msgs

        elif 'Turbo-fast' == prediction_method:
            if query.annotations['ambiguous']:
                raise exceptions.AmbiguousQuerySequenceException()

            nh = sha1()
            nh.update(str(sorted(pred_method_params.get(prediction_method, {}).items())).encode())
            nh_str = nh.hexdigest()

            structures_t, exec_time = turbofold_fast(
                all_seqs=all_hits_list,
                seqs2predict=seqs2predict_list,
                query=query,
                cpu=threads,
                n=pred_method_params.get(prediction_method, {}).get('max_seqs_in_prediction', 3),
                turbofold_params=pred_method_params.get(prediction_method, {}).get('TurboFold', {}),
                len_diff=pred_method_params.get(prediction_method, {}).get('query_max_len_diff', query_max_len_diff),
                pkey=prediction_method,
                sha1val=nh_str,
            )

            structures = BA_support.rebuild_structures_output_from_pred(
                seqs2predict_list,
                structures_t
            )

            del structures_t
            return structures, exec_time, []

    except exceptions.NoHomologousSequenceException:
        msg = nonhomseqwarn(prediction_method)
        return None, None, [msg]
    except exceptions.AmbiguousQuerySequenceException:
        msgfail = "Query sequence contains ambiguous characters. Can't use {}.".format(prediction_method)
        ml.warning(msgfail)
        return None, None, [msgfail]
    except exceptions.SubprocessException as e:
        msg = "{} can't be used. Error message follows: {} \n{}".format(
            prediction_method,
            str(e),
            e.errors
        )
        ml.error(msg)
        return None, None, [str(e)]
    except Exception as e:
        ml.error("{} can't be used. Error message follows: \n{}.".format(
            prediction_method, str(e))
        )
        return None, None, [str(e)]

    assert False, "Should not reach here (bad prediction method name)."
