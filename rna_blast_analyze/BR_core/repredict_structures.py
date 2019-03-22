import logging
import os
import pickle
import re
import gzip
import dill
import numpy as np
import pandas
import json
from tempfile import mkstemp
from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.BA_support as BA_support
from rna_blast_analyze.BR_core.BA_methods import add_loc_to_description
from rna_blast_analyze.BR_core.BA_support import NoHomologousSequenceException, filter_ambiguous_seqs_from_list, \
    AmbiguousQuerySequenceException
from rna_blast_analyze.BR_core.centroid_homfold import me_centroid_homfold, centroid_homfold_fast
from rna_blast_analyze.BR_core.cmalign import RfamInfo, run_cmfetch, get_cm_model, extract_ref_from_cm, run_cmemit
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.infer_homology import infer_hits_cm
from rna_blast_analyze.BR_core.predict_structures import alifold_refold_prediction, rcoffee_refold_prediction, \
    decouple_homologs_alifold_refold_prediction, rnafold_wrap_for_predict, subopt_fold_query,\
    subopt_fold_alifold, cmmodel_rnafold_c, rfam_subopt_pred
from rna_blast_analyze.BR_core.predict_structures import find_nc_and_remove, check_lonely_bp, IUPACmapping
from rna_blast_analyze.BR_core.turbofold import turbofold_fast, turbofold_with_homologous
from rna_blast_analyze.BR_core.output.htmloutput import write_html_output
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecompute2dict
from rna_blast_analyze.BR_core.filter_blast import filter_by_eval, filter_by_bits
from rna_blast_analyze.BR_core.config import CONFIG

ml = logging.getLogger(__name__)

safe_prediction_method = [
    'rnafold',
    'centroid_homfold',
    'TurboFold_fast',
    'rfam_rnafoldc',
]


def wrapped_ending_with_prediction(
        args_inner, analyzed_hits, all_hits_fasta, query, pred_method=None, method_params=None, used_cm_file=None
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
    if pred_method is None:
        pred_method = args_inner.prediction_method

    if isinstance(pred_method, str):
        pred_method = (pred_method,)

    if method_params is None:
        method_params = args_inner.pred_params

    # ======= filter if needed =======
    # do the filtering based on e-val or bitscore
    # homologous hits still gets used for prediction

    all_hits = analyzed_hits.hits.copy()

    if args_inner.filter_by_eval is not None:
        hits2predict = filter_by_eval(analyzed_hits.hits, *args_inner.filter_by_eval)
        analyzed_hits.hits = hits2predict
    elif args_inner.filter_by_bitscore is not None:
        hits2predict = filter_by_bits(analyzed_hits.hits, *args_inner.filter_by_bitscore)
        analyzed_hits.hits = hits2predict

    new_structures, exe_time = repredict_structures_for_homol_seqs(
        analyzed_hits.query,
        analyzed_hits.hits,
        args_inner.threads,
        prediction_method=pred_method,
        pred_method_params=method_params,
        all_hits=all_hits,
        use_cm_file=used_cm_file,
        )

    if 'default' not in pred_method:
        for i, hit in enumerate(analyzed_hits.hits):
            for key in new_structures.keys():
                assert str(hit.extension.seq) == str(new_structures[key][i].seq)
                hit.extension.annotations['sss'] += [key]

                # expects "predicted" in annotations - for now, if not given, default is True, as not all prediction
                #  methods implement "predicted" in their output
                if new_structures[key][i].annotations.get('predicted', True):
                    hit.extension.letter_annotations[key] = new_structures[key][i].letter_annotations['ss0']

    else:
        # default in pred method
        pass

    if 'default' not in pred_method:
        for hit in analyzed_hits.hits:
            del hit.extension.letter_annotations['ss0']
            hit.extension.annotations['sss'].remove('ss0')
    else:
        for i, hit in enumerate(analyzed_hits.hits):
            # assert str(hit.extension.seq) == str(new_structures[key][i].seq)
            hit.extension.annotations['sss'].remove('ss0')
            hit.extension.annotations['sss'] += ['default']
            hit.extension.letter_annotations['default'] = hit.extension.letter_annotations.pop(
                'ss0'
            )

    # remove uid from file descriptor
    add_loc_to_description(analyzed_hits)

    # write html if requested
    if args_inner.html:
        ml.info('Writing html to {}.'.format(args_inner.html))
        with open(args_inner.html, 'w') as h:
            h.write(write_html_output(analyzed_hits))

    # write csv file if requested
    if args_inner.csv:
        ml.info('Writing csv to {}.'.format(args_inner.csv))
        analyzed_hits.to_csv(args_inner.csv)

    # replace with json
    if args_inner.json:
        ml.info('Writing json to {}.'.format(args_inner.json))
        j_obj = json.dumps(blastsearchrecompute2dict(analyzed_hits), indent=2)
        if getattr(args_inner, 'zip_json', False):
            with open(args_inner.json + '.gz', 'wb') as ff:
                ff.write(
                    gzip.compress(j_obj.encode())
                )
        else:
            with open(args_inner.json, 'w') as ff:
                ff.write(j_obj)

    if args_inner.pandas_dump:
        ml.info('Writing pandas pickle to {}.'.format(args_inner.pandas_dump))
        pandas.to_pickle(analyzed_hits.pandas, args_inner.pandas_dump)

    if args_inner.dump:
        ml.info('Writing dump files base: {}.'.format(args_inner.dump))
        with open(args_inner.dump, 'wb') as pp:
            pickle.dump(analyzed_hits, pp, pickle.HIGHEST_PROTOCOL)

        with open(args_inner.dump + '.time_dump', 'wb') as pp:
            pickle.dump(exe_time, pp, pickle.HIGHEST_PROTOCOL)

        with open(args_inner.dump + '.dill', 'wb') as pp:
            dill.dump(analyzed_hits, pp, pickle.HIGHEST_PROTOCOL)


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

    multiple at minimum (2) sequences are needed for profile alignment for some alignmers
    so this function always return two or more sequences or raises exception

    :param sim_threshold_percent:   seq similarity threshold for homology exclusion
    :param all_hits:                list of hits
    :param query:                   blast query
    :param cmscore_tr:              threshold for homology inclusion in bits
    :param cm_threshold_percent:    threshold for homology inclusion in % of query bits
    :param check_unambiguous:       bool wherther to check unambiguous seqs
    :param len_diff:                threshold for exclusion of hits lq=len(query) lq - diff*lq < len(seq) < lq + diff*lq
    :return:
    """
    ml.debug(fname())
    # i need to leave query, even if with umbiguos basepairs in
    # because it is used as an reference during distance computation and subsequence selection,
    # however i dont need to have all homologous seqs there
    if check_unambiguous:
        all_hits = filter_ambiguous_seqs_from_list(all_hits)

    dist_table, homologous_seqs = _trusted_hits_selection_wrapper(
        all_hits,
        query,
        cmscore_tr,
        cm_threshold_percent,
        len_diff_=len_diff
    )

    if dist_table.size == 0:
        raise NoHomologousSequenceException

    to_include = BA_support.select_sequences_from_similarity_rec(
        dist_table,
        sim_threshold_percent=sim_threshold_percent
    )
    nr_homolog_hits = [homologous_seqs[i] for i in to_include]

    # final checking of nr homologs
    # if sequence is filtered here, it is ambiguos basepair in query
    # removing it is fine if multiple homologous sequences are present
    # the problem will arise when only 1 homologous sequence will remain
    # if we added sequence in previous step, raise exception, else behave like in prev step
    # what if trusted hit is only one?
    if len(nr_homolog_hits) < 2 and not check_unambiguous:
        # warn('Only one sequence is unique under defined sim_threshold_percent (including query)')
        ml.warning(
            'Only one sequence is unique under defined sim_threshold_percent (including query)\n'
            'Adding the most disimilar homologous sequence to non redundant sequences list'
        )
        # dis_hom_index = dist_table.index.get_loc(dist_table[0].idxmin())
        dis_hom_index = dist_table[:, 0].argmin()
        nr_homolog_hits.append(SeqRecord(homologous_seqs[dis_hom_index].seq, id='dummy_seq_01'))
        del dis_hom_index

    elif len(nr_homolog_hits) < 2 and check_unambiguous:
        if len(filter_ambiguous_seqs_from_list(nr_homolog_hits)) == 0:
            # this mean query contain ambiguos bases
            raise NoHomologousSequenceException
        else:
            ml.warning(
                'Only one sequence is unique under defined sim_threshold_percent (including query)\n'
                'Adding the most disimilar homologous sequence to non redundant sequences list'
            )
            # dis_hom_index = dist_table.index.get_loc(dist_table[0].idxmin())
            dis_hom_index = dist_table[:, 0].argmin()
            nr_homolog_hits.append(SeqRecord(homologous_seqs[dis_hom_index].seq, id='dummy_seq_01'))
            del dis_hom_index
        homologous_seqs = filter_ambiguous_seqs_from_list(homologous_seqs)

    elif len(nr_homolog_hits) >= 2 and not check_unambiguous:
        pass

    elif len(nr_homolog_hits) > 2 and check_unambiguous:
        nr_homolog_hits = filter_ambiguous_seqs_from_list(nr_homolog_hits)
        homologous_seqs = filter_ambiguous_seqs_from_list(homologous_seqs)

    elif len(nr_homolog_hits) == 2 and check_unambiguous:
        homologous_seqs = filter_ambiguous_seqs_from_list(homologous_seqs)
        if len(filter_ambiguous_seqs_from_list(nr_homolog_hits)) == 1:
            # this mean that query contains ambiguous base
            raise NoHomologousSequenceException

    else:
        raise Exception()

    fd_h, nr_homo_hits_file = mkstemp(prefix='rba_', suffix='_58', dir=CONFIG.tmpdir)
    with os.fdopen(fd_h, 'w') as f:
        BA_support.write_fasta_from_list_of_seqrecords(f, nr_homolog_hits)

    return nr_homo_hits_file, homologous_seqs


def create_nr_homolog_hits_file_MSA_unsafe(sim_threshold_percent=None, all_hits=None, query=None, cmscore_tr=0.0,
                                           cm_threshold_percent=None, len_diff=0.1):
    """
    create non redundant homologous hits file
    """
    ml.debug(fname())
    dist_table, homologous_seqs = _trusted_hits_selection_wrapper(
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
        BA_support.write_fasta_from_list_of_seqrecords(f, nr_homolog_hits)

    return nr_homo_hits_file, homologous_seqs


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
        print('No sequences from BLAST output infered homologous for structure prediction')
        return np.empty(0), [query_]

    # add query to trusted sequences
    trusted_seqs_query = [query_] + trusted_seqs_

    # make nr list of sequences -> faster alignment
    # better selection
    nr_trusted_seqs_query = BA_support.non_redundant_seqs(trusted_seqs_query)

    # check if the homologous sequence is not exact match as query
    #  (ie taking non redundant set would be only one sequence)
    if len(nr_trusted_seqs_query) == 1:
        print('All sequences infered homologous are exact same as query.')
        return np.empty(0), [query_]

    # select only sequences in some predifined length range to query
    # this is needed for longish ncRNAs
    #   tolerate 10 % length difference?
    ref_len = len(query_)
    nr_len_selected_trusted = [
        seq for seq in nr_trusted_seqs_query if ref_len * (1 - len_diff_) < len(seq) < ref_len * (1 + len_diff_)
    ]

    # this is to control if only one sequence remained after filtering for length difference
    if len(nr_len_selected_trusted) == 1:
        print(
            'No homologous sequence satisfy the length difference condition ({}: {}-{})'.format(
                len_diff_,
                ref_len * (1 - len_diff_),
                ref_len * (1 + len_diff_)
            )
        )
        return np.empty(0), [query_]

    # sanitize seq names (muscle has issues with too long names)
    san_hom_seqs, san_dict = BA_support.sanitize_fasta_names_in_seqrec_list(nr_len_selected_trusted)

    c_fd, trusted_sequence_file_ = mkstemp(prefix='rba_', suffix='_60', dir=CONFIG.tmpdir)
    with os.fdopen(c_fd, 'w') as f:
        BA_support.write_fasta_from_list_of_seqrecords(f, san_hom_seqs)

    align_file = BA_support.run_muscle(trusted_sequence_file_, reorder=True)
    alig = AlignIO.read(align_file, format='clustal')
    distance_calc = DistanceCalculator(model='identity')
    dist_mat = distance_calc.get_distance(alig)
    # rebuild index from sanitized
    orig_index = [san_dict[i] for i in dist_mat.names]
    dist_mat_pd = pandas.DataFrame.from_records(dist_mat.matrix, index=orig_index)
    dist_table_ = (1 - dist_mat_pd.values) * 100

    os.remove(align_file)
    os.remove(trusted_sequence_file_)
    return dist_table_, trusted_seqs_query


def nonhomseqwarn(method_name):
    msg = 'No sequence was infered homologous, need at least 1 for {} type of prediction\n' \
          'Try one of following prediction methods: {}.'.format(
              method_name,
              ', '.join(safe_prediction_method)
          )
    ml.warning(msg)


def annotate_ambiguos_bases(seqlist):
    iupac = IUPACmapping()
    reg = re.compile("[^" + "^".join(iupac.unambiguous) + "]+", re.IGNORECASE)
    for seq in seqlist:
        m = re.search(reg, str(seq.seq))
        if m:
            msg = "Ambiguous base dected in {}, violating base {}, pos {}".format(
                seq.id,
                m.group(),
                m.start()
            )
            ml.warning(msg)
            seq.annotations['ambiguous'] = True
            if 'msgs' not in seq.annotations:
                seq.annotations['msgs'] = []
            seq.annotations['msgs'].append(msg)
        else:
            seq.annotations['ambiguous'] = False

    return seqlist


def repredict_structures_for_homol_seqs(
        query, seqs2predict,
        threads=None,
        prediction_method=tuple('alifold_refold',),
        pred_method_params=None,
        all_hits=None,
        use_cm_file=None,
):
    """Run RNA structure prediction based on chosen method and parameters.
    """
    msg = 'Entering structure prediction..'
    ml.info(msg)
    print(msg)

    default_sim_tr_perc = 90
    default_score_tr = 0.0
    query_max_len_diff = 0.1

    fd, seqs2predict_fasta = mkstemp(prefix='rba_', suffix='_83', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as fah:
        for hit in seqs2predict:
            if len(hit.extension.seq) == 0:
                continue
            fah.write('>{}\n{}\n'.format(
                hit.extension.id,
                str(hit.extension.seq))
            )

    # output are structures for all_hits
    all_hits_list = annotate_ambiguos_bases([i.extension for i in all_hits])
    query = annotate_ambiguos_bases([query])[0]
    seqs2predict_list = annotate_ambiguos_bases([i.extension for i in seqs2predict])

    if not isinstance(pred_method_params, dict):
        raise Exception('prediction method parameters must be python dict')

    structures = dict()
    exec_time = dict()

    if 'default' in prediction_method:
        # do nothing
        pass

    if 'rfam_rnafoldc' in prediction_method:
        pkey = 'rfam_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)

        method_parameters = pred_method_params.get(pkey, {})

        if use_cm_file is None:
            fd, temp_query_file = mkstemp(prefix='rba_', suffix='_61', dir=CONFIG.tmpdir)
            with os.fdopen(fd, 'w') as f:
                f.write('>{}\n{}\n'.format(query.id, str(query.seq)))

            best_model = get_cm_model(temp_query_file, params=method_parameters, threads=threads)

            rfam = RfamInfo()
            single_cm_file = run_cmfetch(rfam.file_path, best_model)
            os.remove(temp_query_file)
        else:
            single_cm_file = use_cm_file

        structures[pkey], exec_time[pkey] = cmmodel_rnafold_c(
            seqs2predict_fasta,
            single_cm_file,
            threads=threads,
            params=pred_method_params.get(pkey, None)
        )

        if use_cm_file is None:
            os.remove(single_cm_file)
        del pkey

    if 'rfam_centroid_homfold' in prediction_method:
        pkey = 'rfam_centroid_homfold'
        print('Running: {}...'.format(pkey))
        ml.info(pkey)
        # run cmscan if needed
        # run cmfetch
        # run cmemit -> homologous seqs
        # run centroid_homfold

        method_parameters = pred_method_params.get(pkey, {})

        if use_cm_file is None:
            fd, temp_query_file = mkstemp(prefix='rba_', suffix='_61', dir=CONFIG.tmpdir)
            with os.fdopen(fd, 'w') as f:
                f.write('>{}\n{}\n'.format(query.id, str(query.seq)))

            best_model = get_cm_model(temp_query_file, params=method_parameters, threads=threads)

            rfam = RfamInfo()
            single_cm_file = run_cmfetch(rfam.file_path, best_model)
            os.remove(temp_query_file)
        else:
            single_cm_file = use_cm_file

        cep = method_parameters.get('cmemit', '')
        if '-u' not in cep:
            cep += ' -u'
        if '-N' not in cep:
            cep += ' -N {}'.format(method_parameters.get('n_seqs', 10))

        hf_file = run_cmemit(single_cm_file, params=cep)

        structures[pkey], exec_time[pkey] = me_centroid_homfold(seqs2predict_fasta, hf_file, params=method_parameters)

        if use_cm_file is None:
            os.remove(single_cm_file)
        os.remove(hf_file)
        del pkey

    if 'rfam_subopt' in prediction_method:
        pkey = 'rfam_subopt'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        if use_cm_file is None:
            fd, temp_query_file = mkstemp(prefix='rba_', suffix='_62', dir=CONFIG.tmpdir)
            with os.fdopen(fd, 'w') as f:
                f.write('>{}\n{}\n'.format(query.id, str(query.seq)))

            if pkey in pred_method_params and pred_method_params[pkey]:
                best_model = get_cm_model(temp_query_file, params=pred_method_params[pkey], threads=threads)
            else:
                best_model = get_cm_model(temp_query_file, threads=threads)

            rfam = RfamInfo()
            single_cm_file = run_cmfetch(rfam.file_path, best_model)
            os.remove(temp_query_file)
        else:
            single_cm_file = use_cm_file

        ref_structure = extract_ref_from_cm(single_cm_file)

        structures[pkey], exec_time[pkey] = rfam_subopt_pred(
            seqs2predict_fasta,
            ref_structure,
            params=pred_method_params.get(pkey, None)
        )

        if use_cm_file is None:
            os.remove(single_cm_file)

        del pkey

    if 'rnafold' in prediction_method:
        pkey = 'rnafold'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        structures[pkey], exec_time[pkey] = rnafold_wrap_for_predict(
            seqs2predict_fasta,
            params=pred_method_params.get(pkey, {}).get('RNAfold', '')
        )
        del pkey

    if 'subopt_fold_query' in prediction_method:
        pkey = 'subopt_fold_query'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)

        a, qf = mkstemp(prefix='rba_', suffix='_55', dir=CONFIG.tmpdir)
        with os.fdopen(a, 'w') as fd:
            fd.write('>query\n{}\n'.format(str(query.seq)))

        structures[pkey], exec_time[pkey] = subopt_fold_query(
            seqs2predict_fasta,
            qf,
            params=pred_method_params.get(pkey, None)
        )

        os.remove(qf)

        del pkey

    if 'subopt_fold_clustal_alifold' in prediction_method:
        pkey = 'subopt_fold_clustal_alifold'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file

            f, homologous_sequence_file = mkstemp(prefix='rba_', suffix='_64', dir=CONFIG.tmpdir)
            with os.fdopen(f, 'w') as fh:
                BA_support.write_fasta_from_list_of_seqrecords(fh, homologous_seqs)

            structures[pkey], exec_time[pkey] = subopt_fold_alifold(
                seqs2predict_fasta,
                homologous_sequence_file,
                aligner='clustalo',
                params=method_parameters,
                threads=threads
            )

            os.remove(homologous_sequence_file)
            del homologous_sequence_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'subopt_fold_muscle_alifold' in prediction_method:
        pkey = 'subopt_fold_muscle_alifold'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file

            f, homologous_sequence_file = mkstemp(prefix='rba_', suffix='_65', dir=CONFIG.tmpdir)
            with os.fdopen(f, 'w') as fh:
                BA_support.write_fasta_from_list_of_seqrecords(fh, homologous_seqs)

            structures[pkey], exec_time[pkey] = subopt_fold_alifold(
                seqs2predict_fasta,
                homologous_sequence_file,
                aligner='muscle',
                params=method_parameters,
                threads=threads,
            )

            os.remove(homologous_sequence_file)
            del homologous_sequence_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'clustalo_alifold_refold_rnafoldc' in prediction_method:
        pkey = 'clustalo_alifold_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )

            structures[pkey], exec_time[pkey] = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='clustalo'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'muscle_alifold_refold_rnafoldc' in prediction_method:
        pkey = 'muscle_alifold_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='muscle'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'clustalo_alifold_unpaired_conserved_refold_rnafoldc' in prediction_method:
        pkey = 'clustalo_alifold_unpaired_conserved_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='clustalo'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'muscle_alifold_unpaired_conserved_refold_rnafoldc' in prediction_method:
        pkey = 'muscle_alifold_unpaired_conserved_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = alifold_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                msa_alg='muscle'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'dh_rcoffee_alifold_refold_rnafoldc' in prediction_method:
        pkey = 'dh_rcoffee_alifold_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = decouple_homologs_alifold_refold_prediction(
                nr_homo_hits_file,
                homologous_seqs,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                align='rcoffee'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'dh_rcoffee_alifold_unpaired_conserved_rnafoldc' in prediction_method:
        pkey = 'dh_rcoffee_alifold_unpaired_conserved_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = decouple_homologs_alifold_refold_prediction(
                nr_homo_hits_file,
                homologous_seqs,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                align='rcoffee',
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'dh_clustal_alifold_refold_rnafoldc' in prediction_method:
        pkey = 'dh_clustal_alifold_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = decouple_homologs_alifold_refold_prediction(
                nr_homo_hits_file,
                homologous_seqs,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters,
                align='clustalo'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'dh_clustal_alifold_unpaired_conserved_rnafoldc' in prediction_method:
        pkey = 'dh_clustal_alifold_unpaired_conserved_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, homologous_seqs = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = decouple_homologs_alifold_refold_prediction(
                nr_homo_hits_file,
                homologous_seqs,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters,
                align='clustalo'
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
            del homologous_seqs
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'centroid_homfold' in prediction_method:
        pkey = 'centroid_homfold'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)

        method_parameters = pred_method_params.get(pkey, {})

        nr_homo_hits_file, _ = create_nr_homolog_hits_file_MSA_unsafe(
            all_hits=all_hits_list,
            query=query,
            sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
            cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
            cm_threshold_percent=method_parameters.get('cmscore_percent', None),
            len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
        )

        raw_structures, exec_time[pkey] = me_centroid_homfold(
            seqs2predict_fasta, nr_homo_hits_file,
            params=method_parameters
        )

        # check noncanonical
        if pkey in pred_method_params and pred_method_params[pkey]:
            allow_nc = pred_method_params[pkey].get('allow_noncanonical', False)
            allow_lp = pred_method_params[pkey].get('allow_lonely_pairs', False)
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

        structures[pkey] = raw_structures
        os.remove(nr_homo_hits_file)
        del nr_homo_hits_file
        del pkey

    if 'centroid_homfold_fast' in prediction_method:
        pkey = 'centroid_homfold_fast'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)

        method_parameters = pred_method_params.get(pkey, {})

        raw_structures, exec_time[pkey] = centroid_homfold_fast(
            all_seqs=all_hits_list,
            query=query,
            all_seqs_fasta=seqs2predict_fasta,
            n=method_parameters.get('max_seqs_in_prediction', 10),
            centroid_homfold_params=method_parameters,
            len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff)
        )

        # check noncanonical
        if pkey in pred_method_params and pred_method_params[pkey]:
            allow_nc = pred_method_params[pkey].get('allow_noncanonical', False)
            allow_lp = pred_method_params[pkey].get('allow_lonely_pairs', False)
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

        structures[pkey] = raw_structures
        del pkey

    if 'TurboFold' in prediction_method:
        pkey = 'TurboFold'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        # set arbitrary sim_threshold_percent to 100, because we want to remove only identical sequences from prediction
        #  with Trurbofold. The structure of redundant sequences will be set according to the one in prediction

        try:
            all_hits_filtered = filter_ambiguous_seqs_from_list(all_hits_list)
            seqs2predict_filtered = filter_ambiguous_seqs_from_list(seqs2predict_list)

            if query.annotations['ambiguous']:
                msgfail = "Query sequence contains ambiguous characters. Can't use Turbofold."
                ml.error(msgfail)
                raise AmbiguousQuerySequenceException(msgfail)

            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_homolog_hits_file_MSA_unsafe(
                all_hits=all_hits_filtered,
                query=query,
                sim_threshold_percent=100,
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
                )

            with open(nr_homo_hits_file, 'r') as nrf:
                nr_homo_hits = [seq for seq in SeqIO.parse(nrf, format='fasta')]

            structures_t, exec_time[pkey] = turbofold_with_homologous(
                all_sequences=seqs2predict_filtered,
                nr_homologous=nr_homo_hits,
                params=method_parameters.get('TurboFold', {}),
                n=method_parameters.get('max_seqs_in_prediction', 3),
                cpu=threads
            )

            structures[pkey] = BA_support.rebuild_structures_output_from_pred(
                seqs2predict_list,
                structures_t
            )

            os.remove(nr_homo_hits_file)
            del structures_t
            del nr_homo_hits
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)

        except AmbiguousQuerySequenceException:
            pass

        finally:
            del pkey

    if 'TurboFold_fast' in prediction_method:
        pkey = 'TurboFold_fast'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)

        try:
            if query.annotations['ambiguous']:
                msgfail = "Query sequence contains ambiguous characters. Can't use Turbofold."
                ml.error(msgfail)
                raise AmbiguousQuerySequenceException(msgfail)

            structures_t, exec_time[pkey] = turbofold_fast(
                all_seqs=all_hits_list,
                seqs2predict=seqs2predict_list,
                query=query,
                cpu=threads,
                n=pred_method_params[pkey].get('max_seqs_in_prediction', 3),
                turbofold_params=pred_method_params.get(pkey, {}).get('TurboFold', {}),
                len_diff=pred_method_params[pkey].get('query_max_len_diff', query_max_len_diff)
            )

            structures[pkey] = BA_support.rebuild_structures_output_from_pred(
                seqs2predict_list,
                structures_t
            )

            del structures_t
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)

        except AmbiguousQuerySequenceException:
            pass

        finally:
            del pkey

    if 'rcoffee_alifold_refold_rnafoldc' in prediction_method:
        pkey = 'rcoffee_alifold_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = rcoffee_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='refold_rnafoldc',
                threads=threads,
                params=method_parameters
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    if 'rcoffee_alifold_unpaired_conserved_refold_rnafoldc' in prediction_method:
        pkey = 'rcoffee_alifold_unpaired_conserved_refold_rnafoldc'
        print('Runing: {}...'.format(pkey))
        ml.info(pkey)
        try:
            method_parameters = pred_method_params.get(pkey, {})

            nr_homo_hits_file, _ = create_nr_trusted_hits_file_MSA_safe(
                all_hits=all_hits_list,
                query=query,
                sim_threshold_percent=method_parameters.get('pred_sim_threshold', default_sim_tr_perc),
                cmscore_tr=method_parameters.get('cmscore_tr', default_score_tr),
                cm_threshold_percent=method_parameters.get('cmscore_percent', None),
                len_diff=method_parameters.get('query_max_len_diff', query_max_len_diff),
            )
            structures[pkey], exec_time[pkey] = rcoffee_refold_prediction(
                nr_homo_hits_file,
                seqs2predict_fasta,
                refold='conserved_ss_rnafoldc',
                threads=threads,
                params=method_parameters
            )

            os.remove(nr_homo_hits_file)
            del nr_homo_hits_file
        except NoHomologousSequenceException:
            nonhomseqwarn(pkey)
        finally:
            del pkey

    return structures, exec_time
