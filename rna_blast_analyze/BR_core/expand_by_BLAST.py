import sys
from copy import deepcopy
import logging

import rna_blast_analyze.BR_core.BA_support as BA_support
import rna_blast_analyze.BR_core.extend_hits
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core import exceptions
import rna_blast_analyze.BR_core.luncher

ml = logging.getLogger('rboAnalyzer')


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

        # handle the case when e is not < 0, and indexing is not from the end
        if not e < 0:
            nseq = seq[s:]
        else:
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

    bl = ns.annotations['blast'][1]

    if bl.sbjct_start < bl.sbjct_end:
        bls = bl.sbjct_start - bl.query_start + 1
        ble = bl.sbjct_end + (query_len - bl.query_end)
    elif bl.sbjct_end < bl.sbjct_start:
        bls = bl.sbjct_end - (query_len - bl.query_end)
        ble = bl.sbjct_start + bl.query_start - 1
    else:
        raise exceptions.UnknownStrand("Can't determine HSP strand (sbjct_start appears equal to sbjct_end)")

    # if whole subject sequence too short, this assertion will fail
    if tss or tse or tes or tee:
        msg = 'STATUS: Skipping sequence check ({}) - subject sequence too short.'.format(ns.id)
        ml.info(msg)
    else:
        assert len(ns.seq) == abs(bls - ble) + 1
        assert bls == ns.annotations['extended_start']
        assert ble == ns.annotations['extended_end']

    hit.best_start, hit.best_end = compute_true_location_se(hit, query_len)

    return hit


def extend_simple_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model):
    # the extra here is given "pro forma" the sequence is extended exactly by lenghts of unaligned portions of query
    if args_inner.db_type == "blastdb":
        shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits(
            all_short,
            args_inner.blast_db,
            len(query),
            extra=0,
            blast_regexp=args_inner.blast_regexp,
            skip_missing=args_inner.skip_missing,
            msgs=analyzed_hits.msgs,
        )
    elif args_inner.db_type in ["fasta", "gb", "server"]:
        shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits_from_fasta(
            all_short,
            args_inner.blast_db,
            len(query),
            extra=0,
            blast_regexp=args_inner.blast_regexp,
            skip_missing=args_inner.skip_missing,
            msgs=analyzed_hits.msgs,
            format=args_inner.db_type,
        )
    else:
        raise exceptions.IncorrectDatabaseChoice()

    # check, if blast hits are non - overlapping, if so, add the overlapping hit info to the longer hit
    # reflect this in user output
    # shorts_expanded = merge_blast_hits(shorts_expanded)

    shorts_expanded = trim_before(shorts_expanded)

    shorts_expanded = BA_support.rc_hits_2_rna(shorts_expanded)

    query_seq = query.seq.transcribe()

    # blast only extension
    for exp_hit in shorts_expanded:
        try:
            _out = create_blast_only_report_object(exp_hit, len(query_seq))
            analyzed_hits.hits.append(_out)
        except AssertionError as e:
            exp_hit.annotations['msgs'] += [str(e)]
            analyzed_hits.hits_failed.append(exp_hit)
        except exceptions.UnknownStrand as e:
            exp_hit.annotations['msgs'] += [str(e)]
            analyzed_hits.hits_failed.append(exp_hit)
        except Exception as e:
            ml.error("Unexpected error when extending with 'simple'.")
            exp_hit.annotations['msgs'] += [str(e)]

    if len(analyzed_hits.hits) == 0:
        ml.error(
            "Extension failed for all sequences. Please see the error message. You can also try '--mode locarna'."
        )
        sys.exit(1)

    # assign Locarna score to None as it is not directly accessible from mlocarna
    for hit in analyzed_hits.hits:
        hit.extension.annotations['score'] = None

    # this part predicts homology - it is not truly part of repredict
    homology_prediction, homol_seqs, cm_file_rfam_user = infer_homology(
        analyzed_hits=analyzed_hits, args=args_inner, cm_model_file=ih_model, multi_query=multi_query,
        iteration=iteration
    )
    for hit, pred in zip(analyzed_hits.hits, homology_prediction):
        hit.hpred = pred
    return analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user
