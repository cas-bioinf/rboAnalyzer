import os
import re
from multiprocessing import Pool
from subprocess import call
from tempfile import mkstemp, TemporaryFile
import logging
import shlex

from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.BA_support as BA_support
import rna_blast_analyze.BR_core.extend_hits
from rna_blast_analyze.BR_core.alifold4all import compute_refold
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.infer_homology import infer_homology
from rna_blast_analyze.BR_core.locarna_clustal_like_2stockholm import parse_locarna_alignment
from rna_blast_analyze.BR_core.stockholm_parser import trim_alignment_by_sequence
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core import exceptions

ml = logging.getLogger('rboAnalyzer')


def write_clustal_like_file_with_anchors(fid, seq_name, cseq, my_anchors):
    fid.write('CLUSTAL W --- Clustal like format for locarna\n')
    fid.write('{} {}\n'.format(seq_name, cseq))
    for anch in my_anchors:
        fid.write('{} {}\n'.format(anch[0], anch[1]))


def locarna_worker(pack):
    ml.debug(fname())
    one_expanded_hit, query_seq, locarna_params, anchor_length = pack

    locarna_file1 = locarna_file2 = loc_out_file = None

    try:
        # read the aligned segment and use it as anchors for locarna
        # run locarna in local mode and put the query sequence with the extended sequence with the blast aligned
        # segment as anchor
        blast_entry = one_expanded_hit.annotations['blast'][1]

        anchors = LocarnaAnchor(
            to_rna(blast_entry.query),
            blast_entry.match,
            to_rna(blast_entry.sbjct),
            anchor_length=anchor_length
        )
        # extracted temp is my query

        # access the locarna aligner directly
        fd1, locarna_file1 = mkstemp(prefix='rba_', suffix='_20', dir=CONFIG.tmpdir)
        with os.fdopen(fd1, 'w') as fp_locarna_file_1:
            ql1, ql2 = anchors.anchor_whole_seq(str(query_seq), 'query')
            write_clustal_like_file_with_anchors(fp_locarna_file_1,
                                                 'query',
                                                 str(query_seq),
                                                 (
                                                     ('#A1', ql1.split()[0]),
                                                     ('#A2', ql2.split()[0])
                                                 ))

        fd2, locarna_file2 = mkstemp(prefix='rba_', suffix='_21', dir=CONFIG.tmpdir)
        with os.fdopen(fd2, 'w') as fp_locarna_file_2:
            sl1, sl2 = anchors.anchor_whole_seq(str(one_expanded_hit.seq), 'subject')
            write_clustal_like_file_with_anchors(fp_locarna_file_2,
                                                 one_expanded_hit.id,
                                                 str(one_expanded_hit.seq),
                                                 (
                                                     ('#A1', sl1.split()[0]),
                                                     ('#A2', sl2.split()[0])
                                                 ))

        loc_out_file = run_locarna(
            locarna_file1,
            locarna_file2,
            locarna_params
        )

        # read locarna alignment
        with open(loc_out_file, 'r') as f:
            locarna_alig = parse_locarna_alignment(f)

        if len(locarna_alig) != 2:
            raise exceptions.SubseqMatchError('There must be 2 sequences in Locarna alignment.')

        loc_rep = create_report_object_from_locarna(one_expanded_hit, locarna_alig)

        return loc_rep
    except exceptions.LocarnaException as e:
        one_expanded_hit.annotations['msgs'] = [str(e), e.errors]
        empty_hit = BA_support.Subsequences(one_expanded_hit)
        return empty_hit
    except (exceptions.SubseqMatchError, exceptions.ParsingError) as e:
        one_expanded_hit.annotations['msgs'] = [str(e)]
        empty_hit = BA_support.Subsequences(one_expanded_hit)
        return empty_hit
    except (TypeError, AttributeError, FileNotFoundError) as e:
        one_expanded_hit.annotations['msgs'] = [str(e)]
        empty_hit = BA_support.Subsequences(one_expanded_hit)
        return empty_hit
    finally:
        for f in [locarna_file1, locarna_file2, loc_out_file]:
            if f is not None:
                BA_support.remove_one_file_with_try(f)


def extend_locarna_core(analyzed_hits, query, args_inner, all_short, multi_query, iteration, ih_model):
    # expand hits according to query + 10 nucleotides +-
    if args_inner.db_type == "blastdb":
        shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits(
            all_short,
            args_inner.blast_db,
            len(query),
            extra=args_inner.subseq_window_locarna,
            blast_regexp=args_inner.blast_regexp,
            skip_missing=args_inner.skip_missing,
            msgs=analyzed_hits.msgs,
        )
    elif args_inner.db_type in ["fasta", "gb", "server"]:
        shorts_expanded, _ = rna_blast_analyze.BR_core.extend_hits.expand_hits_from_fasta(
            all_short,
            args_inner.blast_db,
            len(query),
            extra=args_inner.subseq_window_locarna,
            blast_regexp=args_inner.blast_regexp,
            skip_missing=args_inner.skip_missing,
            msgs=analyzed_hits.msgs,
            format=args_inner.db_type,
        )
    else:
        raise exceptions.IncorrectDatabaseChoice()

    shorts_expanded = BA_support.rc_hits_2_rna(shorts_expanded)

    query_seq = query.seq.transcribe()

    # compute alignment here

    if args_inner.threads == 1:
        result = []
        for oeh in shorts_expanded:
            result.append(
                locarna_worker(
                    (
                        oeh,
                        query_seq,
                        args_inner.locarna_params,
                        args_inner.locarna_anchor_length
                    )
                )
            )
    else:
        pack = []
        for oeh in shorts_expanded:
            pack.append(
                (
                    oeh,
                    query_seq,
                    args_inner.locarna_params,
                    args_inner.locarna_anchor_length
                )
            )
        pool = Pool(processes=args_inner.threads)
        result = pool.map(locarna_worker, pack)
        pool.close()

    for res in result:
        if res.extension is None:
            analyzed_hits.hits_failed.append(res)
        else:
            analyzed_hits.hits.append(res)

    # this part predicts homology - it is not truly part of repredict
    homology_prediction, homol_seqs, cm_file_rfam_user = infer_homology(
        analyzed_hits=analyzed_hits, args=args_inner, cm_model_file=ih_model, multi_query=multi_query,
        iteration=iteration
    )
    # add homology prediction to the data
    for hit, pred in zip(analyzed_hits.hits, homology_prediction):
        hit.hpred = pred
    return analyzed_hits, homology_prediction, homol_seqs, cm_file_rfam_user


def _select_refold_structure(refold_structures_, exp_hit_id):
    for seq_refold in refold_structures_:
        # seq_refold can be different due too pipe through CLUSTAL
        if seq_refold.id in exp_hit_id:
            return seq_refold


def create_report_object_from_locarna(exp_hit, locarna_alig):
    """
    create object which will be appended to BlastSearchRecompute class
    This needs to be Subsequences object

    :param exp_hit:
    :param locarna_alig:
    :return:
    """
    ml.debug(fname())
    # chop alignment by seq
    query_ind = [i for i, j in enumerate(locarna_alig) if j.id == 'query']
    if len(query_ind) != 1:
        raise exceptions.SubseqMatchError('Got multiple hits with id "query" in the Locarna alignment.')
    trimmed_locarna_alig = trim_alignment_by_sequence(
        locarna_alig,
        str(locarna_alig[query_ind[0]].seq),
        structure_annotation='SS_cons'
    )

    aligned_subsequence = BA_support.select_analyzed_aligned_hit(trimmed_locarna_alig, exp_hit.id)

    # add annotations from exp hit
    aligned_subsequence.annotations = exp_hit.annotations
    aligned_subsequence.name = exp_hit.name

    # also add annotations from locarna, mainly score
    aligned_subsequence.annotations.update(locarna_alig.annotations)

    # get the structure
    # by refold
    refold_structures = refold_stockholm(trimmed_locarna_alig, trimmed_locarna_alig.column_annotations['SS_cons'])

    # select refold structure for my seq
    seq_refold_structure = _select_refold_structure(refold_structures, exp_hit.id)

    aligned_subsequence.letter_annotations['ss0'] = seq_refold_structure.letter_annotations['ss0']
    aligned_subsequence.annotations['sss'] = ['ss0']

    # prepare seq_record for subsequences
    aligned_subsequence.description = ''
    hit = BA_support.Subsequences(exp_hit)

    hit.extension = aligned_subsequence

    # find the matching sequence
    pos_match = re.search(str(aligned_subsequence.seq), str(exp_hit.seq), flags=re.IGNORECASE)
    if not pos_match:
        raise exceptions.SubseqMatchError(
            'Aligned portion of subject sequence in Locarna alignment was not found in parent sequence.'
        )

    hit.best_start, hit.best_end = compute_true_location_locarna(hit, pos_match)

    return hit


def compute_true_location_locarna(hit, match):
    start, end = match.span()
    if hit.source.annotations['blast'][1].strand == 1:
        if hit.source.annotations['trimmed_ss']:
            s = start + 1
        else:
            s = hit.source.annotations['super_start'] + start

        e = end - start + s + 1
    else:
        if hit.source.annotations['trimmed_ss']:
            s = start + 1
        else:
            s = hit.source.annotations['super_start'] + start

        e = end - start + s - 1

    return s, e


def to_rna(seq):
    rna_seq = re.sub('t', 'u', seq, flags=re.IGNORECASE)
    return rna_seq.upper()


def run_locarna(query_file, subject_file, locarna_params):
    """
    possible settings:
    struct_local : seq_local
    1 : 1   - breaks sequences and joins them !!!! DO NOT USE
    1 : 0   - tries to align globally, end to end - lot of gaps
    0 : 1   - wanted alignment - unaligned ends of input sequences are trimmed in output alignment
    0 : 0   - wanted alignment - returns whole alignment including unaligned ends

    :param query_file:
    :param subject_file:
    :param locarna_params:
    :return:
    """
    ml.debug('running locarna')
    if not os.path.isfile(query_file):
        raise FileNotFoundError('Provided file {} was not found'.format(query_file))
    if not os.path.isfile(subject_file):
        raise FileNotFoundError('Provided file {} was not found'.format(subject_file))

    cmd = '{} {} {} {} > {}'.format(
        shlex.quote('{}locarna'.format(CONFIG.locarna_path)),
        ' '.join([shlex.quote(i) for i in shlex.split(locarna_params)]),
        shlex.quote(query_file),
        shlex.quote(subject_file),
        shlex.quote(subject_file + '.loc_out'),
    )
    ml.debug(cmd)
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = call(cmd, shell=True, stdout=tmp, stderr=tmp)
        if r:
            msgfail = 'Call to locarna failed for files in1:{} in2:{} out:{}'.format(
                query_file, subject_file, subject_file + '.loc_out'
            )
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.LocarnaException(msgfail, tmp.read())

        return subject_file + '.loc_out'


def write_locarna_anchors_with_min_length(match_line, min_anchor_length=1):
    ml.debug(fname())
    h1 = []
    h2 = []
    pa = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    part_desig = 0

    for match in re.finditer(r'\|+', match_line, flags=re.IGNORECASE):
        if len(match.group()) < min_anchor_length:
            # skip the iterations below minimum length
            continue

        for i in range(match.span()[0] - len(h1)):
            h1.append('.')
            h2.append('.')

        c = 0
        part_desig += 1
        for _ in match.group():
            if c == 9:
                c = 0
                part_desig += 1
            c += 1

            h1.append(pa[part_desig])
            h2.append(str(c))

    for i in range(len(match_line) - len(h1)):
        h1.append('.')
        h2.append('.')

    anchor_l1 = ''.join(h1)
    anchor_l2 = ''.join(h2)
    return anchor_l1, anchor_l2


def squeeze_locarna_anchors_to_aligned_seq(aligned_seq, anchor_line1, anchor_line2):
    ml.debug(fname())
    out_seq = []
    out_al1 = []
    out_al2 = []
    for seq_pos, al1, al2 in zip(aligned_seq, anchor_line1, anchor_line2):
        if seq_pos == '-':
            continue
        out_seq.append(seq_pos)
        out_al1.append(al1)
        out_al2.append(al2)

    squeezed_seq = ''.join(out_seq)
    squeezed_al1 = ''.join(out_al1)
    squeezed_al2 = ''.join(out_al2)

    return squeezed_seq, squeezed_al1, squeezed_al2


class LocarnaAnchor(object):
    """
    while initiating LocarnaAnchor object U can specify minimal anchor length to be used
    If default (-1) is kept, then minimal anchor length for succesfull usage for locarna is infered
    and the number is returned in anchor_length parameter
    """
    def __init__(self, query, match, subject, anchor_length=-1):
        self.match = match
        self.query = query
        self.subject = subject
        # self.anchor_l1, self.anchor_l2 = write_locarna_anchors(self.match)
        # compute anchor length

        self.anchor_length = anchor_length
        if anchor_length < 0:
            while True:
                self.anchor_l1, self.anchor_l2 = write_locarna_anchors_with_min_length(self.match, self.anchor_length)
                if '[' in self.anchor_l1:
                    self.anchor_length += 1
                else:
                    break
        else:
            self.anchor_l1, self.anchor_l2 = write_locarna_anchors_with_min_length(self.match, self.anchor_length)

        assert len(self.anchor_l1) == len(self.anchor_l2) == len(self.query) == len(self.subject)

        if anchor_length < 0:
            print('inferred anchor length {}'.format(self.anchor_length))

        self.squeezed_query, self.q_al1, self.q_al2 = squeeze_locarna_anchors_to_aligned_seq(self.query,
                                                                                             self.anchor_l1,
                                                                                             self.anchor_l2)
        self.squeezed_subject, self.s_al1, self.s_al2 = squeeze_locarna_anchors_to_aligned_seq(self.subject,
                                                                                               self.anchor_l1,
                                                                                               self.anchor_l2)

    def anchor_whole_seq(self, seq, seq_line):
        """
        write anchors to whole sequence
        seq input must be appropriate sequence for the seq_line chosen
        :param seq: str or SeqRecord object
        :param seq_line: str - 'query' | 'subject'
        :return:
        """
        if not isinstance(seq, str) | isinstance(seq, SeqRecord):
            raise TypeError("types 'str' and 'SeqRecord' accepted, but {} given.".format(type(seq)))

        if not (seq_line == 'query') | (seq_line == 'subject'):
            raise AttributeError("Expected 'query' or 'subject', but {} given.".format(seq_line))

        if isinstance(seq, SeqRecord):
            sequence = str(seq.seq)
        else:
            sequence = seq

        # find matching sequence in seq
        if seq_line == 'query':
            match_o = re.search(self.squeezed_query, sequence, flags=re.IGNORECASE)
            if not match_o:
                raise exceptions.SubseqMatchError(
                    'Could not find match between provided query sequence and QUERY sequence from BLAST output.'
                )

            start = ''.join(['.' for i in range(match_o.span()[0])])
            end = ''.join(['.' for i in range(match_o.span()[1], len(sequence))])
            l1 = start + self.q_al1 + end + ' #1'
            l2 = start + self.q_al2 + end + ' #2'
        else:
            match_o = re.search(self.squeezed_subject, sequence, flags=re.IGNORECASE)
            if not match_o:
                raise exceptions.SubseqMatchError(
                    'Could not find match between subject sequence from the database and subject sequence from BLAST '
                    'output. This can be DATABASE problem. Try to recreate the database and/or check if the accessions'
                    ' in the BLAST output correspond with accessions in the provided database (i.e.'
                    ' the sequences are the same ones).'
                )

            start = ''.join(['.' for i in range(match_o.span()[0])])
            end = ''.join(['.' for i in range(match_o.span()[1], len(sequence))])
            l1 = start + self.s_al1 + end + ' #1'
            l2 = start + self.s_al2 + end + ' #2'

        return l1, l2


def refold_stockholm(stockholm_alig, consensus_structure):
    """
    compute refold.pl from Vienna RNA package
    :param stockholm_alig:
    :param consensus_structure:
    :return:
    """
    ml.debug(fname())
    # convert to clustal alignment
    fd, clust_tempfile = mkstemp(prefix='rba_', suffix='_23', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as f:
        stockholm_alig.write_clustal(f)

    # write fake alifold output with given consensus structure
    fd, alif_fake_file = mkstemp(prefix='rba_', suffix='_24', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as f:
        # the consensus sequence in alifold file is really not used for anything
        f.write('A'*len(consensus_structure) + '\n')
        f.write(consensus_structure + '\n')

    # compute refold
    # refold_path = locate_refold()
    refold_constrained_file = compute_refold(clust_tempfile, alif_fake_file)

    parsed_seqs = []
    with open(refold_constrained_file, 'r') as f:
        # read the file
        for seq in BA_support.parse_seq_str(f):
            parsed_seqs.append(seq)

    # cleanup
    BA_support.remove_files_with_try([clust_tempfile, alif_fake_file, refold_constrained_file])

    return parsed_seqs
