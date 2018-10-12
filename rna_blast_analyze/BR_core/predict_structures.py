import copy
import itertools
import os
import re
import time
import unicodedata
import logging
from shutil import rmtree
from subprocess import call
from tempfile import mkstemp, mkdtemp

import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.BA_support import devprint, remove_files_with_try, NoHomologousSequenceException
from rna_blast_analyze.BR_core.BA_support import read_seq_str, write_fasta_from_list_of_seqrecords,\
    sanitize_fasta_names_in_seqrec_list, desanitize_fasta_names_in_seqrec_list, run_hybrid_ss_min,\
    run_muscle, ct2db
from rna_blast_analyze.BR_core.RNAshapes import rapid_shapes_list, read_rapid_shapes_list_outfile, structures2shape
from rna_blast_analyze.BR_core.alifold4all import compute_refold, compute_clustalo_clasic, compute_alifold,\
    compute_constrained_prediction
from rna_blast_analyze.BR_core.cmalign import build_stockholm_from_clustal_alig, extract_ref_structure_fromRFAM_CM, run_cmalign_on_fasta, cm_strucutre2br, \
    get_cm_model
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.db2shape import nesting
from rna_blast_analyze.BR_core.infer_homology import _alignment_column_conservation
from rna_blast_analyze.BR_core.par_distance_no_RNAlib import distances_one_thread as rnadistance_one_thread
from rna_blast_analyze.BR_core.stockholm_parser import read_st, trim_cmalign_sequence_by_refseq_one_seq
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)


def _parse_first_record_only(file):
    with open(file, 'r') as f:
        for seq in SeqIO.parse(f, format='fasta'):
            return seq


def _repair_consensus_structure_by_maping(consensus_str, mapping, expected_str_len, gap_char=49):
    # provided lengths must be always the same
    # add . in added regions
    new_consensus = chr(gap_char) * expected_str_len
    # add_origin
    for i, ((j, k), (l, m)) in enumerate(mapping):
        new_consensus = new_consensus[:l] + consensus_str[j:k] + new_consensus[m:]

    if len(new_consensus) != expected_str_len:
        raise Exception('mapping consensus structure to alignment failed')

    return new_consensus


def _map_alignment_columns_from_profile_match(original_match, new_match):
    # asumption 1
    #  new match will have more gaps
    # asumption 2
    #  main introduced gaps are preserved in alignment
    # reality
    #  alignment is organized differently
    # t-coffe has similar functionality but alignment to profile is also shifted
    # different approach - generate mapping by column - it should help to
    # look into cmalign optional output, they deel there with stuff like this
    oseq = str(original_match.seq).upper()
    nseq = str(new_match.seq).upper()

    s = 0
    mapping = []
    # find_seqs
    for match in re.finditer('[A-Za-z]+', nseq):
        om = re.search(match.group(), oseq[s:])
        mapping.append([[i + s for i in om.span()], list(match.span())])
        s += om.span()[1]

    return mapping


def _transform_score_to_pplike_line(score):
    ms = max(score)
    avg_score = [i / ms for i in score]
    l = []
    for s in avg_score:
        if s <= 0.05:
            l.append('0')
        elif 0.05 < s <= 0.15:
            l.append('1')
        elif 0.15 < s <= 0.25:
            l.append('2')
        elif 0.25 < s <= 0.35:
            l.append('3')
        elif 0.35 < s <= 0.45:
            l.append('4')
        elif 0.45 < s <= 0.55:
            l.append('5')
        elif 0.55 < s <= 0.65:
            l.append('6')
        elif 0.65 < s <= 0.75:
            l.append('7')
        elif 0.75 < s <= 0.85:
            l.append('8')
        elif 0.85 < s <= 0.95:
            l.append('9')
        elif 0.95 < s <= 1:
            l.append('*')
        else:
            raise Exception('score error')
    return ''.join(l)


def _refold_with_unpaired_conservation(stockholm_msa, repred_tr='8', conseq_conserved=1):
    """
    get trusted msa with consensus and develop a conservation scoring (based on my column score?)
    output the scoring and output gappless individual sequences together with a scoring and a consensus sequence
    take that and use the scoring
    as a score the PP line in cmalign can be probably used

    # can i use scoring as an SHAPE reactivity constrain somehow?
    # probably yee but on what basis? The SHAPE reactivity is reactivity for nucleotides in very specific conformations
    # (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4337229/)
    # It would be better to add the probabilities directly as energy contribution in the RNAfold_compound,
    #  but it would need some serious benchmark and error check
    :param stockholm_msa:
    :return:
    """

    st_msa = read_st(stockholm_msa)

    # identify consensus unpaired nucleotides
    consensus_conservation_score = _alignment_column_conservation(st_msa, gap_chars='-')

    # transform it to PPlike annotation
    pplike_line = _transform_score_to_pplike_line(consensus_conservation_score)
    st_msa.column_annotations['pp_line'] = pplike_line

    uni_structure = encode_structure_unicode(st_msa.column_annotations['SS_cons'], gap_mark=49)

    # add consensus to each record
    for aligned_rec in st_msa:
        aligned_rec.letter_annotations['cons_uni'] = uni_structure
        aligned_rec.letter_annotations['pp_line'] = pplike_line

    # dealign and degap sequence, its score and consensus structure
    # this can be probably easily done with stockholm_align class, as i've prepared such function there
    out = []
    for unaligned_rec in st_msa.get_unalined_seqs(keep_letter_ann=True):

        # infer constraints before repairing the structure ?
        # because when repairing, new - to predict -  are tagged with '.'
        # infer_constraints

        unaligned_rec.letter_annotations['cons_uni'] = repair_structure_any_variant(
            unaligned_rec.letter_annotations['cons_uni'],
            gap_mark=49,
            rep_mark=48
        )

        inf_c = _constr4conserved_unaligned_from_pplike_score(
            unaligned_rec.letter_annotations['cons_uni'],
            unaligned_rec.letter_annotations['pp_line'],
            repred_tr,
            gap_mark=49
        )

        # replace text occurences of lenght less then specified value
        patt = 'x{' + str(conseq_conserved) + ',}'
        inf_c = re.sub('#', 'x', re.sub('x', '.', re.sub(patt, repl, inf_c)))

        unaligned_rec.letter_annotations['cons_db'] = decode_structure_unicode(
            unaligned_rec.letter_annotations['cons_uni']
        )

        unaligned_rec.letter_annotations['constraints'] = inf_c

        # infer constraints and predict
        structure_pred = rnafoldc(unaligned_rec, constraints_id='constraints')

        unaligned_rec.letter_annotations['ss0'] = structure_pred
        out.append(unaligned_rec)
    return out


def repl(m):
    return '#' * len(m.group())


def _constr4conserved_unaligned_from_pplike_score(consensus, pp_line, tr='8', gap_mark=49):
    """
    infer constraints based on pp_line for conserved unaligned residues
    need dot_bracket structure notation
    :param consensus:
    :param pp_line:
    :param tr: threshold from pp_line
    :return:
    """
    if isinstance(tr, int):
        tr = str(tr)
    alltr = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '*']
    if tr not in alltr:
        raise Exception('conservation level:{} not recognized, valid values are '.format(tr) + ' '.join(alltr))

    gm = chr(gap_mark)
    trin = alltr[alltr.index(tr):]
    constraints = []
    for c, pp in zip(consensus, pp_line):
        if c == gm:
            if pp in trin:
                constraints.append('x')
            else:
                constraints.append('.')
        else:
            constraints.append('.')

    return ''.join(constraints)


def rnafoldc(seqr, constraints_id='cons'):
    """
    predict mfe structure with rnafoldc
    :param seqr:
    :param constraints_id:
    :return:
    """
    new_str = _foldme_rnafoldc_noRNA(str(seqr.seq), seqr.letter_annotations[constraints_id])
    return new_str


def _foldme_rnafoldc_noRNA(seq, constraints):
    """
    mfe rnafold function without RNA lib
    :param seq:
    :param constraints:
    :return:
    """

    fd, tempfile = mkstemp()
    with os.fdopen(fd, 'w') as fh:
        fh.write('>seq01\n{}\n{}\n'.format(
            seq,
            constraints
        ))
    structure, _ = rnafold_prediction(tempfile, params='-C')
    os.remove(tempfile)
    return structure[0].letter_annotations['ss0']


def repair_structure_any_variant(structure, gap_mark=49, rep_mark=48):
    """
    needs special structure encoding
    """
    gm = chr(gap_mark)
    ns = []
    for n in structure:
        c = structure.count(n)
        if c == 1 & c != gm:
            ns.append(chr(rep_mark))
        elif c == 2:
            ns.append(n)
        else:
            ns.append(n)
    return ''.join(ns)


def encode_structure_unicode(structure, br=('(', ')'), gap_mark=49, warn=True):
    """
    encode structure to unicode characters so every basepair is unique so it is possible to know which where sliced
    in the alignment
    :param structure:
    :param br: tuple of chars dennoting a pairing bases
    :param gap_mark: char number dennoting gap, default 49 (1)
    :param warn: bool wherther warning should be raised if gapmark was set nondefault -> this must be reflected in
    calls to related functions
    :return:
    """
    if (gap_mark != 49) & warn:
        print('Warning: gap mark was set to {}, it is mandatory, that it is provided'
              'in calls to decode_structure_unicode and repair_structure_any_variant'.format(gap_mark))

    _, structure_nest = nesting(structure, br=br, gap_mark=gap_mark, initial_val=gap_mark + 1)

    uv = sorted(set(structure_nest) - {gap_mark})

    for i in uv:
        if structure_nest.count(i) <= 2:
            continue

        idx = [j for j, x in enumerate(structure_nest) if x == i][2:]
        # max_structure = max(structure_nest[:idx[0]])
        max_structure = max(structure_nest)
        shift = 1
        for e, f in grouped(idx, 2):
            structure_nest[e] = max_structure + shift
            structure_nest[f] = max_structure + shift
            shift += 1

    # thisi is to ensure that all chars returned are printable
    mut_str = copy.deepcopy(structure_nest)
    mx = max(structure_nest)
    for i, j in enumerate(structure_nest):
        if unicodedata.category(chr(j)) == 'Cc':
            # get possitions of matching nest
            pos = [q for q, o in enumerate(structure_nest) if o == j]

            # choose shift
            shift = j + mx
            while shift in mut_str:
                shift += 1

            mut_str[pos[0]] = shift
            mut_str[pos[1]] = shift

    return ''.join([chr(i) for i in mut_str])


def decode_structure_unicode(structure, br=('(', ')'), gap_char='.', gap_mark=49, rep_mark=48):
    """
    decode structure encoded with encode_structure_unicode
    :param structure: unicode encoded structure
    :param br: one pair of output brackets e.g. '(',')'
    :param gap_char: output gap character
    :param gap_mark: gap mark used for encoding, default 49
    :param rep_mark: used when repairing structure - default 48
    :return:
    """
    original_gap_mark = chr(gap_mark)
    original_rep_mark = chr(rep_mark)
    encoded_str = [i for i in structure]
    un = set(structure)
    out = encoded_str.copy()
    for n in un:
        idxs = [i for i, j in enumerate(encoded_str) if j == n]
        if n in (original_gap_mark, original_rep_mark):
            for idx in idxs:
                out[idx] = gap_char
        else:
            if len(idxs) == 2:
                out[idxs[0]] = br[0]
                out[idxs[1]] = br[1]
            else:
                raise RuntimeError('decoding structure failed, the "{}" key '
                                   'is note defined as a gap_mark or rep_mark'.format(n))

    return ''.join(out)


def grouped(iterable, n):
    """
    return n-tuples from iterable
    :param iterable:
    :param n: n - grouping val
    :return:
    """
    # print return grouped from list
    # adapted from http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
    return zip(*[iter(iterable)]*n)


def timeit_decorator(func):
    def wrapper(*arg, **kwargs):
        t1 = time.time()
        res = func(*arg, **kwargs)
        t2 = time.time()
        # print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        exec_time = t2 - t1
        return res, exec_time
    return wrapper


def tryit_decorator(func):
    def wrapper(*arg, **kwargs):
        # try:
        res = func(*arg, **kwargs)
        return res
        # except:
        #     print('problem in prediction {}'.format(str(func)))
        #     raise
    return wrapper


def sanitize_fasta_file(infasta, used_dict=None):
    fd, out_path = mkstemp()
    with open(infasta, 'r') as inh, os.fdopen(fd, 'w') as outh:
        k = [i for i in SeqIO.parse(inh, format='fasta')]
        san_seqs, san_dict = sanitize_fasta_names_in_seqrec_list(k, used_dict=used_dict)

        SeqIO.write(san_seqs, outh, format='fasta')
        return out_path, san_dict


@timeit_decorator
@tryit_decorator
def alifold_refold_prediction(nr_homologs_hits_fasta, all_hits_fasta, refold='refold', threads=None,
                              params=None, msa_alg='clustalo'):
    """
    return predicted structures for all hits based on provided sequence homologs
    ! beware, clustal mixes order of sequences in profile alignment, correct for it
    possible param keys: "clustal", "alifold", "clustalo_profile", "repred_unpaired_tr"
    """
    ml.debug(fname())
    nr_path, san_dict = sanitize_fasta_file(nr_homologs_hits_fasta)
    all_path, san_dict = sanitize_fasta_file(all_hits_fasta, used_dict=san_dict)

    if params is None:
        params = dict()

    ref_pred = ['refold', 'refold_rnafoldc', 'conserved_ss_rnafoldc']
    if refold not in ref_pred:
        raise Exception('refold procedure not recognized: {}, possible vaues are {}'.format(refold,
                                                                                            ' '.join(ref_pred)))

    cl_file = _aligner_block(nr_path, params, msa_alg, threads)

    # cannot rely on that, the order of a cl_file would be the same as the orded of the nr_homolog_hits_file

    if 'alifold' in params:
        ali_file = compute_alifold(cl_file, alifold_params=params['alifold'])
    else:
        ali_file = compute_alifold(cl_file)

    consensus_record = read_seq_str(ali_file)[0]

    clustalo_profile_params = '--outfmt clustal '
    if 'clustalo_profile' in params:
        clustalo_profile_params += params['clustalo_profile']
    if threads:
        clustalo_profile_params += ' --threads {}'.format(threads)
    realign_file = run_clustal_profile2seqs_align(cl_file, all_path, clustalo_params=clustalo_profile_params)
    realign_alig = AlignIO.read(realign_file, format='clustal')

    # slice alignment ( get seqname from nr_homolog_hits_file, find it in the realign and slice the whole segment off
    #  take care that the id may be the same and it must be checked for multiple occurence

    first_nr_record = _parse_first_record_only(nr_path)

    realign_allseq_possition = [i for i, seq in enumerate(realign_alig) if seq.id == first_nr_record.id]

    new_alig_for_refold = realign_alig[:realign_allseq_possition[-1]]
    old_alig_in_new = realign_alig[realign_allseq_possition[-1]:]

    orig_alignment = AlignIO.read(cl_file, format='clustal')

    first_original_alignment_record = orig_alignment[0]

    match_original_seq_in_new_alig = [i for i in old_alig_in_new if i.id == first_original_alignment_record.id][0]

    mapping = _map_alignment_columns_from_profile_match(first_original_alignment_record,
                                                        match_original_seq_in_new_alig)

    # map and repair structure when mapping is unbiguous
    cs_encode = encode_structure_unicode(consensus_record.letter_annotations['ss0'])
    new_consensus_structure_encoded = _repair_consensus_structure_by_maping(cs_encode,
                                                                            mapping,
                                                                            len(match_original_seq_in_new_alig.seq),
                                                                            gap_char=49)
    new_consensus_structure_repaired = repair_structure_any_variant(new_consensus_structure_encoded)

    new_consensus_structure = decode_structure_unicode(new_consensus_structure_repaired)

    new_consensus_sequence = _repair_consensus_structure_by_maping(str(consensus_record.seq),
                                                                   mapping,
                                                                   len(match_original_seq_in_new_alig.seq),
                                                                   gap_char=ord('_'))

    # write new consensus to a file
    a_fd, new_alifold_consensus_file = mkstemp()
    with os.fdopen(a_fd, 'w') as f:
        f.write(new_consensus_sequence + '\n')
        f.write(new_consensus_structure + '\n')

    # write sliced alignment to a file
    sa_fd, sliced_alignment_file = mkstemp()
    with os.fdopen(sa_fd, 'w') as f:
        AlignIO.write(new_alig_for_refold, f, 'clustal')

    # now process the file, and map alignment to consensus structure (how)
    #  take 1 sequence from original alignment and map it to new sequence alignment ?
    #  mozna lepsi obracene, protoze se gapy mohou pridat v novem alignmentu spis nez by zmizely

    if refold in ['refold', 'refold_rnafoldc']:
        refold_file = compute_refold(sliced_alignment_file, new_alifold_consensus_file)

        if refold == 'refold_rnafoldc':
            structures_rnafold_c_file = compute_constrained_prediction(refold_file)
            seq_str = read_seq_str(structures_rnafold_c_file)
            os.remove(structures_rnafold_c_file)
        else:
            seq_str = read_seq_str(refold_file)

        os.remove(refold_file)

    else:
        st_alig_file = build_stockholm_from_clustal_alig(sliced_alignment_file, new_alifold_consensus_file)
        if 'repred_unpaired_tr' in params:
            repred_tr = str(params['repred_unpaired_tr'])
        else:
            repred_tr = '9'
        if 'conseq_conserved' in params:
            conseq_conserved = params['conseq_conserved']
        else:
            conseq_conserved = 1
        seq_str = _refold_with_unpaired_conservation(st_alig_file,
                                                     repred_tr=repred_tr,
                                                     conseq_conserved=conseq_conserved)
        os.remove(st_alig_file)

    structures_out = desanitize_fasta_names_in_seqrec_list(seq_str, san_dict)

    os.remove(nr_path)
    os.remove(all_path)
    os.remove(sliced_alignment_file)
    os.remove(new_alifold_consensus_file)
    os.remove(cl_file)
    os.remove(ali_file)
    os.remove(realign_file)

    return structures_out


@timeit_decorator
@tryit_decorator
def tcoffee_rcoffee_refold_prediction(nr_homolog_hits_file, all_hits_fasta, refold='refold', threads=None, params=None):
    """
    return predicted structures for all hits based on provided sequence homologs
    possible param keys: "clustal", "alifold", "clustalo_profile", "repred_unpaired_tr"
    """
    ml.debug(fname())
    if params is None:
        params = dict()

    ref_pred = ['refold', 'refold_rnafoldc', 'conserved_ss_rnafoldc']
    if refold not in ref_pred:
        raise Exception('refold procedure not recognized: {}, possible vaues are {}'.format(refold,
                                                                                            ' '.join(ref_pred)))
    # sanitize names input to a new file, sanitize all hits fasta also.
    fid, nr_sanitized_fasta_file = mkstemp()
    fid2, all_hits_fasta_sanitized = mkstemp()

    with open(nr_homolog_hits_file, 'r') as f, os.fdopen(fid, 'w') as nr_fid:
        nr_seqs2san = [nrseq for nrseq in SeqIO.parse(f, 'fasta')]
        nr_san_seqs, san_seqs_dict = sanitize_fasta_names_in_seqrec_list(nr_seqs2san)
        SeqIO.write(nr_san_seqs, nr_fid, 'fasta')

    with open(all_hits_fasta, 'r') as ff, os.fdopen(fid2, 'w') as all_fid:
        all_seqs2san = [aseq for aseq in SeqIO.parse(ff, 'fasta')]
        all_san_seq_hits, san_seqs_dict = sanitize_fasta_names_in_seqrec_list(all_seqs2san, used_dict=san_seqs_dict)
        SeqIO.write(all_san_seq_hits, all_fid, 'fasta')

    # ================================================================================================
    # run tcoffee initial alignment of homologous seqs for consensus prediction

    tcoffee_rcoffe_p = '-quiet '
    if 'tcoffee_rcoffee' in params:
        tcoffee_rcoffe_p += params['tcoffee_rcoffee']

    cl_file = run_tcoffee(nr_sanitized_fasta_file,
                          tcoffee_rcoffee_params=tcoffee_rcoffe_p,
                          threads=threads)

    # ================================================================================================
    # prdict consensus structure
    if 'alifold' in params:
        ali_file = compute_alifold(cl_file, alifold_params=params['alifold'])
    else:
        ali_file = compute_alifold(cl_file)

    consensus_record = read_seq_str(ali_file)[0]

    # ================================================================================================
    # align all sequence to sequence profile
    tcoffe_profile_params = ' -quiet'
    if 'tcoffee_profile' in params:
        tcoffe_profile_params = params['tcoffee_profile']

    realign_file = run_tcoffee_profile_aling(
        all_hits_fasta_sanitized,
        'rcoffee',
        cl_file,
        threads=threads,
        tcoffee_rcoffee_params=tcoffe_profile_params
    )

    # ===============================================================================================
    # read and desanitize output
    realign_alig = AlignIO.read(realign_file, format='clustal')
    des_alig_seqs = desanitize_fasta_names_in_seqrec_list([i.upper() for i in realign_alig], san_seqs_dict)
    realign_alig = AlignIO.MultipleSeqAlignment(des_alig_seqs)

    # ===========================================
    # process the alignment

    # slice alignment ( get seqname from nr_homolog_hits_file, find it in the realign and slice the whole segment off
    #  take care that the id may be the same and it must be checked for multiple occurence

    first_nr_record = _parse_first_record_only(nr_homolog_hits_file)

    realign_allseq_possition = [i for i, seq in enumerate(realign_alig) if seq.id == first_nr_record.id]

    new_alig_for_refold = realign_alig[:realign_allseq_possition[-1]]
    old_alig_in_new = realign_alig[realign_allseq_possition[-1]:]

    orig_alignment = AlignIO.read(cl_file, format='clustal')

    orig_alignment = AlignIO.MultipleSeqAlignment(desanitize_fasta_names_in_seqrec_list(
        [i.upper() for i in orig_alignment], used_dict=san_seqs_dict))

    first_original_alignment_record = orig_alignment[0]

    match_original_seq_in_new_alig = [i for i in old_alig_in_new if i.id == first_original_alignment_record.id][0]

    mapping = _map_alignment_columns_from_profile_match(first_original_alignment_record,
                                                        match_original_seq_in_new_alig)

    # map and repair structure when mapping is unbiguous
    cs_encode = encode_structure_unicode(consensus_record.letter_annotations['ss0'])
    new_consensus_structure_encoded = _repair_consensus_structure_by_maping(cs_encode,
                                                                            mapping,
                                                                            len(match_original_seq_in_new_alig.seq),
                                                                            gap_char=49)
    new_consensus_structure_repaired = repair_structure_any_variant(new_consensus_structure_encoded)

    new_consensus_structure = decode_structure_unicode(new_consensus_structure_repaired)

    new_consensus_sequence = _repair_consensus_structure_by_maping(str(consensus_record.seq),
                                                                   mapping,
                                                                   len(match_original_seq_in_new_alig.seq),
                                                                   gap_char=ord('_'))

    # write new consensus to a file
    a_fd, new_alifold_consensus_file = mkstemp()
    with os.fdopen(a_fd, 'w') as f:
        f.write(new_consensus_sequence + '\n')
        f.write(new_consensus_structure + '\n')

    new_alig_for_refold_san, used_dict_2 = sanitize_fasta_names_in_seqrec_list(new_alig_for_refold)

    # write sliced alignment to a file
    sa_fd, sliced_alignment_file = mkstemp()
    with os.fdopen(sa_fd, 'w') as f:
        AlignIO.write(
            AlignIO.MultipleSeqAlignment(new_alig_for_refold_san),
            f, 'clustal'
        )

    if refold in ['refold', 'refold_rnafoldc']:
        refold_file = compute_refold(sliced_alignment_file, new_alifold_consensus_file)

        if refold == 'refold_rnafoldc':
            structures_rnafold_c_file = compute_constrained_prediction(refold_file)
            seq_str = read_seq_str(structures_rnafold_c_file)
            os.remove(structures_rnafold_c_file)
        else:
            seq_str = read_seq_str(refold_file)

        os.remove(refold_file)
        os.remove(new_alifold_consensus_file)
        os.remove(sliced_alignment_file)
    else:
        st_alig_file = build_stockholm_from_clustal_alig(sliced_alignment_file, new_alifold_consensus_file)
        if 'repred_unpaired_tr' in params:
            repred_tr = str(params['repred_unpaired_tr'])
        else:
            repred_tr = '9'
        if 'conseq_conserved' in params:
            conseq_conserved = params['conseq_conserved']
        else:
            conseq_conserved = 1
        seq_str = _refold_with_unpaired_conservation(st_alig_file,
                                                     repred_tr=repred_tr,
                                                     conseq_conserved=conseq_conserved)
        os.remove(st_alig_file)

    structures_out = desanitize_fasta_names_in_seqrec_list(seq_str, used_dict_2)

    os.remove(cl_file)
    os.remove(ali_file)
    os.remove(realign_file)
    os.remove(nr_sanitized_fasta_file)
    os.remove(all_hits_fasta_sanitized)

    return structures_out


@timeit_decorator
@tryit_decorator
def decouple_homologs_alifold_refold_prediction(nr_homolog_hits_file, homologous_seqs, all_hits_fasta, refold='refold',
                                                threads=None, align='tcoffee',
                                                params=None):
    """
    return predicted structures for all hits based on provided sequence homologs
    instead of alignment of profile with all sequences make a profile, predict the homologous sequences

    which variant?
    a) predict from non redundant - the non redundant are there to ensure that the alifold is not much skewed,
        but there is a risc, that realignment won't be possible
        I could compute the profile, then just realign the homologous ones and the non-homologs separately
    b) remove only same seqeunces
    c) get some consensus prediction algorithm which can weight the position and sequences while predicting
        the consensus structure

    result => make consensus structure and profile, then align all homologs and predict structure and align the
            non homologs ang predict for them
    possible param keys: "clustal", "alifold", "clustalo_profile", "repred_unpaired_tr"
    """
    ml.debug(fname())
    if params is None:
        params = dict()

    ref_pred = ['refold', 'refold_rnafoldc', 'conserved_ss_rnafoldc']
    if refold not in ref_pred:
        raise Exception('refold procedure not recognized: {}, possible vaues are {}'.format(refold,
                                                                                            ' '.join(ref_pred)))

    clustalo_params = '--outfmt clustal '
    if 'clustalo' in params:
        clustalo_params += params['clustalo']
    if threads:
        clustalo_params += ' --threads={}'.format(threads)
    cl_file = compute_clustalo_clasic(nr_homolog_hits_file, clustalo_params=clustalo_params)

    if 'alifold' in params:
        ali_file = compute_alifold(cl_file, alifold_params=params['alifold'])
    else:
        ali_file = compute_alifold(cl_file)

    consensus_record = read_seq_str(ali_file)[0]

    # too here it is the same
    # split all_seqs to homologous and nonhomologous
    hom_seq_db = [str(hseq.seq) for hseq in homologous_seqs]

    remove_nh_duplicit = False

    nonhom_seqs = []
    all_seqs_db = []
    with open(all_hits_fasta, 'r') as f:
        for seq in SeqIO.parse(f, 'fasta'):
            all_seqs_db.append(seq)
            if str(seq.seq) not in hom_seq_db:
                nonhom_seqs.append(seq)
        ah_fd, all_h_seq_fasta = mkstemp()
        an_fd, all_n_seq_fasta = mkstemp()

        if len(nonhom_seqs) == 1 and align != 'tcoffee':
            print('Warning: Only 1 sequence which is considered non-homologous. \n'
                  '         Duplicating record for succesfull execution of profile alignment. '
                  'Added duplicid will be removed\n'
                  '         This may cause inaccurate structure prediction.')
            duplicid = copy.deepcopy(nonhom_seqs[0])
            duplicid.id = 'duplicid' + duplicid.id
            nonhom_seqs.append(duplicid)
            remove_nh_duplicit = True

        san_homologous_seqs, san_dict = sanitize_fasta_names_in_seqrec_list(homologous_seqs)
        san_nonhom_seqs, san_dict = sanitize_fasta_names_in_seqrec_list(nonhom_seqs, used_dict=san_dict)

        with os.fdopen(ah_fd, 'w') as ah, os.fdopen(an_fd, 'w') as an:
            write_fasta_from_list_of_seqrecords(ah, san_homologous_seqs)
            write_fasta_from_list_of_seqrecords(an, san_nonhom_seqs)

    # now predict for homologs and non homologs separately
    # clustal and mafft cannot align only one sequence to profile, which can be the case - need to use tcoffee
    # skipp the computations, if none sequence is in the part

    if len(homologous_seqs) > 0:
        hom_structures = _profile_alig_prediction(params, threads, cl_file, all_h_seq_fasta,
                                                  nr_homolog_hits_file,
                                                  consensus_record,
                                                  refold,
                                                  align)
    else:
        print('skipping profile alignment of homologous - no sequences in the list', flush=True)
        hom_structures = []

    if len(nonhom_seqs) > 0:
        nonhom_structures = _profile_alig_prediction(params, threads, cl_file, all_n_seq_fasta,
                                                     nr_homolog_hits_file,
                                                     consensus_record,
                                                     refold,
                                                     align)
    else:
        print('skipping profile alignment of nonhomologous - no sequences in the list', flush=True)
        nonhom_structures = []

    # desanitize sequence list
    hom_structures = desanitize_fasta_names_in_seqrec_list(hom_structures, san_dict)
    nonhom_structures = desanitize_fasta_names_in_seqrec_list(nonhom_structures, san_dict)

    if remove_nh_duplicit:
        for i, structure in enumerate(nonhom_structures):
            if structure.id.startswith('duplicid'):
                del nonhom_seqs[i]
                del nonhom_structures[i]
        # check if done correctly
        assert len(nonhom_structures) == 1

    # build list in correct order
    seq_str = [1] * len(all_seqs_db)
    all_seqs_id = [seq.id for seq in all_seqs_db]

    # homologous seqs can contain query sequence
    for seq in hom_structures + nonhom_structures:
        # p_i = [i for i, n in enumerate(all_seqs_id) if n.startswith(seq.id)]
        p_i = [i for i, n in enumerate(all_seqs_id) if seq.id in n]
        if len(p_i) != 0:
            seq_str[p_i[0]] = seq

    # check output
    if not all([isinstance(i, SeqRecord) for i in seq_str]):
        raise Exception('structure output was not reconstructed correctly')

    os.remove(cl_file)
    os.remove(ali_file)
    os.remove(all_h_seq_fasta)
    os.remove(all_n_seq_fasta)

    return seq_str


@timeit_decorator
@tryit_decorator
def cmmodel_rnafold_c(allhits_fasta, cmmodel_file, threads=None, params=None):
    ml.debug(fname())
    if params is None:
        params = dict()

    allhits_fasta, san_dict = sanitize_fasta_file(allhits_fasta)

    cmalign_params = ''
    if threads:
        cmalign_params += '--cpu {}'.format(threads)

    if 'cmalign' in params and params['cmalign']:
        cmalign_params += ' ' + params['cmalign']

    if '--notrunc' not in cmalign_params:
        cmalign_params += ' --notrunc'

    alig_file = run_cmalign_on_fasta(allhits_fasta, cmmodel_file, cmalign_params=cmalign_params)

    # multiple sequence cm align
    # split by sequence, then run the rest
    cm_alig = read_st(alig_file)

    structures = []
    for single_alig in trim_cmalign_sequence_by_refseq_one_seq(cm_alig, rs='SS_cons', convert2uppercase=True):
        out_alig, trimmed_seq = trim_and_repair_single_cm_alignment(single_alig)
        conserved_structure_pairs = find_nc_and_remove(str(trimmed_seq.seq), trimmed_seq.letter_annotations['dec_str'])
        trimmed_seq.letter_annotations['constrains'] = conserved_structure_pairs

        # constraint prediction
        # write constraint file
        fd, temp_constraint_file = mkstemp()
        with os.fdopen(fd, 'w') as tmpf:
            tmpf.write('>{}\n{}\n{}\n'.format(
                trimmed_seq.id,
                str(trimmed_seq.seq),
                trimmed_seq.letter_annotations['constrains']
            ))

        single_structure, _ = rnafold_prediction(temp_constraint_file, params='-C')

        os.remove(temp_constraint_file)

        # trimmed_seq.letter_annotations['final'] = single_structure[0].letter_annotations['ss0']
        structures.append(single_structure[0])

    str_out = desanitize_fasta_names_in_seqrec_list(structures, san_dict)

    return str_out


def trim_and_repair_single_cm_alignment(cmalig):
    """
    trim and repair single! cm alignment
    build specificaly for run_rsearch script
    :param cmalig:
    :return:
    """
    ml.debug(fname())
    if len(cmalig) != 1:
        raise Exception('This function is build specificaly for alignment of single sequence to cm model.'
                        'Using it for more sequences would require addition of cm msa trimming module.')

    # recode structure to br
    br_model_structure = cm_strucutre2br(cmalig.column_annotations['SS_cons'])
    encoded_whole_structure = encode_structure_unicode(br_model_structure)

    # add consensus structure to letter_ann
    cmalig[0].letter_annotations['SS_cons'] = cmalig.column_annotations['SS_cons']
    cmalig[0].letter_annotations['en_cons'] = encoded_whole_structure

    # trimming is based on presence of - (dashes) in target sequence
    # essentialy, they are mapped to consensus structure
    unaligned_seq = [i for i in cmalig.get_unalined_seqs(keep_letter_ann=True)][0]

    repaired_structure = repair_structure_any_variant(unaligned_seq.letter_annotations['en_cons'])
    unaligned_seq.letter_annotations['rep_str'] = repaired_structure

    dec_structure = decode_structure_unicode(repaired_structure, gap_char='.')
    unaligned_seq.letter_annotations['dec_str'] = dec_structure

    return cmalig, unaligned_seq


def find_nc_and_remove(sequence, structure, allowed_bp=('AT', 'GC', 'GU', 'AU'), mismatch_char='.'):
    """
    find non canonical base pairs (i.e. basepairs which are not in provided list "allowed_bp")
    :param sequence:
    :param structure:
    :param allowed_bp: tuple of base pairs which are allowed to form
    :param mismatch_char:
    :return:
    """
    # build all possible set of allowed basepairs
    bps = set(list(allowed_bp) + [i[::-1] for i in allowed_bp])

    # encoded structure places 2 unique chars for every
    en_str = encode_structure_unicode(structure, gap_mark=49)

    mutable_str = list(structure)

    for e in set(en_str) - set(chr(49)):
        assert en_str.count(e) == 2
        match_pos = [i for i, j in enumerate(en_str) if j == e]

        pair = ''.join([sequence[match_pos[0]], sequence[match_pos[1]]])
        if pair in bps:
            # allowed
            pass
        else:
            # not allowed
            mutable_str[match_pos[0]] = mismatch_char
            mutable_str[match_pos[1]] = mismatch_char

    return ''.join(mutable_str)


def _profile_alig_prediction(params, usethreads, cl_file, seqs_to_align_and_predict_fasta_file,
                             homolog_profile_fasta_file,
                             consensus_record, refold,
                             align):
    ml.debug(fname())
    if align == 'clustalo':
        # clustal cannot align only one sequence to profile
        # mafft also fails for this case
        # tcoffee can, also in mode rcoffee --- using it
        # beware, that the added sequence is first one

        clustalo_profile_params = '--outfmt clustal '
        if 'clustalo_profile' in params:
            clustalo_profile_params += params['clustalo_profile']
        if usethreads:
            clustalo_profile_params += ' --threads {}'.format(usethreads)
        realign_file = run_clustal_profile2seqs_align(cl_file, seqs_to_align_and_predict_fasta_file,
                                                      clustalo_params=clustalo_profile_params)
    elif align == 'tcoffee':
        tcoffe_profile_params = ' -quiet'
        if 'tcoffee_profile' in params:
            tcoffe_profile_params = params['tcoffee_profile']

        realign_file = run_tcoffee_profile_aling(seqs_to_align_and_predict_fasta_file,
                                                 'rcoffee',
                                                 cl_file,
                                                 threads=usethreads,
                                                 tcoffee_rcoffee_params=tcoffe_profile_params)

    else:
        raise AttributeError('"{}" is not recognized aligner, known are "tcoffee" and "clustalo"')

    realign_alig = AlignIO.read(realign_file, format='clustal')

    # ensure that sequences for refold are in uppercase
    ap = []
    for sline in realign_alig:
        ap.append(sline.upper())
        realign_alig = AlignIO.MultipleSeqAlignment(ap)
    del ap

    # slice alignment ( get seqname from nr_homolog_hits_file, find it in the realign and slice the whole segment off
    #  take care that the id may be the same and it must be checked for multiple occurence

    first_nr_record = _parse_first_record_only(homolog_profile_fasta_file)

    # need to run tcoffee line sanitize because name could be changed.
    if align == 'tcoffee':
        realign_allseq_possition = [i for i, seq in enumerate(realign_alig) if seq.id == tcoffee_line_sanitizer(first_nr_record.id)]
    else:
        realign_allseq_possition = [i for i, seq in enumerate(realign_alig) if seq.id == first_nr_record.id]

    new_alig_for_refold = realign_alig[:realign_allseq_possition[-1]]
    old_alig_in_new = realign_alig[realign_allseq_possition[-1]:]

    orig_alignment = AlignIO.read(cl_file, format='clustal')

    first_original_alignment_record = orig_alignment[0]

    if align == 'tcoffee':
        match_original_seq_in_new_alig = [i for i in old_alig_in_new if i.id == tcoffee_line_sanitizer(first_original_alignment_record.id)][0]
    else:
        match_original_seq_in_new_alig = [i for i in old_alig_in_new if i.id == first_original_alignment_record.id][0]

    mapping = _map_alignment_columns_from_profile_match(first_original_alignment_record,
                                                        match_original_seq_in_new_alig)

    # map and repair structure when mapping is unbiguous
    cs_encode = encode_structure_unicode(consensus_record.letter_annotations['ss0'])
    new_consensus_structure_encoded = _repair_consensus_structure_by_maping(cs_encode,
                                                                            mapping,
                                                                            len(match_original_seq_in_new_alig.seq),
                                                                            gap_char=49)
    new_consensus_structure_repaired = repair_structure_any_variant(new_consensus_structure_encoded)

    new_consensus_structure = decode_structure_unicode(new_consensus_structure_repaired)

    new_consensus_sequence = _repair_consensus_structure_by_maping(str(consensus_record.seq),
                                                                   mapping,
                                                                   len(match_original_seq_in_new_alig.seq),
                                                                   gap_char=ord('_'))

    # write new consensus to a file
    a_fd, new_alifold_consensus_file = mkstemp()
    with os.fdopen(a_fd, 'w') as f:
        f.write(new_consensus_sequence + '\n')
        f.write(new_consensus_structure + '\n')

    # write sliced alignment to a file
    sa_fd, sliced_alignment_file = mkstemp()
    with os.fdopen(sa_fd, 'w') as f:
        AlignIO.write(new_alig_for_refold, f, 'clustal')

    # now process the file, and map alignment to consensus structure (how)
    #  take 1 sequence from original alignment and map it to new sequence alignment ?
    #  mozna lepsi obracene, protoze se gapy mohou pridat v novem alignmentu spis nez by zmizely

    if refold in ['refold', 'refold_rnafoldc']:
        refold_file = compute_refold(sliced_alignment_file, new_alifold_consensus_file)

        if refold == 'refold_rnafoldc':
            structures_rnafold_c_file = compute_constrained_prediction(refold_file)
            seq_str = read_seq_str(structures_rnafold_c_file)
            os.remove(structures_rnafold_c_file)
        else:
            seq_str = read_seq_str(refold_file)

        os.remove(refold_file)
        os.remove(new_alifold_consensus_file)
        os.remove(sliced_alignment_file)
    else:
        st_alig_file = build_stockholm_from_clustal_alig(sliced_alignment_file, new_alifold_consensus_file)
        if 'repred_unpaired_tr' in params:
            repred_tr = str(params['repred_unpaired_tr'])
        else:
            repred_tr = '9'
        if 'conseq_conserved' in params:
            conseq_conserved = params['conseq_conserved']
        else:
            conseq_conserved = 1
        seq_str = _refold_with_unpaired_conservation(st_alig_file,
                                                     repred_tr=repred_tr,
                                                     conseq_conserved=conseq_conserved)
        os.remove(st_alig_file)
    os.remove(realign_file)
    return seq_str


def _aligner_block(nr_homolog_hits_file, params, msa_alg, threads=None):
    """
    returns alignment file in clustal format
    :param nr_homolog_hits_file:
    :param params:
    :param msa_alg:
    :param threads: int
    :return:
    """
    ml.debug(fname())
    if msa_alg == 'clustalo':
        if params and ('clustalo' in params) and params['clustalo']:
            clustal_params = '--outfmt=clustal {}'.format(params['clustalo'])
        else:
            clustal_params = '--outfmt=clustal'
        if threads:
            clustal_params += ' --threads={}'.format(threads)
        alig_file = compute_clustalo_clasic(nr_homolog_hits_file, clustalo_params=clustal_params)

    elif msa_alg == 'muscle':
        if params and ('muscle' in params) and params['muscle']:
            alig_file = run_muscle(nr_homolog_hits_file, muscle_params=params['muscle'], reorder=True)
        else:
            alig_file = run_muscle(nr_homolog_hits_file, reorder=True)

    elif msa_alg == 'rcoffee':
        # sanitize nr hits for rcoffee
        nr_hits_list = [i for i in SeqIO.parse(nr_homolog_hits_file, format='fasta')]
        san_nr_hits, san_dict = sanitize_fasta_names_in_seqrec_list(nr_hits_list)

        nrfd, nr_sanitized_fasta_file = mkstemp()
        with os.fdopen(nrfd, 'w') as f:
            SeqIO.write(san_nr_hits, f, format='fasta')

        tcoffee_rcoffe_p = '-quiet '
        if params and ('rcoffee' in params) and params['rcoffee']:
            tcoffee_rcoffe_p += params['rcoffee']

        print('runing rcofee')
        alig_file = run_tcoffee(nr_sanitized_fasta_file,
                                mode='rcoffee',
                                tcoffee_rcoffee_params=tcoffee_rcoffe_p,
                                threads=threads)

        aligned_sequences = AlignIO.read(alig_file, format='clustal')
        des_as = desanitize_fasta_names_in_seqrec_list(aligned_sequences, san_dict)

        alignment = AlignIO.MultipleSeqAlignment(des_as)

        os.remove(alig_file)
        del alig_file

        afd, alig_file = mkstemp()
        with os.fdopen(afd, 'w') as f:
            AlignIO.write(alignment, f, format='clustal')

        del des_as
        del nrfd
        del san_nr_hits
        del nr_hits_list
        os.remove(nr_sanitized_fasta_file)
        del nr_sanitized_fasta_file

    else:
        print('invalig MSA alg chosen {}, valid are clustalo, muscle, rcoffee'.format(msa_alg))
        raise AttributeError()

    return alig_file


def remove_sharp_hairpins(structure):
    """
    input is an RNA structure in dot bracket format
    it removes the sharp hairpins "(.)" and "(..)" should they be present in the structure

    Intended as a preprocess step after extraction of a template structure from cm model file and before
        the shape extraction analysis.
        The RNAshapes raises error when structure with those aforementioned hairpin substructure are encountered.

    :param structure:
    :return:
    """

    intermediate = re.sub('\(\.\)', '...', structure)
    final = re.sub('\(\.\.\)', '\.\.\.\.', intermediate)
    return final


@tryit_decorator
@timeit_decorator
def rfam_subopt_pred(all_sequence_fasta, query_file, params=None, threads=None):
    ml.debug(fname())
    if params is None:
        params = dict()

    best_model = get_cm_model(query_file, params=params, threads=threads)
    cm_ref_str = extract_ref_structure_fromRFAM_CM(best_model)

    if params and ('mfold' in params) and params['mfold']:
        assert isinstance(params['mfold'], (tuple, list)) and 3 == len(params['mfold'])
        subs = run_hybrid_ss_min(all_sequence_fasta, mfold=params['mfold'])
    else:
        subs = run_hybrid_ss_min(all_sequence_fasta)

    new_structures = []

    # now compute rna distance score
    for seq in subs:
        str2compare = []
        key_list = []
        if seq.annotations['predicted']:
            for key in seq.annotations['sss']:
                str2compare.append((seq.letter_annotations[key], cm_ref_str))
                key_list.append(key)
            rnadist_score = rnadistance_one_thread(str2compare)

            # select best ie lowes score
            mindisti = rnadist_score.index(min(rnadist_score))
            new_structures.append(SeqRecord(seq.seq,
                                            id=seq.id,
                                            annotations={'sss': ['ss0']},
                                            letter_annotations={'ss0': seq.letter_annotations[key_list[mindisti]]}
                                            )
                                  )
        else:
            new_structures.append(seq)

    return new_structures


@tryit_decorator
@timeit_decorator
def cmscan_rapidshapes(all_sequence_fasta, query_file, params=None, threads=None):
    """
    need to pass --allowLP 1 parameter to all Rapidshapes and RNAshapes calls because in rfam reference there can be
    lonely pairs
    """
    ml.debug(fname())
    # todo rapidshapes in parallel
    if params is None:
        params = dict()

    best_model = get_cm_model(query_file, params=params, threads=threads)

    # select best hiting model
    # todo what if not suficient hit is not reached

    cm_ref_str = extract_ref_structure_fromRFAM_CM(best_model)

    # need to precheck extracted structure
    # (.) and (..) are not allowed by RNAshapes
    san_cm_ref_str = remove_sharp_hairpins(cm_ref_str)
    if cm_ref_str != san_cm_ref_str:
        print('reference structure for the {} contains sharp hairpins, which are forbiden in RNAshapes, '
              'removal of this hairpins is essential for sucessful shape extraction and further structure'
              ' computation'.format(best_model))
        print('reference structure:' + cm_ref_str)
        print('repaired ref struct:' + san_cm_ref_str)

    # cm_shape = db2shape(cm_ref_str)
    # or to leverage possibility of selecting shape level
    fd, temp_str_file = mkstemp()
    with os.fdopen(fd, 'w') as f:
        f.write(san_cm_ref_str + '\n')

    if params and ('shape_level' in params) and params['shape_level']:
        sh_level = params['shape_level']
    else:
        sh_level = 5
    out_shape_list = structures2shape(infile=temp_str_file, shape_level=sh_level, params='--allowLP 1')

    out_shape = ' '.join([i.shape for i in out_shape_list])

    os.remove(temp_str_file)
    del temp_str_file

    # same as below to prevent unnecesary computation
    all_seqs_list = [s for s in SeqIO.parse(all_sequence_fasta, format='fasta')]

    uniq_all_seqs_dict = uniquify_seq_list(all_seqs_list)
    uniq_all_structures_dict = dict()

    for i, seq_str in enumerate(uniq_all_seqs_dict.keys()):
        # check if sequence contains noncanonical bases
        if re.search('[^ACGTU]', seq_str, flags=re.IGNORECASE):
            print('running rapid shapes ambiguos prediction #{}'.format(i), flush=True)
            print('for {}'.format(' '.join(uniq_all_seqs_dict[seq_str])))

            if params and ('rapidshapes' in params) and params['rapidshapes']:
                a_structure = rapidshapes_ambigous_prediction(seq_str, out_shape, params.get('rapidshapes') + '--allowLP 1', shape_level=sh_level)
            else:
                a_structure = rapidshapes_ambigous_prediction(seq_str, out_shape, '--allowLP 1', shape_level=sh_level)
            uniq_all_structures_dict[seq_str] = a_structure

            del a_structure

        else:
            devprint('rapidshapes - normal sequence #{}'.format(i), flush=True)
            tp, temp_fasta = mkstemp()
            with os.fdopen(tp, 'w') as f:
                f.write('>seq{}\n{}\n'.format(i, seq_str))

            if params and ('rapidshapes' in params) and params['rapidshapes']:
                rs_file = rapid_shapes_list(temp_fasta, out_shape,
                                            rapidshapes_params=params.get('rapidshapes') + ' --allowLP 1',
                                            shapeLevel=sh_level)
            else:
                rs_file = rapid_shapes_list(temp_fasta, out_shape,
                                            rapidshapes_params='--allowLP 1',
                                            shape_level=sh_level)
            try:
                r_structure = read_rapid_shapes_list_outfile(rs_file)
            except IndexError as e:
                # catch error when structure is not predicted by rapidshapes despite being fed correct parameters and
                # no error is reported whatsoever
                r_structure = [
                    SeqRecord(
                        Seq(seq_str),
                        letter_annotations={out_shape: 'X'*len(seq_str)},
                    )
                ]
                print(e)
                print('Rapidshapes prediction failed - no structure predicted.')

            uniq_all_structures_dict[seq_str] = r_structure[0].letter_annotations[out_shape]

            del tp
            os.remove(temp_fasta)
            os.remove(rs_file)

    # rebuild seq list
    for seq in all_seqs_list:
        seq.annotations['sss'] = ['ss0']
        if 'X' in uniq_all_structures_dict[str(seq.seq)] and all(['X' == i for i in uniq_all_structures_dict[str(seq.seq)]]):
            seq.letter_annotations['ss0'] = '.'*len(seq)
            seq.annotations['predicted'] = False
        else:
            seq.letter_annotations['ss0'] = uniq_all_structures_dict[str(seq.seq)]
            seq.annotations['predicted'] = True

    return all_seqs_list


@timeit_decorator
@tryit_decorator
def msa_alifold_rapidshapes(all_sequence_fasta, nr_homolog_hits_file, params=None, threads=True, msa_alg=''):
    ml.debug(fname())
    # rapidshapes can only work with canonical sequence encoding
    #  ie ACGUT - IUPAC ambigous encoding is no supported
    # but automaticly retrieved sequences frm a database has this encoding
    # todo rapidshapes in parallel
    if params is None:
        params = dict()

    print('aligner block')
    alig_file = _aligner_block(nr_homolog_hits_file, params, msa_alg, threads=threads)

    print('compute alifold')
    if params and ('alifold' in params) and params['alifold']:
        alif_file = compute_alifold(alig_file, alifold_params=params['alifold'])
    else:
        alif_file = compute_alifold(alig_file)

    # possibly need to decode alifold structure
    alif_str = read_seq_str(alif_file)[0]
    consensus_structure = alif_str.letter_annotations['ss0']

    # get shape from consensus
    # shape = db2shape(consensus_structure)
    # shape_str = ''.join(shape[0])

    if params and ('shape_level' in params) and params['shape_level']:
        sh_level = params['shape_level']
    else:
        sh_level = 5

    san_template_str = remove_sharp_hairpins(consensus_structure)
    # get shape from consensus
    fd, temp_str_file = mkstemp()
    with os.fdopen(fd, 'w') as f:
        f.write(san_template_str + '\n')
    out_shape_list = structures2shape(infile=temp_str_file, shape_level=sh_level, params='--allowLP 1')
    shape_str = ' '.join([i.shape for i in out_shape_list])

    os.remove(temp_str_file)
    del temp_str_file

    # predict sequences with that shape
    # use rapidshapes sequence by sequence detenct sequences with ambiguos encoding and do not run for them
    # or use both variants ?
    all_seqs_list = [s for s in SeqIO.parse(all_sequence_fasta, format='fasta')]

    uniq_all_seqs_dict = uniquify_seq_list(all_seqs_list)
    uniq_all_structures_dict = dict()

    for i, seq_str in enumerate(uniq_all_seqs_dict.keys()):
        # check if sequence contains noncanonical bases
        if re.search('[^ACGTU]', seq_str, flags=re.IGNORECASE):
            devprint('running rapid shapes ambiguos prediction #{}'.format(i), flush=True)
            print('for {}'.format(' '.join(uniq_all_seqs_dict[seq_str])))

            if params and ('rapidshapes' in params) and params['rapidshapes']:
                a_structure = rapidshapes_ambigous_prediction(seq_str, shape_str, params.get('rapidshapes'), shape_level=sh_level)
            else:
                a_structure = rapidshapes_ambigous_prediction(seq_str, shape_str, '', shape_level=sh_level)
            uniq_all_structures_dict[seq_str] = a_structure

            del a_structure

        else:
            devprint('rapidshapes - normal sequence #{}'.format(i), flush=True)
            tp, temp_fasta = mkstemp()
            with os.fdopen(tp, 'w') as f:
                f.write('>seq{}\n{}\n'.format(i, seq_str))

            if params and ('rapidshapes' in params) and params['rapidshapes']:
                rs_file = rapid_shapes_list(temp_fasta, shape_str, rapidshapes_params=params.get('rapidshapes'), shape_level=sh_level)
            else:
                rs_file = rapid_shapes_list(temp_fasta, shape_str, shape_level=sh_level)

            try:
                r_structure = read_rapid_shapes_list_outfile(rs_file)
            except IndexError as e:
                # catch error when structure is not predicted by rapidshapes despite being fed correct parameters and
                # no error is reported whatsoever
                r_structure = [
                    SeqRecord(
                        Seq(seq_str),
                        letter_annotations={shape_str: 'X'*len(seq_str)},
                    )
                ]
                print(e)
                print('Rapidshapes prediction failed - no structure predicted.')

            uniq_all_structures_dict[seq_str] = r_structure[0].letter_annotations[shape_str]

            del tp
            os.remove(temp_fasta)
            os.remove(rs_file)

    # rebuild seq list
    for seq in all_seqs_list:
        seq.annotations['sss'] = ['ss0']
        if 'X' in uniq_all_structures_dict[str(seq.seq)] and all(['X' == i for i in uniq_all_structures_dict[str(seq.seq)]]):
            seq.letter_annotations['ss0'] = '.'*len(seq)
            seq.annotations['predicted'] = False
        else:
            seq.letter_annotations['ss0'] = uniq_all_structures_dict[str(seq.seq)]
            seq.annotations['predicted'] = True

    os.remove(alig_file)
    os.remove(alif_file)
    return all_seqs_list


def rapidshapes_ambigous_prediction(sequence_string, shape_str, rs_params, shape_level):
    ml.debug(fname())
    # get all variants
    all_seq_vars = expand_ambiguous_RNA(sequence_string)

    structures = []
    # dum_str = ''.join(['.'] * len(all_seq_vars[0]))
    for i, var_seq in enumerate(all_seq_vars):
        tf, tempfile = mkstemp()
        with os.fdopen(tf, 'w') as f:
            f.write('>seq{}\n{}\n'.format(i, var_seq))

        # todo define exceptions
        try:
            rs_file = rapid_shapes_list(tempfile, shape_str, rs_params, shape_level=shape_level)
            rs_str = read_rapid_shapes_list_outfile(rs_file)
            os.remove(rs_file)
            structures.append(rs_str[0].letter_annotations[shape_str])
        except:
            # add empty structure
            structures.append(None)

        os.remove(tempfile)

    # count empty structures
    e_count = structures.count(None)
    if e_count != 0:
        print('rapidshapes ambiguos prediction - {} form {} were not successfully predicted'.format(
            e_count, len(structures)
        ))

    folded = [fs for fs in structures if fs is not None]

    # expect that all structures should be same
    if (len(folded) != 0) and all([True for i in folded if i == folded[0]]):
        print('rapidshapes ambigous - all structures same: success')

        structure = folded[0]
    else:

        if len(folded) == 0:
            print('rapidshapes ambiguos - prediction failed')
            structure = 'X' * len(sequence_string)
        else:
            print('rapidshapes ambigous - not all structure same: selecting structure most similar to others')
            comb_folded = itertools.combinations(folded, 2)
            comb_dist = rnadistance_one_thread(comb_folded)

            # reformat output distance
            od = np.zeros((len(folded), len(folded)))
            for i in range(len(folded)):
                for j in range(len(folded)):
                    if i == j:
                        continue
                    od[i, j] = comb_dist[i + j]

            sum_dists = od.sum()
            mini = sum_dists.index(min(sum_dists))
            structure = folded[mini]

    return structure


def expand_ambiguous_RNA(rnaseq):
    """
    expand possibilities in ambigously encoded RNA sequence
    with IUPAC code
    :return:
    """

    iupac = IUPACmapping()

    def spawn_branch(seqstr):
        ns = []
        for i, b, in enumerate(seqstr):
            if b in iupac.ambiguous:
                rs = ''.join(ns)
                bb = []
                for j in iupac.rnamapping[b]:
                    bb.append(rs + j + seqstr[i+1:])
                return bb
            else:
                ns.append(b)
        return []

    seqs = spawn_branch(rnaseq)
    while len(seqs) > 0:
        new_seqs = []
        for seq in seqs:
            new_seqs += spawn_branch(seq)

        if len(new_seqs) == 0:
            # comple enumeration in previous step
            break
        seqs = new_seqs

    return seqs


class IUPACmapping(object):
    """
    holder for the IUPAC ambigues information
    """
    def __init__(self):
        self.allowed_bases = tuple('ACGTURYSWKMBDHVN')
        self.ambiguous = tuple('RYSWKMBDHVN')
        self.unambiguous = tuple(set(self.allowed_bases) - set(self.ambiguous))
        self.mapping = {
            # 'A': ('A',),
            # 'C': ('C',),
            # 'G': ('G',),
            # 'T': ('T',),
            # 'U': ('U',),
            'R': ('A', 'G'),
            'Y': ('C', 'T'),
            'S': ('G', 'C'),
            'W': ('A', 'T'),
            'K': ('G', 'T'),
            'M': ('A', 'C'),
            'B': ('C', 'G', 'T'),
            'D': ('A', 'G', 'T'),
            'H': ('A', 'C', 'T'),
            'V': ('A', 'C', 'G'),
            'N': ('A', 'C', 'G', 'T')
        }
        self.rnamapping = self._RNAze()

    def _RNAze(self):
        """
        change the 'T' to 'U' in mapping
        :return:
        """
        rnamapping = dict()
        for key in self.mapping.keys():
            if 'T' in self.mapping[key]:
                rnamapping[key] = tuple([i for i in self.mapping[key] if i != 'T'] + ['U'])
            else:
                rnamapping[key] = self.mapping[key]
        return rnamapping


def uniquify_seq_list(seq_list):
    """
    return dictionary with sequences as keys and sequence ids as values in list
    :param seq_list:
    :return:
    """
    out_dict = dict()
    for seq in seq_list:
        seq_str = str(seq.seq)
        if seq_str in out_dict:
            out_dict[seq_str].append(seq.id)
        else:
            out_dict[seq_str] = []
            out_dict[seq_str].append(seq.id)
    return out_dict


def run_tcoffee_profile_aling(fasta_file, mode, profile_file, threads=None, outfile=None, tcoffee_rcoffee_params=''):
    """
    only aligner which can handle only one sequence input
    :return:
    """
    ml.debug(fname())
    # sanitize profile file for unalowed
    with open(profile_file, 'r') as pin, open(profile_file + '.san', 'w') as pout:
        for line in pin:
            pout.write(tcoffee_line_sanitizer(line))

    tr_params = '-profile {}'.format(profile_file + '.san') + tcoffee_rcoffee_params

    of = run_tcoffee(fasta_file,
                     mode=mode,
                     threads=threads,
                     outfile=outfile,
                     tcoffee_rcoffee_params=tr_params,
                     )
    return of


def tcoffee_line_sanitizer(line):
    return line.replace('|', '_').replace('@', '_').replace(':', '_').replace('%', '_').replace('/', '_')


def run_tcoffee(fasta_file, mode='rcoffee', threads=None, outfile=None, tcoffee_rcoffee_params=''):
    """
    runs t_coffee -mode rcoffee
    accepts t_coffee legal params
    :param fasta_file: fasta file input
    :param threads: how many threads to use default = all available
    :param outfile: output file
    :param tcoffee_rcoffee_params: legal params for t_coffe
    :param mode:
    :return: path to the outfile
    """
    ml.info('Running t-coffee.')
    ml.debug(fname())
    fd, fp = mkstemp()
    os.close(fd)

    exe_path = fp + 'dir'
    os.mkdir(exe_path)

    os.remove(fp)

    wd = os.getcwd()
    os.chdir(exe_path)
    FNULL = open(os.devnull, 'w')
    try:
        if not outfile:
            of, outfile = mkstemp()
            os.close(of)

        tfd, tcoffee_temp_files = mkstemp()
        os.close(tfd)

        if threads:
            cmd = '{}t_coffee {} -mode {} -n_core={} -outfile {} {} > {}'.format(
                CONFIG.tcoffee_path,
                fasta_file,
                mode,
                threads,
                outfile,
                tcoffee_rcoffee_params,
                tcoffee_temp_files
            )
        else:
            cmd = '{}t_coffee {} -mode {} -outfile {} {} > {}'.format(
                CONFIG.tcoffee_path,
                fasta_file,
                mode,
                outfile,
                tcoffee_rcoffee_params,
                tcoffee_temp_files
            )
        ml.debug(cmd)

        if ml.getEffectiveLevel() == 10:
            r = call(cmd, shell=True)
        else:
            r = call(cmd, shell=True, stderr=FNULL)

        # clean after t-coffee
        try:
            rmtree(exe_path)
        except OSError:
            ml.warning('cleanup after tcoffee unsuccessfull - Directory: {} not found.'.format(exe_path))

        remove_files_with_try(
            (tcoffee_temp_files, outfile + '.html'),
            msg_template='cleanup after tcoffee unsuccessfull -'
        )

        if r:
            msgfail = 'call to t_coffee -mode {} failed'.format(mode)
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)

    finally:
        FNULL.close()
        os.chdir(wd)

    return outfile


@timeit_decorator
@tryit_decorator
def rnafold_prediction(fasta2predict, params=''):
    ml.debug(fname())
    a, b = mkstemp()
    os.close(a)
    cmd = '{}RNAfold --noPS {} < {} > {}'.format(
        CONFIG.viennarna_path,
        params,
        fasta2predict,
        b
    )
    r = call(cmd, shell=True)

    if r:
        msgfail = 'call to rnafold failed, please check if rnafold is in path'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)

    structures = read_seq_str(b)
    os.remove(b)
    return structures


@timeit_decorator
@tryit_decorator
def subopt_fold_query(all_fasta_hits_file, query, params=None):
    """
    use folded query sequence as a reference,
    fold all sequences by Unafold, then select structure most similar to query

    accepted parameters:
    "rnafold"
    "mfold"

    :return:
    """
    ml.debug(fname())
    if params is None:
        params = dict()
    # get single query structure
    a, qf = mkstemp()
    with os.fdopen(a, 'w') as fd:
        fd.write('>query\n{}\n'.format(str(query.seq)))

    if params and ('RNAfold' in params) and params['RNAfold']:
        query_structure, _ = rnafold_prediction(qf, params=params['RNAfold'])
    else:
        query_structure, _ = rnafold_prediction(qf)

    if params and ('mfold' in params) and params['mfold']:
        assert isinstance(params['mfold'], (tuple, list)) and 3 == len(params['mfold'])
        subs = run_hybrid_ss_min(all_fasta_hits_file, mfold=params['mfold'])
    else:
        subs = run_hybrid_ss_min(all_fasta_hits_file)

    new_structures = []

    qs_string = query_structure[0].letter_annotations[query_structure[0].annotations['sss'][0]]
    # now compute rna distance score
    for seq in subs:
        str2compare = []
        key_list = []
        if seq.annotations['predicted']:
            for key in seq.annotations['sss']:
                str2compare.append((seq.letter_annotations[key], qs_string))
                key_list.append(key)
            rnadist_score = rnadistance_one_thread(str2compare)

            # select best ie lowes score
            mindisti = rnadist_score.index(min(rnadist_score))
            new_structures.append(SeqRecord(seq.seq,
                                            id=seq.id,
                                            annotations={'sss': ['ss0']},
                                            letter_annotations={'ss0': seq.letter_annotations[key_list[mindisti]]}
                                            ))
        else:
            new_structures.append(seq)

    os.remove(qf)

    return new_structures


@timeit_decorator
@tryit_decorator
def subopt_fold_alifold(all_fasta_hits_file, homologs_file, aligner='muscle', params=None, threads=None):
    """
    run clustal/muscle on selected homologs file
    :return:
    """
    ml.debug(fname())
    if params is None:
        params = dict()
    # run aligner
    # =================================================================================================================
    if 'clustalo' == aligner:
        if params and ('clustalo' in params) and params['clustalo']:
            clustal_params = ' --outfmt=clustal {}'.format(params['clustalo'])
        else:
            clustal_params = ' --outfmt=clustal'
        if threads:
            clustal_params += ' --threads={}'.format(threads)
        alig_file = compute_clustalo_clasic(homologs_file, clustalo_params=clustal_params)

    elif 'muscle' == aligner:
        if params and ('muscle' in params) and params['muscle']:
            alig_file = run_muscle(homologs_file, muscle_params=params['muscle'], reorder=False)
        else:
            alig_file = run_muscle(homologs_file, reorder=False)

    else:
        raise KeyError('provided key ({}) not recognized - avalible: "clustalo" "muscle"'.format(aligner))

    # run consensus prediction
    # =================================================================================================================
    if params and ('alifold' in params) and params['alifold']:
        alif_file = compute_alifold(alig_file, alifold_params=params['alifold'])
    else:
        alif_file = compute_alifold(alig_file)

    # possibly need to decode alifold structure
    alif_str = read_seq_str(alif_file)[0]
    consensus_structure = alif_str.letter_annotations['ss0']

    if params and ('mfold' in params) and params['mfold']:
        assert isinstance(params['mfold'], (tuple, list)) and 3 == len(params['mfold'])
        subs = run_hybrid_ss_min(all_fasta_hits_file, mfold=params['mfold'])
    else:
        subs = run_hybrid_ss_min(all_fasta_hits_file)

    new_structures = []

    # now compute rna distance score
    for seq in subs:
        str2compare = []
        key_list = []
        if seq.annotations['predicted']:
            for key in seq.annotations['sss']:
                str2compare.append((seq.letter_annotations[key], consensus_structure))
                key_list.append(key)
            rnadist_score = rnadistance_one_thread(str2compare)

            # select best ie lowes score
            mindisti = rnadist_score.index(min(rnadist_score))
            new_structures.append(SeqRecord(seq.seq,
                                            id=seq.id,
                                            annotations={'sss': ['ss0']},
                                            letter_annotations={'ss0': seq.letter_annotations[key_list[mindisti]]}
                                            ))
        else:
            new_structures.append(seq)
    os.remove(alif_file)
    os.remove(alig_file)
    return new_structures


def run_clustal_profile2seqs_align(msa_file, fasta_seq_file, clustalo_params='', outfile=None):
    """
    run clustal align MSA to seqs
    aligned columns in input MSA file are preserved and only new sequences are aligned and together they form new
     alignment
    :param msa_file: msa file (works with stockholm)
    :param fasta_seq_file: file with sequences to be aligned (format can be enforced with --infmt in clustalo_params)
    :param clustalo_params: params as accepted by clustalo
    :param outfile: outfile path, if not provided, tempfile will be created with output
    :return: outfile MSA path
    """
    ml.info('Runing clustalo profile.')
    ml.debug(fname())

    def _try_rescue(profile_file):
        # beware AlignIO truncates sequence names so they become non-unique, then clustalo also fails
        ml.warning('Trying rescue for profile alignment if profile ha no gaps clustalo things, '
                'that sequences are no aligned. Appendig trailing gap to overcome the issue.')
        a = AlignIO.read(profile_file, format='clustal')
        s = [SeqRecord(Seq(str(i.seq) + '-'), id=i.id) for i in a]
        fa = AlignIO.MultipleSeqAlignment(s)

        fd, temp = mkstemp()
        with os.fdopen(fd, 'w') as fh:
            # AlignIO.write(fa, fh, format='clustal')
            AlignIO.write(fa, fh, format='fasta')
        return temp

    if outfile:
        clustalo_file = outfile
    else:
        c_fd, clustalo_file = mkstemp()
        os.close(c_fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = '{}clustalo {} --force -i {} --profile1 {} -o {}'.format(
            CONFIG.clustal_path,
            clustalo_params,
            fasta_seq_file,
            msa_file,
            clustalo_file
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd, shell=True)
        else:
            r = call(cmd, shell=True, stdout=FNULL, stderr=FNULL)

        if r:
            ml.warning('Profile align failed.')

            # Initiate rescue attempt
            rewriten_msa = _try_rescue(msa_file)
            cmd2 = '{}clustalo {} --force -i {} --profile1 {} -o {}'.format(
                CONFIG.clustal_path,
                clustalo_params,
                fasta_seq_file,
                rewriten_msa,
                clustalo_file
            )
            ml.debug(cmd2)
            r2 = call(cmd2, shell=True)

            os.remove(rewriten_msa)

            if r2 != 0:
                msgfail = 'call to clustalo profile to sequences failed'
                ml.error(msgfail)
                ml.error(cmd)
                ml.error(cmd2)
                raise ChildProcessError(msgfail + ' ' + cmd)
    return clustalo_file


@timeit_decorator
@tryit_decorator
def turbofold_conservative_prediction(all_sequence_fasta, homologous_file, params=None):
    """
    specific turbofold wrapper
    tasks:
        write conf file
        read and convert output ct structures to db
        handle the files used

    Turbofold instalation
    classic turbofold installation is needed
        RNAstructure/exe must be in PATH

        then enviroment variable must be set, if not, then the RNAstucture (Turbofold) will not work
        export DATAPATH=[directory in which RNAstructure resides]/RNAstructure/data_tables/
    """

    ml.debug(fname())
    env = os.environ.copy()
    if 'DATAPATH' not in env:
        env['DATAPATH'] = CONFIG.rnastructure_datapath

    def _run_turbofold(turbofold_conf_file):
        """
        input is properly configured input file
        exepath is path for turbofold to be able to execute the function will cd to the provided directory and cd out
         after execution - usually the directory data_tables in RNAstructure package suffice
        all paths in conf file must be given in full
        """
        ml.info('Run Turbofold.')
        ml.debug(fname())
        with open(os.devnull, 'w') as FNULL:
            cmd = '{}TurboFold {}'.format(
                CONFIG.turbofold_path,
                turbofold_conf_file
            )
            ml.debug(cmd)
            if ml.getEffectiveLevel() == 10:
                r = call(cmd, shell=True, env=env)
            else:
                r = call(cmd, shell=True, stdout=FNULL, env=env)

            if r:
                msgfail = 'Call to turbofold failed, cmd below:'
                ml.error(msgfail)
                ml.error(cmd)
                raise ChildProcessError(msgfail)

    if (params is not None) and ('turbofold_mode' in params):
        turbo_mode = params['turbofold_mode']
    else:
        turbo_mode = 'MEA'

    tmpdir = mkdtemp()

    # read all sequences
    o = open(all_sequence_fasta, 'r')
    all_seqs = [i for i in SeqIO.parse(o, format='fasta')]
    o.close()

    h = open(homologous_file, 'r')
    homologs = [i for i in SeqIO.parse(h, format='fasta')]
    h.close()

    # write sequences by one in fasta file format
    hom_paths = []
    for i, hseq in enumerate(homologs):
        fastaname = os.path.join(tmpdir, 'hseq{}.fasta'.format(i))
        with open(fastaname, 'w') as f:
            f.write('>{}\n{}\n'.format(
                'hseq{}'.format(i),
                str(hseq.seq)
            ))
        hom_paths.append(fastaname)

    hom_out_paths = [os.path.join(tmpdir, 'hseq{}.ct'.format(i)) for i in range(len(homologs))]

    # ==================================================
    # prediction phase
    #  write predicted seq to a file
    #  write conf file
    #  run prediction
    new_structures = []
    for i, aseq in enumerate(all_seqs):
        fastaname = os.path.join(tmpdir, 'qa_seq{}.fasta'.format(i))
        o_name = os.path.join(tmpdir, 'qa_seq{}_out.ct'.format(i))
        # write the query fasta file
        with open(fastaname, 'w') as f:
            f.write('>{}\n{}\n'.format(
                'qaseq{}'.format(i),
                str(aseq.seq)
            ))

        # write the configuration file
        conf_file = os.path.join(tmpdir, 'conf_file{}.conf'.format(i))
        with open(conf_file, 'w') as c:
            c.write('Mode = {}\n'.format(turbo_mode))

            in_files = '{' + ';'.join(hom_paths) + ';' + fastaname + '}'
            out_files = '{' + ';'.join(hom_out_paths) + ';' + o_name + '}'

            c.write('InSeq = {}\n'.format(in_files))
            c.write('OutCT = {}\n'.format(out_files))

        _run_turbofold(conf_file)

        # convert ct structure 2 db
        o = open(o_name, 'r')
        seq_with_pred_str = ct2db(o, energy_txt='ENERGY')
        o.close()

        new_structures.append(
            SeqRecord(aseq.seq,
                      id=aseq.id,
                      annotations={'sss': ['ss0']},
                      letter_annotations={'ss0': seq_with_pred_str[0].letter_annotations['ss0']}
                      )
        )

    rmtree(tmpdir)
    return new_structures


def write_turbofold_confile(input_sequences, turbofold_params=None, cpus=None, outdir=None):
    ml.debug(fname())
    if outdir:
        tmpdir = outdir
    else:
        tmpdir = mkdtemp()

    params2write = {
        'Mode': 'MEA',
        'OutAln': os.path.join(tmpdir, 'Output.aln')
    }
    if turbofold_params:
        params2write.update(turbofold_params)

    # write sequences by one in fasta file format
    seq_paths = []
    for i, hseq in enumerate(input_sequences):
        fastaname = os.path.join(tmpdir, 'hseq{}.fasta'.format(i))
        with open(fastaname, 'w') as f:
            f.write('>{}\n{}\n'.format(
                'hseq{}'.format(i),
                str(hseq.seq)
            ))
        seq_paths.append(fastaname)

    hom_out_paths = [os.path.join(tmpdir, 'hseq{}.ct'.format(i)) for i in range(len(input_sequences))]

    conf_file = os.path.join(tmpdir, 'conf_file.conf')
    with open(conf_file, 'w') as c:
        for key in params2write.keys():
            c.write('{} = {}\n'.format(key, params2write[key]))

        if cpus and isinstance(cpus, int):
            c.write('Processors = {}'.format(cpus))

        in_files = '{' + ';'.join(seq_paths) + '}'
        out_files = '{' + ';'.join(hom_out_paths) + '}'

        c.write('InSeq = {}\n'.format(in_files))
        c.write('OutCT = {}\n'.format(out_files))

    return tmpdir, conf_file, hom_out_paths


@timeit_decorator
@tryit_decorator
def turbofold_only_homologous(all_sequences, nr_homologous, params):
    """
    Trubofold mode is MEA by default
    :param all_sequences:
    :param nr_homologous:
    :param params: dict with key: value structure for params to be used in turbofold
    :return:
    """
    # ==================================================
    # prediction phase
    # 1) homologous sequences
    #  write predicted seq to a file
    #  write conf file
    #  run prediction
    # 2) nonhomologous sequences
    #  predict with RNAfold or something similar, or predict them together with nr homologous sequence, but use only
    #   structures for nonhomologs from this step
    ml.debug(fname())

    nr_homologous_set = {str(seq.seq) for seq in nr_homologous}

    if len(nr_homologous_set) < 2:
        raise NoHomologousSequenceException

    nonhom_seqs = [seq for seq in all_sequences if str(seq.seq) not in nr_homologous_set]

    hom_seqs_structures = run_turbofold(nr_homologous, params=params)
    hom_seqs_dict = {str(seq.seq): seq.letter_annotations['ss0'] for seq in hom_seqs_structures}

    # fold nonhomologous seqs together with homologous
    nonhom_seqs_structures = run_turbofold(nonhom_seqs + nr_homologous, params=params)
    nonhom_seqs_dict = {str(seq.seq): seq.letter_annotations['ss0'] for seq in nonhom_seqs_structures}

    # merge the dicts together
    # we need to merge in this particular order to ensure that structures for homologous sequences from folding with
    #  nonhomologous sequences are replaced by structures from folding only with homologous seqs

    nonhom_seqs_dict.update(hom_seqs_dict)

    # check if the update run correctly
    assert all([True if hom_seqs_dict[hs] == nonhom_seqs_dict[hs] else False for hs in hom_seqs_dict.keys()])

    predicted = []
    for seq in all_sequences:
        seq.letter_annotations['ss0'] = nonhom_seqs_dict[str(seq.seq)]
        predicted.append(seq)

    return predicted


def run_turbofold(sequences, params):
    ml.info('Running Turbofold.')
    ml.debug(fname())

    env = os.environ.copy()
    if 'DATAPATH' not in env:
        env['DATAPATH'] = CONFIG.rnastructure_datapath

    tmpdir, con_file, output_structure_files = write_turbofold_confile(
        input_sequences=sequences,
        turbofold_params=params,
    )

    # run without prediction progress reporting output
    with open(os.devnull, 'w') as FNULL:
        cmd = [
            '{}TurboFold'.format(CONFIG.turbofold_path),
            con_file
        ]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd, env=env)
        else:
            r = call(cmd, stdout=FNULL, env=env)

        if r:
            msgfail = 'Call to turbofold failed, cmd below:'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)

        # now convert ct files produced by TurboFold
        new_structures = []
        for out_str_file, orig_seq in zip(output_structure_files, sequences):
            o = open(out_str_file, 'r')
            seq_with_pred_str = ct2db(o, energy_txt='ENERGY')
            o.close()

            assert str(seq_with_pred_str[0].seq) == str(orig_seq.seq)

            new_structures.append(
                SeqRecord(
                    orig_seq.seq,
                    id=orig_seq.id,
                    annotations={'sss': ['ss0']},
                    letter_annotations={'ss0': seq_with_pred_str[0].letter_annotations['ss0']}
                )
            )

        rmtree(tmpdir)
        return new_structures


def check_lonely_bp(structure, gap_char='.'):
    """
    check lonely bp in classic dot bracket structure notation
    """
    match = re.search('\.\(\.|\.\)\.', structure)
    if not match:
        return structure
    # print('lp found')
    s = match.start() + 1
    gapmark = 49
    encoded = encode_structure_unicode(structure, gap_mark=gapmark)
    violating = encoded[s]
    repaired = encoded.replace(violating, chr(gapmark))

    repaired_structure = decode_structure_unicode(repaired, gap_char=gap_char, gap_mark=gapmark)
    # go into another round
    return check_lonely_bp(repaired_structure)


if __name__ == '__main__':
    print('encode structure unicode example')
    sample_structure = '(((((((((....((((((((..(((..((((..((....))..))))..)))...)))))).))(((((((.....(((((.(((....))))))))....)))))))))))))))).'
    encoding_example = encode_structure_unicode(sample_structure)
    print(sample_structure)
    print(encoding_example)


