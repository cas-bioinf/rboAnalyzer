import copy
import logging
import os
import re
import unicodedata
from subprocess import call
from tempfile import mkstemp
import shlex
import multiprocessing

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.BA_support import read_seq_str, sanitize_fasta_names_in_seqrec_list, desanitize_fasta_names_in_seqrec_list, run_muscle
from rna_blast_analyze.BR_core.hybrid_ss_min import run_hybrid_ss_min
from rna_blast_analyze.BR_core.alifold4all import compute_refold, compute_clustalo_clasic, compute_alifold
from rna_blast_analyze.BR_core.cmalign import build_stockholm_from_clustal_alig, run_cmalign_on_fasta, cm_strucutre2br
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.db2shape import nesting
from rna_blast_analyze.BR_core.decorators import timeit_decorator
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.infer_homology import alignment_column_conservation
from rna_blast_analyze.BR_core.par_distance import compute_distances
from rna_blast_analyze.BR_core.stockholm_parser import read_st, trim_cmalign_sequence_by_refseq_one_seq

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

    :param stockholm_msa:
    :return:
    """

    st_msa = read_st(stockholm_msa)

    # identify consensus unpaired nucleotides
    consensus_conservation_score = alignment_column_conservation(st_msa, gap_chars='-')

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

    fd, tempfile = mkstemp(prefix='rba_', suffix='_31', dir=CONFIG.tmpdir)
    with os.fdopen(fd, 'w') as fh:
        fh.write('>seq01\n{}\n{}\n'.format(
            str(seqr.seq),
            seqr.letter_annotations[constraints_id]
        ))
    structure = rnafold_prediction(tempfile, params='-C')
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

    # this is to ensure that all chars returned are printable
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


def sanitize_fasta_file(infasta, used_dict=None):
    fd, out_path = mkstemp(prefix='rba_', suffix='_32', dir=CONFIG.tmpdir)
    with open(infasta, 'r') as inh, os.fdopen(fd, 'w') as outh:
        k = [i for i in SeqIO.parse(inh, format='fasta')]
        san_seqs, san_dict = sanitize_fasta_names_in_seqrec_list(k, used_dict=used_dict)

        SeqIO.write(san_seqs, outh, format='fasta')
        return out_path, san_dict


@timeit_decorator
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
        raise Exception('refold procedure not recognized: {}, possible values are {}'.format(
            refold, ' '.join(ref_pred))
        )

    cl_file = _aligner_block(nr_path, params, msa_alg, threads)

    # cannot rely on that, the order of a cl_file would be the same as the order of the nr_homolog_hits_file
    ali_file = compute_alifold(cl_file, alifold_params=params.get('alifold', ''))

    consensus_record = read_seq_str(ali_file)[0]

    clustalo_profile_params = '--outfmt clustal '
    clustalo_profile_params += params.get('clustalo_profile', '')
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
    new_consensus_structure_encoded = _repair_consensus_structure_by_maping(
        cs_encode,
        mapping,
        len(match_original_seq_in_new_alig.seq),
        gap_char=49
    )
    new_consensus_structure_repaired = repair_structure_any_variant(new_consensus_structure_encoded)

    new_consensus_structure = decode_structure_unicode(new_consensus_structure_repaired)

    new_consensus_sequence = _repair_consensus_structure_by_maping(
        str(consensus_record.seq),
        mapping,
        len(match_original_seq_in_new_alig.seq),
        gap_char=ord('_')
    )

    # write new consensus to a file
    a_fd, new_alifold_consensus_file = mkstemp(prefix='rba_', suffix='_33', dir=CONFIG.tmpdir)
    with os.fdopen(a_fd, 'w') as f:
        f.write(new_consensus_sequence + '\n')
        f.write(new_consensus_structure + '\n')

    # write sliced alignment to a file
    sa_fd, sliced_alignment_file = mkstemp(prefix='rba_', suffix='_34', dir=CONFIG.tmpdir)
    with os.fdopen(sa_fd, 'w') as f:
        AlignIO.write(new_alig_for_refold, f, 'clustal')

    # now process the file, and map alignment to consensus structure
    if refold in ['refold', 'refold_rnafoldc']:
        refold_file = compute_refold(sliced_alignment_file, new_alifold_consensus_file)

        if refold == 'refold_rnafoldc':
            rnafold_parameters = params.get('RNAfold', '')
            if '-C' not in rnafold_parameters:
                rnafold_parameters += ' -C'

            seq_str = rnafold_prediction(refold_file, params=rnafold_parameters)

        else:
            seq_str = read_seq_str(refold_file)

        os.remove(refold_file)

    else:
        st_alig_file = build_stockholm_from_clustal_alig(sliced_alignment_file, new_alifold_consensus_file)
        repred_tr = str(params.get('repred_unpaired_tr', '9'))
        conseq_conserved = params.get('conseq_conserved', 1)

        seq_str = _refold_with_unpaired_conservation(
            st_alig_file,
            repred_tr=repred_tr,
            conseq_conserved=conseq_conserved
        )
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
def cmmodel_rnafold_c(allhits_fasta, cmmodel_file, threads=None, params=None):
    ml.debug(fname())
    if params is None:
        params = dict()

    allhits_fasta_file, san_dict = sanitize_fasta_file(allhits_fasta)

    cmalign_params = ''
    if threads:
        cmalign_params += '--cpu {}'.format(threads)

    if 'cmalign' in params and params['cmalign']:
        cmalign_params += ' ' + params['cmalign']

    if '--notrunc' not in cmalign_params:
        cmalign_params += ' --notrunc'

    # rnafold params
    rnafold_params = params.get('RNAfold', '-C')
    assert isinstance(rnafold_params, str)
    if '-C' not in rnafold_params:
        # some parameters given but -C not present
        rnafold_params += ' -C'

    alig_file = run_cmalign_on_fasta(allhits_fasta_file, cmmodel_file, cmalign_params=cmalign_params)
    os.remove(allhits_fasta_file)
    # multiple sequence cm align
    # split by sequence, then run the rest
    cm_alig = read_st(alig_file)
    os.remove(alig_file)

    structures = []
    for single_alig in trim_cmalign_sequence_by_refseq_one_seq(cm_alig, rs='SS_cons', convert2uppercase=True):
        out_alig, trimmed_seq = trim_and_repair_single_cm_alignment(single_alig)
        conserved_structure_pairs = find_nc_and_remove(str(trimmed_seq.seq), trimmed_seq.letter_annotations['dec_str'])
        trimmed_seq.letter_annotations['constrains'] = conserved_structure_pairs

        # constraint prediction
        # write constraint file
        fd, temp_constraint_file = mkstemp(prefix='rba_', suffix='_41', dir=CONFIG.tmpdir)
        with os.fdopen(fd, 'w') as tmpf:
            tmpf.write('>{}\n{}\n{}\n'.format(
                trimmed_seq.id,
                str(trimmed_seq.seq),
                trimmed_seq.letter_annotations['constrains']
            ))

        single_structure = rnafold_prediction(temp_constraint_file, params=rnafold_params)

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
        clustal_params = '--outfmt=clustal --force'
        clustal_params += params.get('clustalo', '')

        if threads:
            clustal_params += ' --threads={}'.format(threads)
        alig_file = compute_clustalo_clasic(nr_homolog_hits_file, clustalo_params=clustal_params)

    elif msa_alg == 'muscle':
        if params and ('muscle' in params) and params['muscle']:
            alig_file = run_muscle(nr_homolog_hits_file, muscle_params=params['muscle'], reorder=True)
        else:
            alig_file = run_muscle(nr_homolog_hits_file, reorder=True)

    else:
        print('invalig MSA alg chosen {}, valid are "clustalo" and "muscle"'.format(msa_alg))
        raise AttributeError()

    return alig_file


@timeit_decorator
def rfam_subopt_pred(all_sequence_fasta, cm_ref_str, params=None, threads=1):
    ml.debug(fname())
    if params is None:
        params = dict()

    if params and ('mfold' in params) and params['mfold']:
        assert isinstance(params['mfold'], (tuple, list)) and 3 == len(params['mfold'])
        subs = run_hybrid_ss_min(all_sequence_fasta, mfold=params['mfold'], threads=threads)
    else:
        subs = run_hybrid_ss_min(all_sequence_fasta, threads=threads)

    # now compute rna distance score
    if threads == 1:
        new_structures = []
        for seq in subs:
            new_structures.append(_helper_subopt(seq, cm_ref_str))
    else:
        with multiprocessing.Pool(processes=threads) as pool:
            tuples = [(seq, cm_ref_str) for seq in subs]
            new_structures = pool.starmap(_helper_subopt, tuples)

    return new_structures


@timeit_decorator
def rnafold_wrap_for_predict(*args, **kwargs):
    return rnafold_prediction(*args, **kwargs)


def rnafold_prediction(fasta2predict, params=''):
    ml.debug(fname())
    fd, structure_output_file = mkstemp(prefix='rba_', suffix='_54', dir=CONFIG.tmpdir)
    os.close(fd)

    structure_output_file = rnafold(fasta2predict, structure_output_file, params)

    structures = read_seq_str(structure_output_file)
    os.remove(structure_output_file)
    return structures


def rnafold(fastafile, outfile, parameters=''):
    """Run RNAfold on commandline

    issue --noPS to prevent plotting of postscript structure
    """
    if '--noPS' not in parameters:
        parameters += ' --noPS'
    cmd = '{} {} < {} > {}'.format(
        shlex.quote('{}RNAfold'.format(CONFIG.viennarna_path)),
        ' '.join([shlex.quote(i) for i in parameters.split()]),
        shlex.quote(fastafile),
        shlex.quote(outfile)
    )
    r = call(cmd, shell=True)

    if r:
        msgfail = 'call to rnafold failed, please check if rnafold is in path'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)
    return outfile


@timeit_decorator
def subopt_fold_query(all_fasta_hits_file, query_file, params=None, threads=1):
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
    query_structure = rnafold_prediction(query_file, params.get('RNAfold', ''))

    if params and ('mfold' in params) and params['mfold']:
        assert isinstance(params['mfold'], (tuple, list)) and 3 == len(params['mfold'])
        subs = run_hybrid_ss_min(all_fasta_hits_file, mfold=params['mfold'], threads=threads)
    else:
        subs = run_hybrid_ss_min(all_fasta_hits_file, threads=threads)

    qs_string = query_structure[0].letter_annotations[query_structure[0].annotations['sss'][0]]

    # now compute rna distance score
    if threads == 1:
        new_structures = []
        for seq in subs:
            new_structures.append(_helper_subopt(seq, qs_string))
    else:
        with multiprocessing.Pool(processes=threads) as pool:
            tuples = [(seq, qs_string) for seq in subs]
            new_structures = pool.starmap(_helper_subopt, tuples)

    return new_structures


@timeit_decorator
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
        clustal_params = ' --outfmt=clustal --force'
        clustal_params += params.get('clustalo', '')

        if threads:
            clustal_params += ' --threads={}'.format(threads)
        alig_file = compute_clustalo_clasic(homologs_file, clustalo_params=clustal_params)

    elif 'muscle' == aligner:
        alig_file = run_muscle(homologs_file, muscle_params=params.get('muscle', ''), reorder=False)

    else:
        raise KeyError('provided key ({}) not recognized - avalible: "clustalo" "muscle"'.format(aligner))

    # run consensus prediction
    # =================================================================================================================
    alif_file = compute_alifold(alig_file, alifold_params=params.get('alifold', ''))

    # possibly need to decode alifold structure
    alif_str = read_seq_str(alif_file)[0]
    consensus_structure = alif_str.letter_annotations['ss0']

    subs = run_hybrid_ss_min(all_fasta_hits_file, mfold=params.get('mfold', (10, 2, 20)), threads=threads)

    # now compute rna distance score
    if threads == 1:
        new_structures = []
        for seq in subs:
            new_structures.append(_helper_subopt(seq, consensus_structure))
    else:
        with multiprocessing.Pool(processes=threads) as pool:
            tuples = [(seq, consensus_structure) for seq in subs]
            new_structures = pool.starmap(_helper_subopt, tuples)

    os.remove(alif_file)
    os.remove(alig_file)
    return new_structures


def _helper_subopt(seq, consensus_structure):
    if seq.annotations['predicted']:
        str2compare = []
        key_list = []
        for key in seq.annotations['sss']:
            str2compare.append((seq.letter_annotations[key], consensus_structure))
            key_list.append(key)
        rnadist_score = compute_distances(str2compare)

        # select best ie lowes score
        mindisti = rnadist_score.index(min(rnadist_score))
        return SeqRecord(
                seq.seq,
                id=seq.id,
                annotations={'sss': ['ss0']},
                letter_annotations={'ss0': seq.letter_annotations[key_list[mindisti]]}
            )
    else:
        return seq


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
        ml.warning(
            'Trying rescue for profile alignment if profile has no gaps, sequences appears not aligned. '
            'Appending trailing gap to overcome the issue.'
        )
        a = AlignIO.read(profile_file, format='clustal')
        s = [SeqRecord(Seq(str(i.seq) + '-'), id=i.id) for i in a]
        fa = AlignIO.MultipleSeqAlignment(s)

        fd, temp = mkstemp(prefix='rba_', suffix='_56', dir=CONFIG.tmpdir)
        with os.fdopen(fd, 'w') as fh:
            AlignIO.write(fa, fh, format='fasta')
        return temp

    if outfile:
        clustalo_file = outfile
    else:
        c_fd, clustalo_file = mkstemp(prefix='rba_', suffix='_57', dir=CONFIG.tmpdir)
        os.close(c_fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = [
            '{}clustalo'.format(CONFIG.clustal_path),
            '--force',
            '-i', fasta_seq_file,
            '--profile1', msa_file,
            '-o', clustalo_file
        ]
        if clustalo_params != '':
            cmd += clustalo_params.split()

        # cmd = '{}clustalo {} --force -i {} --profile1 {} -o {}'.format(
        #     CONFIG.clustal_path,
        #     clustalo_params,
        #     fasta_seq_file,
        #     msa_file,
        #     clustalo_file
        # )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL, stderr=FNULL)

        if r:
            ml.warning('Profile align failed.')

            # Initiate rescue attempt
            rewriten_msa = _try_rescue(msa_file)
            cmd2 = [
                '{}clustalo'.format(CONFIG.clustal_path),
                '--force',
                '-i', fasta_seq_file,
                '--profile1', rewriten_msa,
                '-o', clustalo_file
            ]
            if clustalo_params:
                cmd2 += clustalo_params.split()
            # cmd2 = '{}clustalo {} --force -i {} --profile1 {} -o {}'.format(
            #     CONFIG.clustal_path,
            #     clustalo_params,
            #     fasta_seq_file,
            #     rewriten_msa,
            #     clustalo_file
            # )
            ml.debug(cmd2)
            r2 = call(cmd2)

            os.remove(rewriten_msa)

            if r2 != 0:
                msgfail = 'call to clustalo profile to sequences failed'
                ml.error(msgfail)
                ml.error(cmd)
                ml.error(cmd2)
                raise ChildProcessError(msgfail + ' ' + cmd)
    return clustalo_file