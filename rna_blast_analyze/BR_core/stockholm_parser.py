import re
from warnings import warn

from rna_blast_analyze.BR_core.db2shape import nesting
from Bio.AlignIO import read as clustal_read
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.stockholm_alig import StockholmFeatureStock, StockholmAlig


def parse_single_stockholm_alig(file_handle, features=None):
    """
    universal parser for stockholm files
    parser cannot parse white spaces in letter annotations
    :param file_handle:
    :param features: enables user to pass custom parsing parameters
    :return:
    """

    if not hasattr(file_handle, 'readline'):
        raise AttributeError('File handle does not have readline attribute')

    # reference by wikipedia https://en.wikipedia.org/wiki/Stockholm_format
    # #=GF <feature> <Generic per-File annotation, free text>
    # #=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
    # #=GS <seqname> <feature> <Generic per-Sequence annotation, free text>
    # #=GR <seqname> <feature> <Generic per-Residue annotation, exactly 1 char per residue>
    # aligned sequences are not prepended by any means

    if not features:
        features = StockholmFeatureStock()

    sequence_dict = {}
    seq_letter_ann = {}
    alignment_annotations = {}
    alig_letter_annotations = {}
    seq_annotation = {}
    seq_order = []

    # read annotations, but must recognize parts which belongs to aligned sequences
    txt = file_handle.readline()
    if 'STOCKHOLM' not in txt:
        raise Exception('File does not appear to be STOCKHOLM file (first line contain "STOCKHOLM")')
    txt = file_handle.readline()

    while txt:
        if txt.strip() == '':
            txt = file_handle.readline()
            continue

        if txt.strip() == '//':
            # yield or return here
            # // marks end of alignment but does not necessarily mean EOF
            return sequence_dict, seq_letter_ann, alignment_annotations, \
                   alig_letter_annotations, seq_annotation, seq_order

        if txt[:4] == '#=GF':
            # annotation
            # multi-line annotations with same call symbol will be merged together
            # #=GF <feature> <Generic per-File annotation, free text>
            tline = txt[4:].split()
            alig_spec = tline[0]
            ann_text = ' '.join(tline[1:])

            if alig_spec in features.GF:
                if alig_spec in alignment_annotations:
                    alignment_annotations[alig_spec] += ann_text.strip()
                else:
                    alignment_annotations[alig_spec] = ann_text.strip()
            else:
                warn('Unrecognized annotation {}, omitting:\n{}'.format(alig_spec, txt))

        elif txt[:4] == '#=GS':
            # per sequence annotation
            # #=GS <seqname> <feature> <Generic per-Sequence annotation, free text>

            tline = txt[4:].split()

            seqid = tline[0]
            alig_spec = tline[1]
            seq_l = ' '.join(tline[2:])

            if seqid not in seq_annotation:
                seq_annotation[seqid] = {}

            if alig_spec in features.GS:
                if alig_spec in seq_annotation[seqid]:
                    seq_annotation[seqid][alig_spec] += seq_l
                else:
                    seq_annotation[seqid][alig_spec] = seq_l
            else:
                warn('Unrecognized sequence specific annotation {}, omitting:\n{}'.format(tline[1], txt))

        elif txt[:4] == '#=GC':
            # generic per column annotation
            # #=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
            tline = txt[4:].split()
            alig_spec = tline[0]
            letter_ann = tline[1]

            if alig_spec in features.GC:
                if alig_spec in alig_letter_annotations:
                    alig_letter_annotations[alig_spec] += letter_ann.strip()
                else:
                    alig_letter_annotations[alig_spec] = letter_ann.strip()
            else:
                warn('Unrecognized annotation {}, omitting {}'.format(alig_spec, txt))

        elif txt[:4] == '#=GR':
            # per residue annotation - sequence specific
            # there can be multiple per residue annotations
            tline = txt[4:].split()
            if len(tline) != 3:
                print('Except: expected letter annotation line in format:'
                      '#=GR seq_id ann_specifier letter_annotation')
                raise Exception('do not understand this alignment line\n{}'.format(txt))
            seqid = tline[0]
            alig_spec = tline[1]
            seq_l = tline[2]

            if seqid not in seq_letter_ann:
                seq_letter_ann[seqid] = {}

            if alig_spec in features.GR:
                if alig_spec in seq_letter_ann[seqid]:
                    seq_letter_ann[seqid][alig_spec] += seq_l
                else:
                    seq_letter_ann[seqid][alig_spec] = seq_l
            else:
                warn('Unrecognized letter_annotation {}, omitting {}'.format(tline[1], txt))

        else:
            # txt is sequence | sequence must be determined before its annotation
            # split on space - first is the id, then the seq alignment
            tline = txt.split()
            if len(tline) != 2:
                raise Exception('do not understand this alignment line\n{}'.format(txt))
            seqid = tline[0]
            seq_l = tline[1]
            if seqid in sequence_dict:
                sequence_dict[seqid] += seq_l
            else:
                sequence_dict[seqid] = seq_l
                seq_order.append(seqid)

        # read next line
        txt = file_handle.readline()

    # if this part is reached, raise error
    # if this is EOF all would be empty
    if sequence_dict or seq_letter_ann or alignment_annotations or alig_letter_annotations or seq_annotation:
        raise Exception('Unexpected EOF, alignemnt broken, termination signal not found: missing: "//"')


def process_parsed_data(sequences, seq_letter_ann, alignment_annotations,
                        alig_letter_annotations, seq_annotation, order):

    final_alig = StockholmAlig()
    final_alig.annotations = alignment_annotations
    final_alig.column_annotations = alig_letter_annotations

    for seqid in order:
        # construct Seq object
        if seqid in seq_letter_ann:
            appendable_letter_annotations = seq_letter_ann[seqid]
        else:
            appendable_letter_annotations = {}
        if seqid in seq_annotation:
            appendable_annotations = seq_annotation[seqid]
        else:
            appendable_annotations = {}
        seq_alig_line = SeqRecord(Seq(sequences[seqid], SingleLetterAlphabet),
                                  id=seqid,
                                  letter_annotations=appendable_letter_annotations,
                                  annotations=appendable_annotations)
        final_alig.append(seq_alig_line)

    return final_alig


def stockholm_parse(file_handle, features=None):
    """
    return iterator, there can be multiple alignment object in alignment file
    :param file_handle:
    :param features:
    :return:
    """
    a, b, c, d, e, f = parse_single_stockholm_alig(file_handle, features=features)
    one_alig = process_parsed_data(a, b, c, d, e, f)
    yield one_alig


def stockholm_read(file_handle, features=None):
    """
    returns StockholmAlig object, only one alignment from file will be returned
    If multiple alignment are present, they will be dropped
    :param file_handle:
    :param features:
    :return:
    """
    a, b, c, d, e, f = parse_single_stockholm_alig(file_handle, features=features)
    one_alig = process_parsed_data(a, b, c, d, e, f)
    return one_alig


def read_st(st):
    with open(st, 'r') as f:
        return stockholm_read(f)


def _clustalize_alignment(alig):
    new_alig = StockholmAlig()
    for al in alig:
        ns = SeqRecord(Seq(re.sub('[.~]', '-', str(al.seq))),
                       id=al.id)
        new_alig.append(ns)
    return new_alig


def clustal2stockholm(file):
    with open(file, 'r') as f:
        cl = clustal_read(f, 'clustal')
        st = StockholmAlig()
        st._records = list(cl)
        return st


def trim_cmalign_sequence_by_refseq(alignment, rf='RF', rs='RS', lt='PP', ref_name='query'):
    """
    this function takes cmalign like input and returns list of alignments

    return instance of stockholm alig? need to preprocess the output further to complete seqs

    :param alignment:
    :param rf: key for reference sequnence in stockhlolm cmalign input alignment
    :param rs: key for reference structure in stockholm cmalign input alignment
    :param lt: key for per sequence aligned probabilities
    :param ref_name: name of reference sequence used in result alignment output
    :return:
    """
    if rf not in alignment.column_annotations.keys():
        raise KeyError('Provided key "{}" not found in column_annotations'.format(rf))
    if rs not in alignment.column_annotations.keys():
        raise KeyError('Provided key "{}" not found in column_annotations'.format(rs))

    for aligned_seq in alignment:

        if lt not in aligned_seq.letter_annotations.keys():
            raise KeyError('Provided key "{}" not found in column_annotations'.format(lt))

        seq, pp, refseq, refstr = hmmlike_alig_line2uninserted_seqs(str(aligned_seq.seq),
                                                                    aligned_seq.letter_annotations[lt],
                                                                    alignment.column_annotations[rf],
                                                                    alignment.column_annotations[rs])

        new_rec = SeqRecord(Seq(seq),
                            id=aligned_seq.id,
                            description=aligned_seq.description,
                            letter_annotations={lt: pp})

        new_ref = SeqRecord(Seq(refseq),
                            id=ref_name,
                            description='CM model reference sequence')

        # rehandle this so the output is alignment
        single_alig = StockholmAlig()
        single_alig.append(new_ref)
        single_alig.append(new_rec)
        single_alig.column_annotations['SS_cons'] = refstr
        yield single_alig


def trim_cmalign_sequence_by_refseq_one_seq(alignment, rf='RF', rs='RS', lt='PP', convert2uppercase=False):
    """
    this function takes cmalign like input and returns list of alignments

    return instance of stockholm alig? need to preprocess the output further to complete seqs

    :param alignment:
    :param rf: key for reference sequnence in stockhlolm cmalign input alignment
    :param rs: key for reference structure in stockholm cmalign input alignment
    :param lt: key for per sequence aligned probabilities
    :return:
    """
    if rf not in alignment.column_annotations.keys():
        raise KeyError('Provided key "{}" not found in column_annotations'.format(rf))
    if rs not in alignment.column_annotations.keys():
        raise KeyError('Provided key "{}" not found in column_annotations'.format(rs))

    for aligned_seq in alignment:

        if lt not in aligned_seq.letter_annotations.keys():
            raise KeyError('Provided key "{}" not found in column_annotations'.format(lt))

        seq, pp, refseq, refstr = hmmlike_alig_line2uninserted_seqs(str(aligned_seq.seq),
                                                                    aligned_seq.letter_annotations[lt],
                                                                    alignment.column_annotations[rf],
                                                                    alignment.column_annotations[rs])

        if convert2uppercase:
            seq = seq.upper()
            refseq = refseq.upper()

        new_rec = SeqRecord(Seq(seq),
                            id=aligned_seq.id,
                            description=aligned_seq.description,
                            letter_annotations={lt: pp})

        # new_ref = SeqRecord(Seq(refseq),
        #                     id=ref_name,
        #                     description='CM model reference sequence')

        # rehandle this so the output is alignment
        single_alig = StockholmAlig()
        # single_alig.append(new_ref)
        single_alig.append(new_rec)
        single_alig.column_annotations['SS_cons'] = refstr
        single_alig.column_annotations['RF'] = refseq
        yield single_alig


def hmmlike_alig_line2uninserted_seqs(*args, gapchars=('.', '~')):
    """
    function reads the sequence string, score string, refseq, string, and refseq structure from cmalign output
    the '.' here means MSA insertion, which is taken out from all 4 input files
    the dashes are kept - they point to places where there is some sequence missing either in sequence or refseq
    no inference is done here
    :param args: sequences (strings) from cmalign to be searched through, passed as separate arguments
    :param gapchars: chracters that are considered gaps - ie if they are present in passed sequences (strings), that
                        column will be removed
    :return: sequences (strings) without excessive alignment columns
    """
    assert len({len(i) for i in args}) == 1
    # process
    joined_list = []
    for i in zip(*args):
        if all([any([j == ch for ch in gapchars]) for j in i]):
            continue
        else:
            joined_list += i

    # reformat
    res_seq = []
    for i in range(len(args)):
        res_seq.append(''.join(joined_list[i::len(args)]))

    return res_seq


def trim_alignment_by_sequence(alignment, aligned_query_seq, structure_annotation=''):
    """
    for locarna seq global - and presumably other aligners there is a need for infering where sequences starts
    and ends
    ie trimming the sequence by another one

    for locarna global, that is taking the query and trim by its aligned portions
    - issue with this is that poorly aligned portions at the ends of alignment have great impact on the start
    and end position of the sequence and also the result sequence cannot be shorter then trimming sequence

    structure trim only works for structures in dot bracket format
    :param aligned_query_seq:
    :return:
    """
    # get seq_inds
    start_match = re.search('^[-.]+', aligned_query_seq)
    end_match = re.search('[-.]+$', aligned_query_seq)
    if start_match:
        start_seq = start_match.span()[1]
    else:
        start_seq = 0

    if end_match:
        end_seq = end_match.span()[0]
    else:
        end_seq = len(aligned_query_seq)

    # compute also for structure
    if structure_annotation != '':
        ta = alignment.column_annotations[structure_annotation][start_seq:end_seq]

        repaired_structure = repair_trimmed_structure(ta)

        ra = alignment.slice_alig(start_seq, end_seq)

        ra.column_annotations[structure_annotation] = repaired_structure
        return ra
    else:
        return alignment.slice_alig(start_seq, end_seq)


def repair_trimmed_structure(ta):
    """
    repair structure - only for use with locarna
    :param ta:
    :return:
    """
    position, c = missing_bracket(ta)
    if position < len(ta) - 1 and c != 1:
        # missing on start
        ta = ta[:position] + '.' + ta[position + 1:]
        # print('structure were repaired at pos {}'.format(position))

    if position == len(ta) - 1 and c != 1:
        # missing on end
        # bud it can be different than
        # tb = ta[::-1]
        # new_pos = tb.index(')')
        # rs = tb[:new_pos] + '.' + tb[new_pos + 1:]
        # ta = rs[::-1]

        # do not count gaps
        _, str_ann = nesting(ta, gap_mark='G')
        set_of_pairs = set(str_ann)
        set_of_pairs.discard('G')

        # for i in sorted(set_of_pairs):
        #     if str_ann.count(i) % 2 == 1:
        #         miss_nest = i
        #         break

        miss_nest = _return_from_sorted_pairs(set_of_pairs, str_ann)

        # missing on the end
        r_ann = str_ann[::-1]
        r_ta = ta[::-1]
        mi = r_ann.index(miss_nest)
        qq = r_ta[:mi] + '.' + r_ta[mi + 1:]
        ta = qq[::-1]
        # print('structure were repaired at pos {}'.format(len(ta) - mi))

    _, c = missing_bracket(ta)
    if c != 1:
        ta = repair_trimmed_structure(ta)
    return ta


def _return_from_sorted_pairs(_set_of_pairs, _str_ann):
    for i in sorted(_set_of_pairs):
        if _str_ann.count(i) % 2 == 1:
            return i
    else:
        raise StopIteration


def missing_bracket(seq, br=('(', ')')):
    """br must be nesting brackets of the kind you want to analyze"""
    c = 1
    i = 0
    for i, char in enumerate(seq):
        if char == br[0]:
            c += 1
        elif char == br[1]:
            c -= 1
        if c == 0:
            return i, c
    return i, c


if __name__ == '__main__':

    F = StockholmFeatureStock()
    F.add_custom_parser_tags('GC', {'AB': 'asdkajsd'})

    with open('stockholm_test.sto', 'r') as f:
        stockholm_parse(f, features=F)