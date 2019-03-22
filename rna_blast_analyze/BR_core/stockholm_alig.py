import re
from Bio.AlignIO.ClustalIO import ClustalWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class StockholmFeatureStock(object):
    def __init__(self):
        # GF
        feature_annotations = {'AC': 'Accession number',
                               'ID': 'Identification',
                               'DE': 'Definition',
                               'AU': 'Author',
                               'SE': 'Source of seed',
                               'SS': 'Source of structure',
                               'BM': 'Build method',
                               'SM': 'Search method',
                               'GA': 'Gathering threshold',
                               'TC': 'Trusted Cutoff',
                               'NC': 'Noise Cutoff',
                               'TP': 'Type',
                               'SQ': 'Sequence',    # number of sequences in the alignment
                               'DC': 'Database Comment',
                               'DR': 'Database Reference',
                               'RC': 'Reference Comment',
                               'RN': 'Reference Number',
                               'RM': 'Reference Medline',
                               'RT': 'Reference Title',
                               'RA': 'Reference Author',
                               'RL': 'Reference Location',
                               'PI': 'Previous identifier',
                               'KW': 'Keywords',
                               'CC': 'Comment',
                               'NE': 'Pfam accession',
                               'NL': 'Location',
                               'WK': 'Wikipedia link',
                               'CL': 'Clan',
                               'MB': 'Membership',
                               'NH':  'New Hampshire',
                               'TN':  'Tree ID',
                               'FR': 'False discovery Rate',
                               'CB': 'Calibration method'
                               }

        # GS
        seq_feature = {'AC': 'ACcession number',
                       'DE': 'DEscription',
                       'DR': 'Database Reference',
                       'OS': 'Organism',
                       'OC': 'Organism Classification',
                       'LO': 'Look'}

        # GR
        residue_annotation = {'SS': 'Secondary Structure',
                              'SA': 'Surface Accessibility',
                              'TM': 'TransMembrane',
                              'PP': 'Posterior Probability',
                              'LI': 'LIgand binding',
                              'AS': 'Active Site',
                              'pAS': 'AS - Pfam predicted',
                              'sAS': 'AS - from SwissProt',
                              'IN': 'INtron (in or after)',
                              'tWW': 'WC/WC',
                              'cWH': 'WC/Hoogsteen cis',
                              'cWS': 'WC/SugarEdge cis',
                              'tWS': 'WC/SugarEdge trans',
                              }
        # GC
        per_column_ann = {'RF': 'ReFerence annotations',
                          'MM': 'Model Mask',
                          'SS_cons': 'Consensus structure',
                          }

        self.GC = per_column_ann
        self.GR = residue_annotation
        self.GS = seq_feature
        self.GF = feature_annotations

    def add_custom_parser_tags(self, type, tags_dict):
        """
        adds custom fields to stockholm parser
        there are 4 types of tags namely GC, GR, GS, and GF
        you must specify the type for dict addition as they have different meaning and are parser differently
        see wikipedia https://en.wikipedia.org/wiki/Stockholm_format for reference

        the line must follow the general format of the tagged line
        # #=GF <feature> <Generic per-File annotation, free text>
        # #=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
        # #=GS <seqname> <feature> <Generic per-Sequence annotation, free text>
        # #=GR <seqname> <feature> <Generic per-Residue annotation, exactly 1 char per residue>

        :param type: str: 'GC' or 'GR' or 'GS' or 'GF'
        :param tags_dict: dict of the custom tags see example
        :return:

        example:
        my_dict = {'AB': 'bbaabbabbabbab'}
        feature_obj = StockholmFeatureStock()
        feature_obj.add_custom_parser_tags('GC', my_dict)

        will add ability to recognize 'AB' tag in #=GC line
        the line must follow the general format of the tagged line

        also this wile fail spectacularly if input requirements are not met
        """
        self.__getattribute__(type).update(tags_dict)


class StockholmAlig(object):
    """
    this object is lot like a list of SeqRecord objects
    """
    def __init__(self):
        super(StockholmAlig, self).__init__()
        self.annotations = {}
        self.column_annotations = {}
        self._records = []

    def append(self, p_object):
        if not isinstance(p_object, SeqRecord):
            raise TypeError('Expected type SeqRecord, but given {}'.format(type(p_object)))
        # super(StockholmAlig, self).append(p_object)
        self._records.append(p_object)

    def __iter__(self):
        """
        iterate over records
        :return:
        """
        return iter(self._records)

    def __len__(self):
        """
        return number of sequences in alignment
        :return:
        """
        return len(self._records)

    def __str__(self):
        return '\n'.join([str(i.seq) for i in self._records])

    def __repr__(self):
        return str(self)

    def get_unalined_seqs(self, keep_letter_ann=False, ref_seq=0):
        """
        return iterator
        if keep_letter_ann is True, then letter annotations will be included
        and gap characters will be removed from it
        :return:
        """
        for aligned_seq in self.__iter__():

            seq = str(aligned_seq.seq)
            allins = re.finditer('[^-.]+', seq)
            seq_without_gaps = ''

            for match in allins:
                m_start = match.span()[0]
                m_end = match.span()[1]
                seq_without_gaps += seq[m_start:m_end]

            if keep_letter_ann:
                seq_letter_annotatitons = dict.fromkeys(aligned_seq.letter_annotations.keys(), '')
                allins = re.finditer('[^-.]+', seq)
                for match in allins:
                    m_start = match.span()[0]
                    m_end = match.span()[1]
                    for key in seq_letter_annotatitons.keys():
                        seq_letter_annotatitons[key] += aligned_seq.letter_annotations[key][m_start:m_end]

                gapless = SeqRecord(Seq(seq_without_gaps, aligned_seq.seq.alphabet),
                                    id=aligned_seq.id,
                                    letter_annotations=seq_letter_annotatitons,
                                    annotations=aligned_seq.annotations)
            else:
                gapless = SeqRecord(Seq(seq_without_gaps, aligned_seq.seq.alphabet),
                                    id=aligned_seq.id)

            yield gapless

    def slice_alig(self, slice_start, slice_end):

        ret_alig = StockholmAlig()
        for record in self.__iter__():
            ret_alig.append(record[slice_start:slice_end])

        return ret_alig

    def __getitem__(self, index):
        """
        Access part of th alignment in style of Bio.Align class
        :param index:
        :return:
        """

        if isinstance(index, int):
            return self._records[index]
        elif isinstance(index, slice):
            sub_align = StockholmAlig()
            sub_align._records = self._records[index]
            sub_align.column_annotations = self.column_annotations
            sub_align.annotations = self.annotations
            return sub_align
        elif len(index) != 2:
            raise TypeError('Invalid index type')

        row_index, col_index = index
        if isinstance(row_index, int):
            return self._records[row_index][col_index]
        elif isinstance(col_index, int):
            return ''.join(rec[col_index] for rec in self._records[row_index])
        else:
            sub_align = StockholmAlig()
            sub_align._records = [rec[col_index] for rec in self._records[row_index]]
            # for rec in self._records[row_index]:
            #     nr = rec[col_index]
            #     for la in rec.letter_annotations.keys():
            #         nr.letter_annotations[la] = rec.letter_annotations[la][col_index]
            sub_col_ann = dict()
            for key in self.column_annotations.keys():
                sub_col_ann[key] = self.column_annotations[key][col_index]
            sub_align.column_annotations = sub_col_ann
            sub_align.annotations = self.annotations
            return sub_align

    def __add__(self, other):
        """
        bluntly define add on st_alig for alig concatation
        preserve id from left alig only

        consider this to be part of new class, because of id throaway

        :param other:
        :return:
        """
        if not isinstance(other, StockholmAlig):
            raise TypeError('expected other to be type {} but type {} given'.format(type(self), type(other)))

        ns = StockholmAlig()
        if len(self.column_annotations) != 0:
            if self.column_annotations == other.column_annotations:
                for key in self.column_annotations.keys():
                    ns.column_annotations[key] = self.column_annotations[key] + other.column_annotations[key]
            else:
                print('cannot join column annotations as they have different names')

        for l, r in zip(self, other):
            n = SeqRecord(Seq(str(l.seq) + str(r.seq)),
                          id=l.id)
            if len(l.letter_annotations) != 0:
                if l.letter_annotations.keys() == r.letter_annotations.keys():
                    for key in zip(l.letter_annotations.keys()):
                        n.letter_annotations[key] = l.letter_annotations[key] + r.letter_annotations[key]
                else:
                    print('cannot join letter annotations as they have different names')

            ns.append(n)

        return ns

    def ladd_gap(self, gapchar='-'):
        na = StockholmAlig()
        if len(self.column_annotations) != 0:
            for key in self.column_annotations.keys():
                na.column_annotations[key] = gapchar + self.column_annotations[key]

        for seq in self:
            ns = SeqRecord(Seq(gapchar + str(seq.seq)),
                           id=seq.id)
            if len(seq.letter_annotations) != 0:
                for key in seq.letter_annotations.keys():
                    ns.letter_annotations[key] = gapchar + seq.letter_annotations[key]
            na.append(ns)
        return na

    def vertcat(self, other):
        """
        this will vertically concat the alignment
        CAUTION this operation will drop any
        :param other:
        :return:
        """
        if not isinstance(other, StockholmAlig):
            raise TypeError('expected type {} but type {} given'.format(type(self), type(other)))
        assert len(self._records[0]) == len(other._records[0])

        ns = StockholmAlig()
        ns._records = self._records + other._records
        return ns

    def write_clustal(self, f):
        """
        write Clustal alignment to a file
        code from Bio.AlingIO.ClustalIO is used
        :param f:
        :return:
        """
        # transform possible unclustal gap characters to dashes
        clustalized_alig = _clustalize_alignment(self)
        writer = ClustalWriter(f)
        writer.write_alignment(clustalized_alig)

    def write_stockholm(self, f):
        """
        write stockholm alignment file
        :param f:
        :return:
        """

        f.write('# STOCKHOLM 1.0\n')
        if hasattr(self, 'version'):
            f.write('#=GF FMT written by {}\n'.format(self.version))
        if hasattr(self, 'name'):
            f.write('#=GF ID {}'.format(self.name))
        f.write(self._GF())
        f.write(self._GS())
        f.write(self._seq_and_GR())
        f.write(self._GC())
        f.write('//\n')

    def write_fasta(self, f):
        for seq in self._records:
            f.write('>{}\n{}\n'.format(seq.id, str(seq.seq)))

    def _GF(self):
        # annotations
        gf = ''
        for key in self.annotations.keys():
            gf += '#=GF ' + key + ' ' + ' ' + self.annotations[key] + '\n'

        return gf

    def _GS(self):
        # sequence_annotations
        gs = ''
        for seq in self.__iter__():
            for key in seq.annotations.keys():
                gs += '#=GS ' + seq.id + ' ' + key + ' ' + seq.annotations[key] + '\n'

        return gs

    def _GC(self):
        # column_annotations
        gc = ''
        for key in self.column_annotations.keys():
            gc += '#=GC ' + key + ' ' + self.column_annotations[key] + '\n'

        return gc

    def _seq_and_GR(self):
        # sequence and per sequence residue annotations
        sgr = ''
        for seq in self.__iter__():
            sgr += seq.id + ' ' + str(seq.seq) + '\n'
            for key in seq.letter_annotations.keys():
                sgr += '#=GR ' + seq.id + ' ' + key + ' ' + seq.letter_annotations[key] + '\n'

        return sgr

    def get_alignment_length(self):
        """
        get maximum length of aligned sequnces a specified in description of this function
         - all aligned sequences should be same length
        :return: maximum length of aligned sequences
        """
        return max([len(align.seq) for align in self.__iter__()])


def _clustalize_alignment(alig):
    new_alig = StockholmAlig()
    for al in alig:
        ns = SeqRecord(Seq(re.sub('[.~]', '-', str(al.seq))),
                       id=al.id)
        new_alig.append(ns)
    return new_alig