# NOTE this is derived code from Biopython
#   We add only a few modifications needed for custom parser

# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Record classes to hold BLAST output.

Classes:
Blast              Holds all the information from a blast search.
PSIBlast           Holds all the information from a psi-blast search.

Header             Holds information from the header.
Description        Holds information about one hit description.
Alignment          Holds information about one alignment hit.
HSP                Holds information about one HSP.
MultipleAlignment  Holds information about a multiple alignment.
DatabaseReport     Holds information from the database report.
Parameters         Holds information from the parameters.

"""
# XXX finish printable BLAST output

# These are imported in original file but not used
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment

# These are declared in original file, however, they are not modified
#  if they are needed without modification, import directly from Bio.Blast.Record
# from Bio.Blast.Record import Description
# from Bio.Blast.Record import Alignment
# from Bio.Blast.Record import MultipleAlignment
# from Bio.Blast.Record import Round
# from Bio.Blast.Record import PSIBlast

from Bio.Blast.Record import DatabaseReport
from Bio.Blast.Record import Parameters


class Header(object):
    """Saves information from a blast header.

    Members:
    application         The name of the BLAST flavor that generated this data.
    version             Version of blast used.
    date                Date this data was generated.
    reference           Reference for blast.

    query               Name of query sequence.
    query_letters       Number of letters in the query sequence.  (int)

    database            Name of the database.
    database_sequences  Number of sequences in the database.  (int)
    database_letters    Number of letters in the database.  (int)

    rid                 NCBI request id

    """
    def __init__(self):
        self.application = ''
        self.version = ''
        self.date = ''
        self.reference = ''
        self.rid = None
        self.query = ''
        self.query_letters = None

        self.database = ''
        self.database_sequences = None
        self.database_letters = None


class HSP(object):
    """Stores information about one hsp in an alignment hit.

    Members:
        - score           BLAST score of hit.  (float)
        - bits            Number of bits for that score.  (float)
        - expect          Expect value.  (float)
        - num_alignments  Number of alignments for same subject.  (int)
        - identities      Number of identities (int) if using the XML parser.
          Tuple of number of identities/total aligned (int, int)
          if using the (obsolete) plain text parser.
        - positives       Number of positives (int) if using the XML parser.
          Tuple of number of positives/total aligned (int, int)
          if using the (obsolete) plain text parser.
        - gaps            Number of gaps (int) if using the XML parser.
          Tuple of number of gaps/total aligned (int, int) if
          using the (obsolete) plain text parser.
        - align_length    Length of the alignment. (int)
        - strand          Tuple of (query, target) strand.
        - frame           Tuple of 1 or 2 frame shifts, depending on the flavor.

        - query           The query sequence.
        - query_start     The start residue for the query sequence.  (1-based)
        - query_end       The end residue for the query sequence.  (1-based)
        - match           The match sequence.
        - sbjct           The sbjct sequence.
        - sbjct_start     The start residue for the sbjct sequence.  (1-based)
        - sbjct_end       The end residue for the sbjct sequence.  (1-based)

    Not all flavors of BLAST return values for every attribute::

                  score     expect     identities   positives    strand  frame
        BLASTP     X          X            X            X
        BLASTN     X          X            X            X          X
        BLASTX     X          X            X            X                  X
        TBLASTN    X          X            X            X                  X
        TBLASTX    X          X            X            X                 X/X

    Note: for BLASTX, the query sequence is shown as a protein sequence,
    but the numbering is based on the nucleotides.  Thus, the numbering
    is 3x larger than the number of amino acid residues.  A similar effect
    can be seen for the sbjct sequence in TBLASTN, and for both sequences
    in TBLASTX.

    Also, for negative frames, the sequence numbering starts from
    query_start and counts down.

    """
    def __init__(self):
        self.score = None
        self.bits = None
        self.expect = None
        self.num_alignments = None
        self.identities = (None, None)
        self.positives = (None, None)
        self.gaps = (None, None)
        self.align_length = None
        self.strand = (None, None)
        self.frame = ()

        self.query = ''
        self.query_start = None
        self.query_end = None
        self.match = ''
        self.sbjct = ''
        self.sbjct_start = None
        self.sbjct_end = None
        self.features = None

    def __str__(self):
        lines = ["Score %i (%i bits), expectation %0.1e, alignment length %i"
                 % (self.score, self.bits, self.expect, self.align_length)]
        if self.align_length < 50:
            lines.append("Query:%s %s %s" % (str(self.query_start).rjust(8),
                                       str(self.query),
                                       str(self.query_end)))
            lines.append("               %s"
                         % (str(self.match)))
            lines.append("Sbjct:%s %s %s" % (str(self.sbjct_start).rjust(8),
                                       str(self.sbjct),
                                       str(self.sbjct_end)))
        else:
            lines.append("Query:%s %s...%s %s"
                         % (str(self.query_start).rjust(8),
                            str(self.query)[:45],
                            str(self.query)[-3:],
                            str(self.query_end)))
            lines.append("               %s...%s"
                         % (str(self.match)[:45],
                            str(self.match)[-3:]))
            lines.append("Sbjct:%s %s...%s %s"
                         % (str(self.sbjct_start).rjust(8),
                            str(self.sbjct)[:45],
                            str(self.sbjct)[-3:],
                            str(self.sbjct_end)))
        return "\n".join(lines)


# TODO - Add a friendly __str__ method to BLAST results
class Blast(Header, DatabaseReport, Parameters):
    """Saves the results from a blast search.

    Members:
    descriptions        A list of Description objects.
    alignments          A list of Alignment objects.
    multiple_alignment  A MultipleAlignment object.
    + members inherited from base classes

    """
    def __init__(self):
        Header.__init__(self)
        DatabaseReport.__init__(self)
        Parameters.__init__(self)
        self.descriptions = []
        self.alignments = []
        self.multiple_alignment = None

    def update(self, commonparseddata):
        """
        :param commonparseddata: dict
        :return:
        """
        for key in commonparseddata.keys():
            if hasattr(self, key):
                setattr(self, key, commonparseddata[key])
