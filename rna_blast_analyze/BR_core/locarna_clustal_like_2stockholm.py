from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rna_blast_analyze.BR_core.exceptions import ParsingError
from rna_blast_analyze.BR_core.stockholm_alig import StockholmAlig


def parse_locarna_alignment(fid_file):

    txt = fid_file.readline()
    if 'Score' not in txt:
        raise ParsingError('Unrecognized format - expected "Score as a first argument"')
    st = StockholmAlig()
    st.annotations['score'] = int(txt.split()[-1])

    namespace = {}
    txt = fid_file.readline()

    c = 1
    while txt:
        if not txt.strip():
            txt = fid_file.readline()
            c += 1
            continue

        l = txt.split()
        if len(l) != 2:
            raise ParsingError('Unrecognized format (unexpected white-space encountered on line {})'.format(c))
        if l[0] in namespace.keys():
            namespace[l[0]] += l[1]
        else:
            namespace[l[0]] = l[1]

        txt = fid_file.readline()
        c += 1

    known_attr = ['#S',
                  '#FS',
                  '#A1',
                  '#A2',
                  'alifold',
                  ]
    for seqn in namespace.keys():
        if seqn not in known_attr:
            new_seq = SeqRecord(Seq(namespace[seqn]),
                                id=seqn)
            st.append(new_seq)
        else:
            if seqn == 'alifold':
                st.column_annotations['SS_cons'] = namespace[seqn]
            else:
                st.column_annotations[seqn] = namespace[seqn]

    return st