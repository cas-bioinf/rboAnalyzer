import binascii
import pickle
import time
import copy
from argparse import Namespace

from Bio.Blast import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

import rna_blast_analyze.BR_core.BA_methods
from rna_blast_analyze.BR_core.BA_support import Subsequences


def seq2dict(seq):
    if not isinstance(seq, Seq):
        raise TypeError('accepts only Bio.Seq')
    out = {
        'seq': str(seq),
        'alphabet': binascii.hexlify(pickle.dumps(seq.alphabet)).decode()
    }
    return out


def seqfromdict(indict):
    return Seq(
        indict['seq'],
        alphabet=pickle.loads(
            binascii.unhexlify(
                indict['alphabet']
            )
        )
    )


def seqrecord2dict(seqrec):
    """
    using pickle for features (usually not used)
     Bio.SeqFeature.SeqFeature class with Bio.SeqFeature.FeatureLocation should be used
     - implement direct conversion if needed
    :param seqrec:
    :return:
    """
    if seqrec is None:
        return None
    if not isinstance(seqrec, SeqRecord):
        raise TypeError('accepts only Bio.SeqRecord but "{}" provided'.format(type(seqrec)))
    out = {
        'Bio.Seq': seq2dict(seqrec.seq),
        'annotations': annotations_items2dict(seqrec.annotations),
        'dbxrefs': seqrec.dbxrefs,
        'description': seqrec.description,
        'features': [binascii.hexlify(pickle.dumps(ff)).decode() for ff in seqrec.features],
        'id': seqrec.id,
        'name': seqrec.name,
        'letter_annotations': seqrec.letter_annotations,
    }
    return out


def seqrecordfromdict(indict):
    if indict is None:
        return None
    return SeqRecord(
        seqfromdict(indict['Bio.Seq']),
        annotations=annotations_items_from_dict(indict['annotations']),
        dbxrefs=indict['dbxrefs'],
        description=indict['description'],
        features=[pickle.loads(binascii.unhexlify(ff)) for ff in indict['features']],
        id=indict['id'],
        name=indict['name'],
        letter_annotations=indict['letter_annotations'],
    )


class NativeDict(dict):
    """
        Helper class to ensure that only native types are in the dicts produced by
        :func:`to_dict() <pandas.DataFrame.to_dict>`

        .. note::

            Needed until `#21256 <https://github.com/pandas-dev/pandas/issues/21256>`_ is resolved.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(((k, self.convert_if_needed(v)) for row in args for k, v in row), **kwargs)

    @staticmethod
    def convert_if_needed(value):
        """
            Converts `value` to native python type.

            .. warning::

                Only :class:`Timestamp <pandas.Timestamp>` and numpy :class:`dtypes <numpy.dtype>` are converted.
        """
        if pd.isnull(value):
            return None
        if isinstance(value, pd.Timestamp):
            return value.to_pydatetime()
        if hasattr(value, 'dtype'):
            mapper = {'i': int, 'u': int, 'f': float}
            _type = mapper.get(value.dtype.kind, lambda x: x)
            return _type(value)
        return value


def annotations_items2dict(item):
    out = item.copy()
    if 'blast' in item:
        out['blast'] = (item['blast'][0], hsptodict(item['blast'][1]))
    if 'cmstat' in item:
        out['cmstat'] = item['cmstat'].to_dict(into=NativeDict)

    return out


def annotations_items_from_dict(indict):
    if 'blast' in indict:
        indict['blast'] = (indict['blast'][0], hspfromdict(indict['blast'][1]))
    if 'cmstat' in indict:
        if isinstance(indict['cmstat'], str):
            indict['cmstat'] = pd.read_json(indict['cmstat'], typ='series')
    return indict


def subsequences2dict(subsequences):
    if not isinstance(subsequences, Subsequences):
        raise TypeError('accepts only BA_support.Subsequences')
    out = {
        'extension': seqrecord2dict(subsequences.extension),
        'source': seqrecord2dict(subsequences.source),
        'query_name': subsequences.query_name,
        'best_start': subsequences.best_start,
        'best_end': subsequences.best_end,
        'templates': subsequences.templates,
    }
    return out


def subsequencesfromdict(indict):
    out = Subsequences(
        seqrecordfromdict(indict['source'])
    )
    if 'extension' in indict:
        out.extension = seqrecordfromdict(indict['extension'])
    elif 'subs' in indict and 'ret_keys' in indict:
        # this is for backward compatibility
        out.extension = seqrecordfromdict(indict['subs'][indict['ret_keys'][0]])
    else:
        raise KeyError
    out.query_name = indict['query_name']
    out.best_start = indict['best_start']
    out.best_end = indict['best_end']
    out.templates = indict['templates']
    return out


def hitlist2dict(hl):
    if not isinstance(hl, rna_blast_analyze.BR_core.BA_methods.HitList):
        raise TypeError('Accepts only BA_methods.HitList. {} given.'.format(type(hl)))
    return [subsequences2dict(i) for i in hl]
    # out = {}
    # for i in range(len(hl)):
    #     out[str(i)] = subsequences2dict(hl[i])
    # return out


def hitlistfromdict(indict):
    out = rna_blast_analyze.BR_core.BA_methods.HitList()
    if isinstance(indict, list):
        for i in indict:
            out.append(subsequencesfromdict(i))
        return out
    elif isinstance(indict, dict):
        # for backward compatibility
        for i in range(len(indict)):
            out.append(subsequencesfromdict(indict[str(i)]))
        return out
    else:
        raise TypeError


def blastsearchrecompute2dict(bsr):
    if not isinstance(bsr, rna_blast_analyze.BR_core.BA_methods.BlastSearchRecompute):
        raise TypeError('Accepts only BA_methods.BlastSearchRecompute. {} given'.format(type(bsr)))
    out = {
        'hits': hitlist2dict(bsr.hits),
        'query': seqrecord2dict(bsr.query),
        'creation': bsr.creation,
        '_runstat': bsr._runstat,
        'args': vars(bsr.args),
        'date_of_run': list(bsr.date_of_run)
    }
    return out


def blastsearchrecomputefromdict(indict):
    out = rna_blast_analyze.BR_core.BA_methods.BlastSearchRecompute()
    out.hits = hitlistfromdict(indict['hits'])
    out.query = seqrecordfromdict(indict['query'])
    out.creation = indict['creation']
    out._runstat = indict['_runstat']
    out.args = Namespace(**indict['args'])
    out.date_of_run = time.struct_time(indict['date_of_run'])
    return out


def headertodict(header):
    if not isinstance(header, Record.Header):
        raise TypeError('accepts only BlastRecord.Header')
    return vars(header)


def headerfromdict(indict):
    return universalfromdict(indict, Record.Header)


def hsptodict(hsp):
    if not isinstance(hsp, Record.HSP):
        raise TypeError('accepts only BlastRecord.HSP')
    return vars(hsp)


def hspfromdict(indict):
    return universalfromdict(indict, Record.HSP)


def blastalignmenttodict(blastalig):
    if not isinstance(blastalig, Record.Alignment):
        raise TypeError('accepts only Bio.Blast.Record.Alignment')
    out = vars(copy.copy(blastalig))
    out['hsps'] = [hsptodict(h) for h in blastalig.hsps]
    return out


def blastalignmentfromdict(indict):
    out = Record.Alignment()
    keys = set(indict.keys())
    keys.remove('hsps')
    for key in keys:
        setattr(out, key, indict[key])
    out.hsps = [hspfromdict(h) for h in indict['hsps']]
    return out


def blasttodict(blast):
    if not isinstance(blast, Record.Blast):
        raise TypeError('accepts only BlastRecord.Blast')

    out = dict()
    for key in vars(blast).keys():
        if key == 'descriptions':
            out[key] = [vars(d) for d in blast.descriptions]
        elif key == 'alignments':
            out[key] = [blastalignmenttodict(balig) for balig in blast.alignments]
        elif key == 'multiple_alignments':
            out[key] = []
        else:
            out[key] = getattr(blast, key)
    return out


def blastfromdict(indict):
    out = Record.Blast()
    for key in indict.keys():
        if key == 'descriptions':
            out.descriptions = [universalfromdict(d, Record.Description) for d in indict['descriptions']]
        elif key == 'alignments':
            out.alignments = [blastalignmentfromdict(balig) for balig in indict['alignments']]
        else:
            setattr(out, key, indict[key])
    return out


def universalfromdict(inputdict, outclasspointer):
    outputobject = outclasspointer()
    for key in inputdict.keys():
        setattr(outputobject, key, inputdict[key])

    return outputobject
