import operator
import logging
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)

OPERATIONS = {
    '>': operator.gt,
    '<': operator.lt,
    '=': operator.eq,
    '>=': operator.ge,
    '<=': operator.le,
}


def filter_by_eval(blast_hitlist, getter, relation, eval):
    ml.debug(fname())
    filtered = []
    for h in blast_hitlist:
        if OPERATIONS[relation](getter(h).expect, eval):
            filtered.append(h)
    return filtered


def filter_by_bits(blast_hitlist, getter, relation, bits):
    ml.debug(fname())
    filtered = []
    for h in blast_hitlist:
        if OPERATIONS[relation](getter(h).bits, bits):
            filtered.append(h)
    return filtered


