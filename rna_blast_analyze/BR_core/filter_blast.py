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


def filter_by_eval(blast_hitlist, relation, eval):
    ml.debug(fname())
    filtered = []
    for h in blast_hitlist:
        if OPERATIONS[relation](h[1].expect, eval):
            filtered.append(h)
    return filtered


def filter_by_bits(blast_hitlist, relation, bits):
    ml.debug(fname())
    filtered = []
    for h in blast_hitlist:
        if OPERATIONS[relation](h[1].bits, bits):
            filtered.append(h)
    return filtered


