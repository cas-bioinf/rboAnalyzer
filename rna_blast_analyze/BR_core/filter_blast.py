import operator
import logging
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger('rboAnalyzer')

OPERATIONS = {
    '>': operator.gt,
    '<': operator.lt,
    '=': operator.eq,
    '>=': operator.ge,
    '<=': operator.le,
}


def filter_by_eval(blast_hitlist, getter, filter_conditions):
    ml.debug(fname())

    result = blast_hitlist
    for relation, condition in filter_conditions:
        result = [h for h in result if OPERATIONS[relation](getter(h).expect, condition)]
    return result


def filter_by_bits(blast_hitlist, getter, filter_conditions):
    ml.debug(fname())

    result = blast_hitlist
    for relation, condition in filter_conditions:
        result = [h for h in result if OPERATIONS[relation](getter(h).bits, condition)]
    return result


