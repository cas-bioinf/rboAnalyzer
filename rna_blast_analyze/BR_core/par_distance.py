import logging

ml = logging.getLogger('rboAnalyzer')

try:
    import RNA
    from rna_blast_analyze.BR_core.par_distance_efective import *
    ml.info('Importing RNAdistance functions from RNAlib')
except ImportError:
    ml.info('Import of RNAdistance from RNAlib failed. Using commandline version.')
    from rna_blast_analyze.BR_core.par_distance_no_RNAlib import *
    from rna_blast_analyze.BR_core.viennaRNA import RNAfold, wrap_RNAfold
