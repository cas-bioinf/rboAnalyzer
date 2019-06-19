import os
import re
import logging
import rna_blast_analyze.BR_core.tools_versions as tools
from rna_blast_analyze.BR_core.filter_blast import OPERATIONS

ml = logging.getLogger(__name__)

SPECIAL_CHARS = "`!@#$^&*(){}|[]\\;\'\",<>?"


def check_file(f):
    if isinstance(f, str):
        if not shell_safe(f):
            raise ValueError
        return os.path.isfile(f)
    else:
        raise TypeError("Expected 'str' got '{}' instead.".format(type(f)))


def check_int(x):
    if not isinstance(x, int):
        raise TypeError("Expected 'int' got '{}' instead.".format(type(x)))


def shell_safe(x:str):
    for c in SPECIAL_CHARS:
        if c in x:
            ml.error("String '{}' can't contain '{}'".format(x, c))
            return False
    return True


def validate_args(args):
    if not check_file(args.blast_in):
        ml.error("Provided BLAST output file NOT found.")
        return False

    if not check_file(args.blast_query):
        ml.error("Provided query file NOT found.")
        return False

    if args.cm_file is not None and not check_file(args.cm_file):
        ml.error("Provided CM file NOT found.")
        return False

    if args.config_file is not None and not check_file(args.config_file):
        ml.error("Provided config file NOT found.")
        return False

    if args.b_type not in ['guess', 'xml', 'plain']:
        ml.error("Wrong provided BLAST type.")
        return False

    dbtype = ["blastdb", "fasta", "gb", "server"]
    if args.db_type not in dbtype:
        ml.error("'db_type' can only be {}".format(dbtype))
        return False

    if args.db_type == "blastdb":
        if not isinstance(args.blast_db, str):
            raise TypeError

        ext_try = ['nsd', 'nhr', 'nog', 'nsi', 'nin']
        if not os.path.isfile(args.blast_db + '.nal'):
            for ext in ext_try:
                if not os.path.isfile(args.blast_db + '.' + ext):
                    ml.error('Expected file : "' + args.blast_db + '.' + ext + '" NOT found.')
                    return False
    else:
        if not os.path.isdir(args.blast_db):
            return False

    if not isinstance(args.blast_regexp, str):
        raise TypeError
    else:
        try:
            re.compile(args.blast_regexp)
        except re.error:
            ml.error("Provided BLAST regexp is invalid.")
            return False

    if not isinstance(args.keep_all, bool):
        raise TypeError

    if not isinstance(args.zip_json, bool):
        raise TypeError

    if not isinstance(args.show_gene_browser, bool):
        raise TypeError

    if not isinstance(args.skip_missing, bool):
        raise TypeError

    if not isinstance(args.dev_pred, bool):
        raise TypeError

    if not isinstance(args.use_rfam, bool):
        raise TypeError

    if not isinstance(args.download_rfam, bool):
        raise TypeError

    if not isinstance(args.enable_overwrite, bool):
        raise TypeError

    check_int(args.subseq_window_locarna)
    if args.subseq_window_locarna <= 0:
        ml.error("Argument 'subseq_window_locarna' can't be negative.")
        return False

    check_int(args.subseq_window_simple_ext)
    if args.subseq_window_simple_ext <= 0:
        ml.error("Argument 'subseq_window_simple_ext' can't be negative.")
        return False

    check_int(args.locarna_anchor_length)
    if args.locarna_anchor_length <= 1:
        ml.error("Argument 'locarna_anchor_length' must be > 1.")
        return False

    if args.threads is not None:
        check_int(args.threads)
        if args.threads < 1:
            ml.error("If argument 'threads' given, it must be > 1 (was {}).".format(args.threads))
            return False

    if not shell_safe(args.locarna_params):
        raise ValueError

    def precheck_file_exist(f):
        if f is not None and check_file(f):
            return True

    if precheck_file_exist(args.dump) and not args.enable_overwrite:
        ml.error("Refusing to overwrite 'dump' {}.".format(args.dump))
        return False

    if precheck_file_exist(args.pandas_dump) and not args.enable_overwrite:
        ml.error("Refusing to overwrite 'pandas_dump' {}.".format(args.pandas_dump))
        return False

    if precheck_file_exist(args.csv) and not args.enable_overwrite:
        ml.error("Refusing to overwrite 'csv' {}.".format(args.csv))
        return False

    if precheck_file_exist(args.json) and not args.enable_overwrite:
        ml.error("Refusing to overwrite 'json' {}.".format(args.json))
        return False

    if precheck_file_exist(args.html) and not args.enable_overwrite:
        ml.error("Refusing to overwrite 'html' {}.".format(args.html))
        return False

    assert isinstance(args.logmsgs, list)

    if any(set(args.prediction_method) - tools.prediction_methods):
        ml.error("Required prediction method(s) not avalible: {}.".format(
            set(args.prediction_method) - tools.prediction_methods)
        )

    # handle the filter_by
    if args.filter_by_eval is None:
        pass
    elif not isinstance(args.filter_by_eval, (list, tuple)):
        raise TypeError
    elif len(args.filter_by_eval) != 2:
        raise ValueError
    else:
        d, f = args.filter_by_eval
        if d not in OPERATIONS:
            ml.error("Relation descriptor {} not allowed".format(d))
            return False

        if not isinstance(f, float):
            raise TypeError

    if args.filter_by_bitscore is None:
        pass
    elif not isinstance(args.filter_by_bitscore, (list, tuple)):
        raise TypeError
    elif len(args.filter_by_bitscore) != 2:
        raise ValueError
    else:
        d, f = args.filter_by_bitscore
        if d not in OPERATIONS:
            ml.error("Relation descriptor {} not allowed".format(d))
            return False

        if not isinstance(f, float):
            raise TypeError

    if args.logfile is not None and not shell_safe(args.logfile):
        return False

    if args.repredict_file is not None and not shell_safe(args.repredict_file):
        return False

    return True