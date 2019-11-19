import os
import re
import sys
import logging
import rna_blast_analyze.BR_core.tools_versions as tools
from rna_blast_analyze.BR_core.filter_blast import OPERATIONS
from Bio import SeqIO
from rna_blast_analyze.BR_core.BA_support import IUPACmapping
from hashlib import sha1

ml = logging.getLogger('rboAnalyzer')

SPECIAL_CHARS = "`!@#$^&*(){}|[]\\;\'\",<>?"


def check_file(f):
    if isinstance(f, str):
        if not shell_safe(f):
            ml.error("'{}' does not appear to be shell safe.".format(f))
            return False
        return os.path.isfile(f)
    else:
        ml.error("Expected 'str' got '{}' instead.".format(type(f)))
        return False


def check_int(x):
    if not isinstance(x, int):
        ml.error("Expected 'int' got '{}' instead.".format(type(x)))
        return False
    else:
        return True


def shell_safe(x:str):
    for c in SPECIAL_CHARS:
        if c in x:
            ml.error("String '{}' can't contain '{}'".format(x, c))
            return False
    return True


def validate_args(args):
    if not check_file(args.blast_in):
        ml.error("Provided BLAST output file NOT found. ({})".format(args.blast_in))
        return False

    if not check_file(args.blast_query):
        ml.error("Provided query file NOT found. ({})".format(args.blast_query))
        return False

    if args.cm_file is not None and not check_file(args.cm_file):
        ml.error("Provided CM file NOT found. ({})".format(args.cm_file))
        return False

    if args.config_file is not None and not check_file(args.config_file):
        ml.error("Provided config file NOT found. ({})".format(args.config_file))
        return False

    if args.b_type not in ['guess', 'xml', 'plain']:
        ml.error("Wrong provided BLAST type.")
        return False

    # check if file appears to be correct fasta file
    check_fasta(args.blast_query)

    dbtype = ["blastdb", "fasta", "gb", "server", "entrez"]
    if args.db_type not in dbtype:
        ml.error("'db_type' can only be {}".format(dbtype))
        return False

    if args.db_type == "blastdb":
        if not isinstance(args.blast_db, str):
            ml.error('Blastdb must be PATH (string).')
            return False

        ext_try = ['nsd', 'nhr', 'nog', 'nsi', 'nin']
        if not os.path.isfile(args.blast_db + '.nal'):
            for ext in ext_try:
                if not os.path.isfile(args.blast_db + '.' + ext):
                    ml.error('Expected file : "' + args.blast_db + '.' + ext + '" NOT found.')
                    return False
    elif args.db_type == 'entrez':
        if not re.fullmatch(r'[^@]+@[^@ ]+\.[^@ ]+', args.entrez):
            ml.error('The provided email appears to be invalid.')
            return False
    else:
        if not os.path.isdir(args.blast_db):
            return False

    msg = 'Provided BLAST regexp must be valid regular expression (https://docs.python.org/3.6/library/re.html).'
    if not isinstance(args.blast_regexp, str):
        ml.error(msg)
        return False
    else:
        try:
            re.compile(args.blast_regexp)
        except re.error:
            ml.error(msg)
            return False

    for arg in ['zip_json', 'skip_missing', 'dev_pred', 'use_rfam', 'enable_overwrite']:
        if not isinstance(getattr(args, arg), bool):
            ml.error("parameter '{}' should be 'bool' not {}".format(arg, type(args.get(arg))))
            return False

    if not any([args.json, args.html, args.csv, args.pandas_dump, args.dump]):
        ml.error(
            "It appears that no output file was requested."
            " Please provide --html and/or --json and/or --csv argument(s)."
        )
        return False

    if not check_int(args.subseq_window_locarna):
        return False
    if args.subseq_window_locarna <= 0:
        ml.error("Argument 'subseq_window_locarna' can't be negative or zero.")
        return False

    if not check_int(args.locarna_anchor_length):
        return False
    if args.locarna_anchor_length <= 1:
        ml.error("Argument 'locarna_anchor_length' must be > 1.")
        return False

    if args.threads is not None:
        if not check_int(args.threads):
            return False
        if args.threads < 1:
            ml.error("If argument 'threads' is given, it must be > 1 (was {}).".format(args.threads))
            return False

    if not shell_safe(args.locarna_params):
        ml.error('Provided arguments for locarna appears not to be shell safe.')
        return False

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

    if any(set(args.prediction_method) - tools.prediction_methods):
        ml.error("Required prediction method(s) not avalible: {}.".format(
            set(args.prediction_method) - tools.prediction_methods)
        )

    # check the filter_by
    if not _check_filter(args.filter_by_eval, 'filter_by_eval'):
        return False

    if not _check_filter(args.filter_by_bitscore, 'filter_by_bitscore'):
        return False

    if args.logfile is not None and not shell_safe(args.logfile):
        return False

    if args.repredict_file is not None and not shell_safe(args.repredict_file):
        return False

    return True


def _check_filter(fv, name):
    if fv is None:
        pass
    elif not isinstance(fv, (list, tuple)):
        ml.error("'{}' should be list or tuple length 2.".format(name))
        return False
    elif any(len(i) != 2 for i in fv):
        ml.error("'{}' should consist of list or tuple length 2.".format(name))
        return False
    else:
        for i in fv:
            d, f = i
            if d not in OPERATIONS:
                ml.error("Relation descriptor {} not allowed".format(d))
                return False

            if not isinstance(f, float):
                ml.error("Filtering Value should be number (float).")
                return False
    return True


def check_fasta(file):
    status = False
    try:
        iupac = IUPACmapping()
        am = set(iupac.ambiguous) | set(i.lower() for i in iupac.ambiguous)
        un = set(iupac.unambiguous) | set(i.lower() for i in iupac.unambiguous)
        with open(file, 'r') as f:
            c = 0
            for s in SeqIO.parse(f, format='fasta'):
                c += 1
                # check if fasta is rna/dna
                seq_seq = set(str(s.seq))
                if len(seq_seq - un) == 0:
                    # sequence is OK
                    pass
                elif len(seq_seq - am) == 0:
                    # sequence contains ambiguous bases
                    pass
                else:
                    # sequence contains not allowed characters
                    ml.error('Sequence {} contains disallowed characters. Please check the sequence.'.format(s.id))
                    status = True

                # check fasta id
                # we allow A-Za-z0-9_|[]- apart form control ">"
                if not re.fullmatch(r'[A-Za-z0-9_\|\[\]\-\.]*', s.id):
                    ml.error(
                        'Sequence id contains disallowed characters. '
                        'Please check the sequence id. '
                        'We allow following characters: "A-Za-z0-9_|[]-".')
                    status = True

            if c == 0:
                raise Exception(
                    'The provided FASTA file appears to be empty. '
                    '(No sequence in FASTA format (starting with ">") detected.)'
                )

        if status:
            raise Exception('One or more sequences contains characters not allowed in DNA/RNA FASTA file.')
    except Exception as e:
        ml.error('There was an error with provided FASTA file "{}":\n{}'.format(file, e))
        sys.exit(1)


def verify_query_blast(blast, query):
    """
    verify if query from fasta file matches query sequence described in BLAST.
    :param blast:
    :param query:
    :return:
    """
    if not (query.id == blast.query or query.description == blast.query):
        ml.warning(
            'Provided query id ({}) do not match query id in BLAST output ({})'.format(
                query.id,
                blast.query
            )
        )
    if len(query) != blast.query_length:
        ml.error(
            'Provided query lenght ({}: {}) do not match BLAST query length ({}: {}).\n'
            'Please provide correct query sequence.'.format(
                query.id,
                len(query),
                blast.query,
                blast.query_length
            )
        )
        sys.exit(1)
    return


def check_blast(p_blast):
    # since we dont want to check the BLAST program
    #  we just going to check if all accessions are found and if all sbjct sequences in BLAST appears to be nucleotide
    iupac = IUPACmapping()
    allowed = set(iupac.allowed_bases) | set(i.lower() for i in iupac.allowed_bases) | set('-')
    for search in p_blast:
        for aln in search.alignments:
            for hsp in aln.hsps:
                c = set(hsp.sbjct) - allowed
                if c:
                    ml.error(
                        "Unexpected character encountered in BLAST output.\n"
                        "The character(s) was '{}' in {}.\n"
                        "Please check that nucleotide database was used.".format(
                        c, aln.hit_id
                    ))
                    sys.exit(1)


def compute_args_hash(args, rfam_file=None):
    change_prone = [
        'subseq_window_locarna',
        'locarna_params',
        'locarna_anchor_length',
        'mode',
        'skip_missing',
        'blast_regexp',
        'use_rfam',
    ]

    files = [
        getattr(args, argname) for argname in ['blast_in', 'blast_query', 'cm_file', 'config_file', ]
    ] + [rfam_file]

    hh = sha1()
    for argname in change_prone:
        hh.update(str(getattr(args, argname)).encode())

    for fn in files:
        try:
            if fn is None:
                continue
            with open(fn, 'rb') as f:
                while True:
                    cdata = f.read(4096)
                    if not cdata:
                        break
                    hh.update(cdata)
        except Exception:
            print(fn)
            ml.error('Failed to generate hash from given parameters.')
            sys.exit(1)

    nh = hh.hexdigest()
    return nh