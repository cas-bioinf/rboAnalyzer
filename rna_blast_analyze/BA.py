#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
import argcomplete
import argparse
import sys
import json
import os

# this must precede the CONFIG
from rna_blast_analyze.BR_core.tools_versions import method_required_tools, prediction_methods, pred_params
from rna_blast_analyze.BR_core.parse_accession import accession_regex

cfg_name = '--config_file'
download_name = '--download_rfam'

with open(os.path.join(os.path.dirname(__file__), 'VERSION'), 'r') as o:
    version = o.read().strip()


class ParseFilter(argparse.Action):
    def __init__(self, *args, **kwargs):
        self.ops1 = {'>', '<', '='}
        self.ops2 = {'>=', '<='}
        super().__init__(*args, **kwargs)

    def __call__(self, parser, args, values, option_string=None):
        vv = values.split(',')
        if len(vv) in {1, 2}:
            setattr(
                args,
                self.dest,
                [self._return_parsed_val(v) for v in vv]
            )
        else:
            raise argparse.ArgumentError(self, 'Incorrect filtering argument.')

    def _return_parsed_val(self, argval):
            vv = argval.strip()
            if vv[:2] in self.ops2:
                try:
                    val = float(vv[2:])
                    return vv[:2], val
                except ValueError:
                    raise argparse.ArgumentError(self, "Filtering value must be a number, not '{}'.".format(vv[2:]))
            elif vv[0] in self.ops1:
                try:
                    val = float(vv[1:])
                    return vv[0], val
                except ValueError:
                    raise argparse.ArgumentError(self, "Filtering value must be a number, not '{}'.".format(vv[1:]))
            else:
                raise argparse.ArgumentError(self, 'Do not understand operator. Allowed are: >, <, =, <= and >=.')


def str2bool(v):
    # from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def f_parser():
    """
    Input parser
    :return: args structure
    """
    parser = argparse.ArgumentParser(
        description='rboAnalyzer - tool for analyzing BLAST search output for RNA sequences',
    )
    input_group = parser.add_argument_group(
        'INPUT'
    )
    output_group = parser.add_argument_group(
        'OUTPUT'
    )
    parameters_group = parser.add_argument_group(
        'PARAMETERS'
    )
    misc_group = parser.add_argument_group(
        'MISC'
    )
    input_group.add_argument(
        '-in',
        '--blast_in',
        type=str,
        required=True,
        metavar='PATH',
        help='BLAST output file with hits to analyze.'
    )
    input_group.add_argument(
        '-q',
        '--blast_query',
        default=None,
        required=True,
        type=str,
        metavar='PATH',
        help='The Blast query fasta file.'
    )
    blastdb_group = parser.add_mutually_exclusive_group(required=True)
    blastdb_group.add_argument(
        '-db',
        '--blast_db',
        metavar='path',
        type=str,
        help='Provide path to blast database, '
             'that is the complete path with blast db name without any extension '
             '(*.nin, nsd, nog, nsi, nhr, nsq, nal).'
    )
    blastdb_group.add_argument(
        '--entrez',
        type=str,
        default=None,
        help='EMAIL - '
             'Indicate that you want to use NCBI Entrez service to download required regions of sequences at runtime. '
             'To comply with NCBI service rules you are required to provide valid email address '
             'at which the NCBI staff could contact you if they need to.'
    )
    misc_group.add_argument(
        '--db_type',
        choices=['blastdb', 'fasta', 'gb', 'server'],
        default='blastdb',
        help="Type of a database provided. "
             "If 'fasta' or 'gb' then --blast_db must be directory "
             "containing files with names in accession.version format. "
             "Example '/home/my-best-db/' with files like 'CP000001.1'. "
    )
    misc_group.add_argument(
        '--b_type',
        default='guess',
        choices=['guess', 'xml', 'plain']
    )
    misc_group.add_argument(
        '--blast_regexp',
        type=str,
        default=accession_regex,
        help='Provide python valid regular expression which capture the index key to blastdb'
        ' (usualy the accession.version number).'
    )
    parameters_group.add_argument(
        '--mode',
        type=str,
        default='locarna',
        choices=['simple', 'locarna', 'meta'],
        help=(
            'Choose mode of hit elongation: '
            'simple (extend by unaligned parts of query) '
            'locarna (run locarna algorithm - uses secondary structure for better alignment) '
            'meta (uses both methods and chooses the alignment which has better RSEARCH score).'
        )
    )
    parameters_group.add_argument(
        '--turbo_fast_preset',
        default=False,
        action='store_true',
        help=(
            "Act's as parameter preset for Turbo-fast setting the max_seqs_in_prediction to 2."
            " This means that only the query sequence is used as homologous sequence."
            " It is useful if analyzing very distant BLAST HITs."
        )
    )
    parameters_group.add_argument(
        '--centroid_fast_preset',
        default=False,
        action='store_true',
        help=(
            "Parameter preset for centroid-fast. Set's the max_seqs_in_prediction to 1."
            " This means that only the query sequence is used as homologous sequence for prediction."
            " It is useful if analyzing very distant BLAST HITs."
        )
    )
    output_group.add_argument(
        '--html',
        metavar='PATH',
        type=str,
        help='Output html file with secondary structure pictures and other useful stuff.'
    )
    misc_group.add_argument(
        '--threads',
        default=None,
        type=int,
        metavar='N',
        help='Number of threads to use (default = N of logical cores detected).'
    )
    output_group.add_argument(
        '--csv',
        default=None,
        help='Output in csv table, infered sequence and structure present.'
    )
    output_group.add_argument(
        '--json',
        type=str,
        metavar='PATH',
        default=None,
        help='Dump all stored data to JSON (developer only - it is possible to convert to all other output formats).'
    )
    mu2 = misc_group.add_mutually_exclusive_group()
    mu2.add_argument(
        '--cm_file',
        default=None,
        type=str,
        metavar='CM_file',
        help='Provided covariance model will be used for homology inference instead of RSEARCH model.'
    )
    mu2.add_argument(
        '--use_rfam',
        action='store_true',
        default=False,
        help='Search in rfam database for covariance model to infer homology with instead of RSEARCH model.'
    )
    misc_group.add_argument(
        download_name,
        action='store_true',
        default=False,
        help='Retrieve RFAM covariance models database. Will download only if new version avalible.'
    )
    misc_group.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + version
    )
    parameters_group.add_argument(
        cfg_name,   # --config_file
        type=str,
        default=None,
        metavar='PATH',
        help='Provide config file if tools and data are in non-default paths.'
    )
    parameters_group.add_argument(
        '-pm',
        '--prediction_method',
        nargs='*',
        type=str,
        metavar='prediction_method_name',
        default=['TurboFold', 'rfam-Rc', 'rnafold'],
        choices=prediction_methods | {'all'},
        help=(
            'Prediction method to use. Multiple prediction methods are allowed. '
            'Possible values: {}'.format(' '.join(prediction_methods))
        )
    )
    parameters_group.add_argument(
        '--pm_param_file',
        type=str,
        metavar='PATH',
        default=None,
        help=(
            "Path to file with parameters for prediction methods in JSON. "
            "Prediction methods not declared within provided file are used with default values. "
            "File is in json format. Default values (also example how to provide parameters) are stored in "
            "'[install location]/rna_blast_analyze/BR_core/prediction_parameters.json'"
        )
    )
    misc_group.add_argument(
        '--logfile',
        type=str,
        default=None,
        metavar='logfile',
        help='Path to where logfile should be written.'
    )
    parameters_group.add_argument(
        '--subseq_window_locarna',
        type=int,
        default=30,
        help='N of nucleotides to add to expected start/end of sequence before realignement. '
             'The unaligned nucleotides are not included in reported sequence.'
    )
    parameters_group.add_argument(
        '--locarna_params',
        type=str,
        default='--struct-local=0 --sequ-local=0 --free-endgaps=++++',
        help=argparse.SUPPRESS
        # help="Parameters for locarna execution. Used when 'mode' is 'locarna' or 'meta'."
    )
    parameters_group.add_argument(
        '--locarna_anchor_length',
        type=int,
        default=7,
        help='Minimal number of adjacent matching bases in BLAST hit to create an anchor for Locarna.'
    )
    parser.add_argument(
        '--repredict_file',
        type=str,
        metavar='PATH',
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        '--dev_pred',
        action='store_true',
        default=False,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        '--dump',
        type=str,
        metavar='PATH',
        default=None,
        help=argparse.SUPPRESS,
        # help='if given, result data will be dump in python pickle, the datastructure can change',
    )
    parser.add_argument(
        '--pandas_dump',
        type=str,
        metavar='PATH',
        default=None,
        help=argparse.SUPPRESS,
        # help='same data as with --csv but in binary (pandas pickle) format'
    )
    # parser.add_argument(
    #     '--no_backup',
    #     action='store_true',
    #     default=False,
    #     help="rboAnalyzer will NOT store the backup data (the [BLAST FILE].r-[HASH] files).",
    # )
    # parser.add_argument(
    #     '--clean',
    #     action='store_true',
    #     default=False,
    #     help="rboAnalyzer will not use the backup data (the [BLAST FILE].r-[HASH] files) even if present."
    # )
    parser.add_argument(
        '--show_gene_browser',
        default=True,
        type=str2bool,
        # help='option to hide gene browser for debugging output web page'
        help=argparse.SUPPRESS,
    )
    mu = misc_group.add_mutually_exclusive_group()
    mu.add_argument(
        '--filter_by_eval',
        default=None,
        action=ParseFilter,
        help='Filter the input blast by E-value. Only hits following the rule will be kept. '
             'Example ">10e-10" will keep only hits with eval greater then 10e-10. '
             'Interval can be specified with "," e.g. ">10e-100, <10e-1". '
             'The homologous sequences used with certain prediction methods are taken from all hits '
             '(regardless of the filtering).'
    )
    mu.add_argument(
        '--filter_by_bitscore',
        default=None,
        action=ParseFilter,
        help='Filter the input blast by bit score. Only hits following the rule will be kept. '
             'Example "<20" will keep only hits with bit score less then 20. '
             'Interval can be specified with "," e.g. ">30, <45". '
             'The homologous sequences used with certain prediction methods are taken from all hits '
             '(regardless of the filtering).'
    )
    misc_group.add_argument(
        '-v', '--verbose',
        dest='verbose',
        action='count',
        default=0,
        help='output verbosity -> most detailed -vv (lot of output)'
    )
    misc_group.add_argument(
        '--skip_missing',
        action='store_true',
        default=False,
        help='If given, the missing records in given blast database will be skipped.'
             ' This may alter the results of bit_score computation (for homology prediction)'
             ' and secondary structure prediction for several methods.'
    )
    misc_group.add_argument(
        '--zip_json',
        action="store_true",
        default=False,
        help=argparse.SUPPRESS
    )
    misc_group.add_argument(
        '--enable_overwrite',
        action='store_true',
        default=False,
        # help="Output will overwrite existing file. (default False)"
        help=argparse.SUPPRESS
    )

    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    args.command = sys.argv

    # handle prediction params input file
    # fallback to default
    default_param_file = os.path.abspath(
            os.path.dirname(__file__) + os.sep + os.path.join('BR_core', 'prediction_parameters.json')
        )
    params = dict()
    with open(default_param_file, 'r') as ff:
        default_params = json.load(ff)
        params.update(default_params)

    if args.pm_param_file:
        with open(args.pm_param_file, 'r') as ff:
            provided_params = json.load(ff)
            params.update(provided_params)

    if args.turbo_fast_preset:
        print('Using Turbo-fast preset')
        params['Turbo-fast']['max_seqs_in_prediction'] = 2

    if args.centroid_fast_preset:
        print('Using centroid-fast preset')
        params['centroid-fast']['max_seqs_in_prediction'] = 1

    if 'all' in args.prediction_method:
        args.prediction_method = list(prediction_methods)

    if args.entrez is not None:
        args.db_type = 'entrez'

    # check if prediction method params are valid
    check_params(params)

    args.pred_params = params
    return args


def check_if_rfam_needed(inargs):
    """
    check if we need to have rfam defined
    :param inargs:
    :return:
    """
    if inargs.use_rfam:
        return True
    elif any([True for pmname in inargs.prediction_method if 'rfam' in pmname]):
        return True
    else:
        return False


def main():
    try:
        # outer envelope for the script
        # ========= perform argument parsing =========
        if download_name in sys.argv and not ('-q' in sys.argv or '--blast_query' in sys.argv):
            # run download rfam here
            # do not run, if given with normal run request
            from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
            from rna_blast_analyze.BR_core import cmalign
            if cfg_name in sys.argv:
                CONFIG.override(tools_paths(config_file=sys.argv[sys.argv.index(cfg_name) + 1]))
            cmalign.download_cmmodels_file()
            # rfam database downloaded
            sys.exit(0)

        args = f_parser()

        _ = lunch_with_args(args)

        # if we reach here, exit with 0
        sys.exit(0)
    except Exception as e:
        print('Something went wrong.')
        try:
            import traceback
            print(
                'The error traceback is written to rboAnalyzer.log . '
                'Please send it along with the query file and BLAST input to the developers.'
            )

            with open('rboAnalyzer.log', 'w') as fd:
                fd.write(str(e))
                fd.write(traceback.format_exc())
        except Exception:
            pass
        sys.exit(1)


def lunch_with_args(args):
    # ========= imports ==========
    # move slow imports here, so the argcomplete would be fast
    import logging
    from rna_blast_analyze.BR_core import BA_verify
    from rna_blast_analyze.BR_core import cmalign
    from rna_blast_analyze.BR_core.config import tools_paths, CONFIG
    from rna_blast_analyze.BR_core.validate_args import validate_args, compute_args_hash
    from rna_blast_analyze.BR_core.luncher import lunch_computation
    from rna_blast_analyze.BR_core.convert_classes import blastsearchrecompute2dict
    from rna_blast_analyze.BR_core.cmalign import RfamInfo

    logger = logging.getLogger('rboAnalyzer')

    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    # set logger level
    logger.setLevel(max(3 - args.verbose, 1) * 10)

    logger.debug('parsed arguments: {}'.format(args))

    # create logging file if requested
    if args.logfile:
        fh = logging.FileHandler(args.logfile)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    start_msg = 'STATUS: starting rboAnalyzer...'
    print(start_msg)
    logger.info(start_msg)

    logger.info('BLAST file: {}'.format(args.blast_in))
    logger.info('Query file: {}'.format(args.blast_query))
    logger.info('BLAST db:   {}'.format(args.blast_db))
    if args.config_file:
        logger.info('configfile: {}'.format(args.config_file))

    # ========= load optional cfg file =========
    CONFIG.override(tools_paths(config_file=args.config_file))

    # ========= check rfam =========
    if not args.download_rfam and not cmalign.check_rfam_present():
        msgfail = 'RFAM models file is not present in specified path. ' \
                  'Please enable rfam download (--download_rfam) or provide prepared directory.'
        logger.error(msgfail)
        sys.exit(1)

    if args.download_rfam:
        cmalign.download_cmmodels_file()

    # ========= check if tools needed for requested methods are installed =========
    BA_verify.check_necessery_tools(methods=args.prediction_method + [args.mode])

    # ========= check if parameters make sense =========
    if not validate_args(args):
        print("There was an error with provided arguments. Please see the error message.")
        sys.exit(1)

    # ========= compute args hash =========
    rfam = RfamInfo()
    hashstring = compute_args_hash(args, os.path.join(rfam.rfam_dir, rfam.gzname))
    setattr(args, 'sha1', hashstring)

    # ========= run =========
    blast_fn = os.path.basename(args.blast_in) + '.r-' + hashstring[:10]
    blast_dir = os.path.dirname(args.blast_in)
    if blast_dir == '':
        blast_dir = os.getcwd()
    potential_matches = [f for f in os.listdir(blast_dir) if f == blast_fn]
    if len(potential_matches) == 0:
        with open(args.blast_in + '.r-' + args.sha1[:10], 'w') as f:
            json.dump(None, f)

    _, results = lunch_computation(args)
    with open(os.path.join(blast_dir, blast_fn), 'w') as f:
        json.dump([blastsearchrecompute2dict(r) for r in results], f, indent=2)

    return results


def check_params(par: dict):
    """
    Verify if provided parameters are valid.
    :param par:
    :return:
    """
    known_methods = set(method_required_tools.keys())
    for provided_key in par.keys():
        if provided_key not in known_methods:
            raise ValueError("Method '{}' does not exist.".format(provided_key))
        for pk in par[provided_key]:
            if pk not in pred_params:
                raise ValueError("Parameter '{}' for method '{}' does not exist.".format(pk, provided_key))


if __name__ == '__main__':
    main()


