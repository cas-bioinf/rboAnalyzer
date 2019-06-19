import os
import re
from io import StringIO
from subprocess import call, check_output, STDOUT
from tempfile import mkstemp
import logging
import gzip
import shutil

from Bio import AlignIO
import pandas as pd

from rna_blast_analyze.BR_core.BA_support import parse_seq_str
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.stockholm_parser import stockholm_read
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)
# this file holds files needed for running and parsing infernal tools


def run_cmscan(fastafile, cmmodels_file=None, params=None, outfile=None, threads=None):
    """
    run cmscan program for finding suitable known CM model for a sequence
    :return:
    """
    ml.info('Runing cmscan.')
    ml.debug(fname())
    rfam = RfamInfo()

    if outfile:
        out = outfile
    else:
        fd, out = mkstemp(prefix='rba_', suffix='_10', dir=CONFIG.tmpdir)
        os.close(fd)

    if threads:
        params += '--cpu {}'.format(threads)

    if cmmodels_file:
        cm_file = cmmodels_file
    else:
        cm_file = os.path.join(rfam.rfam_dir, rfam.rfam_file_name)

    with open(os.devnull, 'w') as FNULL:

        # build commandline
        cmd = ['{}cmscan'.format(CONFIG.infernal_path)]
        if params != '':
            cmd += params.split()
        cmd += [
            '--tblout', out,
            cm_file, fastafile
        ]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'cmscan failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)

    return out


def run_cmfetch(cmfile, modelid, outfile=None):
    """

    :param cmfile:
    :param modelid:
    :return:
    """
    ml.info('Runing cmfetch.')
    ml.debug(fname())
    if outfile:
        out = outfile
    else:
        fd, out = mkstemp(prefix='rba_', suffix='_11', dir=CONFIG.tmpdir)
        os.close(fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = [
            '{}cmfetch'.format(CONFIG.infernal_path),
            '-o', out,
            cmfile,
            modelid
        ]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'call to cmfetch failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)
        return out


def run_cmemit(model, params='', out_file=None):
    """

    :param model:
    :param params:
    :return:
    """
    ml.info('Run cmemit.')
    ml.debug(fname())
    if out_file:
        out = out_file
    else:
        fd, out = mkstemp(prefix='rba_', suffix='_12', dir=CONFIG.tmpdir)
        os.close(fd)

    with open(os.devnull, 'w') as FNULL:
        # build commandline
        cmd = ['{}cmemit'.format(CONFIG.infernal_path)]
        if params != '':
            cmd += params.split()
        cmd += ['-o', out, model]

        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'call to cmemit failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)
    return out


def extract_ref_structure_fromRFAM_CM(model_name):
    """
    Extract reference structure encoded in covariance model to dot bracket notation.
    :param model_name: model name in cm file
    :return: string
    """
    ml.debug(fname())
    rfam = RfamInfo()

    single_cm_file = run_cmfetch(rfam.file_path, model_name)

    ref_structure = extract_ref_from_cm(single_cm_file)
    os.remove(single_cm_file)
    return ref_structure


def extract_ref_from_cm(cm_file):
    ml.debug(fname())
    single_alig_file = run_cmemit(cm_file, params='-a -N 1')
    o = open(single_alig_file, 'r')
    salig = stockholm_read(o)
    o.close()

    os.remove(single_alig_file)

    if len(salig) != 1:
        raise Exception('cmemit file does not have only one record in (not including reference)')

    # recode structure, return it
    ss = salig.column_annotations['SS_cons']
    # inserts = str(salig[0].seq)
    inserts = str(salig.column_annotations['RF'])

    gapchars = '.~'

    structure_list = []
    for i, j in zip(ss, inserts):
        if i in gapchars:
            continue
        structure_list.append(i)

    # recode structure list
    structure = cm_strucutre2br(''.join(structure_list))
    return structure


def cm_strucutre2br(seq):
    """
    Converts cm structure from stockholm notation to dot bracket.
    replaces all bracketing chars with coresponding round brackets
    replaces all chars for unpaired sequence with dot
    :param seq:
    :return:
    """
    ns = re.sub('[-:,_~]', '.', seq)
    ns = re.sub('[<({\[]', '(', ns)
    ns = re.sub('[>)}\]]', ')', ns)
    return ns


class RfamInfo(object):
    def __init__(self):
        self.rfam_dir = CONFIG.rfam_dir
        self.rfam_file_name = 'Rfam.cm'
        self.url = CONFIG.rfam_url
        self.gzname = 'Rfam.cm.gz'
        self.file_path = os.path.join(self.rfam_dir, self.rfam_file_name)


def check_rfam_present():
    """
    Check if RFAM file is present and converted to binary format required by cmscan
     by running program cmpres.
    If present but not converted, conversion is attempted.

    :return: bool
    """
    ml.debug(fname())
    rfam = RfamInfo()
    cm_present = os.path.isfile(rfam.file_path)
    if cm_present:
        if not check_if_cmpress_processed():
            run_cmpress(rfam.file_path)
        return True
    else:
        return False


def check_if_cmpress_processed():
    """
    Check if we have cmpressed file in rfam directory
    This is defined by presence of binary files [i1f, i1i, i1m, i1p]
    :return:
    """
    ml.debug(fname())
    rfam = RfamInfo()
    files = os.listdir(rfam.rfam_dir)
    for suff in ['.i1f', '.i1i', '.i1m', '.i1p']:
        if rfam.rfam_file_name + suff not in files:
            return False
    return True


def download_cmmodels_file(path=None, url=None):
    """
    downloads cm model from rfam database
    default retrieve url is: 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz'
    :param path:
    :param url:
    :return:
    """
    print('Running CM download from RFAM.')
    ml.debug(fname())
    rfam = RfamInfo()
    if path is None:
        path = rfam.rfam_dir
    if url is None:
        url = rfam.url

    if not os.path.exists(path):
        os.makedirs(path)

    cmd = ['wget', '-N', '-P', path, url]
    ml.debug(cmd)

    ml.info('Downloading RFAM database (aprox 300Mb). This may take a while...')
    r = check_output(cmd, stderr=STDOUT)

    if not r:
        msgfail = 'call to wget failed'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)

    if 'Remote file no newer than local file ‘Rfam.cm.gz’ -- not retrieving.' in r.decode():
        # do not download
        ml.info('No new data. Nothing to do.')
    else:
        # unzip using build in gzip
        with gzip.open(os.path.join(path, rfam.gzname), 'rb') as fin:
            with open(os.path.join(path, rfam.rfam_file_name), 'wb') as fout:
                shutil.copyfileobj(fin, fout)

        # run cmpress to create binary files needed to run cmscan
        run_cmpress(os.path.join(path, rfam.rfam_file_name))

    return os.path.join(path, rfam.rfam_file_name)


def run_cmpress(file2process):
    ml.info('Running cmpress.')
    ml.debug(fname())
    with open(os.devnull, 'w') as FNULL:
        cmd = ['{}cmpress'.format(CONFIG.infernal_path), '-F', file2process]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL, stderr=FNULL)

        if r:
            msgfail = 'call to cmpress failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)


def run_cmbuild(cmbuild_input_file, cmbuild_params=''):
    """
    run cmbuild procedure
    input must be MSA in stockholm format with secondary structure prediction

    note: consider what to do if only one sequence is available

    :param cmbuild_input_file: Stockholm or selex alignment file
    :param cmbuild_params: additional params to cmbuild
    :return:
    """
    ml.info('Runing cmbuild.')
    ml.debug(fname())
    cm_fd, cm_file = mkstemp(prefix='rba_', suffix='_13', dir=CONFIG.tmpdir)
    os.close(cm_fd)

    FNULL = open(os.devnull, 'w')
    try:
        cmd = ['{}cmbuild'.format(CONFIG.infernal_path), '-F']
        if cmbuild_params != '':
            cmd += cmbuild_params.split()
        cmd += [cm_file, cmbuild_input_file]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'Call to cmbuild failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(
                '{} - files: {} {} \nparams: {}'.format(
                    msgfail,
                    cmbuild_input_file,
                    cm_file,
                    cmbuild_params
                )
            )
    finally:
        FNULL.close()

    return cm_file


def run_cmalign_on_fasta(fasta_file, model_file, cmalign_params='--notrunc', alig_format='stockholm'):
    """
    run cmalign program with provided CM model file
    :param fasta_file: input fasta to be aligned to cm model
    :param model_file: file containing one or more cm models
    :param cmalign_params: parameter of the search
    :return:
    """
    ml.info('Runing cmaling.')
    ml.debug(fname())
    cma_fd, cma_file = mkstemp(prefix='rba_', suffix='_14', dir=CONFIG.tmpdir)
    os.close(cma_fd)
    FNULL = open(os.devnull, 'w')
    try:
        cmd = [
            '{}cmalign'.format(CONFIG.infernal_path),
            '--informat', 'fasta',
            '--outformat', alig_format,
        ]
        if cmalign_params != '':
            cmd += cmalign_params.split()
        cmd += ['-o', cma_file, model_file, fasta_file]

        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'call to cmalign failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError('{} for files:\nmodel file:{}fasta file:{}\noutput file:{}'.format(
                msgfail, model_file, fasta_file, cma_file))
    finally:
        FNULL.close()

    return cma_file


def build_stockholm_from_clustal_alig(clustal_file, alif_file):
    """
    build stockholm alignment
    :return:
    """
    ml.debug(fname())
    with open(clustal_file, 'r') as cf, open(alif_file, 'r') as af:
        # write stockholm align to buffer and read it with my parser
        clust = AlignIO.read(cf, format='clustal')
        temp = StringIO()
        AlignIO.write(clust, temp, format='stockholm')
        temp.seek(0)
        st_alig = stockholm_read(temp)

        # parse alifold output and add structure to stockholm alignment
        for i, alif in enumerate(parse_seq_str(af)):
            alifold_structure = alif.letter_annotations['ss0']
            st_alig.column_annotations['SS_cons'] = alifold_structure
            if i == 0:
                break

        st_fd, st_file = mkstemp(prefix='rba_', suffix='_15', dir=CONFIG.tmpdir)
        with os.fdopen(st_fd, 'w') as sf:
            st_alig.write_stockholm(sf)

            return st_file


def read_cmalign_sfile(f):
    """
    reads cmalign optional sfile output
    beware that index here is counted from 1, rather then zero
    :param f: file handle or file path
    :return:
    """
    # there is an issue that if the table has more then 10000 entries, the header is repeated
    ml.debug(fname())
    sfile = pd.read_csv(
        f,
        skiprows=4,
        header=None,
        names=(
            'seq_name', 'length', 'cm_from', 'cm_to', 'trunc', 'bit_sc',
            'avg_pp', 'band_calc', 'alignment', 'total', 'mem'
        ),
        sep='\s+',
        comment='#'
    )

    return sfile


def parse_cmalign_infernal_table(tbl):
    # first row => names but with spaces
    # second row => guide line
    ml.debug(fname())
    expected_names = ['#target_name',
                      'accession',
                      'query_name',
                      'accession',
                      'mdl',
                      'mdl_from',
                      'mdl_to',
                      'seq_from',
                      'seq_to',
                      'strand',
                      'trunc',
                      'pass',
                      'gc',
                      'bias',
                      'score',
                      'E-value',
                      'inc',
                      'description_of_target']

    output_names = ['target_name',
                    'accession_seq',
                    'query_name',
                    'accession_mdl',
                    'mdl',
                    'mld_from',
                    'mld_to',
                    'seq_from',
                    'seq_to',
                    'strand',
                    'trunc',
                    'pass',
                    'gc',
                    'bias',
                    'score',
                    'E-value',
                    'inc',
                    'description_of_target']

    name_line = tbl.readline()
    lc_line = tbl.readline()
    lc_loc = [i for i in re.finditer('#?-+', lc_line)]
    span = [list(i.span()) for i in lc_loc]
    span[-1][1] = None

    out = dict()
    for i, (s, e) in enumerate(span):
        name = re.sub('\s+', '_', name_line[s:e].strip())
        if name != expected_names[i]:
            raise Exception('unknown column name in table {}'
                            '\nUpdate expected_names and output_names variables.'.format(name))
        out[output_names[i]] = []

    txt = tbl.readline()
    while txt:
        if txt[0] == '#':
            txt = tbl.readline()
            continue

        for i, (s, e) in enumerate(span):
            out[output_names[i]].append(txt[s:e].strip())
        txt = tbl.readline()

    pd_out = pd.DataFrame.from_dict(out)
    pd_out['E-value'] = pd_out['E-value'].astype('float')
    pd_out['score'] = pd_out['score'].astype('float')
    return pd_out


def get_cm_model(query_file, params=None, threads=None):
    ml.debug(fname())
    cmscan_data = get_cm_model_table(query_file, params, threads)
    best_model_row = select_best_matching_model_from_cmscan(cmscan_data)

    best_model = best_model_row['target_name']

    ml.info('Best matching model: {}'.format(best_model))
    return best_model


def get_cm_model_table(query_file, params=None, threads=None):
    ml.debug(fname())
    if params is None:
        params = dict()

    cmscan_params = '-g '
    if params and ('cmscan' in params) and params['cmscan']:
        cmscan_params += params['cmscan']
    out_table = run_cmscan(query_file, params=cmscan_params, threads=threads)
    f = open(out_table, 'r')
    cmscan_data = parse_cmalign_infernal_table(f)
    f.close()
    os.remove(out_table)
    return cmscan_data


def select_best_matching_model_from_cmscan(cmscan_data):
    ei = cmscan_data['E-value'].idxmin()
    si = cmscan_data['score'].idxmax()
    assert ei == si

    best_model = cmscan_data.loc[ei].to_dict()
    return best_model
