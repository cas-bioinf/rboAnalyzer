import os
import re
from io import StringIO
from subprocess import call, check_output
from tempfile import mkstemp

from Bio import AlignIO
import pandas as pd

from rna_blast_analyze.BR_core.BA_support import parse_seq_str
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.stockholm_parser import stockholm_read


# this file holds files needed for running and parsing infernal tools


def run_cmscan(fastafile, cmmodels_file=None, params=None, outfile=None, threads=None):
    """
    run cmscan program for finding suitable known CM model for a sequence
    :return:
    """

    rfam = RfamInfo()

    if outfile:
        out = outfile
    else:
        fd, out = mkstemp()
        os.close(fd)

    if threads:
        params += '--cpu {}'.format(threads)

    if cmmodels_file:
        cm_file = cmmodels_file
    else:
        cm_file = os.path.join(rfam.rfam_dir, rfam.rfam_file_name)

    with open(os.devnull, 'w') as FNULL:

        cmd = '{}cmscan {} --tblout {} {} {}'.format(
            CONFIG.infernal_path,
            params,
            out,
            cm_file,
            fastafile
        )

        r = call(
            cmd,
            shell=True,
            stdout=FNULL
        )

        if r:
            print('cmd: {}'.format(cmd))
            raise ChildProcessError('cmscan failed')

    return out


def run_cmfetch(cmfile, modelid, outfile=None):
    """

    :param cmfile:
    :param modelid:
    :return:
    """
    if outfile:
        out = outfile
    else:
        fd, out = mkstemp()
        os.close(fd)

    r = call(
        '{}cmfetch -o {} {} {}'.format(
            CONFIG.infernal_path,
            out,
            cmfile,
            modelid
        ),
        shell=True
    )
    if r:
        raise ChildProcessError('call to cmfetch failed')
    return out


def run_cmemit(model, params='', out_file=None):
    """

    :param model:
    :param params:
    :return:
    """
    if out_file:
        out = out_file
    else:
        fd, out = mkstemp()
        os.close(fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = '{}cmemit {} -o {} {}'.format(
            CONFIG.infernal_path,
            params,
            out,
            model
        )
        r = call(
            cmd,
            shell=True,
            stdout=FNULL
        )
        if r:
            print('call to cmemit failed')
            print('cmd: {}'.format(cmd))
            raise ChildProcessError('cmemit failed')
    return out


def extract_ref_structure_fromRFAM_CM(model_name):
    """
    Extract reference structure encoded in covariance model to dot bracket notation.
    :param model_name: model name in cm file
    :return: string
    """
    rfam = RfamInfo()

    single_cm_file = run_cmfetch(rfam.file_path, model_name)

    ref_structure = extract_ref_from_cm(single_cm_file)
    os.remove(single_cm_file)
    return ref_structure


def extract_ref_from_cm(cm_file):
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
    rfam = RfamInfo()
    if not path:
        path = rfam.rfam_dir
    if not url:
        url = rfam.url

    if not os.path.exists(path):
        os.makedirs(path)

    filename = rfam.gzname
    try:
        r = check_output('wget -N -P {} {}'.format(path, url), shell=True)

        if 'Remote file no newer than local file ‘Rfam.cm.gz’ -- not retrieving.' in r:
            # do not download
            pass
        else:
            r = call(
                'gunzip -k -f {}'.format(
                    os.path.join(
                        path,
                        filename
                    )
                ),
                shell=True
            )
            if r:
                print('call to gunzip failed')
                print('cmd: {}'.format('gunzip -k -f {}'.format(os.path.join(path, filename))))
                raise ChildProcessError('call to tar -xzf failed')

            # run cmpress to create binary files needed to run cmscan
            run_cmpress(os.path.join(path, filename))

    except ChildProcessError:
        print('call to wget failed')
        print('cmd: {}'.format('wget -N -P {} {}'.format(path, url)))
        raise

    return os.path.join(path, filename)


def run_cmpress(file2process):
    r = call(
        '{}cmpress -F {}'.format(
            CONFIG.infernal_path,
            file2process
        ), shell=True
    )
    if r:
        print('call to cmpress failed')
        print('cmd: {}'.format('{}cmpress -F {}'.format(CONFIG.infernal_path, file2process)))
        raise ChildProcessError('call to cmpress failed')


def run_cmbuild(cmbuild_input_file, cmbuild_params=''):
    """
    run cmbuild procedure
    input must be MSA in stockholm format with secondary structure prediction

    note: consider what to do if only one sequence is available

    :param cmbuild_input_file: Stockholm or selex alignment file
    :param cmbuild_params: additional params to cmbuild
    :return:
    """
    cm_fd, cm_file = mkstemp()
    os.close(cm_fd)
    print('run_cmbuild for input:{}'.format(cmbuild_input_file))
    FNULL = open(os.devnull, 'w')
    try:
        r = call(
            '{}cmbuild -F {} {} {}'.format(
                CONFIG.infernal_path,
                cmbuild_params,
                cm_file,
                cmbuild_input_file
            ),
            shell=True,
            stdout=FNULL
        )

        if r:
            # print('if oneThread parameter is used, the numactl program must be installed ()')
            raise ChildProcessError(
                'Call to cmbuild failed - files: {} {} \nparams: {}'.format(
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
    cma_fd, cma_file = mkstemp()
    os.close(cma_fd)
    FNULL = open(os.devnull, 'w')
    print('run_cmalign_on_fasta for model_file:{} fasta_file:{} output_file:{}'.format(model_file,
                                                                                       fasta_file,
                                                                                       cma_file))
    try:
        r = call(
            '{}cmalign --informat fasta --outformat {} {} -o {} {} {}'.format(
                CONFIG.infernal_path,
                alig_format,
                cmalign_params,
                cma_file,
                model_file,
                fasta_file,
            ),
            shell=True,
            stdout=FNULL
        )
        if r:
            raise ChildProcessError('call to cmalign failed for files:\n'
                                    'model file:{}\n'
                                    'fasta file:{}\n'
                                    'output file:{}'.format(model_file, fasta_file, cma_file))
    finally:
        FNULL.close()

    return cma_file


def build_stockholm_from_clustal_alig(clustal_file, alif_file):
    """
    build stockholm alignment
    :return:
    """
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

        st_fd, st_file = mkstemp()
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

    sfile = pd.read_table(
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
    # todo revrite this to be fully featured pd.table read not this hacky function returning strings
    # first row => names but with spaces
    # second row => guide line

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
    if params is None:
        params = dict()

    cmscan_params = '-g '
    if params and ('cmscan' in params) and params['cmscan']:
        cmscan_params += params['cmscan']
    out_table = run_cmscan(query_file, params=cmscan_params, threads=threads)
    f = open(out_table, 'r')
    cmscan_data = parse_cmalign_infernal_table(f)
    f.close()

    ei = cmscan_data['E-value'].idxmin()
    si = cmscan_data['score'].idxmax()
    assert ei == si

    best_model = cmscan_data['target_name'][ei]

    os.remove(out_table)

    return best_model
