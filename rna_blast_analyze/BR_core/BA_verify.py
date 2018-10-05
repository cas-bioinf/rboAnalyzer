import logging
import os
import re
import sys
from subprocess import check_output, STDOUT, CalledProcessError

from rna_blast_analyze.BR_core import cmalign
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.tools_versions import blast_minimal_version, locarna_minimal_version, \
    infernal_minimal_version, vrna_minimal_version, clustalo_minimal_version, muscle_minimal_version, \
    tcoffee_rcoffee_minimal_version, centroid_homfold_minimal_version, turbofold_minimal_version, mfold_minimal_version, \
    pred_method_required_tools

ml = logging.getLogger(__name__)


def verify_query_blast(blast, query):
    """
    verify if query from fasta file matches query sequence described in BLAST.
    :param blast:
    :param query:
    :return:
    """
    if not (query.id == blast.query or query.description == blast.query):
        ml.warn(
            'Provided query id ({}) do not match query id in BLAST output ({})'.format(
                query.id,
                blast.query
            )
        )
    if len(query) != blast.query_length:
        ml.warn(
            'Provided query lenght ({}: {}) do not match BLAST query length ({}: {}).'.format(
                query.id,
                len(query),
                blast.query,
                blast.query_length
            )
        )
    return


def verify_blastdbcmd(minimal_version):
    """verify if blastdbcmd is present in supported version
    """
    msgversion = 'blastcmd not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}blastcmd could not be located (not in PATH)'.format(CONFIG.blast_path)
    msgsuccess = 'blastcmd is installed in required version'
    try:
        a = check_output(
            [
                '{}blastdbcmd'.format(CONFIG.blast_path),
                '-version'
            ]
        )
        a = a.decode(encoding='utf-8')
        r = re.search('(?<=blastdbcmd: )[0-9.]+', a)
        if r:
            ver = r.group().split('.')
            ver = [int(i) for i in ver]
            for v, minv in zip(ver, minimal_version):
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_locarna(minimal_version):
    msgversion = 'Locarna is not installed in required version, required version is {}.{}.{}'.format(*locarna_minimal_version)
    msgpath = '{}LocARNA could not be located (not in PATH).'.format(CONFIG.locarna_path)
    msgsuccess = 'Locarna is installed in required version'
    try:
        a = check_output(
            [
                '{}locarna'.format(CONFIG.locarna_path),
                '--version'
            ]
        )
        a = a.decode()
        if a.startswith('LocARNA'):
            r = re.finditer('[0-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_infernal(program, minimal_version):
    msgversion = '{} (part of INFERNAL) is not installed in required version, required is {}.{}.{}'.format(program, *minimal_version)
    msgpath = '{}{} could not be located (not in PATH).'.format(CONFIG.infernal_path, program)
    msgsuccess = '{} is installed in required version'.format(program)
    try:
        a = check_output(
            [
                CONFIG.infernal_path + program,
                '-h'
            ]
        )
        a = a.decode()
        if a.startswith('# {}'.format(program)):
            r = re.search('(?<=# INFERNAL )[0-9.]+', a)
            ver = r.group().split('.')
            ver = [int(i) for i in ver]
            for v, minv in zip(ver, minimal_version):
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        raise


def verify_viennarna_program(program, minimal_version):
    msgversion = '{} is not installed in required version, required version is {}.{}.{}'.format(program, *minimal_version)
    msgpath = '{}{} could not be located'.format(CONFIG.viennarna_path, program)
    msgsuccess = '{} is installed in required version'.format(program)
    try:
        a = check_output(
            [
                CONFIG.viennarna_path + program,
                '--version'
            ]
        )
        a = a.decode()
        if a.startswith(program):
            ver = a.split()[1].split('.')
            for v, minv in zip(ver, minimal_version):
                if int(v) > minv:
                    ml.info(msgsuccess)
                    return True
                elif int(v) < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_viennarna(programs, vrna_minv):
    installed = set()
    for prog in programs:
        if verify_viennarna_program(prog, vrna_minv):
            installed.add(prog)
    return installed


def verify_vrna_refold():
    msgversion = 'refold.pl is not installed in required version.'
    msgpath = '{}refold.pl could not be located (not in PATH).'.format(CONFIG.refold_path)
    msgsuccess = 'refold.pl is instaled in required version'
    try:
        try:
            a = check_output(
                [
                    CONFIG.refold_path + 'refold.pl',
                    '-h'
                ],
                stderr=STDOUT
            )
        except CalledProcessError as e:
            a = e.output

        a = a.decode()
        if a.strip().startswith('refold.pl [-t threshold] myseqs.aln alidot.ps | RNAfold -C'):
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_clustalo(minimal_version):
    msgversion = 'clustalo is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}clustalo could not be located (not in PATH).'.format(CONFIG.clustal_path)
    msgsuccess = 'clustalo is installed in required version'
    try:
        a = check_output(
            [
                '{}clustalo'.format(CONFIG.clustal_path),
                '--version'
            ]
        )
        a = a.decode()
        if a:
            r = re.finditer('[0-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_muscle(minimal_version):
    msgversion = 'muscle is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}muscle could not be located (not in PATH).'.format(CONFIG.muscle_path)
    msgsuccess = 'muslce is installed in required version'
    try:
        a = check_output(
            [
                '{}muscle'.format(CONFIG.muscle_path),
                '-version'
            ]
        )
        a = a.decode()
        if a.startswith('MUSCLE'):
            r = re.finditer('[0-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_tcoffee(minimal_version):
    # do not check the revision
    msgversion = 't_coffee is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}t_coffee could not be located (not in PATH).'.format(CONFIG.tcoffee_path)
    msgsuccess = 't_coffee is installed in required version'
    try:
        a = check_output(
            [
                '{}t_coffee'.format(CONFIG.tcoffee_path),
                '-version'
            ]
        )
        a = a.decode()
        if re.search('PROGRAM: T-COFFEE', a):
            im = re.search('(?<=version_)[0-9.]+', a, flags=re.IGNORECASE)
            if im:
                b = im.group()
                r = re.finditer('[0-9]+', b)
                for match, minv in zip(r, minimal_version):
                    v = int(match.group())
                    if v > minv:
                        ml.info(msgsuccess)
                        return True
                    elif v < minv:
                        ml.warn(msgversion)
                        return False
                ml.info(msgsuccess)
                return True
            else:
                ml.warn(msgversion)
                return False
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_rcoffee(minimal_version):
    msgversion = 't_coffee -mode rcoffee is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}t_coffee -mode rcoffee could not be located (not in PATH).'.format(CONFIG.muscle_path)
    msgsuccess = 't_coffee -mode rcoffee is installed in required version'

    if not verify_tcoffee(minimal_version):
        ml.warn(
            't_coffee is not installed in required version, required is {}.{}'.format(
                *minimal_version
            )
        )
        return False
    try:
        try:
            a = check_output(
                [
                    '{}t_coffee'.format(CONFIG.muscle_path),
                    '-mode',
                    'rcoffee'
                ],
                stderr=STDOUT
            )
        except CalledProcessError as e:
            a = e.output

        a = a.decode()
        t1 = re.search('-- ERROR: You have not provided any sequence', a)
        t2 = re.search('-- COM:(.)*t_coffee -mode rcoffee', a)

        n1 = re.search("ERROR: special_mode rcoffee is unknown \[FATAL:T-COFFEE\]", a)

        if t1 and t2 and not n1:
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False

    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_centroid_homfold(minimal_version):
    msgversion = 'centroid_homfold is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}centroid_homfold could not be located (not in PATH).'.format(CONFIG.centriod_path)
    msgsuccess = 'centroid_homfold is installed in required version'
    # double try because help returns exit status 1
    try:
        try:
            a = check_output(
                [
                    '{}centroid_homfold'.format(CONFIG.centriod_path),
                    '-h'
                ],
                stderr=STDOUT
            )
        except CalledProcessError as e:
            a = e.output

        a = a.decode()
        b = a.split()
        if b[0] == 'CentroidHomfold':
            r = re.finditer('[0-9]+', b[1])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_turbofold(minimal_version):
    msgversion = 'TurboFold is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}TurboFold could not be located (not in PATH).'.format(CONFIG.turbofold_path)
    msgsuccess = 'TruboFold is installed in required version'
    try:
        try:
            a = check_output(
                [
                    '{}TurboFold'.format(CONFIG.turbofold_path),
                    '-v'
                ],
                stderr=STDOUT
            )
        except CalledProcessError as e:
            a = e.output

        a = a.decode()
        b = a.split()
        if b[0] == 'TurboFold:':
            r = re.finditer('[0-9]+', b[2])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def verify_turbofold_datapath():
    """
    verify that datapath enviroment variable is set (it is not set by default when installing turbofold from conda)
    :return:
    """
    try:
        dp = os.environ['DATAPATH']
        return True
    except KeyError:
        if CONFIG.rnastructure_datapath is None:
            return False
        else:
            return True


def verify_mfold(minimal_version):
    msgversion = 'hybrid-ss-min (UNAfold) is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}hybrid-ss-min could not be located (not in PATH).'.format(CONFIG.mfold_path)
    msgsuccess = 'hybrid-ss-min (UNAfold) is installed in required version'
    try:
        a = check_output(
            [
                '{}hybrid-ss-min'.format(CONFIG.mfold_path),
                '-V'
            ]
        )
        a = a.decode()
        b = a.split()
        if b[0] == 'hybrid-ss-min':
            r = re.finditer('[0-9]+', b[2])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    ml.info(msgsuccess)
                    return True
                elif v < minv:
                    ml.warn(msgversion)
                    return False
            ml.info(msgsuccess)
            return True
        else:
            ml.warn(msgversion)
            return False
    except FileNotFoundError:
        ml.warn(msgpath)
        return False


def check_3rd_party_tools():
    installed = set()
    if verify_blastdbcmd(blast_minimal_version):
        installed.add('blastdbcmd')

    if verify_locarna(locarna_minimal_version):
        installed.add('locarna')

    if verify_infernal('cmbuild', infernal_minimal_version):
        installed.add('infernal')

    if verify_infernal('cmalign', infernal_minimal_version):
        installed.add('infernal')

    installed |= verify_viennarna(['RNAfold', 'RNAplot', 'RNAdistance'], vrna_minimal_version)

    return installed


def check_3rd_party_data():
    installed = set()
    if cmalign.check_rfam_present():
        installed.add('rfam')
    return installed


def check_3rd_party_prediction_tools():
    # clustalo
    # alifold
    # tcoffee rcoffee
    # refold
    # centroid_homfold
    # TurboFold
    # mfold (hybrid-ss-min)

    installed = set()

    if verify_viennarna(['RNAalifold',], vrna_minimal_version):
        installed.add('RNAalifold')

    if verify_vrna_refold():
        installed.add('refold.pl')

    if verify_clustalo(clustalo_minimal_version):
        installed.add('clustalo')

    if verify_rcoffee(tcoffee_rcoffee_minimal_version):
        installed.add('tcoffee')
        installed.add('rcoffee')

    if verify_centroid_homfold(centroid_homfold_minimal_version):
        installed.add('centroid_homfold')

    if verify_turbofold(turbofold_minimal_version):
        installed.add('turbofold')

    if verify_mfold(mfold_minimal_version):
        installed.add('mfold')

    if verify_muscle(muscle_minimal_version):
        installed.add('muscle')

    return installed


def check_necessery_tools(methods):
    avalible_tools = check_3rd_party_tools()
    avalible_tools |= check_3rd_party_data()
    avalible_tools |= check_3rd_party_prediction_tools()

    for met in methods:
        needed = pred_method_required_tools[met] - avalible_tools

        if 'refold.pl' in needed:
            msgfail = 'Please add the refold.pl to PATH or add the path to refold.pl to configfile.'
            is_conda = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
            if is_conda:
                ml.info('Trying to find refold.pl in "CONDA_ROOT/share"')
                out = check_output('find {}/share -type f -name refold.pl'.format(sys.prefix), shell=True)
                op = out.decode().strip()
                if op != '':
                    op = os.path.dirname(op.split('/n')[0])
                    ml.info('Inferred refold.pl in {}'.format(op))
                    ml.info('writing configuration to {}'.format(CONFIG.conf_file))
                    CONFIG.tool_paths['refold'] = op
                    if 'TOOL_PATHS' not in CONFIG.config_obj:
                        CONFIG.config_obj['TOOL_PATHS'] = {}
                    CONFIG.config_obj['TOOL_PATHS']['refold'] = op
                    with open(CONFIG.conf_file, 'w') as updated_cfh:
                        CONFIG.config_obj.write(updated_cfh)
                    avalible_tools.add('refold.pl')
                    needed.remove('refold.pl')
                else:
                    ml.error(msgfail)
                    raise EnvironmentError(msgfail)
            else:
                ml.error(msgfail)

        if needed:
            msgfail = 'Missing {} (needed for {}).'.format(' '.join(needed), met)
            ml.error(msgfail)
            raise EnvironmentError(msgfail)

    if 'TurboFold' in methods or 'TurboFold_conservative' in methods:
        msgfail = 'Please provide DATAPATH for Turbofold from RNAstructure package. Either as DATAPATH ' \
                  'environment variable or as rnastructure_DATAPATH entry in config DATA section.'
        if not verify_turbofold_datapath():
            ml.warn('The turbofold is installed but the DATAPATH environment variable is not set nor present in config.txt.')
            is_conda = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
            if is_conda:
                ml.info('Trying to find required data in "CONDA_ROOT/share"')
                out = check_output('find {}/share -type d -name data_tables'.format(sys.prefix), shell=True)
                op = out.decode().strip()
                if op != '':
                    op = op.split('/n')[0]
                    ml.info('Inferred datapath in {}'.format(op))
                    ml.info('writing configuration to {}'.format(CONFIG.conf_file))
                    CONFIG.data_paths['rnastructure_datapath'] = op
                    if 'DATA' not in CONFIG.config_obj:
                        CONFIG.config_obj['DATA'] = {}
                    CONFIG.config_obj['DATA']['rnastructure_datapath'] = op
                    with open(CONFIG.conf_file, 'w') as updated_cfh:
                        CONFIG.config_obj.write(updated_cfh)
                else:
                    ml.error(msgfail)
                    raise EnvironmentError(msgfail)
            else:
                ml.error(msgfail)
                raise EnvironmentError(msgfail)

if __name__ == '__main__':
    print(check_3rd_party_tools())
    print(check_3rd_party_prediction_tools())
