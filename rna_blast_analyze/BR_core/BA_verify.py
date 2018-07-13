import re
import warnings
from subprocess import check_output, STDOUT, CalledProcessError

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core import cmalign

blast_minimal_version = [2, 2, 28]
locarna_minimal_version = [1, 9, 2]
infernal_minimal_version = [1, 1, 2]
vrna_minimal_version = [2, 3, 5]
clustalo_minimal_version = [1, 2, 4]
muscle_minimal_version = [3, 8, 31]
tcoffee_rcoffee_minimal_version = [10, 0]
centroid_homfold_minimal_version = [0, 0, 15]
turbofold_minimal_version = [6, 0]
mfold_minimal_version = [3, 6]

# always needed: clustal, muscle, rnafold

pred_method_required_tools = {
    # ===== All possible prediction methods =====
    #  list only methods which are not required by default
    'rfam_rnafoldc': {'rfam',},
    'rfam_subopt': {'mfold',},
    'rfam_rapidshapes': {'rfam', 'rapidshapes',},
    'clustalo_alifold_rapidshapes': {'RNAalifold', 'rapidshapes',},
    'muscle_alifold_rapidshapes': {'RNAalifold', 'rapidshapes',},
    'rcoffee_alifold_rapidshapes': {'rcoffee', 'RNAalifold', 'rapidshapes',},
    'alifold_refold': {'RNAalifold', 'refold.pl',},
    'muscle_alifold_refold': {'RNAalifold', 'refold.pl',},
    'rnafold': set(),
    'subopt_fold_query': {'mfold', 'RNAdistance',},
    'subopt_fold_clustal_alifold': {'mfold', 'RNAalifold', 'RNAdistance',},
    'subopt_fold_muscle_alifold': {'mfold', 'RNAalifold', 'RNAdistance',},
    'alifold_refold_rnafold_c': {'RNAalifold', 'refold.pl',},
    'muscle_alifold_refold_rnafold_c': {'RNAalifold', 'refold.pl',},
    'alifold_unpaired_conserved_refold': {'RNAalifold', 'refold.pl',},
    'muscle_alifold_unpaired_conserved_refold': {'RNAalifold', 'refold.pl',},
    'dh_tcoffee_alifold_refold': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_tcoffee_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_tcoffee_alifold_conserved_ss_rnafoldc': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_clustal_alifold_refold': {'RNAalifold', 'refold.pl',},
    'dh_clustal_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl',},
    'dh_clustal_alifold_conserved_ss_rnafoldc': {'RNAalifold', 'refold.pl',},
    'pairwise_centroid_homfold': {'centroid_homfold',},
    'TurboFold_conservative': {'turbofold',},
    'TurboFold': {'turbofold'},
    'tcoffee_rcoffee_alifold_refold': {'rcoffee', 'RNAalifold', 'refold.pl',},
    'tcoffee_rcoffee_alifold_refold_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl',},
    'tcoffee_rcoffee_alifold_conserved_ss_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl',},
}


def verify_query_blast(blast, query):
    """
    verify if query from fasta file matches query sequence described in BLAST.
    :param blast:
    :param query:
    :return:
    """
    if not (query.id == blast.query or query.description == blast.query):
        warnings.warn(
            'Provided query id ({}) do not match query id in BLAST output ({})'.format(
                query.id,
                blast.query
            )
        )
    if len(query) != blast.query_length:
        warnings.warn(
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
    try:
        a = check_output(
            [
                '{}blastdbcmd'.format(CONFIG.blast_path),
                '-version'
            ]
        )
        a = a.decode(encoding='utf-8')
        r = re.search('(?<=blastdbcmd: )[1-9.]+', a)
        if r:
            ver = r.group().split('.')
            ver = [int(i) for i in ver]
            for v, minv in zip(ver, minimal_version):
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_locarna(minimal_version):
    msgversion = 'Locarna is not installed in required version, required version is {}.{}.{}'.format(*locarna_minimal_version)
    msgpath = '{}LocARNA could not be located (not in PATH).'.format(CONFIG.locarna_path)
    try:
        a = check_output(
            [
                '{}locarna'.format(CONFIG.locarna_path),
                '--version'
            ]
        )
        a = a.decode()
        if a.startswith('LocARNA'):
            r = re.finditer('[1-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_infernal(program, minimal_version):
    msgversion = '{} (part of INFERNAL) is not installed in required version, required is {}.{}.{}'.format(program, *minimal_version)
    msgpath = '{}{} could not be located (not in PATH).'.format(CONFIG.infernal_path, program)
    try:
        a = check_output(
            [
                CONFIG.infernal_path + program,
                '-h'
            ]
        )
        a = a.decode()
        if a.startswith('# {}'.format(program)):
            r = re.search('(?<=# INFERNAL )[1-9.]+', a)
            ver = r.group().split('.')
            ver = [int(i) for i in ver]
            for v, minv in zip(ver, minimal_version):
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        raise


def verify_viennarna_program(program, minimal_version):
    msgversion = '{} is not installed in required version, required version is {}.{}.{}'.format(program, *minimal_version)
    msgpath = '{}{} could not be located'.format(CONFIG.viennarna_path, program)
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
                    return True
                elif int(v) < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
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
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_clustalo(minimal_version):
    msgversion = 'clustalo is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}clustalo could not be located (not in PATH).'.format(CONFIG.clustal_path)
    try:
        a = check_output(
            [
                '{}clustalo'.format(CONFIG.clustal_path),
                '--version'
            ]
        )
        a = a.decode()
        if a:
            r = re.finditer('[1-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_muscle(minimal_version):
    msgversion = 'muscle is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}muscle could not be located (not in PATH).'.format(CONFIG.muscle_path)
    try:
        a = check_output(
            [
                '{}muscle'.format(CONFIG.muscle_path),
                '-version'
            ]
        )
        a = a.decode()
        if a.startswith('MUSCLE'):
            r = re.finditer('[1-9]+', a)
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_tcoffee(minimal_version):
    # do not check the revision
    msgversion = 't_coffee is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}t_coffee could not be located (not in PATH).'.format(CONFIG.tcoffee_path)
    try:
        a = check_output(
            [
                '{}t_coffee'.format(CONFIG.tcoffee_path),
                '-version'
            ]
        )
        a = a.decode()
        if a.startswith('PROGRAM: T-COFFEE'):
            im = re.search('(?<=version_)[0-9.]+', a, flags=re.IGNORECASE)
            if im:
                b = im.group()
                r = re.finditer('[0-9]+', b)
                for match, minv in zip(r, minimal_version):
                    v = int(match.group())
                    if v > minv:
                        return True
                    elif v < minv:
                        print(msgversion)
                        return False
                return True
            else:
                print(msgversion)
                return False
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_rcoffee(minimal_version):
    msgversion = 't_coffee -mode rcoffee is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}t_coffee -mode rcoffee could not be located (not in PATH).'.format(CONFIG.muscle_path)

    if not verify_tcoffee(minimal_version):
        print(
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
        t2 = re.search('-- COM: t_coffee -mode rcoffee', a)

        n1 = re.search("ERROR: special_mode rcoffee is unknown \[FATAL:T-COFFEE\]", a)

        if t1 and t2 and not n1:
            return True
        else:
            print(msgversion)
            return False

    except FileNotFoundError:
        print(msgpath)
        return False


def verify_centroid_homfold(minimal_version):
    msgversion = 'centroid_homfold is not installed in required version, required version is {}.{}.{}'.format(*minimal_version)
    msgpath = '{}centroid_homfold could not be located (not in PATH).'.format(CONFIG.centriod_path)
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
            r = re.finditer('[1-9]+', b[1])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_turbofold(minimal_version):
    msgversion = 'TurboFold is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}TurboFold could not be located (not in PATH).'.format(CONFIG.turbofold_path)
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
            r = re.finditer('[1-9]+', b[2])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
        return False


def verify_mfold(minimal_version):
    msgversion = 'mfold is not installed in required version, required version is {}.{}'.format(*minimal_version)
    msgpath = '{}mfold could not be located (not in PATH).'.format(CONFIG.mfold_path)
    try:
        a = check_output(
            [
                '{}mfold'.format(CONFIG.mfold_path),
                '-v'
            ]
        )
        a = a.decode()
        b = a.split()
        if b[0] == 'mfold':
            r = re.finditer('[1-9]+', b[2])
            for match, minv in zip(r, minimal_version):
                v = int(match.group())
                if v > minv:
                    return True
                elif v < minv:
                    print(msgversion)
                    return False
            return True
        else:
            print(msgversion)
            return False
    except FileNotFoundError:
        print(msgpath)
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
        if needed:
            raise EnvironmentError('Missing {} (needed for {}).'.format(' '.join(needed), met))

if __name__ == '__main__':
    print(check_3rd_party_tools())
    print(check_3rd_party_prediction_tools())
