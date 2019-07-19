import logging
import os
from subprocess import call
from tempfile import mkstemp
import shlex

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger('rboAnalyzer')


def compute_clustalo_clasic(file, clustalo_params=''):
    ml.info('Running clustalo.')
    ml.debug(fname())

    fd, out_path = mkstemp(prefix='rba_', suffix='_01', dir=CONFIG.tmpdir)
    os.close(fd)
    FNULL = open(os.devnull, 'w')

    try:
        parlist = clustalo_params.split()
        cmd = [
            '{}clustalo'.format(CONFIG.clustal_path),
            ] + parlist + [
            '-i', file,
            '-o', out_path
        ]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd)
        else:
            r = call(cmd, stdout=FNULL)

        if r:
            msgfail = 'call to clustalo failed for files: in: {}, out: {}'.format(file, out_path)
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)
    finally:
        FNULL.close()

    return out_path


def compute_alifold(msa_file, alifold_params=''):
    ml.info('Running RNAalifold.')
    ml.debug(fname())
    fd, out_path = mkstemp(prefix='rba_', suffix='_02', dir=CONFIG.tmpdir)
    os.close(fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = '{} --noPS -f C {} < {} > {}'.format(
            shlex.quote('{}RNAalifold'.format(CONFIG.viennarna_path)),
            ' '.join([shlex.quote(i) for i in shlex.split(alifold_params)]),
            shlex.quote(msa_file),
            shlex.quote(out_path)
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd, shell=True)
        else:
            r = call(cmd, shell=True, stderr=FNULL)

        if r:
            msgfail = 'RNAalifold failed for files: in: {}, out {}'.format(msa_file, out_path)
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)

        return out_path


def compute_refold(alig_file, cons_file):
    """
    runs refold program
    :param alig_file: MSA alignment file in clustal format
    :param cons_file: file with consensus structure in alifold format
    :return:
    """
    ml.debug(fname())
    fd, out_path = mkstemp(prefix='rba_', suffix='_03', dir=CONFIG.tmpdir)
    os.close(fd)
    cmd = '{} {} {} > {}'.format(
        shlex.quote('{}refold.pl'.format(CONFIG.refold_path)),
        shlex.quote(alig_file),
        shlex.quote(cons_file),
        shlex.quote(out_path)
    )
    ml.debug(cmd)
    k = call(cmd, shell=True)
    if k:
        msgfail = 'refold.pl failed for files in: {}, {}, out: {}'.format(alig_file, cons_file, out_path)
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)

    return out_path
