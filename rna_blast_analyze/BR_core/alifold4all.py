import os
import logging
from subprocess import call
from tempfile import mkstemp

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)


def compute_clustalo_clasic(file, clustalo_params=''):
    ml.info('Running clustalo.')
    ml.debug(fname())

    fd, out_path = mkstemp(prefix='rba_', suffix='_01')
    os.close(fd)
    FNULL = open(os.devnull, 'w')

    try:
        cmd = '{}clustalo {} -i {} > {}'.format(
            CONFIG.clustal_path,
            clustalo_params,
            file,
            out_path
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd, shell=True)
        else:
            r = call(cmd, shell=True, stdout=FNULL)

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
    fd, out_path = mkstemp(prefix='rba_', suffix='_02')
    os.close(fd)

    with open(os.devnull, 'w') as FNULL:
        cmd = '{}RNAalifold --noPS -f C {} < {} > {}'.format(
            CONFIG.viennarna_path,
            alifold_params,
            msa_file,
            out_path
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
    fd, out_path = mkstemp(prefix='rba_', suffix='_03')
    os.close(fd)
    cmd = '{}refold.pl {} {} > {}'.format(
        CONFIG.refold_path,
        alig_file,
        cons_file,
        out_path
    )
    ml.debug(cmd)
    k = call(cmd, shell=True)
    if k:
        msgfail = 'refold.pl failed for files in: {}, {}, out: {}'.format(alig_file, cons_file, out_path)
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)

    return out_path


def compute_constrained_prediction(constrained_file):
    ml.debug(fname())
    fd, outfile = mkstemp(prefix='rba_', suffix='_04')
    os.close(fd)
    cmd = '{}RNAfold -C --noPS < {} > {}'.format(
        CONFIG.viennarna_path,
        constrained_file,
        outfile
    )
    ml.debug(cmd)
    r = call(cmd, shell=True)
    if r:
        msgfail = 'call to RNAfold -C failed'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)
    return outfile
