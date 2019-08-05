import logging
import os
from subprocess import call
from tempfile import mkstemp, TemporaryFile
import shlex

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core import exceptions

ml = logging.getLogger('rboAnalyzer')


def compute_clustalo_clasic(file, clustalo_params=''):
    ml.info('Running clustalo.')
    ml.debug(fname())

    fd, out_path = mkstemp(prefix='rba_', suffix='_01', dir=CONFIG.tmpdir)
    os.close(fd)

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        parlist = clustalo_params.split()
        cmd = [
            '{}clustalo'.format(CONFIG.clustal_path),
            ] + parlist + [
            '-i', file,
            '-o', out_path
        ]
        ml.debug(cmd)
        r = call(cmd, stdout=tmp, stderr=tmp)

        if r:
            msgfail = 'Call to clustalo failed.'
            tmp.seek(0)
            ml.error(msgfail)
            ml.error(cmd)
            raise exceptions.ClustaloException(msgfail, tmp.read())

    return out_path


def compute_alifold(msa_file, alifold_params=''):
    ml.info('Running RNAalifold.')
    ml.debug(fname())
    fd, out_path = mkstemp(prefix='rba_', suffix='_02', dir=CONFIG.tmpdir)
    os.close(fd)

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        cmd = '{} --noPS -f C {} < {} > {}'.format(
            shlex.quote('{}RNAalifold'.format(CONFIG.viennarna_path)),
            ' '.join([shlex.quote(i) for i in shlex.split(alifold_params)]),
            shlex.quote(msa_file),
            shlex.quote(out_path)
        )
        ml.debug(cmd)
        r = call(cmd, shell=True, stderr=tmp, stdout=tmp)

        if r:
            msgfail = 'RNAalifold failed.'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.RNAalifoldException(msgfail, tmp.read())

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
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = call(cmd, shell=True, stdout=tmp, stderr=tmp)
        if r:
            msgfail = 'Call to refold.pl failed.'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.RefoldException(msgfail, tmp.read())

        return out_path
