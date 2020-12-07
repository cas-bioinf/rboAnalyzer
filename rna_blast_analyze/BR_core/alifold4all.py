import logging
import os
import subprocess
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
        r = subprocess.call(cmd, stdout=tmp, stderr=tmp)

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

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp, os.fdopen(fd, 'w') as output, open(msa_file, 'r') as inp:
        cmd = [
            '{}RNAalifold'.format(CONFIG.viennarna_path),
            '--noPS',
            '-f', 'C',
        ] + shlex.split(alifold_params)
        ml.debug(cmd)

        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=output,
            stderr=tmp,
            universal_newlines=True
        )
        p.communicate(input=inp.read())

        if p.returncode:
            msgfail = 'RNAalifold failed.'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.RNAalifoldException(msgfail, tmp.read())

        return out_path


def compute_refold(alig_file, cons_file, timeout=None):
    """
    runs refold program
    :param alig_file: MSA alignment file in clustal format
    :param cons_file: file with consensus structure in alifold format
    :return:
    """
    ml.debug(fname())
    fd, out_path = mkstemp(prefix='rba_', suffix='_03', dir=CONFIG.tmpdir)
    cmd = [
        '{}refold.pl'.format(CONFIG.refold_path),
        alig_file,
        cons_file
    ]
    ml.debug(cmd)
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp, os.fdopen(fd, 'w') as output:
        with subprocess.Popen(cmd, stdout=output, stderr=tmp) as p:
            try:
                p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()
                raise
            if p.returncode:
                msgfail = 'Call to refold.pl failed.'
                ml.error(msgfail)
                tmp.seek(0)
                raise exceptions.RefoldException(msgfail, tmp.read())

        return out_path
