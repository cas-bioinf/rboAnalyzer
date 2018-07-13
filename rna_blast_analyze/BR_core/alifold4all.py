import os
from subprocess import call
from tempfile import mkstemp

from rna_blast_analyze.BR_core.config import CONFIG


def compute_clustalo_clasic(file, clustalo_params=''):

    fd, out_path = mkstemp()
    os.close(fd)
    FNULL = open(os.devnull, 'w')

    print('compute_clustalo_clasic infile {} outfile {}'.format(file, out_path))

    try:
        k = call(
            '{}clustalo {} -i {} > {}'.format(
                CONFIG.clustal_path,
                clustalo_params,
                file,
                out_path
            ),
            shell=True,
            stdout=FNULL
        )
        if k:
            raise ChildProcessError('call to clustalo failed for files: in: {}, out: {}'.format(file, out_path))
    finally:
        FNULL.close()

    return out_path


def compute_alifold(msa_file, alifold_params=''):

    fd, out_path = mkstemp()
    os.close(fd)

    k = call(
        '{}RNAalifold --noPS -f C {} < {} > {}'.format(
            CONFIG.viennarna_path,
            alifold_params,
            msa_file,
            out_path
        ),
        shell=True
    )
    if k:
        raise ChildProcessError('RNAalifold failed for files: in: {}, out {}'.format(msa_file, out_path))

    return out_path


def compute_refold(alig_file, cons_file):
    """
    runs refold program
    :param alig_file: MSA alignment file in clustal format
    :param cons_file: file with consensus structure in alifold format
    :return:
    """
    fd, out_path = mkstemp()
    os.close(fd)
    k = call(
        '{}refold.pl {} {} > {}'.format(
            CONFIG.refold_path,
            alig_file,
            cons_file,
            out_path
        ),
        shell=True
    )
    if k:
        raise ChildProcessError('refold.pl failed for files in: {}, {}, out: {}'.format(alig_file, cons_file, out_path))

    return out_path


def compute_constrained_prediction(constrained_file):
    # constrained prediction je funkcni
    fd, outfile = mkstemp()
    os.close(fd)
    r = call(
        '{}RNAfold -C --noPS < {} > {}'.format(
            CONFIG.viennarna_path,
            constrained_file,
            outfile
        ),
        shell=True
    )
    if r:
        raise ChildProcessError('call to RNAfold -C failed')
    return outfile
