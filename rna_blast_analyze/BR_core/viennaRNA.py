import os
import shlex
import logging
import shutil
import subprocess
from subprocess import call, check_output
from tempfile import TemporaryFile, gettempdir
from multiprocessing import Pool

from rna_blast_analyze.BR_core import exceptions
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.BA_support import generate_random_name

ml = logging.getLogger('rboAnalyzer')


def rnafold_fasta(fastafile, outfile, parameters='', timeout=None):
    """Run RNAfold on commandline

    issue --noPS to prevent plotting of postscript structure
    """
    if '--noPS' not in parameters:
        parameters += ' --noPS'

    cmd = ['{}RNAfold'.format(CONFIG.viennarna_path),] + shlex.split(parameters)

    ml.debug(cmd)

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp, open(fastafile, 'r') as inp, open(outfile, 'w') as output:
        with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=output, stderr=tmp, universal_newlines=True) as p:
            try:
                p.communicate(input=inp.read(), timeout=timeout)
            except subprocess.TimeoutExpired:
                p.kill()
                p.communicate()
                raise

            if p.returncode:
                msgfail = 'Call to RNAfold failed, please check if RNAfold is in path.'
                ml.error(msgfail)
                tmp.seek(0)
                raise exceptions.RNAfoldException(msgfail, tmp.read())
        return outfile


def run_rnaplot(seq, structure=None, format='svg', outfile=None, timeout=None):
    """
    run rnaplot in desired format
    if seq
    :param seq:
    :param structure:
    :return:
    """
    ml.debug(fname())
    if structure is None:
        sequence = str(seq.seq)
        structure = seq.letter_annotations['ss0']
    else:
        sequence = seq

    assert len(sequence) == len(structure)

    allowed_formats = {'ps', 'svg', 'gml'}
    if format not in allowed_formats:
        raise TypeError('Format can be only from {}.'.format(allowed_formats))

    currdirr = os.getcwd()
    tmpdir = gettempdir()
    os.chdir(tmpdir)

    cmd = ['{}RNAplot'.format(CONFIG.viennarna_path), '--output-format={}'.format(format)]
    ml.debug(cmd)

    rnaname = generate_random_name(10)

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        with subprocess.Popen(cmd, universal_newlines=True, stdin=subprocess.PIPE, stdout=tmp, stderr=tmp) as p:
            try:
                p.communicate(input='>{}\n{}\n{}\n'.format(
                    rnaname,
                    sequence,
                    structure
                ), timeout=timeout)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()
                raise

            if p.returncode:
                msgfail = 'Call to RNAplot failed.'
                ml.error(msgfail)
                os.chdir(currdirr)
                tmp.seek(0)
                details = tmp.read()
                ml.debug(details)
                raise exceptions.RNAplotException(msgfail, details)

        plot_output_file = os.path.join(tmpdir, rnaname + '_ss.' + format)
        os.chdir(currdirr)
        if outfile is None:
            return plot_output_file
        else:
            shutil.move(os.path.join(tmpdir, plot_output_file), outfile)
            return outfile


def RNAfold(sequence, timeout=None):
    ml.debug(fname())
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = check_output(
            [
                '{}RNAfold'.format(CONFIG.viennarna_path),
                '--noPS',
            ],
            input=sequence.encode(),
            stderr=tmp,
            timeout=timeout
        )

        if isinstance(r, Exception):
            msgfail = 'RNAfold failed.'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.RNAfoldException(msgfail, tmp.read())

        # more robust decode
        out_str = r.decode()
        spl = out_str.split('\n')
        seq = spl[0]
        structure = spl[1][:len(seq)]
        energy = float(spl[1][len(seq)+2:-1])
        # seq, structure, energy = r.decode().split()
        # return seq, structure, float(energy[1:-1])
        return structure, energy


def wrap_RNAfold(sequences, threads=1):
    if threads == 1:
        return [RNAfold(seq) for seq in sequences]
    else:
        with Pool(processes=threads) as pool:
            return pool.map(RNAfold, sequences)
