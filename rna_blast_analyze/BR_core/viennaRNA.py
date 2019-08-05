import os
import shlex
import logging
import shutil
from subprocess import call
from tempfile import TemporaryFile, mkstemp, gettempdir

from rna_blast_analyze.BR_core import exceptions
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.BA_support import remove_one_file_with_try

ml = logging.getLogger('rboAnalyzer')


def rnafold(fastafile, outfile, parameters=''):
    """Run RNAfold on commandline

    issue --noPS to prevent plotting of postscript structure
    """
    if '--noPS' not in parameters:
        parameters += ' --noPS'
    cmd = '{} {} < {} > {}'.format(
        shlex.quote('{}RNAfold'.format(CONFIG.viennarna_path)),
        ' '.join([shlex.quote(i) for i in parameters.split()]),
        shlex.quote(fastafile),
        shlex.quote(outfile)
    )
    ml.debug(cmd)
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = call(cmd, shell=True, stderr=tmp, stdout=tmp)

        if r:
            msgfail = 'Call to RNAfold failed, please check if RNAfold is in path.'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.RNAfoldException(msgfail, tmp.read())
        return outfile


def run_rnaplot(seq, structure=None, format='svg', outfile=None):
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

    allowed_formats = {'ps', 'svg', 'gml', 'xrna'}
    if format not in allowed_formats:
        raise TypeError('Format can be only from {}.'.format(allowed_formats))

    fd, tmpfile = mkstemp(prefix='rba_', suffix='_08', dir=CONFIG.tmpdir)

    rnaname = tmpfile.split('/')[-1].split('\\')[-1]

    currdirr = os.getcwd()
    tmpdir = gettempdir()
    os.chdir(tmpdir)

    with os.fdopen(fd, 'w') as fh:
        fh.write('>{}\n{}\n{}\n'.format(
            rnaname,
            sequence,
            structure
        ))

    cmd = '{} --output-format={} < {}'.format(
        shlex.quote('{}RNAplot'.format(CONFIG.viennarna_path)),
        shlex.quote(format),
        shlex.quote(tmpfile)
    )
    ml.debug(cmd)
    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = call(cmd, shell=True, stdout=tmp, stderr=tmp)
        if r:
            msgfail = 'Call to RNAplot failed.'
            ml.error(msgfail)
            os.chdir(currdirr)
            tmp.seek(0)
            details = tmp.read()
            ml.debug(details)
            raise exceptions.RNAplotException(msgfail, details)

        # output file is name of the sequence (in this case name of the file) + "_ss." + chosen format
        remove_one_file_with_try(tmpfile)

        plot_output_file = os.path.join(tmpdir, rnaname + '_ss.' + format)
        os.chdir(currdirr)
        if outfile is None:
            return plot_output_file
        else:
            shutil.move(os.path.join(tmpdir, plot_output_file), outfile)
            return outfile