import os
import logging
import multiprocessing
from subprocess import call
from tempfile import mkstemp, TemporaryFile

from Bio import SeqIO

from rna_blast_analyze.BR_core.BA_support import rebuild_structures_output_from_pred, ct2db, remove_one_file_with_try, remove_files_with_try
from rna_blast_analyze.BR_core import exceptions
from rna_blast_analyze.BR_core.config import CONFIG

ml = logging.getLogger('rboAnalyzer')


def _run_hybrid_ss_min_wrapper(seq, P, W, M):
    fd, tmp_fasta = mkstemp(prefix='rba_', suffix='_06', dir=CONFIG.tmpdir)
    try:
        with os.fdopen(fd, 'w') as fid:
            fid.write('>{}\n{}\n'.format(seq.id, str(seq.seq)))

        predicted_ss = _run_hybrid_ss_min_single(tmp_fasta, P, W, M)
        return predicted_ss[0]

    except exceptions.HybridssminException as e:
        return None

    finally:
        remove_one_file_with_try(tmp_fasta)


def _run_hybrid_ss_min_single(file_path, P, W, M):
    with TemporaryFile(mode='w', encoding='utf-8') as tmp:
        cmd3 = [
            '{}hybrid-ss-min'.format(CONFIG.mfold_path),
            '--suffix=DAT',
            '--NA=RNA',
            '--noisolate',
            '--mfold=' + str(P) + ',' + str(W) + ',' + str(M),
            file_path
        ]
        ml.debug(cmd3)

        rt = call(
            cmd3,
            cwd=''.join(os.path.join(os.path.split(file_path)[:-1])),
            stdout=tmp,
            stderr=tmp
        )

        if rt:
            msgfail = 'Execution of hybrid-ss-min failed'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.HybridssminException(msgfail, tmp.read())

        if not os.path.isfile(file_path + '.ct'):
            msgfail = 'Execution of hybrid-ss-min failed - not output file'
            ml.error(msgfail)
            tmp.seek(0)
            raise exceptions.HybridssminException(msgfail, tmp.read())

        with open(file_path + '.ct', 'r') as sout:
            pred_structures = ct2db(sout)

        remove_files_with_try([
            file_path + '.run',
            file_path + '.plot',
            file_path + '.dG',
            file_path + '.ct',
            file_path + '.ann',
        ])
        return pred_structures


def run_hybrid_ss_min(in_path, mfold=(10, 2, 20), threads=1):
    """
    forbid lonely pairs when calling hybrid-ss-min, rnashapes cannot work with them

    # beware of

    :param in_path: path to file to compute suboptimal structures for
    :param mfold: mfold parameters (P, W, M)
    :returns: list of SeqRecord objects wit suboptimal structures as letter_annotations in them
    """
    ml.info('Runing hybrid-ss-min')
    assert isinstance(mfold, (tuple, list)) and 3 == len(mfold)

    P = mfold[0]
    W = mfold[1]
    M = mfold[2]

    if threads == 1:
        try:
            suboptimals = _run_hybrid_ss_min_single(in_path, P, W, M)
        except exceptions.HybridssminException as e:
            suboptimals = []
            for seq in SeqIO.parse(in_path, format='fasta'):
                suboptimals.append(
                    _run_hybrid_ss_min_wrapper(seq, P, W, M)
                )
    else:

        tuples = [(seq, P, W, M) for seq in SeqIO.parse(in_path, format='fasta')]

        with multiprocessing.Pool(processes=threads) as pool:
            suboptimals = pool.starmap(_run_hybrid_ss_min_wrapper, tuples)

    with open(in_path, 'r') as fin:
        inseqs = [s for s in SeqIO.parse(fin, format='fasta')]
        suboptimals = rebuild_structures_output_from_pred(inseqs, suboptimals, method='hybrid-ss-min')

        return suboptimals