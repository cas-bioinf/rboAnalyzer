import os
import logging
from subprocess import call
from tempfile import mkstemp

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.decorators import timeit_decorator
from rna_blast_analyze.BR_core import BA_support
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.exceptions import AmbiguousQuerySequenceException

ml = logging.getLogger('rboAnalyzer')


def run_centroid_homfold(fasta2predict, fasta_homologous_seqs, centroid_homfold_params='', outfile=None):
    if outfile:
        ch_outfile = outfile
    else:
        ch, ch_outfile = mkstemp(prefix='rba_', suffix='_09', dir=CONFIG.tmpdir)
        os.close(ch)

    # build commandline
    cmd = [
        '{}centroid_homfold'.format(CONFIG.centriod_path),
        '-H', fasta_homologous_seqs,
    ]
    if centroid_homfold_params != '':
        cmd += centroid_homfold_params.split()
    cmd += [
        '-o', ch_outfile,
        fasta2predict
    ]
    r = call(cmd)
    if r:
        raise ChildProcessError(
            'call to run centroid_homfold failed for files'
            ' in:{} homologs:{} outfile:{}'.format(
                fasta2predict,
                fasta_homologous_seqs,
                ch_outfile
            )
        )
    return ch_outfile


@timeit_decorator
def me_centroid_homfold(fasta2predict, fasta_homologous_seqs, params=None):
    """
    run centroid_homefold several times and vary -g parameter, to predict the best possible structure
    :param fasta2predict:
    :param fasta_homologous_seqs:
    :param params:
    :return:
    """

    # first run centroid homefold for several stages of g (-1)
    # find the most stable structure value of g
    # structure of output

    if params is None:
        params = dict()

    ch_params = ''
    if params and ('centroid_homfold' in params) and params['centroid_homfold']:
        ch_params += params['centroid_homfold']

    if '-g ' in ch_params and '-g -1' not in ch_params or '-t ' in params:
        print('We only allow to run centroid homfold wit automatic mode where the structure is predicted with multiple'
              ' weights and then best scoring is selected, thresholding is also forbiden as it implies -g.')
        raise AttributeError('Error centroid homfold not permited to run with "-g" or "-t".')
    ch_params += ' -g -1'

    first_structures = run_centroid_homfold(fasta2predict, fasta_homologous_seqs, centroid_homfold_params=ch_params)
    structures2return = [ch_struc for ch_struc in centroid_homfold_select_best(first_structures)]
    os.remove(first_structures)
    return structures2return


@timeit_decorator
def centroid_homfold_fast(all_seqs, query, all_seqs_fasta, n, centroid_homfold_params, len_diff):
    ml.debug(fname())

    homologous_file = centroid_homfold_fast_prep(all_seqs, query, n, len_diff)

    structures, _ = me_centroid_homfold(all_seqs_fasta, homologous_file, params=centroid_homfold_params)
    os.remove(homologous_file)
    return structures


def centroid_homfold_fast_prep(all_seqs, query, n, len_diff):
    ml.debug(fname())

    assert n >= 1

    if query.annotations['ambiguous']:
        msgfail = "Query sequence contains ambiguous characters. Can't use centroid-fast."
        ml.error(msgfail)
        raise AmbiguousQuerySequenceException(msgfail)

    nr_na_ld = BA_support.sel_seq_simple(all_seqs, query, len_diff)
    nr_na_ld_n = nr_na_ld[:int(n)]

    ch, ch_hom = mkstemp(prefix='rba_', suffix='_74', dir=CONFIG.tmpdir)
    with os.fdopen(ch, 'w') as h:
        BA_support.write_fasta_from_list_of_seqrecords(h, nr_na_ld_n)

    return ch_hom


def centroid_homfold_select_best(first_structures):
    for cen_hom_proposed_structures in _parse_centroid_homefold_output_file(first_structures):
        best_structure_by_e = dict()
        for key in cen_hom_proposed_structures.annotations['sss']:
            cen_pred_params = cen_hom_proposed_structures.annotations[key].strip('()').split('=')
            if len(cen_pred_params) != 4:
                raise ValueError('unexpected number of centroid homfold prediction params')
            best_structure_by_e[round(float(cen_pred_params[-1]), 2)] = key

        best_structure_key = best_structure_by_e[min(best_structure_by_e.keys())]
        best_structure = SeqRecord(cen_hom_proposed_structures.seq,
                                   id=cen_hom_proposed_structures.id)
        # rename selected (best) structure to ss0
        best_structure.letter_annotations['ss0'] = cen_hom_proposed_structures.letter_annotations[best_structure_key]
        best_structure.annotations['sss'] = []
        best_structure.annotations['sss'].append('ss0')
        best_structure.annotations['ss0'] = cen_hom_proposed_structures.annotations[best_structure_key]
        yield best_structure


def _parse_centroid_homefold_output_file(file):
    with open(file, 'r') as f:
        for sr in BA_support.parse_one_rec_in_multiline_structure(f):
            cf = sr.strip().splitlines()

            cfr = SeqRecord(Seq(cf[1]), id=cf[0])
            cfr.annotations['sss'] = []
            for i, ll in enumerate(cf[2:]):
                structure, ann = ll.split()
                cfr.letter_annotations['ss' + str(i)] = structure
                cfr.annotations['ss' + str(i)] = ann
                cfr.annotations['sss'].append('ss' + str(i))

            yield cfr
