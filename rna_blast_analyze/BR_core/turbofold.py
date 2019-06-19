import logging
import multiprocessing
import os
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp

from Bio.SeqRecord import SeqRecord

import rna_blast_analyze.BR_core.exceptions
from rna_blast_analyze.BR_core import BA_support
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.decorators import timeit_decorator
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)


@timeit_decorator
def turbofold_fast(all_seqs, seqs2predict, query, cpu, n, turbofold_params, len_diff):
    # ambiguos sequence cannot be predicted with turbofold
    # do not predict sequence twice
    # seqs2predict must be subset of all_seqs
    # todo add test for this

    ml.debug(fname())

    if query.annotations['ambiguous']:
        msgfail = "Query sequence contains ambiguous characters. Can't use TurboFold_fast."
        ml.error(msgfail)
        raise ValueError(msgfail)

    nr_na_ld = BA_support.sel_seq_simple(all_seqs, query, len_diff)

    return turbofold_ext_nr_fast(seqs2predict, nr_na_ld, cpu, n, turbofold_params)


def turbofold_ext_nr_fast(all_seqs, nrset, cpu, n, turbofold_params):
    # ambiguos sequence cannot be predicted with turbofold
    # do not predict sequence twice
    ml.debug(fname())

    msg_short_list = 'Turbofold: Number of sequences is less then required.'

    if not len(nrset) == len(BA_support.filter_ambiguous_seqs_from_list(nrset)) == len(BA_support.non_redundant_seqs(nrset)):
        msgfail = 'Wrong nr set specification.'
        if len(nrset) != len(BA_support.filter_ambiguous_seqs_from_list(nrset)):
            msgfail += 'nr set contains seq(s) with ambiguous character.'
        if len(nrset) != len(BA_support.non_redundant_seqs(nrset)):
            msgfail += 'nr set contain non redundant sequence(s).'
        ml.error(msgfail)
        raise AssertionError(msgfail)

    list2predict = []
    for seq in all_seqs:
        if seq.annotations.get('ambiguous', False):
            ml.warning('Skipping TurboFold prediction for {} (ambiguous base)'.format(seq.id))
            continue
        seq_set = _prepare_set_n(seq, nrset, n)

        if len(seq_set) < 2:
            ml.error("Turbofold can't be used with less then 2 sequences.")
            raise rna_blast_analyze.BR_core.exceptions.NoHomologousSequenceException

        if len(seq_set) < n:
            ml.warning(msg_short_list)
            seq.annotations['msgs'].append(msg_short_list)
        list2predict.append((seq_set, turbofold_params, seq.id))

    if cpu == 1:
        pred_list = []
        for oneseqset, tpar, _ in list2predict:
            pred_list.append(run_turbofold(oneseqset, tpar))
    else:
        pool = multiprocessing.Pool(processes=cpu)
        pred_list = pool.map(_rt_wrapper, list2predict)
        pool.close()

    return [[o for o in out if o.id == l_in[2]][0] for out, l_in in zip(pred_list, list2predict)]


def _prepare_set_n(seq, nr_seqs, n):
    if len(nr_seqs) - n + 1 < 0:
        return BA_support.non_redundant_seqs([seq] + nr_seqs)

    i = 0
    while i != len(nr_seqs) - n + 1:
        seqs = [seq] + nr_seqs[:n - 1 + i]
        if len(BA_support.non_redundant_seqs(seqs)) == n:
            return seqs
        i += 1

    return BA_support.non_redundant_seqs([seq] + nr_seqs)


def _rt_wrapper(pack):
    return run_turbofold(pack[0], pack[1])


def write_turbofold_confile(input_sequences, turbofold_params=None, cpus=None, outdir=None):
    ml.debug(fname())
    if outdir:
        tmpdir = outdir
    else:
        tmpdir = mkdtemp(prefix='rba_')

    params2write = {
        'Mode': 'MEA',
        'OutAln': os.path.join(tmpdir, 'Output.aln')
    }
    if turbofold_params is not None:
        params2write.update(turbofold_params)

    # write sequences by one in fasta file format
    seq_paths = []
    for i, hseq in enumerate(input_sequences):
        fastaname = os.path.join(tmpdir, 'hseq{}.fasta'.format(i))
        with open(fastaname, 'w') as f:
            f.write('>{}\n{}\n'.format(
                'hseq{}'.format(i),
                str(hseq.seq)
            ))
        seq_paths.append(fastaname)

    hom_out_paths = [os.path.join(tmpdir, 'hseq{}.ct'.format(i)) for i in range(len(input_sequences))]

    conf_file = os.path.join(tmpdir, 'conf_file.conf')
    with open(conf_file, 'w') as c:
        for key in params2write.keys():
            c.write('{} = {}\n'.format(key, params2write[key]))

        if cpus and isinstance(cpus, int):
            c.write('Processors = {}'.format(cpus))

        in_files = '{' + ';'.join(seq_paths) + '}'
        out_files = '{' + ';'.join(hom_out_paths) + '}'

        c.write('InSeq = {}\n'.format(in_files))
        c.write('OutCT = {}\n'.format(out_files))

    return tmpdir, conf_file, hom_out_paths


def run_turbofold(sequences, params):
    ml.info('Running Turbofold.')
    ml.debug(fname())

    env = os.environ.copy()
    if 'DATAPATH' not in env:
        env['DATAPATH'] = CONFIG.rnastructure_datapath

    tmpdir, con_file, output_structure_files = write_turbofold_confile(
        input_sequences=sequences,
        turbofold_params=params,
    )

    # run without prediction progress reporting output
    with open(os.devnull, 'w') as FNULL:
        cmd = [
            '{}TurboFold'.format(CONFIG.turbofold_path),
            con_file
        ]
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = call(cmd, env=env)
        else:
            r = call(cmd, stdout=FNULL, env=env)

        if r:
            msgfail = 'Call to turbofold failed, cmd below:'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)

        # now convert ct files produced by TurboFold
        new_structures = []
        for out_str_file, orig_seq in zip(output_structure_files, sequences):
            o = open(out_str_file, 'r')
            seq_with_pred_str = BA_support.ct2db(o, energy_txt='ENERGY')
            o.close()

            assert str(seq_with_pred_str[0].seq) == str(orig_seq.seq)

            new_structures.append(
                SeqRecord(
                    orig_seq.seq,
                    id=orig_seq.id,
                    annotations={'sss': ['ss0']},
                    letter_annotations={'ss0': seq_with_pred_str[0].letter_annotations['ss0']}
                )
            )

        rmtree(tmpdir)
        return new_structures


@timeit_decorator
def turbofold_with_homologous(all_sequences, nr_homologous, params, n, cpu):
    """
    Trubofold mode is MEA by default
    :param all_sequences:
    :param nr_homologous:
    :param params: dict with key: value structure for params to be used in turbofold
    :return:
    """
    ml.debug(fname())

    msg_short_list = 'Turbofold: Number of sequences is less then required.'

    nr_homologous_set = {str(seq.seq) for seq in nr_homologous}

    if len(nr_homologous_set) < 1:
        raise rna_blast_analyze.BR_core.exceptions.NoHomologousSequenceException

    list2predict = []
    for seq in all_sequences:
        seq_set = _prepare_set_n(seq, nr_homologous, n)

        if len(seq_set) < 2:
            ml.error("Turbofold can't be used with less then 2 sequences.")
            raise rna_blast_analyze.BR_core.exceptions.NoHomologousSequenceException

        if len(seq_set) < n:
            ml.warning(msg_short_list)
            seq.annotations['msgs'].append(msg_short_list)
        list2predict.append((seq_set, params, seq.id))

    if cpu == 1:
        pred_list = []
        for oneseqset, tpar, _ in list2predict:
            pred_list.append(run_turbofold(oneseqset, tpar))
    else:
        pool = multiprocessing.Pool(processes=cpu)
        pred_list = pool.map(_rt_wrapper, list2predict)
        pool.close()

    return [[o for o in out if o.id == l_in[2]][0] for out, l_in in zip(pred_list, list2predict)]
