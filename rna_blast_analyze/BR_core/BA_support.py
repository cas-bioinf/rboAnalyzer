import datetime
import locale
import os
import re
import shutil
from io import StringIO
from random import choice, shuffle
from subprocess import call, STDOUT, check_output, CalledProcessError
from tempfile import mkstemp, gettempdir
from warnings import warn

import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.blast_bio.parser_to_bio_blast import blast_parse_txt as blast_minimal_parser
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.parse_accession import accession_regex

# idiotic matplotlib
locale.setlocale(locale.LC_ALL, 'C')


def write_fasta_from_list_of_seqrecords(f, seq_list):
    for seq in seq_list:
        f.write('>{}\n{}\n'.format(seq.id,
                                   str(seq.seq)))


def run_hybrid_ss_min(in_path, mfold=(10, 2, 20)):
    """
    forbid lonely pairs when calling hybrid-ss-min, rnashapes cannot work with them

    # beware of

    :param in_path: path to file to compute suboptimal structures for
    :param mfold: mfold parameters (P, W, M)
    :returns: list of SeqRecord objects wit suboptimal structures as letter_annotations in them
    """
    P = mfold[0]
    W = mfold[1]
    M = mfold[2]

    FNULL = open(os.devnull, 'w')

    # run_prefix = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(6))
    # hybrid ss operates from call directory, output are ctfiles
    if not os.path.isfile(in_path):
        FNULL.close()
        raise FileNotFoundError('Provided file path not found %s', in_path)
    repeat = 0
    done = False
    while repeat < 5:
        try:
            rt = call(
                [
                    '{}hybrid-ss-min'.format(CONFIG.mfold_path),
                    '--suffix=DAT',
                    '--NA=RNA',
                    '--noisolate',
                    '--mfold=' + str(P) + ',' + str(W) + ',' + str(M),
                    in_path
                ],
                # cwd=''.join(os.path.join(os.path.split(in_path)[:-1])),
                cwd=os.path.dirname(in_path),
                stdout=FNULL,
                stderr=STDOUT
            )

            if rt:
                raise ChildProcessError('Execution of hybrid-ss-min failed')

            # well, there is only one file, with the same name as the input file had
            if not os.path.isfile(in_path + '.ct'):
                raise ChildProcessError('Execution of hybrid-ss-min failed - not output file\n'
                                        'expected {}'.format(in_path + '.ct'))

            with open(in_path + '.ct', 'r') as sout:
                suboptimals = ct2db(sout)

            done = True
            break
        except ChildProcessError:
            # try again
            repeat += 1
            warn('hybrid-ss-min failed, trying again {} times'.format(5-repeat))

    if repeat >= 5 and not done:
        warn('There is an issue with hybrid-ss-min, that it does not work for certain combination of energies,'
             ' window parameters and desired number of structures.'
             'The issue also appears to be related to exact order of sequences in the input file'
             'we well try to resolve this by predicting each sequence separetly')

        repeat = 1
        # get the order of file
        fasta_file = []
        with open(in_path, 'r') as f:
            for seqr in SeqIO.parse(f, 'fasta'):
                fasta_file.append(seqr)

        while repeat < 5:
            try:

                # shuffle input file
                # the shuffle operates directly on the list, it doesnt return anything
                shuff_fasta = fasta_file.copy()
                shuffle(shuff_fasta)

                # write shuffled file and rewrite input path
                fd, retry_path = mkstemp()
                with os.fdopen(fd, 'w') as fp:
                    for r in shuff_fasta:
                        fp.write('>{}\n{}\n'.format(r.id, str(r.seq)))

                rt = call(
                    [
                        '{}hybrid-ss-min'.format(CONFIG.mfold_path),
                        '--suffix=DAT',
                        '--NA=RNA',
                        '--noisolate',
                        '--mfold=' + str(P) + ',' + str(W) + ',' + str(M),
                        retry_path
                    ],
                    cwd=''.join(os.path.join(os.path.split(retry_path)[:-1])),
                    stdout=FNULL,
                    stderr=STDOUT
                )

                if rt:
                    print('exit code of child process: {}'.format(rt))
                    raise ChildProcessError('Execution of hybrid-ss-min failed')

                # well, there is only one file, with the same name as the input file had

                if not os.path.isfile(retry_path + '.ct'):
                    raise ChildProcessError('Execution of hybrid-ss-min failed - not output file\n'
                                            'expected {} to be present'.format(retry_path + '.ct'))

                with open(retry_path + '.ct', 'r') as sout:
                    shuff_suboptimals = ct2db(sout)

                # rearrange entries (maybe unnecessary, but to be sure)
                suboptimals = []
                s_s_names = [i.id for i in shuff_suboptimals]
                for seqr in fasta_file:
                    suboptimals.append(shuff_suboptimals[s_s_names.index(seqr.id)])

                done = True
                break
            except ChildProcessError:
                # try again
                repeat += 1
                # adjust the max amount of produced structures
                warn('hybrid-ss-min failed, trying again {} times'.format(5-repeat))
            finally:
                if retry_path:
                    for ext in ['.run', '.plot', '.dG', '.ct', '.ann', '']:
                        try:
                            rpath = retry_path + ext
                            os.remove(rpath.strip())
                        except FileNotFoundError:
                            warn('cannot remove file: {}, file not found'.format(retry_path))
                        except OSError:
                            warn('cannot remove file: {}, file is directory'.format(retry_path))

        if not done:
            # make last try: predict one by one
            print('Last try to predict structures with hybrid-ss-min')
            suboptimals = []
            for seq in fasta_file:

                fd, retry_path = mkstemp()
                with os.fdopen(fd, 'w') as fid:
                    fid.write('>{}\n{}\n'.format(seq.id, str(seq.seq)))

                rt = call(
                    [
                        '{}hybrid-ss-min'.format(CONFIG.mfold_path),
                        '--suffix=DAT',
                        '--NA=RNA',
                        '--noisolate',
                        '--mfold=' + str(P) + ',' + str(W) + ',' + str(M),
                        retry_path
                    ],
                    cwd=''.join(os.path.join(os.path.split(retry_path)[:-1])),
                    stdout=FNULL,
                    stderr=STDOUT
                )

                if rt:
                    raise ChildProcessError('Execution of hybrid-ss-min failed')

                if not os.path.isfile(retry_path + '.ct'):
                    raise ChildProcessError('Execution of hybrid-ss-min failed - not output file')

                with open(retry_path + '.ct', 'r') as sout:
                    suboptimals.append(ct2db(sout)[0])

                os.remove(retry_path + '.run')
                os.remove(retry_path + '.plot')
                os.remove(retry_path + '.dG')
                os.remove(retry_path + '.ct')
                os.remove(retry_path + '.ann')
                os.remove(retry_path)

            done = True

        if done:
            pass
        else:
            FNULL.close()
            raise Exception('hybrid-ss-failed with unexpected behaviour, please contact developers')

    # remove files that hybrid-ss-min created
    # the source file is removed elsewhere
    os.remove(in_path + '.run')
    os.remove(in_path + '.plot')
    os.remove(in_path + '.dG')
    os.remove(in_path + '.ct')
    os.remove(in_path + '.ann')

    with open(in_path, 'r') as fin:
        inseqs = [s for s in SeqIO.parse(fin, format='fasta')]
        suboptimals = rebuild_structures_output_from_pred(inseqs, suboptimals, method='hybrid-ss-min')

        return suboptimals


def ct2db(f, energy_txt='dG'):
    """read ct Zuker file
    there are two versions, first, in ct file there in principle can be pseudoknots encoded
    and simple db notation cannot handle them
    second, when pseudoknots are present, wee need to use different king of brackets
    also, there is WUSS notation for different types of pairs present
    so cogent package is weird, there are references, that appears out of nowhere, and i don't even know
    if the code has worked in the first place
    :returns list of Seq Record objects,
    each Seq Record object contains suboptimal structures in its 'letter_annotations'
    in dict with keys termed ss1 ss2 ..., the list of keys to letter_annotations is in annotations under the key 'sss'"""
    # if not isinstance(f, 'file'):
    #     print('accepts open handle to a file')
    #     raise FileNotFoundError
    txt = f.readline()
    out = []    # list of seq record objects
    prev_seq = ''
    prev_name = ''

    # sl_re = '\d+\s*(?=dG)'
    sl_re = '\d+\s*(?=' + energy_txt + ')'

    # q_re = '(?<=dG\s=)\s*-?\d+\.?\d*'
    q_re = '(?<=' + energy_txt + '\s=)\s*-?\d+\.?\d*'

    while txt:
        seq_len = int(re.search(sl_re,txt).group().lstrip().rstrip())
        q = re.search(q_re, txt)
        # seq_e = float(q.group().lstrip().rstrip())
        seq_name = txt[q.end():].lstrip().rstrip()

        seq_seq = bytearray('#' * seq_len, 'utf-8')
        seq_str = bytearray('#' * seq_len, 'utf-8')
        txt = f.readline()
        while txt:
            # this block parses one structure
            l = txt.split()
            # adjust pos for python indexing
            pos = int(l[0]) - 1
            # print(str(pos) + ' ' + l[1] + '\n')
            seq_seq[pos] = ord(l[1])

            # i is position of pair which is pos paired with
            # adjust for python indexing
            i = int(l[4]) - 1
            if i == -1:
                seq_str[pos] = 46   # '.'
            elif pos < i:
                seq_str[pos] = 40   # '('
                seq_str[i] = 41     # ')'
            else:
                seq_str[pos] = 41   # ')'
                seq_str[i] = 40     # '('

            txt = f.readline()
            # if txt is EOF the loops break itsefl
            if pos + 1 == seq_len:
                # standard end
                # print('standart_end')
                break
            elif 'dG' in txt:
                # print('didnt_reach_end')
                # didn't reach end
                break
        # add parsed structure
        if seq_seq == prev_seq and seq_name == prev_name:
            # only add new structure
            # get keys to dict
            n = len(curr_record.annotations.get('sss'))
            curr_record.letter_annotations['ss' + str(n)] = seq_str.decode('utf-8')
            curr_record.annotations['sss'].append('ss' + str(n))
            prev_seq = seq_seq
            prev_name = seq_name
            if txt == '':
                out.append(curr_record)

        else:
            if not prev_seq == '':
                out.append(curr_record)
            # here create an new record
            curr_record = SeqRecord(Seq(seq_seq.decode('utf-8'), IUPAC.IUPACUnambiguousRNA), id=seq_name, name=seq_name)
            curr_record.annotations['sss'] = ['ss0']
            curr_record.letter_annotations['ss0'] = seq_str.decode('utf-8')
            prev_seq = seq_seq
            prev_name = seq_name

            if txt == '':
                out.append(curr_record)

    return out


def rc_hits_2_rna(seqs, strand=None):
    """ take list of Seq objects and run reverse complement by the state of strand
    strand var possibilities (list or tuple of +1 / -1 // Plus / Minus // + / -)
    also turn all seqs in rna (if some of them was in dna) (the transcribe method)
    :returns list of Seq objects
    """

    # sanity check
    if strand:
        if len(seqs) != len(strand):
            raise AssertionError('seqs and strand must be same length')
    else:
        strand = []
        for seq in seqs:
            assert 'strand' in seq.annotations
            strand.append(seq.annotations['strand'])

    try:
        strand = tuple(strand)
    except:
        print('Unrecognized type of strand, recommended are list and tuple')
        raise

    # determine type of strand
    if strand[0] == 1 or strand[0] == -1:
        p_type = 1
    elif strand[0] == 'Plus' or strand[0] == 'Minus':
        p_type = 'Plus'
    elif strand[0] == '+' or strand[0] == '-':
        p_type = '-'
    else:
        raise RuntimeError('Unrecognized designation of strand, allowed entries are +1 / -1 // "Plus" / "Minus" // "+" '
                           '/ "-"')

    out = []
    for i, seq in enumerate(seqs):
        if strand[i] == p_type:
            out.append(SeqRecord(seq.seq.transcribe(),
                            id=seq.id + 'fw',
                            name=seq.name + 'fw',
                            annotations=seq.annotations,
                            description=seq.description))
        else:
            out.append(SeqRecord(seq.seq.reverse_complement().transcribe(),
                            id=seq.id + 'rc',
                            name=seq.name + 'rc',
                            annotations=seq.annotations,
                            description=seq.description))
    return out


def expand_hits(hits, blast_db, query_length, extra=0, keep_all=0, blast_regexp=None):
    """takes list of blast.HSP objects as first argument and
    path to local blast database as second argument
    then it uses blastdbcmd from blast+ installation to obtain desired sequence
    Two temporary files are used in this call and are deleted at final stage
    :return list of SeqRecord objects (parsed fasta file)
    """
    ext_try = ['nsd', 'nhr', 'nog', 'nsi', 'nin']
    if not os.path.isfile(blast_db + '.nal'):
        for ext in ext_try:
            if not os.path.isfile(blast_db + '.' + ext):
                raise FileNotFoundError('expected file : "' + blast_db + '.' + ext + '" was not found')

    fd, temp_filename = mkstemp()
    fdb, blast_tempfile = mkstemp()
    os.close(fdb)
    exp_hits = []
    strand = []
    loc = []
    temp_file = os.fdopen(fd, 'w')
    for i, hit in enumerate(hits):
        if hit[1].strand[1] == 'Minus':
            start = hit[1].sbjct_end - hit[1].query_start + 1 - extra
            end = hit[1].sbjct_start + _positive_index(query_length - hit[1].query_end) + extra
            strand.append(-1)
            d = {'query_start': hit[1].sbjct_end, 'query_end': hit[1].sbjct_start,
                 'extended_start': hit[1].sbjct_end - hit[1].query_start + 1,
                 'extended_end': hit[1].sbjct_start + _positive_index(query_length - hit[1].query_end),
                 'strand': -1}
        else:
            # ss - positive(qs) - extra
            start = hit[1].sbjct_start - hit[1].query_start + 1 - extra
            # se + positive(len(q) - qe) + extra
            end = hit[1].sbjct_end + _positive_index(query_length - hit[1].query_end) + extra
            strand.append(1)
            d = {'query_start': hit[1].sbjct_start, 'query_end': hit[1].sbjct_end,
                 'extended_start': hit[1].sbjct_start - hit[1].query_start + 1,
                 'extended_end': hit[1].sbjct_end + _positive_index(query_length - hit[1].query_end),
                 'strand': 1}

        # ====== information about possible trim ======
        # assume ok
        d['trimmed_start'] = False
        d['trimmed_end'] = False

        if start <= 0:
            start = 1
            d['trimmed_start'] = True
        d['super_start'] = start
        d['super_end'] = end
        # repair possible extended start violation
        if d['extended_start'] < 1:
            d['extended_start'] = 1
            d['trimmed_start'] = True
        loc.append(d)

        # add blast record
        d['blast'] = hit

        # get a index name to the blastdb
        if blast_regexp is None:
            blast_regexp = accession_regex
        hname = re.search(blast_regexp, hit[0])
        if not hname:
            print(blast_regexp)
            raise RuntimeError('provided regexp returned no result for %s,'
                               ' please provide regexp valid even for this name', hit[0])
        temp_file.write(hname.group() + ' ' + '-region ' + str(start) + '-' + str(end) + '\n')

    temp_file.close()

    cmd = [
            '{}blastdbcmd'.format(CONFIG.blast_path),
            '-dbtype',
            'nucl',
            '-db',
            str(blast_db),
            '-entry_batch',
            temp_filename,
            '-out',
            blast_tempfile
        ]

    r_val = call(
        cmd
    )
    if r_val:
        ChildProcessError('Unable to run blastdbcmd command, please check its availability')

    # records at minus are not in rc
    with open(blast_tempfile, 'r') as bf:
        for i, parsed_record in enumerate(SeqIO.parse(bf, 'fasta')):
            parsed_record.annotations = loc[i]
            # add uid to ensure that all hits are unique
            parsed_record.id = 'uid:' + str(i) + '|' + parsed_record.id

            # ==================================================================
            # add checking, that hits are retrieved whole for edgecases, when
            # hit is close to begening or end of the database sequence
            # repair possible extended end violation
            # relevant to extended_end and super_end ? or super_end only
            if loc[i]['super_start'] + len(parsed_record.seq) - 1 < loc[i]['super_end']:
                print(
                    'Sequence retrieved from database based on hit {} is shorther then expected.'
                    'This either means database difference or that the hit is located near the end of the database '
                    'and could not be extended as requested. Please note that when examining the results.'
                    ''.format(parsed_record.id)
                )
                parsed_record.annotations['super_end'] = parsed_record.annotations['super_start'] + len(parsed_record.seq) - 1
                parsed_record.annotations['trimmed_end'] = True

                if loc[i]['super_start'] + len(parsed_record.seq) - 1 < loc[i]['extended_end']:
                    print(
                        'Retrieved sequence is shorter even for simple extention by unaligned length of query. '
                        'We use whole avalible sequence as deposited in the database.'
                    )

                    parsed_record.annotations['extended_end'] = parsed_record['super_start'] + len(parsed_record.seq) - 1

            exp_hits.append(parsed_record)
    if not keep_all:
        os.remove(temp_filename)
        os.remove(blast_tempfile)

    if len(exp_hits) != len(hits):
        raise LookupError(
            'Some IDs were not found in provided BLAST database.'
            ' Please update the database or remove the hits.'
        )
    return exp_hits, strand


def _positive_index(o):
    """
    this ensures that indexes stay positive and would not colapse extended sequence
    :param o:
    :return:
    """
    if o < 0:
        o = 0
    return o


def blast_hsps2list(parsed_blast):
    """
    reformat blast record to list of hits
    :param parsed_blast:
    :return:
    """
    hit_list = []
    for ba in parsed_blast.alignments:
        for bh in ba.hsps:
            hit_list.append([ba.hit_id, bh])
    return hit_list


def parse_seq_str(fid):
    """
    parses fasta-like file with sequences nad structures
    can be output of Vienna RNAfold - in this case, energy is ommited
    :param fid: open file handle
    :return:
    """
    if not hasattr(fid, 'readline'):
        raise Exception('cannot read from file')
    c = 0
    while True:
        c += 1
        k = fid.readline()
        if k == '':
            break
        if k == '\n':
            continue
        k = k.rstrip('\n')

        if '>' == k[0]:
            name = k[1:].strip()
            seq = fid.readline().rstrip('\n')
            structure = fid.readline().rstrip('\n').rsplit(' ')[0]
        else:
            name = str(c)
            seq = k
            structure = fid.readline().rstrip('\n').rsplit(' ')[0]
        oi = SeqRecord(Seq(seq, IUPAC.IUPACAmbiguousRNA),
                       id=name,
                       name=name)
        oi.annotations['sss'] = ['ss0']
        oi.letter_annotations['ss0'] = structure
        yield oi


def read_seq_str(fasta_like_file):
    with open(fasta_like_file, 'r') as flf:
        return [rec for rec in parse_seq_str(flf)]


def print_parameters(args, f):
    """
    print all input parameters to a file so it is ascertainable what where the commands which produced the data
    :param args:
    :param f:
    :return:
    """

    A = vars(args)

    for ar in A.keys():
        f.write('{}: {}\n'.format(ar, A[ar]))


def print_time():
    return '{:%Y-%m-%d %H:%M}'.format(datetime.datetime.now())


def prevent_mem_in_fasta_name(func):
    """
    decorator function preventing return of a string which starts with 'Met'
    implemented because of problems with t-coffee rcoffee, which fails on sequences starting with 'Met' in sequence name
    :param func:
    :return:
    """
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if res.startswith('Mem'):
            while res.startswith('Mem'):
                res = func(*args, **kwargs)

        return res
    return wrapper


@prevent_mem_in_fasta_name
def generate_random_name(len, rn=()):
    # generate random names
    l = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'

    s = ''.join([choice(l) for i in range(len)])
    while s in rn:
         s = ''.join([choice(l) for i in range(len)])
    return s


def sanitize_fasta_names_in_seqrec_list(seqrec_list, used_dict=None):
    if used_dict:
        sanitize_dict = used_dict
    else:
        sanitize_dict = {}

    sanitized_list = []
    for rec in seqrec_list:
        uname = generate_random_name(8, list(sanitize_dict.keys()))
        sanitize_dict[uname] = rec.id
        ns = SeqRecord(rec.seq, id=uname)

        if len(rec.letter_annotations) != 0:
            ns.letter_annotations = rec.letter_annotations
        if len(rec.annotations) != 0:
            ns.annotations = rec.annotations

        sanitized_list.append(ns)
    return sanitized_list, sanitize_dict


def desanitize_fasta_names_in_seqrec_list(seqrec_list, used_dict):
    reverted = []
    for rec in seqrec_list:
        ns = SeqRecord(rec.seq, id=used_dict[rec.id])

        if len(rec.letter_annotations) != 0:
            ns.letter_annotations = rec.letter_annotations
        if len(rec.annotations) != 0:
            ns.annotations = rec.annotations

        reverted.append(ns)
    return reverted


def run_mlocarna(file_name, mlocarna_params='--sequ-local 0 --struct-local 0 --free-endgaps --stockholm'):
    """
    possible settings:
    struct_local : seq_local
    1 : 1   - breaks sequences and joins them !!!! DO NOT USE
    1 : 0   - tries to align globally, end to end - lot of gaps
    0 : 1   - wanted alignment - unaligned ends of input sequences are trimmed in output alignment
    0 : 0   - wanted alignment - returns whole alignment including unaligned ends - also not good - tries to align
                ends of sequences together - then i have to trimm alignemnt and check if structure is correct
                which it isn't

    return locarna output directory
    :param file_name:
    :param mlocarna_params: mlocarna_params
    :return:
    """

    FNULL = open(os.devnull, 'w')

    if not os.path.isfile(file_name):
        raise FileNotFoundError('Provided file {} was not found'.format(file_name))
    # if not os.path.isfile(subject_file):
    #     raise FileNotFoundError('Provided file {} was not found'.format(subject_file))

    r = call(
        '{}mlocarna {} --tgtdir {} {}'.format(
            CONFIG.locarna_path,
            mlocarna_params,
            file_name + 'out',
            file_name
        ),
        shell=True, stdout=FNULL
    )

    FNULL.close()

    if r:
        raise ChildProcessError('call to mlocarna failed for files in:{}'.format(file_name))

    return file_name + 'out'


def select_sequences_from_similarity_rec(dist_mat: np.ndarray, sim_threshold_percent=90) -> list:
    """
    :param dist_mat: distmat table, by default obtained from read_clustal_distmat_file, values in percent
    :param sim_threshold_percent: threshold for similarity in percent
    :return:
    """

    # dists = np.triu(dist_mat.as_matrix(), 1)          # removes unwanted similarities
    if dist_mat is None:
        return 0,
    dists = dist_mat.transpose()
    # row, col = where(dists > sim_threshold_percent) # determine where the similarities are
    include = set()
    exclude = set()
    a = np.array(range(len(dists)))
    for i, r in enumerate(dists):
        pr = r[~np.isnan(r)]
        pa = a[~np.isnan(r)]
        if (i in exclude) | (any(pr >= sim_threshold_percent)):
            pu = np.where(pr >= sim_threshold_percent)
            u = pa[pu]
            if i not in exclude:
                include |= {i}
            to_ex = set(u.tolist()) - include
            exclude |= to_ex                         # union operation
        else:
            include |= {i}

    return sorted(include)


def run_muscle(fasta_file, out_file=None, muscle_params='', reorder=True):
    """
    beware, muscle does not keep sequence order and the --stable switch is broken
    :param fasta_file:
    :param out_file:
    :param muscle_params:
    :return:
    """
    if out_file:
        cl_file = out_file
    else:
        cl_fd, cl_file = mkstemp()
        os.close(cl_fd)

    r = call(
        '{}muscle -clwstrict -seqtype rna -out {} -in {} {} -quiet'.format(
            CONFIG.muscle_path,
            cl_file,
            fasta_file,
            muscle_params
        ),
        shell=True
    )
    if r:
        raise ChildProcessError('call to muscle failed')

    if reorder:
        # reorder sequences acording to input file
        with open(fasta_file, 'r') as ff, open(cl_file, 'r+') as oo:
            orig_seqs = [i.id for i in SeqIO.parse(ff, format='fasta')]
            muscle_align = {i.id: i for i in AlignIO.read(oo, format='clustal')}

            # reorder
            reo_alig = []
            for s_name in orig_seqs:
                # muscle cuts names
                reo_alig.append(muscle_align[s_name[:32]])
            alig = AlignIO.MultipleSeqAlignment(reo_alig)
            # write
            oo.seek(0)
            AlignIO.write(alig, oo, format='clustal')
            oo.truncate()

    return cl_file


def run_mafft(fasta_file, mafft_mode='qinsi', out_file=None, mafft_params='', thread=1):
    avalible_modes = {'qinsi', 'linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwns', 'nwnsi', 'profile'}
    if mafft_mode not in avalible_modes:
        raise Exception('MAFFT mode can be only {}, not {}'.format(avalible_modes, mafft_mode))

    if out_file:
        cl_file = out_file
    else:
        cl_fd, cl_file = mkstemp()
        os.close(cl_fd)

    r = call(
        '{}mafft-{} --clustalout {} --thread {} {} > {}'.format(
            CONFIG.mafft_path,
            mafft_mode,
            mafft_params,
            thread,
            fasta_file,
            cl_file
        ),
        shell=True
    )
    if r:
        raise ChildProcessError('call to mafft failed for files {} {}'.format(fasta_file, cl_file))
    return cl_file


def RNAfold(sequence):
    r = check_output(
        [
            '{}RNAfold'.format(CONFIG.viennarna_path),
            '--noPS',
        ],
        input=sequence.encode()
    )

    if isinstance(r, CalledProcessError):
        print('RNAfold failed')
        raise ChildProcessError

    # more robust decode
    out_str = r.decode()
    spl = out_str.split('\n')
    seq = spl[0]
    structure = spl[1][:len(seq)]
    energy = float(spl[1][len(seq)+2:-1])
    # seq, structure, energy = r.decode().split()
    # return seq, structure, float(energy[1:-1])
    return seq, structure, energy


def select_analyzed_aligned_hit(one_alig, exp_hit_id):
    # select analyzed hit, drop query
    for seq_in_alig in one_alig.get_unalined_seqs(keep_letter_ann=True):
        if seq_in_alig.id == exp_hit_id:
            return seq_in_alig


def rebuild_structures_output_from_pred(reference_sequences_list, predicted_structures_list, method=None):
    structuresids = [seq.id for seq in predicted_structures_list]
    structures_list = []
    for seq in reference_sequences_list:
        nr = SeqRecord(
            seq.seq,
            id=seq.id,
            name=seq.name,
            description=seq.description,
            annotations=seq.annotations,
            letter_annotations=seq.letter_annotations
        )
        if seq.id in structuresids:
            n = structuresids.index(seq.id)
            # nr.letter_annotations['ss0'] = predicted_structures_list[n].letter_annotations['ss0']
            nr.letter_annotations.update(predicted_structures_list[n].letter_annotations)
            nr.annotations.update(predicted_structures_list[n].annotations)
            nr.annotations['predicted'] = True
        else:
            if method:
                wmsg = 'method: {} failed to predict structure for seq {}'.format(method, nr.id)
                warn(wmsg)
            nr.annotations['predicted'] = False

        structures_list.append(nr)
        del nr

    return structures_list


def blast_in(blast_in, b='guess'):
    """
    gueses blast format
    :returns list of SearchIO objects, one set of hits per per query sequence per field
    """
    multiq = []
    with open(blast_in, 'r') as f:
        if b == 'guess':
            # gues the format
            l = f.readline()
            f.seek(0,0)     # seek to begining
            if re.search('^BLASTN \d+\.\d+\.\d+',l):
                # blast object prob plaintext
                b_type = 'plain'
            elif re.search('<\?xml version', l):
                # run xml parser
                b_type = 'xml'
            else:
                print('could not guess the blast format, preferred format is NCBI xml')
                raise AssertionError
        else:
            b_type = b

        if b_type == 'plain':
            try:
                for p in blast_minimal_parser(f):
                    multiq.append(p)
            except:
                raise IOError('Failed to parse provided file {} as text blast output'.format(blast_in))
        elif b_type == 'xml':
            try:
                for p in NCBIXML.parse(f):
                    multiq.append(p)
            except:
                raise IOError('Failed to parse provided file {} as xml'.format(blast_in))
        else:
            raise AttributeError('Type not known: allowed types: plain, xml, guess')

        # todo from SearchIO add remaining functionality
    return multiq


def run_rnaplot(seq, structure=None, format='svg', outfile=None):
    """
    run rnaplot in desired format
    if seq
    :param seq:
    :param structure:
    :return:
    """
    if structure is None:
        sequence = str(seq.seq)
        structure = seq.letter_annotations['ss0']
    else:
        sequence = seq

    assert len(sequence) == len(structure)

    allowed_formats = {'ps', 'svg', 'gml', 'xrna'}
    if format not in allowed_formats:
        raise TypeError('Format can be only from {}.'.format(allowed_formats))

    fd, tempfile = mkstemp()

    rnaname = tempfile.split('/')[-1].split('\\')[-1]

    currdirr = os.getcwd()
    tmpdir = gettempdir()
    os.chdir(tmpdir)

    with os.fdopen(fd, 'w') as fh:
        fh.write('>{}\n{}\n{}\n'.format(
            rnaname,
            sequence,
            structure
        ))

    r = call(
        '{}RNAplot --output-format={} < {}'.format(
            CONFIG.viennarna_path,
            format,
            tempfile
        ),
        shell=True
    )
    if r:
        print('Call to RNAplot failed')
        print('cmd: {}'.format('{}RNAplot --output-format={} < {}'.format(CONFIG.viennarna_path, format, tempfile)))
        os.chdir(currdirr)
        raise ChildProcessError

    # output file is name of the sequence (in this case name of the file) + "_ss." + chosen format

    os.remove(tempfile)

    plot_output_file = os.path.join(tmpdir, rnaname + '_ss.' + format)
    os.chdir(currdirr)
    if outfile is None:
        return plot_output_file
    else:
        shutil.move(os.path.join(tmpdir, plot_output_file), outfile)
        return outfile


class Subsequences(object):
    """
    wraper class for subsequences
    attributes
    source = source sequence for subsequences in SeqRecord format
    subs = dict where subsequences should be placed
    """

    def __init__(self, rec):
        if not isinstance(rec, SeqRecord):
            raise Exception('Class can only be invoked with its source sequence in SeqRecord format from Biopython'
                            ' package')
        self.subs = {}
        self.source = rec
        self.ret_keys = None
        self.query_name = ''
        self.best_start = None
        self.best_end = None
        self.templates = {}


def format_blast_hit(hit, linelength=90, delim='_'):
    """
    format blast hit for printing
    :param hit: blast hit hsp data structure
    :param linelength: limit for line length for printing
    :param delim: char to use as a fill (whitespace sometimes gets wrapped and behaves incorrectly)
    :return:
    """
    qs = hit.query_start
    qe = hit.query_end
    ss = hit.sbjct_start
    se = hit.sbjct_end
    ml = max([len(i) for i in [str(qs), str(qe), str(ss), str(se)]])
    ad = linelength - (2*ml + 10 + 4)
    qq = StringIO(hit.query)
    mm = StringIO(hit.match)
    sss = StringIO(hit.sbjct)
    try:
        f_string = ''
        q = qq.read(ad)
        m = mm.read(ad)
        s = sss.read(ad)
        while True:
            l = len(q)
            if l == 0:
                break
            if len(m) < len(q):
                m = '{:{width}}'.format(m, width=len(q))

            qe = qs + l - 1 - q.count('-')
            if hit.sbjct_start < hit.sbjct_end:
                se = ss + l - 1 - s.count('-')
            else:
                se = ss - l + s.count('-') + 1

            f_string = f_string + 'Query' + delim + delim*(ml-len(str(qs))) + str(qs) + \
                       ' ' + q + ' ' + \
                       str(qe) + delim*(ml-len(str(qe))) + '\n' + \
                       delim*(ml + 6) + \
                       ' ' + m + ' ' + \
                       delim*(ml) + '\n' + \
                       'Sbjct' + delim + delim*(ml-len(str(ss))) + str(ss) + \
                       ' ' + s + ' ' + \
                       str(se) + delim*(ml-len(str(se))) + '\n\n'

            qs = qe + 1
            if hit.sbjct_start < hit.sbjct_end:
                ss = se + 1
            else:
                ss = se - 1

            if l < ad:
                break

            q = qq.read(ad)
            m = mm.read(ad)
            s = sss.read(ad)
    finally:
        qq.close()
        mm.close()
        sss.close()
    if not (qe == hit.query_end and se == hit.sbjct_end):
        raise Exception('blast write fail')
    return f_string


def devprint(s, **kwargs):
    print('{} - {}'.format(os.getpid(), s), **kwargs)


def blasthsp2pre(bhsp):
    """
    printing method for blast high scoring pair
    :param bhsp: Bio.Blast.Record.HSP
    :return:
    """
    # todo add printing with merged blast hits (ie overlapping hits or similar)
    eformat='.2E'

    if bhsp.strand[1] is None:
        strand = ''
    else:
        strand = bhsp.strand[1]

    score_template = ''.join(
        (
            " Score = {score:.1f} bits ({bits:.1f}), Expect = {expect:{eformat:}}\n"
            " Identities = {identities:d}/{length:d} ({id_p:d}%), Gaps = {gaps:d}/{length:d} ({gap_p:d}%)\n"
            " Strand = Plus/{mystrand:}\n"
        )
    )
    myf = score_template.format(
        score=bhsp.score,
        bits=bhsp.bits,
        expect=bhsp.expect,
        eformat=eformat,
        identities=bhsp.identities[0],
        gaps=bhsp.gaps[0],
        length=bhsp.align_length,
        id_p=round(bhsp.identities[0]/bhsp.align_length * 100),
        gap_p=round(bhsp.gaps[0]/bhsp.align_length * 100),
        mystrand=strand,
    )

    if hasattr(bhsp, 'features') and bhsp.features:
        if isinstance(bhsp.features, str):
            fk = 'in'
            features_string = bhsp.features
        elif isinstance(bhsp.features, (list, tuple)):
            fk = 'flanking'
            features_string = ''.join(bhsp.features)
        else:
            raise Exception("unknown features type")
        features_template = ''.join(
            (
                " Features {in_or_flank:} this part of subject sequence:\n",
                "{features:}\n",
            )
        )
        features = features_template.format(
            in_or_flank=fk,
            features=features_string,
        )

        myf += '\n' + features

    blast_text_hit_body = format_blast_hit(bhsp, linelength=90, delim=' ')
    return myf + blast_text_hit_body


def remove_one_file_with_try(file, msg_template=''):
    """
    Removes file with call to os.remove, prints messsage if problem instead of raising exception
    :param file: path (str)
    :param msg_template: message to embed with errors
    :return:
    """
    try:
        os.remove(file)
    except FileNotFoundError:
        print('{} File: {} Not Found.'.format(msg_template, file))
    except OSError:
        print('{} File: {} is a directory.'.format(msg_template, file))
    return


def remove_files_with_try(files, msg_template):
    for fi in files:
        remove_one_file_with_try(fi, msg_template)
    return


class NoHomologousSequenceException(Exception):
    pass


def non_redundant_seqs(sequences: list) -> list:
    """
    return list of non redundant sequences
     the first sequence encountered is stored
    :param sequences:
    :return:
    """
    nr_seqs = set()
    out = []
    for seq in sequences:
        str_seq = str(seq.seq)
        if str_seq not in nr_seqs:
            out.append(seq)
            nr_seqs.add(str_seq)

    return out


def parse_named_structure_file(file):
    with open(file, 'r') as f:
        for sr in parse_one_rec_in_multiline_structure(f):
            cf = sr.strip().splitlines()
            cfr = SeqRecord(Seq(cf[1]), id=cf[0])
            cfr.annotations['sss'] = []
            for i, ll in enumerate(cf[2:]):
                structure, str_name = ll.split()
                cfr.letter_annotations[str_name] = structure
                cfr.annotations['sss'].append(str_name)

            yield cfr


def parse_one_rec_in_multiline_structure(fh):
    """
    in fasta-like multiple structures file iterates through records one by one
    :param fh:
    :return:
    """
    r = []
    txt = fh.readline()
    while txt:
        if txt[0] == '>' and len(r) != 0:
            l = ''.join(r)
            r = []
            yield l
        r.append(txt)
        txt = fh.readline()

    if len(r) != 0:
        yield ''.join(r)
