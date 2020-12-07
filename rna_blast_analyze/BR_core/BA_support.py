import locale
import logging
import os
import re
from io import StringIO
from random import choice
from subprocess import call
from tempfile import mkstemp, TemporaryFile
import shlex
import sys

from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.parser_to_bio_blast import blast_parse_txt
from rna_blast_analyze.BR_core import exceptions

# idiotic matplotlib

locale.setlocale(locale.LC_ALL, 'C')

ml = logging.getLogger('rboAnalyzer')


def rc_hits_2_rna(seqs, strand=None):
    """ take list of Seq objects and run reverse complement by the state of strand
    strand var possibilities (list or tuple of +1 / -1 // Plus / Minus // + / -)
    also turn all seqs in rna (if some of them was in dna) (the transcribe method)
    :returns list of Seq objects
    """
    ml.debug(fname())
    # sanity check
    if strand:
        if len(seqs) != len(strand):
            raise AssertionError('seqs and strand must be same length')
    else:
        strand = []
        for seq in seqs:
            assert 'strand' in seq.annotations
            strand.append(seq.annotations['strand'])

    if not isinstance(strand, (list, tuple)):
        raise AssertionError('Unrecognized type of strand, recommended are list and tuple')

    plus_strands = {1, '+', 'Plus'}
    minus_strands = {-1, '-', 'Minus'}

    out = []
    for i, seq in enumerate(seqs):
        if strand[i] in plus_strands:
            out.append(
                SeqRecord(
                    seq.seq.transcribe(),
                    id=seq.id + 'fw',
                    name=seq.name + 'fw',
                    annotations=seq.annotations,
                    description=seq.description
                )
            )
        elif strand[i] in minus_strands:
            out.append(
                SeqRecord(
                    seq.seq.reverse_complement().transcribe(),
                    id=seq.id + 'rc',
                    name=seq.name + 'rc',
                    annotations=seq.annotations,
                    description=seq.description
                )
            )
        else:
            raise AssertionError(
                'Unrecognized designation of strand, allowed entries are +1 / -1 // "Plus" / "Minus" // "+" / "-"'
            )
    return out


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
            name = k[1:].strip().split(' ')
            id = name[0]
            description = ' '.join(name[1:])
            seq = fid.readline().rstrip('\n')
            structure = fid.readline().rstrip('\n').rsplit(' ')[0]
        else:
            id = str(c)
            seq = k
            description = ''
            structure = fid.readline().rstrip('\n').rsplit(' ')[0]
        oi = SeqRecord(
            Seq(seq),
            id=id,
            description=description
        )
        oi.annotations['sss'] = ['ss0']
        oi.letter_annotations['ss0'] = structure
        yield oi


def read_seq_str(fasta_like_file):
    with open(fasta_like_file, 'r') as flf:
        return [rec for rec in parse_seq_str(flf)]


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
def generate_random_name(lenght, rn=()):
    # generate random names
    l = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_'

    s = ''.join([choice(l) for _ in range(lenght)])
    while s in rn:
        s = ''.join([choice(l) for _ in range(lenght)])
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


def run_muscle(fasta_file, out_file=None, muscle_params='', reorder=True):
    """
    beware, muscle does not keep sequence order and the --stable switch is broken
    :param fasta_file:
    :param out_file:
    :param muscle_params:
    :param reorder:
    :return:
    """
    ml.info('Running muscle.')
    ml.debug(fname())
    if out_file:
        cl_file = out_file
    else:
        cl_fd, cl_file = mkstemp(prefix='rba_', suffix='_07', dir=CONFIG.tmpdir)
        os.close(cl_fd)

    cmd = [
        '{}muscle'.format(CONFIG.muscle_path),
        '-clwstrict',
        '-seqtype', 'rna',
        '-out', cl_file,
        '-in', fasta_file,
        '-quiet']
    if muscle_params != '':
        cmd += [
            ' '.join([shlex.quote(i) for i in shlex.split(muscle_params)])
        ]
    ml.debug(cmd)

    with TemporaryFile(mode='w+', encoding='utf-8') as tmp:
        r = call(cmd, stdout=tmp, stderr=tmp)
        if r:
            msgfail = 'Call to muscle failed.'
            ml.error(msgfail)

            tmp.seek(0)
            raise exceptions.MuscleException(msgfail, tmp.read())

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


def select_analyzed_aligned_hit(one_alig, exp_hit_id):
    # select analyzed hit, drop query
    for seq_in_alig in one_alig.get_unalined_seqs(keep_letter_ann=True):
        if seq_in_alig.id == exp_hit_id:
            return seq_in_alig


def rebuild_structures_output_from_pred(reference_sequences_list, predicted_structures_list, method=None):
    ml.debug(fname())
    structuresids = [seq.id for seq in predicted_structures_list if hasattr(seq, 'id')]
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
            nr.letter_annotations.update(predicted_structures_list[n].letter_annotations)
            nr.annotations.update(predicted_structures_list[n].annotations)
            nr.annotations['predicted'] = True
        else:
            if method:
                wmsg = '{} failed to predict structure for seq {}.'.format(method, nr.id)
                ml.warning(wmsg)
            nr.annotations['predicted'] = False

        structures_list.append(nr)
        del nr

    return structures_list


def blast_in(blast_input_file, b='guess'):
    """
    gueses blast format
    :returns list of SearchIO objects, one set of hits per per query sequence per field
    """
    ml.debug(fname())

    try:
        with open(blast_input_file, 'r') as f:
            multiq = blast_in_handle(f, b)
            if multiq is None:
                raise Exception
            return multiq
    except Exception as e:
        msgfail = 'Failed to parse provided file {}'.format(blast_input_file)
        ml.error(msgfail)
        ml.error(str(e))
        sys.exit(1)


def blast_in_handle(handle, b='guess', log=True):
    """
    gueses blast format
    :returns list of SearchIO objects, one set of hits per per query sequence per field
    """
    if log:
        ml.debug(fname())
    multiq = []
    if b == 'guess':
        # gues the format
        l = handle.readline()
        handle.seek(0, 0)     # seek to begining
        if re.search(r'^BLASTN \d+\.\d+\.\d+', l):
            # blast object prob plaintext
            b_type = 'plain'
            if log:
                ml.info('Inferred BLAST format: txt.')
        elif re.search(r'<\?xml version', l):
            # run xml parser
            b_type = 'xml'
            if log:
                ml.info('Inferred BLAST format: xml.')
        else:
            if log:
                ml.error('Could not guess the BLAST format, preferred format is NCBI xml.')
            return None
    else:
        b_type = b

    if b_type == 'plain':
        for p in blast_parse_txt(handle):
            multiq.append(p)

        return multiq

    elif b_type == 'xml':
        for p in NCBIXML.parse(handle):
            multiq.append(p)

        return multiq
    else:
        if log:
            ml.error('BLAST type not known: allowed types: plain, xml, guess')
        return None


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
        self.extension = None
        self.source = rec
        self.query_name = ''
        self.best_start = None
        self.best_end = None
        self.templates = {}


def blasthsp2pre(bhsp):
    """
    printing method for blast high scoring pair
    :param bhsp: Bio.Blast.Record.HSP
    :return:
    """
    eformat = '.2E'

    if bhsp.sbjct_start < bhsp.sbjct_end:
        strand = 'Plus'
    else:
        strand = 'Minus'

    score_template = ''.join(
        (
            " Score = {score:.1f} bits ({bits:.1f}), Expect = {expect:{eformat:}}\n"
            " Identities = {identities:d}/{length:d} ({id_p:d}%), Gaps = {gaps:d}/{length:d} ({gap_p:d}%)\n"
            " Strand = Plus/{mystrand:}\n"
        )
    )
    if isinstance(bhsp.identities, (tuple, list)):
        idd = bhsp.identities[0]
    elif isinstance(bhsp.identities, int):
        idd = bhsp.identities
    else:
        ml.warning('Format HSP: identities in unexpected format')
        idd = 0
    if isinstance(bhsp.gaps, (tuple, list)):
        gap = bhsp.gaps[0]
    elif isinstance(bhsp.gaps, int):
        gap = bhsp.gaps
    else:
        ml.warning('Format HSP: gaps in unexpected format')
        gap = 0

    myf = score_template.format(
        score=bhsp.score,
        bits=bhsp.bits,
        expect=bhsp.expect,
        eformat=eformat,
        identities=idd,
        gaps=gap,
        length=bhsp.align_length,
        id_p=round(idd/bhsp.align_length * 100),
        gap_p=round(gap/bhsp.align_length * 100),
        mystrand=strand,
    )

    try:
        if hasattr(bhsp, 'features') and bhsp.features:
            if isinstance(bhsp.features, str):
                fk = 'in'
                features_string = bhsp.features
            elif isinstance(bhsp.features, (list, tuple)):
                fk = 'flanking'
                features_string = ''.join(bhsp.features)
            else:
                raise exceptions.BlastFormatException("Format HSP: unknown features type")

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
    except exceptions.BlastFormatException as e:
        ml.warning(str(e))

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
        ml.debug('{} File: {} Not Found.'.format(msg_template, file))
    except OSError:
        ml.debug('{} File: {} is a directory.'.format(msg_template, file))
    return


def remove_files_with_try(files, msg_template=''):
    for fi in files:
        remove_one_file_with_try(fi, msg_template)
    return


def non_redundant_seqs(sequences: list) -> list:
    """
    return list of non redundant sequences
     the first sequence encountered is stored
    :param sequences:
    :return:
    """
    ml.debug(fname())
    nr_seqs = set()
    out = []
    for seq in sequences:
        str_seq = str(seq.seq)
        if str_seq not in nr_seqs:
            out.append(seq)
            nr_seqs.add(str_seq)

    return out


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


def filter_ambiguous_seqs_from_list(seqlist):
    return [seq for seq in seqlist if not seq.annotations['ambiguous']]


def filter_by_length_diff(sequences, ref_len, len_diff_):
    return [
        seq for seq in sequences if ref_len * (1 - len_diff_) < len(seq) < ref_len * (1 + len_diff_)
    ]


def sel_seq_simple(all_seqs, query, len_diff):
    """Select nr seqs simply

    Will return non redundant non-ambiguous seqs with cm bit_sc > 0 (must be precomputed)

    """
    q_all = [query] + all_seqs
    q_cm = [s for s in q_all if s.annotations['cmstat'].bit_sc > 0]
    not_amb = filter_ambiguous_seqs_from_list(q_cm)
    nr_na = non_redundant_seqs(not_amb)
    nr_na_ld = filter_by_length_diff(nr_na, len(query), len_diff)
    return nr_na_ld


def annotate_ambiguos_bases(seqlist):
    iupac = IUPACmapping()
    reg = re.compile("[^" + "^".join(iupac.unambiguous) + "]+", re.IGNORECASE)
    for seq in seqlist:
        annotate_ambiguos_base(seq, iupac=iupac, reg=reg)
    return seqlist


def annotate_ambiguos_base(seq, iupac=None, reg=None):
    if 'ambiguous' in seq.annotations:
        return seq
    if iupac is None:
        iupac = IUPACmapping()
    if reg is None:
        reg = re.compile("[^" + "^".join(iupac.unambiguous) + "]+", re.IGNORECASE)
    m = re.search(reg, str(seq.seq))
    if m:
        msg = "Ambiguous base detected in {}, violating base {}, pos {}".format(
            seq.id,
            m.group(),
            m.start()
        )
        ml.warning(msg)
        seq.annotations['ambiguous'] = True
        if 'msgs' not in seq.annotations:
            seq.annotations['msgs'] = []
        seq.annotations['msgs'].append(msg)
    else:
        seq.annotations['ambiguous'] = False

    return seq


class IUPACmapping(object):
    """
    holder for the IUPAC ambigues information
    """
    def __init__(self):
        self.allowed_bases = tuple('ACGTURYSWKMBDHVN')
        self.ambiguous = tuple('RYSWKMBDHVN')
        self.unambiguous = tuple(set(self.allowed_bases) - set(self.ambiguous))
        self.mapping = {
            # 'A': ('A',),
            # 'C': ('C',),
            # 'G': ('G',),
            # 'T': ('T',),
            # 'U': ('U',),
            'R': ('A', 'G'),
            'Y': ('C', 'T'),
            'S': ('G', 'C'),
            'W': ('A', 'T'),
            'K': ('G', 'T'),
            'M': ('A', 'C'),
            'B': ('C', 'G', 'T'),
            'D': ('A', 'G', 'T'),
            'H': ('A', 'C', 'T'),
            'V': ('A', 'C', 'G'),
            'N': ('A', 'C', 'G', 'T')
        }
        self.rnamapping = self._RNAze()

    def _RNAze(self):
        """
        change the 'T' to 'U' in mapping
        :return:
        """
        rnamapping = dict()
        for key in self.mapping.keys():
            if 'T' in self.mapping[key]:
                rnamapping[key] = tuple([i for i in self.mapping[key] if i != 'T'] + ['U'])
            else:
                rnamapping[key] = self.mapping[key]
        return rnamapping


def iter2file_name(orig_file, multi_query, iteration):
    if multi_query:
        if '.' in orig_file:
            parts = orig_file.split('.')
            file = '.'.join(parts[:-1]) + '_({}).{}'.format(iteration, parts[-1])
        else:
            file = '{}_({})'.format(orig_file, iteration)
    else:
        file = orig_file
    return file


def blast_hit_getter_from_hits(x):
    return x[1]


def blast_hit_getter_from_subseq(x):
    return x.source.annotations['blast'][1]


def ct2db(f, energy_txt='dG'):
    """
    read ct Zuker file
        there are two versions,
        first, in ct file there in principle can be pseudoknots encoded and simple db notation cannot handle them
        second, when pseudoknots are present, wee need to use different king of brackets
        also, there is WUSS notation for different types of pairs present
    :returns list of Seq Record objects,
        each Seq Record object contains suboptimal structures in its 'letter_annotations'
        in dict with keys termed ss1 ss2 ...,
        the list of keys to letter_annotations is in annotations under the key 'sss'
    """
    ml.debug(fname())

    txt = f.readline()
    out = []    # list of seq record objects
    prev_seq = ''
    prev_name = ''

    # sl_re = '\d+\s*(?=dG)'
    sl_re = r'\d+\s*(?=' + energy_txt + ')'

    # q_re = '(?<=dG\s=)\s*-?\d+\.?\d*'
    q_re = '(?<=' + energy_txt + r'\s=)\s*-?\d+\.?\d*'

    while txt:
        seq_len = int(re.search(sl_re, txt).group().lstrip().rstrip())
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


def get_hit_n(x):
    return int(re.split('[|:]', x.source.id)[1])


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
    max_l = max([len(i) for i in [str(qs), str(qe), str(ss), str(se)]])
    ad = linelength - (2 * max_l + 10 + 4)
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

            f_string = f_string + 'Query' + delim + delim * (max_l - len(str(qs))) + str(qs) + \
                       ' ' + q + ' ' + \
                       str(qe) + delim * (max_l - len(str(qe))) + '\n' + \
                       delim * (max_l + 6) + \
                       ' ' + m + ' ' + \
                       delim * max_l + '\n' + \
                       'Sbjct' + delim + delim * (max_l - len(str(ss))) + str(ss) + \
                       ' ' + s + ' ' + \
                       str(se) + delim * (max_l - len(str(se))) + '\n\n'

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
        raise exceptions.BlastFormatException('Failed to format BLAST HSP properly.')
    return f_string


def add_loc_to_description(analyzed_hits):
    for hit in analyzed_hits.hits:
        d2a = '{}-{}'.format(hit.best_start, hit.best_end)
        # hit.source.description += d2a
        hit.extension.description += d2a


def match_acc(hit_id, blast_regexp):
    # get a index name to the blastdb
    hname = re.search(blast_regexp, hit_id)

    if hname and 'pdb|' in hit_id:
        # case when sequence is from pdb e.g. gi|1276810701|pdb|5VT0|R
        # but it is available in the database under 5VT0_R
        bdb_accession = hname.group().replace('|', '_')
    elif hname is not None:
        bdb_accession = hname.group()
    else:
        msg = 'Provided regexp returned no result for {}.\n' \
              'Please provide regular expression capturing sequence id.'.format(hit_id)
        raise exceptions.AccessionMatchException(msg)
    return bdb_accession