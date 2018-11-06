import os
import re
import logging
from subprocess import Popen, PIPE
from tempfile import mkstemp

from Bio import SeqIO

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.parse_accession import accession_regex

ml = logging.getLogger(__name__)


def expand_hits(hits, blast_db, query_length, extra=0, keep_all=0, blast_regexp=None, skip_missing=False, msgs=None):
    """takes list of blast.HSP objects as first argument and
    path to local blast database as second argument
    then it uses blastdbcmd from blast+ installation to obtain desired sequence
    Two temporary files are used in this call and are deleted at final stage
    :return list of SeqRecord objects (parsed fasta file)
    """
    ml.info('Retrieving sequence neighborhoods for blast hits.')
    ml.debug(fname())
    ext_try = ['nsd', 'nhr', 'nog', 'nsi', 'nin']
    if not os.path.isfile(blast_db + '.nal'):
        for ext in ext_try:
            if not os.path.isfile(blast_db + '.' + ext):
                raise FileNotFoundError('expected file : "' + blast_db + '.' + ext + '" was not found')

    fd, temp_filename = mkstemp(prefix='rba_', suffix='_25')
    fdb, blast_tempfile = mkstemp(prefix='rba_', suffix='_26')
    os.close(fdb)
    exp_hits = []
    strand = []
    loc = []
    temp_file = os.fdopen(fd, 'w')
    for index, hit in enumerate(hits):
        # +1 here because blastdbcmd counts sequences from 1
        if hit[1].sbjct_end < hit[1].sbjct_start:
            # this is hit to minus strand
            start = hit[1].sbjct_end - _positive_index(query_length - hit[1].query_end) - extra
            end = hit[1].sbjct_start + hit[1].query_start + extra - 1
            strand.append(-1)
            d = {'query_start': hit[1].sbjct_end, 'query_end': hit[1].sbjct_start,
                 'extended_start': hit[1].sbjct_end - _positive_index(query_length - hit[1].query_end),
                 'extended_end': hit[1].sbjct_start + hit[1].query_start - 1,
                 'strand': -1}
        else:
            # this is hit to plus strand
            start = hit[1].sbjct_start - hit[1].query_start + 1 - extra
            end = hit[1].sbjct_end + _positive_index(query_length - hit[1].query_end) + extra
            strand.append(1)
            d = {'query_start': hit[1].sbjct_start, 'query_end': hit[1].sbjct_end,
                 'extended_start': hit[1].sbjct_start - hit[1].query_start + 1,
                 'extended_end': hit[1].sbjct_end + _positive_index(query_length - hit[1].query_end),
                 'strand': 1}

        # ====== information about possible trim ======
        # assume ok
        d['trimmed_ss'] = False
        d['trimmed_se'] = False
        d['trimmed_es'] = False
        d['trimmed_ee'] = False

        d['super_start'] = start
        d['super_end'] = end

        if start < 1:
            start = 1                    # index from which sequence should be retrieved from the db
            d['trimmed_ss'] = True

        # repair possible extended start violation
        if d['extended_start'] < 1:
            d['trimmed_es'] = True

        # add blast record
        d['blast'] = hit

        # get a index name to the blastdb
        if blast_regexp is None:
            blast_regexp = accession_regex
        hname = re.search(blast_regexp, hit[0])
        if not hname:
            print(blast_regexp)
            os.remove(temp_filename)
            os.remove(blast_tempfile)
            raise RuntimeError('provided regexp returned no result for %s,'
                               ' please provide regexp valid even for this name', hit[0])
        bdb_accession = hname.group()
        d['blast'][0] = bdb_accession
        loc.append(d)

        temp_file.write(bdb_accession + ' ' + '-region ' + str(start) + '-' + str(end) + '\n')

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
    ml.debug(cmd)

    try:
        pcall = Popen(
            cmd,
            stdout=PIPE,
            stderr=PIPE,
            universal_newlines=True
        )
        out, errs = pcall.communicate()
    except FileNotFoundError:
        msgfail = 'Unable to run blastdbcmd command, please check its availability.'
        ml.error(msgfail)
        ml.error(cmd)
        os.remove(temp_filename)
        os.remove(blast_tempfile)
        raise ChildProcessError(msgfail)

    # inspect the stdout (the blastdbcmd returns exit code 1 even if only one sequence is missing)

    msgfail = 'Incomplete database. Some sequences not found in database.'
    if ml.getEffectiveLevel() == 10:
        msgfail += ' Details: ' + errs

    if errs and not skip_missing:
        ml.error(msgfail)
        os.remove(temp_filename)
        os.remove(blast_tempfile)
        raise LookupError(msgfail)
    elif errs and skip_missing:
        ml.warning(msgfail)

    requested_ids = {l['blast'][0] for l in loc}
    obtained_ids = set()

    with open(blast_tempfile, 'r') as bf:
        index = 0
        for parsed_record in SeqIO.parse(bf, 'fasta'):
            obtained_ids.add(parsed_record.id)
            index = _get_correct_blast_hit(parsed_record, loc, index, skip=skip_missing)

            parsed_record.annotations = loc[index]
            parsed_record.annotations['msgs'] = []
            # add uid to ensure that all hits are unique
            parsed_record.id = 'uid:' + str(index) + '|' + parsed_record.id

            if loc[index]['trimmed_ss']:
                if loc[index]['super_start'] + len(parsed_record.seq) < loc[index]['super_end'] + loc[index]['super_start']:
                    parsed_record.annotations['trimmed_se'] = True
                    if loc[index]['trimmed_es']:
                        if len(parsed_record.seq) < loc[index]['extended_end'] + loc[index]['extended_start']:
                            parsed_record.annotations['trimmed_ee'] = True
                    else:
                        if len(parsed_record.seq) < loc[index]['extended_end'] - loc[index]['extended_start']:
                            parsed_record.annotations['trimmed_ee'] = True
            else:
                if loc[index]['super_start'] + len(parsed_record.seq) - 1 < loc[index]['super_end']:
                    parsed_record.annotations['trimmed_se'] = True
                    if loc[index]['super_start'] + len(parsed_record.seq) - 1 < loc[index]['extended_end']:
                        parsed_record.annotations['trimmed_ee'] = True

            msgsub = '{}: Sequence cannot be extended sufficiently'.format(parsed_record.id)
            if parsed_record.annotations['trimmed_ss']:
                msgwarn = msgsub + '. Missing {} nt upstream in the genome.'.format(parsed_record.annotations['super_start'])
                parsed_record.annotations['msgs'].append(msgwarn)
                ml.warning(msgwarn)
            if parsed_record.annotations['trimmed_se']:
                msgwarn = msgsub + '. Missing nt downstream in the genome.'.format(parsed_record.id)
                parsed_record.annotations['msgs'].append(msgwarn)
                ml.warning(msgwarn)
            if parsed_record.annotations['trimmed_es']:
                msgwarn = msgsub + ' by unaligned portion of query. THIS IS PROBABLY FRAGMENT!'
                msgwarn += ' Trimmed upstream.'
                parsed_record.annotations['msgs'].append(msgwarn)
                ml.warning(msgwarn)
            if parsed_record.annotations['trimmed_ee']:
                msgwarn = msgsub + ' by unalined portion of query. THIS IS PROBABLY FRAGMENT!'
                msgwarn += ' Trimmed downstream.'
                parsed_record.annotations['msgs'].append(msgwarn)
                ml.warning(msgwarn)

            exp_hits.append(parsed_record)
            index += 1

    if not keep_all:
        os.remove(temp_filename)
        os.remove(blast_tempfile)

    if len(requested_ids - obtained_ids) != 0:
        msgs.append('Incomplete database. Sequences with following ids were not found:')
        for m in requested_ids - obtained_ids:
            msgs.append(m)

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


def _get_correct_blast_hit(db_rec, loc, index, skip):
    if loc[index]['blast'][0] == db_rec.id:
        return index
    elif loc[index]['blast'][0] != db_rec.id and skip:
        msgwarn = 'Sequence {} not found in provided db. Skipping.'.format(loc[index]['blast'][0])
        ml.warning(msgwarn)
        return _get_correct_blast_hit(db_rec, loc, index + 1, skip)
    else:
        msgerror = 'Sequence {} not found in provided db. ' \
                   'Please provide correct database or give "--skip_missing" flag.'.format(
            loc[index]['blast'][0]
        )
        ml.error(msgerror)
        raise LookupError(msgerror)
