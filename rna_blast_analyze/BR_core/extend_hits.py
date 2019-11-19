import os
import logging
from subprocess import Popen, PIPE
from tempfile import mkstemp, TemporaryFile, gettempdir
import sys

from Bio import SeqIO

from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.fname import fname
from rna_blast_analyze.BR_core.BA_support import remove_files_with_try, match_acc, remove_one_file_with_try
from rna_blast_analyze.BR_core.load_from_bgzip import load_genome
from rna_blast_analyze.BR_core import exceptions
from Bio import Entrez
from http.client import HTTPException
from math import floor

ml = logging.getLogger('rboAnalyzer')


def fetch_accession_range(acc: str, start: int, stop: int):
    with Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text', seq_start=start, seq_stop=stop) as h, \
            TemporaryFile(mode='w+') as temp:
        temp.write(h.read())
        temp.seek(0)
        return SeqIO.read(temp, format='fasta')


def expand_hits(hits, blast_db, query_length, extra=0, blast_regexp=None, skip_missing=False, msgs=None):
    """takes list of blast.HSP objects as first argument and
    path to local blast database as second argument
    then it uses blastdbcmd from blast+ installation to obtain desired sequence
    Two temporary files are used in this call and are deleted at final stage
    :return list of SeqRecord objects (parsed fasta file)
    """
    ml.info('Retrieving sequence neighborhoods for blast hits.')
    ml.debug(fname())

    fd, temp_filename = mkstemp(prefix='rba_', suffix='_25', dir=CONFIG.tmpdir)
    fdb, blast_tempfile = mkstemp(prefix='rba_', suffix='_26', dir=CONFIG.tmpdir)
    os.close(fdb)
    exp_hits = []
    strand = []
    loc = []
    try:
        with os.fdopen(fd, 'w') as temp_file:
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

                try:
                    bdb_accession = match_acc(hit[0], blast_regexp)
                except exceptions.AccessionMatchException as e:
                    remove_files_with_try([temp_filename, blast_tempfile])
                    raise e

                d['blast'][0] = bdb_accession
                loc.append(d)

                temp_file.write(bdb_accession + ' ' + '-region ' + str(start) + '-' + str(end) + '\n')
    except RuntimeError as e:
        ml.error(str(e))
        sys.exit(1)

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
        remove_files_with_try([temp_filename, blast_tempfile])
        sys.exit(1)

    # inspect the stdout (the blastdbcmd returns exit code 1 even if only one sequence is missing)
    msgfail = 'Incomplete database. Some sequences not found in database.'
    msgfail += ' Details: ' + errs

    if errs and not skip_missing:
        ml.error(msgfail)
        remove_files_with_try([temp_filename, blast_tempfile])
        sys.exit(1)

    elif errs and skip_missing:
        ml.warning(msgfail)

    requested_ids = {l['blast'][0] for l in loc}
    obtained_ids = set()

    with open(blast_tempfile, 'r') as bf:
        index = 0
        for parsed_record in SeqIO.parse(bf, 'fasta'):
            record_id = parsed_record.id.split(':')[0]

            if parsed_record.description.startswith(parsed_record.id):
                parsed_record.description = parsed_record.description[len(parsed_record.id):].strip()

            obtained_ids.add(record_id)
            index = _get_correct_blast_hit(parsed_record, loc, index, skip=skip_missing)

            parsed_record.annotations = loc[index]
            parsed_record.annotations['msgs'] = []
            # add uid to ensure that all hits are unique
            parsed_record.id = 'uid:' + str(index) + '|' + record_id

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

            msgsub = '{}: Sequence cannot be extended sufficiently'.format(record_id)
            if parsed_record.annotations['trimmed_ss']:
                msgwarn = msgsub + '. Missing {} nt upstream in the genome.'.format(parsed_record.annotations['super_start'])
                parsed_record.annotations['msgs'].append(msgwarn)
                ml.warning(msgwarn)
            if parsed_record.annotations['trimmed_se']:
                msgwarn = msgsub + '. Missing nt downstream in the genome.'.format(record_id)
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

    remove_files_with_try(
        [temp_filename, blast_tempfile]
    )

    if len(requested_ids - obtained_ids) != 0:
        msgs.append('Incomplete database. Sequences with following ids were not found:')
        for m in requested_ids - obtained_ids:
            msgs.append(m)

    return exp_hits, strand


def expand_hits_from_fasta(hits, database, query_length, extra=0, blast_regexp=None, skip_missing=False, msgs=None, format='fasta', entrez_email=None, blast_input_file=None):
    """takes list of blast.HSP objects and return extended sequences
    :return list of SeqRecord objects (parsed fasta file)
    """
    ml.info('Retrieving sequence neighborhoods for blast hits.')
    ml.debug(fname())

    if CONFIG.tmpdir is None:
        temp_entrez_file = os.path.join(
            gettempdir(), os.path.basename(blast_input_file + '.r-temp_entrez')
        )
    else:
        temp_entrez_file = os.path.join(
            CONFIG.tmpdir, os.path.basename(blast_input_file + '.r-temp_entrez')
        )

    if format == 'entrez':
        try:
            known_seq_index = SeqIO.index(temp_entrez_file, format='fasta')
            ml.info("File {} loaded.".format(temp_entrez_file))
            if len(known_seq_index) == 0:
                remove_one_file_with_try(temp_entrez_file)
        except FileNotFoundError:
            # ignore that we don't have that file (usual)
            known_seq_index = {}
        except Exception as e:
            ml.info("Could not load the temporary file {}.".format(temp_entrez_file))
            known_seq_index = {}
    else:
        known_seq_index = {}

    exp_hits = []
    strand = []
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

        try:
            bdb_accession = match_acc(hit[0], blast_regexp)
        except exceptions.AccessionMatchException as e:
            raise e

        d['blast'][0] = bdb_accession

        # read from file
        if format in ['fasta', 'gb']:
            ff = os.path.join(database, bdb_accession)
            if not os.path.isfile(ff):
                if skip_missing:
                    msgwarn = 'Sequence {} not found in provided db. Skipping.'.format(bdb_accession)
                    msgs.append(msgwarn)
                    ml.warning(msgwarn)
                else:
                    msgerror = 'Sequence {} not found in provided db. ' \
                               'Please provide correct database or give "--skip_missing" flag.'.format(
                        bdb_accession
                    )
                    ml.error(msgerror)
                    raise LookupError(msgerror)

            with open(ff, 'r') as handle:
                ext_seq = next(SeqIO.parse(handle, format=format))
                parsed_record = ext_seq[start - 1:end]

        elif format == 'server':
            # only used when server
            parsed_record = load_genome(database, bdb_accession, start - 1, end)
        elif format == 'entrez':
            if index == 0:
                sys.stdout.write('STATUS: Downloading required sequences from NCBI with entrez.\n')
                sys.stdout.write('{:3d}% {}'.format(index, bdb_accession))
            else:
                sys.stdout.write('\r{:3d}% {}'.format(floor(index * 100 / len(hits)), bdb_accession))

            seq_id = '{}:{}-{}'.format(bdb_accession, start, end)
            if seq_id not in known_seq_index:
                try:
                    Entrez.email = entrez_email
                    parsed_record = fetch_accession_range(bdb_accession, start, end)

                    with open(temp_entrez_file, 'a') as tmpf:
                        SeqIO.write([parsed_record], tmpf, format='fasta')

                    if index == len(hits) - 1:
                        sys.stdout.write('\r Done.\n')
                except HTTPException as e:
                    msg = 'HTTP exception encountered: {}' \
                          'Please check your internet connection and availability of NCBI ENTREZ web services.\n ' \
                          'Also check that the requested accession number "{}" is available in NCBI "nucleotide" database.'.format(
                        e, bdb_accession)
                    if skip_missing:
                        ml.warning(msg)
                        continue
                    else:
                        ml.error(msg)
                        sys.exit(1)
                except ValueError:
                    msg = 'Received malformed fasta file. ' \
                          'Please check that requested accession number "{}" is available in NCBI nucleotide database.'.format(
                        bdb_accession)
                    if skip_missing:
                        ml.warning(msg)
                        continue
                    else:
                        ml.error(msg)
                        sys.exit(1)
            else:
                parsed_record = known_seq_index[seq_id]
        else:
            raise NotImplementedError

        record_id = parsed_record.id.split(':')[0]

        if parsed_record.description.startswith(parsed_record.id):
            parsed_record.description = parsed_record.description[len(parsed_record.id):].strip()

        parsed_record.annotations = d
        parsed_record.annotations['msgs'] = []
        # add uid to ensure that all hits are unique
        parsed_record.id = 'uid:' + str(index) + '|' + record_id

        if d['trimmed_ss']:
            if d['super_start'] + len(parsed_record.seq) < d['super_end'] + d['super_start']:
                parsed_record.annotations['trimmed_se'] = True
                if d['trimmed_es']:
                    if len(parsed_record.seq) < d['extended_end'] + d['extended_start']:
                        parsed_record.annotations['trimmed_ee'] = True
                else:
                    if len(parsed_record.seq) < d['extended_end'] - d['extended_start']:
                        parsed_record.annotations['trimmed_ee'] = True
        else:
            if d['super_start'] + len(parsed_record.seq) - 1 < d['super_end']:
                parsed_record.annotations['trimmed_se'] = True
                if d['super_start'] + len(parsed_record.seq) - 1 < d['extended_end']:
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

    # ==== Remove the entrez tempfile =====
    # here we have all the sequences and we can safely delete the tempfile
    if database == 'entrez':
        remove_one_file_with_try(temp_entrez_file)
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
