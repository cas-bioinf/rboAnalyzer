#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
from Bio import Entrez
import argparse
from subprocess import call
import sys
import os
import re
from rna_blast_analyze.BR_core.BA_support import blast_in
from rna_blast_analyze.BR_core.parse_accession import accession_regex
from http.client import HTTPException


def parser():
    p = argparse.ArgumentParser(
        description=(
            'Download whole nucleotide sequences whose accessions are present in BLAST output from NCBI.'
            ' The script will download records in BLAST output from NCBI nucleotide database using entrez.'
            ' And create the BLAST database needed for running rboAnalyzer.'
            ' The intermediate FASTA file (appending ".fasta" to the BLAST database name) is saved as sequences can be added to it.'
            ' If the intermediate FASTA file exists, it will be scanned for fasta header IDs in form of accession.version.'
            ' Required records which present in the intermediate file will not be downloaded.'
            ' New records will be appended to the file. The BLAST db will overwrite old BLAST db of same name.'
        ),
        usage='genomes_from_blast -e YOUR@EMAIL -in BLAST_FILE'
    )
    p.add_argument(
        '-e',
        '--email',
        type=str,
        required=True,
        help='Enter valid email so NCBI staff could contact you if needed.'
    )
    p.add_argument(
        '-in',
        '--blast_in',
        type=str,
        required=True,
        help='BLAST output file (txt or xml) created with NCBI BLAST against one of their databases.'
             ' Accession numbers in BLAST output must correspond to records in NCBI nucl database'
             ' for this script to work.'
    )
    p.add_argument(
        '-o',
        '--out',
        type=str,
        default=None,
        help=(
            'Output BLAST database file. Default "[BLAST IN]+.bdb".'
            ' The intermediate FASTA file name is then --out FILE + ".fasta".'
        )
    )
    p.add_argument(
        '--blast_regexp',
        default=None,
        help='BLAST regexp. Must capture the accession numbers.'
    )
    p.add_argument(
        '--retry',
        default=10,
        type=int,
        help='Number of retries to download required sequences.'
    )
    return p.parse_args()


def fetch_by_accession(aclist):
    with Entrez.efetch(db='nucleotide', id=','.join(aclist), rettype='fasta', retmode='text') as h:
        while True:
            data = h.read(4096)
            if not data:
                break
            yield data


def check_accessions(handle):
    handle.seek(0)
    known_accs = []
    last_pos = None
    while True:
        line = handle.readline()
        if line == '':
            break
        if line[0] == '>':
            last_pos = handle.tell() - len(line)
            known_accs.append(line.split()[0][1:])
    return known_accs, last_pos


def main():
    args = parser()
    Entrez.email = args.email

    accessions = set()
    for b in blast_in(args.blast_in):
        for a in b.alignments:
            if args.blast_regexp is None:
                blast_regexp = accession_regex
            else:
                blast_regexp = args.blast_regexp
            hname = re.search(blast_regexp, a.hit_id)
            if not hname:
                msg = 'Provided regexp returned no result for {}, ' \
                      'please provide regexp valid even for this name'.format(a)
                print(msg)
                print(blast_regexp)
                sys.exit(1)

            bdb_accession = hname.group()
            accessions.add(bdb_accession)

    if args.out:
        od = args.out
        of = args.out + '.fasta'
    else:
        od = args.blast_in + '.bdb'
        of = args.blast_in + '.bdb.fasta'

    sys.stdout.write('Total of {} sequences are required.\n'.format(len(accessions)))

    c = 0

    while True:
        if os.path.exists(of):
            # check completeness
            with open(of, 'a+') as fh:
                # rollback in necessary (incomplete file)
                known_accs_list = _rollback(fh)

                known_accs = set(known_accs_list)
                if c >= args.retry and len(accessions - known_accs) != 0:
                    sys.stdout.write(
                        "Failed to retrieve all sequences. BLAST db will still be created.\n"
                        "You can rerun the command for another try at downloading all required sequences.\n"
                        "Option --skipp_missing might be required to rna_blast_analyze pipeline to run. "
                        "This may be caused by connection problems or by old BLAST search with invalid accessions."
                    )
                    break
                elif len(accessions - known_accs) == 0:
                    # nothing to retrieve & check successful
                    break

                # seek to EOF
                fh.seek(0, 2)
                # remove already retrieved accessions
                accessions -= known_accs
                sys.stdout.write('\n{} sequences to retrieve.\n'.format(len(accessions)))

                retrieve_and_write(fh, accessions)
        else:
            with open(of, 'w+') as h:
                retrieve_and_write(h, accessions)

        c += 1

    # convert to blastdb
    cmd = [
        'makeblastdb',
        '-in', of,
        '-input_type', 'fasta',
        '-dbtype', 'nucl',
        '-out', od,
        '-parse_seqids'
    ]
    r = call(cmd)
    if r:
        print('Call to "makeblastdb" failed. Please check that the "makeblastdb" is in your PATH.')
        sys.exit(1)


def assert_complete_file(handle):
    handle.seek(0, 2)
    cp = handle.tell()
    if cp > 3:
        handle.seek(cp - 3, 0)
        cs = handle.read()
        if cs == "\n\n\n":
            return True
    # return False otherwise
    return False


def _rollback(handle):
    # check completeness
    # if not complete record (indicated by "\n\n\n") truncate the file to last entry
    known_accs_list, last_pos = check_accessions(handle)
    if not assert_complete_file(handle):
        sys.stdout.write(
            "Output file appear to be incomplete. Removing last potentially incomplete record. "
            "It will be downloaded again if necessary."
        )
        if last_pos is None:
            # this appears to contain no records
            handle.seek(0)
        else:
            handle.seek(last_pos, 0)
        handle.truncate()
        handle.write('\n\n\n')
        known_accs_list = known_accs_list[:-1]
    return known_accs_list


def retrieve_and_write(handle, accessions):
    # seek to eof
    handle.seek(0, 2)
    try:
        formater = '\rRetrieving: {:' + str(max([len(i) for i in accessions])) + '}'

        for chunk in fetch_by_accession(list(accessions)):
            if '>' in chunk:
                for acc_l in chunk.split('>')[1:]:
                    toprint = acc_l.split()[0]
                    sys.stdout.write(formater.format(toprint))
            handle.write(chunk)
        # write 3 linebreaks denoting that retrieval was successful
        handle.write('\n\n\n')
        return

    except HTTPException:
        sys.stdout.write(
            'Connection problem encountered. Rolling back last retrieved record as it may be incomplete.\n'
        )
        _rollback(handle)

    except IOError:
        sys.stdout.write(
            'Network problem encountered. Rolling back last retrieved record as it may be incomplete.\n'
        )
        _rollback(handle)

    except KeyboardInterrupt:
        sys.stdout.write(
            'Script interrupted. Trying to rollback to last complete entry.\n'
        )
        _rollback(handle)

    except Exception as e:
        sys.stdout.write(
            'Unexpected error when retrieving sequences. Trying to rollback to last complete entry.\n'
            'Error message: {}'.format(str(e))
        )
        _rollback(handle)


if __name__ == '__main__':
    try:
        main()
        sys.exit(0)
    except Exception as e:
        print('Something went wrong. Please contact the developers.')
        print(str(e))
        sys.exit(1)
