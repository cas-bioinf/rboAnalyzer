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
            'Download whole sequences whose accessions are present in BLAST output from NCBI.'
            ' The script will download records in BLAST output from NCBI nucleotide database using entrez.'
            ' The BLAST database needed for running rboAnalyzer will also be created.'
            ' If the --out FILE exists we assume that it is FASTA file - it will be scanned for fasta header IDs in form of accession.version.'
            ' Required records present in the file will not be downloaded.'
            ' New records will be appended to the file. The BLAST db will overwrite old BLAST db of same name.'
        ),
        usage='genomes_from_blast -email YOUR@EMAIL -blast_in BLAST_FILE'
    )
    p.add_argument(
        '-email',
        type=str,
        required=True,
        help='Enter valid email so NCBI staff could contact you if needed.'
    )
    p.add_argument(
        '-in',
        '--blast_in',
        type=str,
        required=True,
        help='BLAST output file (txt or xml) created with NCBI BLAST agains one of their databases.'
             ' Accession numbers in BLAST output must correspond to records in NCBI nucl database'
             ' for this script to work.'
    )
    p.add_argument(
        '--out',
        type=str,
        default=None,
        help=(
            'Output fasta file. Default "[BLAST IN]+.fasta". The BLAST db name is then --out FILE + ".bdb".'
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


def check_accessions(file):
    with open(file, 'r+') as h:
        known_accs = set()
        while True:
            l = h.readline()
            if l == '':
                break
            if l[0] == '>':
                known_accs.add(l.split()[0][1:])
        return known_accs


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
                print(blast_regexp)
                raise RuntimeError('provided regexp returned no result for %s,'
                                   ' please provide regexp valid even for this name', a)
            bdb_accession = hname.group()
            accessions.add(bdb_accession)

    if args.out:
        of = args.out
    else:
        of = args.blast_in + '.fasta'

    sys.stdout.write('Total of {} sequences are required.\n'.format(len(accessions)))

    c = 0

    while True:
        if os.path.exists(of):
            known_accs = check_accessions(of)
            if c >= args.retry and len(accessions - known_accs) != 0:
                sys.stdout.write(
                    "Failed to retrieve all sequences. BLAST db will still be created. "
                    "Option --skipp_missing might be required to rna_blast_analyze pipeline to run. "
                    "This may be caused by connection problems or by old BLAST search with invalid accessions."
                )
                break
            elif len(accessions - known_accs) == 0:
                # nothing to retrieve & check successful
                break

            with open(of, 'r+') as h:
                # seek to EOF
                h.seek(0, 2)
                # remove already retrieved accessions
                accessions -= known_accs
                sys.stdout.write('\n{} sequences to retrieve.\n'.format(len(accessions)))

                retrieve_and_write(h, accessions)
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
        '-out', of + '.bdb',
        '-parse_seqids'
    ]
    r = call(cmd)
    if r:
        print(cmd)
        raise ChildProcessError('Call to makeblastdb failed.')


def retrieve_and_write(handle, accessions):
    try:
        formater = '\rRetrieving: {:' + str(max([len(i) for i in accessions])) + '}'

        for chunk in fetch_by_accession(list(accessions)):
            if '>' in chunk:
                for acc_l in chunk.split('>')[1:]:
                    toprint = acc_l.split()[0]
                    sys.stdout.write(formater.format(toprint))
            handle.write(chunk)
        return

    except HTTPException:
        sys.stdout.write(
            'Connection problem encountered. Rolling back last retrieved record as it may be incomplete.\n'
        )
        rollback_fasta(handle)

    except IOError:
        sys.stdout.write(
            'Network problem encountered. Rolling back last retrieved record as it may be incomplete.\n'
        )
        rollback_fasta(handle)


def rollback_fasta(handle):
    handle.seek(0, 2)
    curr_size = handle.tell()

    n = 1
    while 1024 * n < curr_size:
        handle.seek(curr_size - 1024 * n)

        fasta_start = []
        while True:
            line = handle.readline()
            if not line:
                break

            if line[0] == '>':
                fasta_start.append(handle.tell() - len(line))

        if any(fasta_start):
            break
        n += 1
    else:
        raise EOFError(
            'The fasta file does not contain any fasta header (line starting with ">"). Cannot rollback.'
        )

    handle.truncate(fasta_start[-1])
    sys.stdout.write('Rollback successful.\n')

if __name__ == '__main__':
    main()
