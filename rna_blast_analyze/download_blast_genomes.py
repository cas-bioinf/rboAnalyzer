#!/usr/bin/env python3
from Bio import Entrez
import argparse
from subprocess import call
import sys
import os
import re
from rna_blast_analyze.BR_core.BA_support import blast_in
from rna_blast_analyze.BR_core.parse_accession import accession_regex


def parser():
    p = argparse.ArgumentParser(
        description='Download sequences from BLAST output. The script will download records in BLAST output from '
                    'NCBI nucl database using entrez. The BLASTdb needed for running rna_blast_analyze pipeline '
                    'will also be created.',
        usage='download_seqs_from_blast -email YOUR@EMAIL -blast_in BLAST_FILE'
    )
    p.add_argument(
        '-email',
        type=str,
        required=True,
        help='Enter valid email so NCBI staff could contact you if needed.'
    )
    p.add_argument(
        '-blast_in',
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
        help='output fasta file'
    )
    p.add_argument(
        '--blast_regexp',
        default=None,
        help='BLAST regexp. Must capture the accession numbers.'
    )
    return p.parse_args()


def fetch_by_accession(aclist):
    with Entrez.efetch(db='nucleotide', id=','.join(aclist), rettype='fasta', retmode='text') as h:
        while True:
            data = h.read(4096)
            if not data:
                break
            yield data


def main():
    args = parser()
    Entrez.email = args.email

    accessions = set()
    for b in blast_in(args.blast_in):
        for a in b.alignments:
            if hasattr(a, 'accession'):
                accessions.add(a.accession)
            else:
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

    if os.path.exists(args.out):
        with open(args.out, 'r+') as h:
            known_accs = []
            p = 0
            while True:
                l = h.readline()
                if l == '':
                    break
                if l[0] == '>':
                    p = h.tell() - len(l)
                    known_accs.append(l.split()[0][1:])
            known_accs.pop()
            # remove already retrieved accessions
            accessions - set(known_accs)
            h.seek(p)
            for chunk in fetch_by_accession(list(accessions)):
                if '>' in chunk:
                    toprint = chunk.split('>')[1].split()[0]
                    sys.stdout.write('\rRetrieving: ' + toprint)
                h.write(chunk)
    else:
        with open(args.out, 'w') as h:
            for chunk in fetch_by_accession(list(accessions)):
                if '>' in chunk:
                    toprint = chunk.split('>')[1].split()[0]
                    sys.stdout.write('\rRetrieving: ' + toprint)
                h.write(chunk)

    # convert to blastdb
    cmd = [
        'makeblastdb',
        '-in', args.out,
        '-input_type', 'fasta',
        '-dbtype', 'nucl',
        '-out', args.out + '.bdb',
        '-parse_seqids'
    ]
    r = call(cmd)
    if r:
        print(cmd)
        raise ChildProcessError('Call to makeblastdb failed.')


if __name__ == '__main__':
    main()
