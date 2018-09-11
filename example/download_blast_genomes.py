#!/usr/bin/env python3
from Bio import Entrez
from Bio.Blast import NCBIXML
import argparse
from subprocess import call
import sys
import os


def parser():
    p = argparse.ArgumentParser(
        description='Download sequences from BLAST output',
    )
    p.add_argument(
        '-email',
        type=str,
        help='Enter valid email so ncbi staff could contact you.'
    )
    p.add_argument(
        '-blast_in',
        type=str,
    )
    p.add_argument(
        '-out',
        type=str,
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
    with open(args.blast_in, 'r') as bh:
        accessions = set()
        for b in NCBIXML.parse(bh):
            for a in b.alignments:
                if hasattr(a, 'accession'):
                    accessions.add(a.accession)
                else:
                    accessions.add(a.hit_id)

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
        r = call(
            ['makeblastdb',
             '-in', args.out,
             '-input_type', 'fasta',
             '-dbtype', 'nucl',
             '-out', args.out + '.bdb',
             '-parse_seqids']
        )


if __name__ == '__main__':
    main()
