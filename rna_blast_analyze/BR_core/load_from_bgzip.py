import os
import pysam
import gzip
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_genome(root_dir, accession, start, end):
    ff = resolver(root_dir, accession) + '.fasta.gz'

    g = pysam.FastaFile(ff)
    extracted = g.fetch(reference=accession, start=start, end=end)
    g.close()

    ace = accession.encode()
    with gzip.open(ff) as gz:
        for line in gz:
            if line[1:].startswith(ace):
                return SeqRecord(Seq(extracted), id=accession, description=line[1:].decode().strip())

    raise RuntimeError("Failed to find accession {} in {}".format(accession, ff))


def resolver(root, accession):
    """Return path for accession

    Path for genome file. Includes subdirectories based on last 4 characters from accession.
    The subdirectories are created if not exists

    :param root: PATH where the genome database should be created
    :param accession: genome pointer based on ACCESSION.VERSION
    :return:
    """

    ac = accession.split('.')[0]
    acc = ac.replace('_', '')
    dp = os.path.join(root, acc[:2], acc[2:4], acc)
    return os.path.join(dp, accession)
