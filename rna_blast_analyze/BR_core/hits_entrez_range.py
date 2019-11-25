from Bio import Entrez


def fetch_accession_range(acc: str, start: int, stop: int):
    with Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text', seq_start=start, seq_stop=stop) as h:
        return h.read()
