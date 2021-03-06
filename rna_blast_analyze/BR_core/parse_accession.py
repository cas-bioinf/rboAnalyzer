import re
# list of prefixes and rules for parsing accession numbers

# ===== UniProt accession =====
# uniprot from: https://www.uniprot.org/help/accession_numbers
# uniprot = ['[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}']

# uniprot regular expression may overlap genebank accession

# ===== GenBank classic accession =====
# from: https://www.ncbi.nlm.nih.gov/Sequin/acc.html
# update: Feb 2019
# we need Nucleotide and WGS
genbank_nucl = [
    r'[A-Z][0-9]{5}\.[0-9]+',
    r'[A-Z]{2}[0-9]{6}\.[0-9]+',
    r'[A-Z]{2}[0-9]{8}\.[0-9]+',
    ]
genbank_prot = [
    r'[A-Z]{3}[0-9]{5}\.[0-9]+',
    r'[A-Z]{3}[0-9]{7}\.[0-9]+',
    ]
genbank_mga = [
    r'[A-Z]{5}[0-9]{7}\.[0-9]+'
]
genbank_wgs = [
    r'[A-Z]{4}[0-9]{8,}\.[0-9]+',
    r'[A-Z]{6}[0-9]{9,}\.[0-9]+',
]

genbank_wgs_scafolds = [
    r'[A-Z]{4}[0-9]{2}S?[0-9]{6,}\.[0-9]+',
    r'[A-Z]{6}[0-9]{2}S?[0-9]{7,}\.[0-9]+',
]

# ===== RefSeq format =====
# refseq from: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
# they have underscore in thrid position

refseq_re = [r"[A-Z]{2}_[0-9A-Z]+\.[0-9]+"]

old = [r'ZP_[0-9]{8}\.[0-9]+', r'NS_[0-9]{6}\.[0-9]+']

pdb_on_steroids = [r"(gi\|\d+\|)?(?(1))((?<=\|pdb\|)[0-9A-Z]{4}[_|][0-9A-Za-z]{1,4}|^[0-9A-Z]{4}[_|][0-9A-Za-z]{1,4})"]

exceptions = [
    '1KPD',  # although it has chain, in the NCBI nt database it is listed without chain
    r'GPS_[0-9]{9}\.[0-9]+',  # from refseq
]

known_acc_formats = genbank_nucl + genbank_wgs + refseq_re + genbank_mga + genbank_prot + genbank_wgs_scafolds + exceptions + old + pdb_on_steroids
accession_regex = '|'.join(known_acc_formats)
compiled_accession_regex = re.compile(accession_regex)


if __name__ == '__main__':
    jxre = genbank_nucl + genbank_wgs + refseq_re + genbank_mga + genbank_prot + genbank_wgs_scafolds + exceptions + old + [r"[0-9A-Z]{4}\|[0-9A-Za-z]{1,4}"]
    jxre_txt = '|'.join(jxre)
    _javascript_xml_re = '(' + jxre_txt + ')(?=\|[A-Z0-9_]*</Hit_id>)'

    tre = genbank_nucl + genbank_wgs + refseq_re + genbank_mga + genbank_prot + genbank_wgs_scafolds + exceptions + old + [r"[0-9A-Z]{4}_[0-9A-Za-z]{1,4}"]
    tre_txt = '|'.join(tre)
    _javascript_txt_re = '^>(' + tre_txt + ')'

    print(_javascript_xml_re.__repr__())
    print(_javascript_txt_re.__repr__())

    import sys
    import gzip
    import re
    import io
    import multiprocessing
    import itertools

    if len(sys.argv) == 1:
        print('Test for known accession with ncbi accession2taxid gzipped files. Missed accession added to parsing_errors file.')

    regex = re.compile(accession_regex)
    e = open('parsing_errors', 'a+')

    def _worker(acc_v):
        if not re.fullmatch(regex, acc_v.decode().strip()):
            e.write(acc_v.decode())

    def _worker_plain(acc_v):
        if not re.fullmatch(regex, acc_v.strip()):
            e.write(acc_v)

    def grouper(n, iterable):
        iterable = iter(iterable)
        return iter(lambda: list(itertools.islice(iterable, n)), [])

    try:
        with multiprocessing.Pool() as pool:
            for file in sys.argv[1:]:
                print(file)
                if file.endswith('.gz'):
                    with gzip.open(file, 'rb') as f:
                        for chunk in grouper(100000, io.BufferedReader(f)):
                            pool.map(_worker, chunk)
                else:
                    with open(file, 'r') as f:
                        for chunk in grouper(100000, f):
                            pool.map(_worker_plain, chunk)
    finally:
        e.close()

