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
# they have unrerscore(underbar)
# AC_	Genomic	Complete genomic molecule, usually alternate assembly
# NC_	Genomic	Complete genomic molecule, usually reference assembly
# NG_	Genomic	Incomplete genomic region
# NT_	Genomic	Contig or scaffold, clone-based or WGSa
# NW_	Genomic	Contig or scaffold, primarily WGSa
# NZ_b	Genomic	Complete genomes and unfinished WGS data
# NM_	mRNA	Protein-coding transcripts (usually curated)
# NR_	RNA	Non-protein-coding transcripts
# XM_c	mRNA	Predicted model protein-coding transcript
# XR_c	RNA	Predicted model non-protein-coding transcript
# AP_	Protein	Annotated on AC_ alternate assembly
# NP_	Protein	Associated with an NM_ or NC_ accession
# YP_c	Protein	Annotated on genomic molecules without an instantiated
# transcript record
# XP_c	Protein	Predicted model, associated with an XM_ accession
# WP_	Protein	Non-redundant across multiple strains and species
refseq = [
    'AC_',
    'NC_',
    'NG_',
    'NT_',
    'NW_',
    'NZ_',
    'NM_',
    'NR_',
    'XM_',
    'XR_',
    'AP_',
    'NP_',
    'YP_',
    'XP_',
    'WP_',
]
refseq_re = [prefix + r"[0-9A-Z]+\.[0-9]+" for prefix in refseq]

pdb = ["[0-9A-Z]{4}[_|][0-9A-Z]{1,2}",]

exceptions = [
    '1KPD',  # although it has chain, in the NCBI nt database it is listed without chain
]

known_acc_formats = genbank_nucl + genbank_wgs + refseq_re + genbank_mga + genbank_prot + pdb + genbank_wgs_scafolds + exceptions
accession_regex = '|'.join(known_acc_formats)


if __name__ == '__main__':
    import sys
    import gzip
    import re
    import io

    if len(sys.argv) == 1:
        print('Test for known accession with ncbi accession2taxid gzipped files. Missed accession added to parsing_errors file.')

    regex = re.compile(accession_regex)

    for file in sys.argv[1:]:
        with gzip.open(file, 'rb') as f, open('parsing_errors', 'a+') as e:
            for line in io.BufferedReader(f):
                acc_v = line.decode().split('\t')[1]
                if not re.search(regex, acc_v):
                    e.write(acc_v + ' {}\n'.format(file))
