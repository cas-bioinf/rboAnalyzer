# list of prefixes and rules for parsing accession numbers

import re

# ===== UniProt accession =====
# uniprot from: https://www.uniprot.org/help/accession_numbers
uniprot = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'

# ===== GenBank classic accession =====
# from: https://www.ncbi.nlm.nih.gov/Sequin/acc.html
#
# -- Basic rules --
# Nucleotide:	1 letter + 5 numerals OR 2 letters + 6 numerals
# Protein:	3 letters + 5 numerals
# WGS:	4 letters + 2 numerals for WGS assembly version + 6-8 numerals
# MGA:	5 letters + 7 numerals
#
# we need Nucleotide and WGS
genbank_nucl = '[A-Z][0-9]{5}\.[0-9]+|[A-Z]{2}[0-9]{6}\.[0-9]+'
genbank_prot = '[A-Z]{3}[0-9]{5}\.[0-9]+'
genbank_mga = '[A-Z]{5}[0-9]{7}\.[0-9]+'
genbank_wgs = '[A-Z]{4}[0-9]{2,9}\.[0-9]+'

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
refseq_re = '|'.join([i + '[0-9]+\.[0-9]+|' + genbank_nucl + '|' + genbank_wgs for i in refseq])

known_acc_formats = [genbank_nucl, genbank_wgs, refseq_re, genbank_mga, uniprot, genbank_prot]
accession_regex = re.compile('|'.join(known_acc_formats))
