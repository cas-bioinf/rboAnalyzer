#!/usr/bin/env bash
echo "running test installed..."
rna_blast_analyze -blast_in test_data/RF00001.blastout -blast_query test_data/RF00001.fasta -blast_db test_data/blastdb/RF00001-art.blastdb --blast_regexp "(?<=\|)[A-Z0-9]*\.?\d*$" --prediction_method rnafold --html test_output.html --json test_output.json --csv test_output.csv --pandas_dump test_output.pandas_dump
