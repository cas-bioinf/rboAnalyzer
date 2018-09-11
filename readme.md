# RNA_BLAST_ANALYZE
Pipeline for analyzing BLAST search output for non-coding RNAs (ncRNAs).

## Short description
Provided with query rna sequence and BLAST output, the pipeline will
 realign each BLAST hit with Locarna to account for more RNA sequence
 variability. The realign step is done to recover potential ncRNA
 sequence at that loci precisely even from low-scoring BLAST hits.

Sequences are analyzed with RSEARCH and its bit-score is used as
 similarity criteria.

Next sequences are predicted with one or more methods and predicted
 structures are merged with blast output to be easy to understand.

## Prequisities
* python >= 3.4
* ncbi-blast+
* locarna
* infernal
* clustalo
* muscle

For prediction:
* t-coffee (with r-coffee)
* viennarna (with refold.pl in PATH)
* centroid_homfold

Optional (some prediction methods are not avalible without these):
* RNAstructure >= 6.0 (TurboFold)
* RapidShapes - RNAshapes
* mfold >= 3.8

Data:
* local copy of the database in which analyzed query were searched
(often this will be ncbi-blast db)

## Installation
Install all prequisities and run
```
python3 setup.py install
```


## Basic Usage
```
rna_blast_analyze -blast_in BLAST_OUTPUT.txt -blast_db USED_DATABASE_PATH -blast_query BLAST_QUERY.fasta
```

## Example
Examples are provided in example directory.
#### Example 1:

Analyzing subset of NCBI blast HITs for 6S ncRNA.
1) Install rna_blast_analyze.
2) Obtain a copy of `example` directory.
3) cd to `example` directory.
4) Now you need to obtain a copy of the BLAST database.
This is best done by using the `download_blast_genomes.py` by calling: 
    ```
    python3 download_blast_genome.py -email YOUR_EMAIL -blast_in 6S_short.xml -out 6S.fasta
    ```
    This command will download all needed genomes and create blast database for you (if `makeblastdb` command is avalible).
    The EMAIL is needed so the NCBI would contact you in case of missuse of their resources.
    The 6S.fasta file is not needed if the blast database was created successfully and you can delete it.
    
    The blast database with name 6S.fasta.bdb was created for you if everything was successful.
    You will need in next step.

5) Now you can run the pipeline itself:
    ```
    rna_blast_analyze -blast_in 6S_short.xml -blast_query 6S_query.fasta -blast_db 6S.fasta.bdb -html 6S_out.html
    ```
6) The output is single html file. On open the NCBI genome viewer will fetch data to render the genomic loci of the hit.

### Notes
This is Beta version. Especially the default parameters for prediction can change.

## References
- ViennaRNA: Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
 ViennaRNA Package 2.0. <https://doi.org/10.1186/1748-7188-6-26>.
 We include copy of refold.pl script with Viennarna license here for convinience. [website][https://www.tbi.univie.ac.at/RNA/]
- Locarna: Sebastian Will, Tejal Joshi, Ivo L. Hofacker, Peter F. Stadler, and Rolf Backofen.
LocARNA-P: Accurate boundary prediction and improved detection of structural RNAs
RNA, 18 no. 5, pp. 900-14, 2012. <https://doi.org/10.1261/rna.029041.111>, [website][http://rna.informatik.uni-freiburg.de/LocARNA/Input.jsp]
- NCBI BLAST: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008)
 "BLAST+: architecture and applications." BMC Bioinformatics 10:421. <https://doi.org/10.1186/1471-2105-10-421>, [website][https://blast.ncbi.nlm.nih.gov/Blast.cgi]
- RSEARCH: Finding Homologs of Single Structured RNA Sequences.
 R. J. Klein, S. R. Eddy. BMC Bioinformatics, 4:44, 2003. <https://doi.org/10.1186/1471-2105-4-44>, [website][http://eddylab.org/software.html#rsearch]
- Infernal: Infernal 1.1: 100-fold Faster RNA Homology Searches.
 E. P. Nawrocki, S. R. Eddy. Bioinformatics, 29:2933-2935, 2013. <https://doi.org/10.1093/bioinformatics/btt509>, [website][http://eddylab.org/infernal/]
- Clustal omega: Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011).
 Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.
  Molecular Systems Biology 7:539. <https://doi.org/10.1038/msb.2011.75>, [website][http://www.clustal.org/omega/]
- MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput
 Nucleic Acids Res. 32(5):1792-1797. <https://doi.org/10.1093/nar/gkh340>, [website][https://www.drive5.com/muscle/]
- T-coffee (r-coffee): Wilm, A., Higgins, D.G., Notredame, C.
 R-Coffee: a method for multiple alignment of non-coding RNA.
  Nucleic Acids Res., 36(9):e52 (2008). <https://doi.org/10.1093/nar/gkn174>, [website][http://www.tcoffee.org/Projects/tcoffee/]
- Centroid homfold: Michiaki Hamada, Koichiro Yamada, Kengo Sato, Martin C. Frith, Kiyoshi Asai;
 CentroidHomfold-LAST: accurate prediction of RNA secondary structure using automatically collected homologous sequences,
 Nucleic Acids Research, Volume 39, Issue suppl_2, 1 July 2011, Pages W100–W106.
  <https://doi.org/10.1093/nar/gkr290>, [github][https://github.com/satoken/centroid-rna-package/]
- TurboFold (RNAstructure): Tan, Z., Fu, Y., Sharma, G., & Mathews, D. H. (2017).
 TurboFold II: RNA structural alignment and secondary structure prediction informed by multiple homologs.
  Nucleic Acids Research. 45: 11570-11581. <https://doi.org/10.1093/nar/gkx815>, [website][http://rna.urmc.rochester.edu/RNAstructure.html]
- UNAFold: Markham N.R., Zuker M. (2008) UNAFold. In: Keith J.M. (eds) Bioinformatics.
 Methods in Molecular Biology™, vol 453. Humana Press. <https://doi.org/10.1007/978-1-60327-429-6_1>, [website][http://unafold.rna.albany.edu/]
- Biopython: Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics.
 Bioinformatics 2009 Jun 1; 25(11) 1422-3. <http://dx.doi.org/10.1093/bioinformatics/btp163>, [website][https://biopython.org/]