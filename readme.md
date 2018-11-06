# RNA_BLAST_ANALYZE
Pipeline for analyzing BLAST search output for non-coding RNAs (ncRNAs).

## Short description
Provided with query RNA sequence and BLAST output, the pipeline will
 realign each BLAST hit with Locarna to account for more RNA sequence
 variability. The realign step is done to recover potential ncRNA
 sequence at that loci precisely even from low-scoring BLAST hits.

Sequences are analyzed with RSEARCH and its bit-score is used as
 similarity criteria.

Next sequences are predicted with one or more methods and predicted
 structures are merged with blast output to be easy to understand.

More information on the pipeline is avalible in the [docs](docs/help.md) section.


## Installation

### Install via Conda
 The easies way to install this pipeline is to use conda. This package is avalible
 from bioconda channel.

 If you don't have conda, install it from [here](https://conda.io/docs/index.html).

 Then open terminal and run
```
# for system-wide installation
conda install -c conda-forge -c bioconda rna_blast_analyze

# install to virtual environment (not tested)
conda create YOUR_VIRTUAL_ENV_NAME
source activate YOUR_VIRTUAL_ENV_NAME
conda install -c conda-forge -c bioconda rna_blast_analyze
```
If virtual environment is used, then you need to activate virtual
 environment before usage


### Install from source

 __Prequisities__
* python >= 3.4
* ncbi-blast+ >= 2.6
* locarna >= 1.9
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
* UNAFold >= 3.8

Clone or download this repository. Go to root folder and run
```
python3 setup.py install
```

The rna_blast_analyze executable should be created.
To test it, restart terminal (close and open new) and run
```
rna_blast_analyze --version
```
Which should return installed version number.

## Preparation
### Shell autocomplete
The pipeline is equipped with argument completion for bash shell.
To enable this feature you need to register the script (more info [here](https://pypi.org/project/argcomplete/)).

To get the autocompletion working run:
```
register-python-argcomplete rna_blast_analyze >> ~/.bashrc
```

### Downloading the BLAST database
For each analysis you need to provide the BLAST database which was used for BLAST search.
The latest databases are rather non-intuitively provided here [NCBI LATEST](ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST).
And this code snippet can be used to obtain and update the database:
```
for the "nt" database
wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/nt*

for other provided databases (insert database name without square brackets)
wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/[database name]*
```

If you do not wish to download whole blastdb you may use prepared script
 `genomes_from_blast`, which downloads only the needed sequences
  (those in the BLAST output) and build the blastdb from them.
```
genomes_from_blast -email YOUR_EMAIL -blast_in BLAST_OUT_FILE -out FASTA_FILE_OUT
```
This command will download all needed genomes and create blast database for you (if `makeblastdb` command is avalible).
The `YOUR_EMAIL` is needed so the NCBI would contact you in case of missuse of their resources.

### Installing UNAFold
Prediction methods using suboptimal structures need UNAFold software to work.
 It is avalible here <http://unafold.rna.albany.edu/>.
 Follow installation instructions. The pipeline uses the `hybrid-ss-min` program.
 Either add it to PATH or add path to it to `config.txt` file.

### Download RFAM database
There are 2 ways:
1) Download `Rfam.cm.gz` from [RFAM CURRENT](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT),
unpack it and add the path to `Rfam.cm` to your `config.txt` file.
2) Use build-in download by issuing the `--download_rfam` flag.
This will download Rfam covariance models to designated directory
(by default this `rna_blast_analyze/3rd_party_source/rfam`).

## Basic Usage
### Help
```
rna_blast_analyze -h
```

### Usage
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
4) Now you need to obtain a copy of the BLAST database with all
 accessions which are in the BLAST output. This is most simply
  done by using the `genomes_from_blast` by calling:
    ```
    genomes_from_blast -email YOUR_EMAIL -blast_in 6S_short.xml -out 6S.fasta
    ```
    The 6S.fasta file is not needed if the blast database was created
    successfully and you can delete it.

    The blast database with name 6S.fasta.bdb was created for you if everything was successful.
    You will need in next step.

5) Now you can run the pipeline itself:
    ```
    if first time (if rfam was not downloaded):
    rna_blast_analyze -blast_in 6S_short.xml -blast_query 6S_query.fasta -blast_db 6S.fasta.bdb -html 6S_out.html --download_rfam

    else:
    rna_blast_analyze -blast_in 6S_short.xml -blast_query 6S_query.fasta -blast_db 6S.fasta.bdb -html 6S_out.html
    ```
6) The output is single html file. On open the NCBI genome viewer will fetch data to render the genomic loci of the hit.

### Solving issues:
1) One or more records not found.

    Solution: Provide correct blast database
    (update current or create special by `genomes_from_blast`).
    Another option is to call pipeline with `--skip_missing` flag.
    This will skip the missing sequences.
    No blast hit to that sequence will be included in pipeline output
    and RSEARCH score and some prediction methods may be influenced
    by that missing sequence.


### Notes
This is Beta version. Especially the default parameters for prediction can change.

## References
- ViennaRNA: Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
 ViennaRNA Package 2.0. <https://doi.org/10.1186/1748-7188-6-26>.
 We include copy of refold.pl script with Viennarna license here for convinience. [website](https://www.tbi.univie.ac.at/RNA/)
- Locarna: Sebastian Will, Tejal Joshi, Ivo L. Hofacker, Peter F. Stadler, and Rolf Backofen.
LocARNA-P: Accurate boundary prediction and improved detection of structural RNAs
RNA, 18 no. 5, pp. 900-14, 2012. <https://doi.org/10.1261/rna.029041.111>, [website](http://rna.informatik.uni-freiburg.de/LocARNA/Input.jsp)
- NCBI BLAST: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008)
 "BLAST+: architecture and applications." BMC Bioinformatics 10:421. <https://doi.org/10.1186/1471-2105-10-421>, [website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- RSEARCH: Finding Homologs of Single Structured RNA Sequences.
 R. J. Klein, S. R. Eddy. BMC Bioinformatics, 4:44, 2003. <https://doi.org/10.1186/1471-2105-4-44>, [website](http://eddylab.org/software.html#rsearch)
- Infernal: Infernal 1.1: 100-fold Faster RNA Homology Searches.
 E. P. Nawrocki, S. R. Eddy. Bioinformatics, 29:2933-2935, 2013. <https://doi.org/10.1093/bioinformatics/btt509>, [website](http://eddylab.org/infernal/)
- Clustal omega: Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011).
 Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega.
  Molecular Systems Biology 7:539. <https://doi.org/10.1038/msb.2011.75>, [website](http://www.clustal.org/omega/)
- MUSCLE: Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput
 Nucleic Acids Res. 32(5):1792-1797. <https://doi.org/10.1093/nar/gkh340>, [website](https://www.drive5.com/muscle/)
- T-coffee (r-coffee): Wilm, A., Higgins, D.G., Notredame, C.
 R-Coffee: a method for multiple alignment of non-coding RNA.
  Nucleic Acids Res., 36(9):e52 (2008). <https://doi.org/10.1093/nar/gkn174>, [website](http://www.tcoffee.org/Projects/tcoffee/)
- Centroid homfold: Michiaki Hamada, Koichiro Yamada, Kengo Sato, Martin C. Frith, Kiyoshi Asai;
 CentroidHomfold-LAST: accurate prediction of RNA secondary structure using automatically collected homologous sequences,
 Nucleic Acids Research, Volume 39, Issue suppl_2, 1 July 2011, Pages W100–W106.
  <https://doi.org/10.1093/nar/gkr290>, [github](https://github.com/satoken/centroid-rna-package/)
- TurboFold (RNAstructure): Tan, Z., Fu, Y., Sharma, G., & Mathews, D. H. (2017).
 TurboFold II: RNA structural alignment and secondary structure prediction informed by multiple homologs.
  Nucleic Acids Research. 45: 11570-11581. <https://doi.org/10.1093/nar/gkx815>, [website](http://rna.urmc.rochester.edu/RNAstructure.html)
- UNAFold: Markham N.R., Zuker M. (2008) UNAFold. In: Keith J.M. (eds) Bioinformatics.
 Methods in Molecular Biology™, vol 453. Humana Press. <https://doi.org/10.1007/978-1-60327-429-6_1>, [website](http://unafold.rna.albany.edu/)
- Biopython: Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics.
 Bioinformatics 2009 Jun 1; 25(11) 1422-3. <http://dx.doi.org/10.1093/bioinformatics/btp163>, [website](https://biopython.org/)