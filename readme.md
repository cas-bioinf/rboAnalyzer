# rboAnalyzer
A tool for analyzing BLAST search output for RNA sequences.

## Short description
rboAnalyzer is a tool for complement BLAST algorithm when searching for query
 sequence that is RNA with secondary structure (which does not have to be known).

The HSPs in BLAST output are often incomplete
 (ie. the alignment in HSP does not cover whole query sequence).
This is major drawback when trying to characterize the potential ncRNA
 indicated by the HSP.

Therefor the rboAnalyzer estimates full-length RNA sequences
 from the incomplete HSPs from the BLAST output and predicts secondary structure
 with one or more methods.
Score for similarity (as proxy to homology) between the estimated full-length sequence
 and query sequence is also computed.
The BLAST output is combined with computed data and presented in form of interactive HTML page.

For this the rboAnalyzer need the BLAST input and output.
That is the query sequence, the BLAST database containing sequences within the output
 and the BLAST output itself.


## Installation

<!---
### Install via Conda
 The easies way to install this pipeline is to use conda. This package is available
 from bioconda channel.

 If you don't have conda, install it. [conda docs](https://conda.io/docs/index.html)

 Then open terminal and run

System wide installation
 The `rboAnalyzer` and `genomes_from_blast` executables and their
  dependencies will be available from terminal.

```shell
conda install -c conda-forge -c bioconda rboAnalyzer
```

Installation to virtual environment
 The `rboAnalyzer` and its dependencies will be available only in shell
  for which the virtual environment was activated. If virtual environment is
  used, then you need to activate virtual environment before usage.

```shell
# create virtual environment
conda create -n YOUR_VIRTUAL_ENV_NAME

# activate it
# new conda version (> 4.4)
conda activate YOUR_VIRTUAL_ENV_NAME

## old conda version
#source activate YOUR_VIRTUAL_ENV_NAME

# run installation
conda install -c conda-forge -c bioconda rboAnalyzer
```
-->
### Install from source

 __Prerequisites__
* python >= 3.4
* ncbi-blast+ >= 2.6
* locarna >= 1.9
* infernal >= 1.1.2
* clustalo
* muscle

For prediction:
* viennarna (with refold.pl in PATH)
* centroid_homfold
* RNAstructure >= 6.0 (TurboFold)

Optional (some prediction methods are not available without):
* UNAFold >= 3.8

Clone or download this repository. Go to root folder and run

```shell
python3 setup.py install
```

The rboAnalyzer executable should be created.
To test it, restart terminal (close and open new) and run

```shell
rboAnalyzer --version
```
Which should return number of the installed version.

## Preparation

<a name="rfamdownload" id="rfamdownload"></a>

### Obtain Rfam database
For correct function the rboAnalyze needs a copy of Rfam database.

There are 2 ways:
1. Run rboAnalyzer with `--download_rfam` flag.
    ```shell
    rboAnalyzer --download_rfam
    ```
This will download Rfam covariance models to default directory
(`[INSTALL_LOCATION]/rna_blast_analyze/3rd_party_source/rfam`).

2. Alternatively download the `Rfam.cm.gz` file from
[Rfam CURRENT](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT),
unpack it and add the path to directory containing `Rfam.cm` file to your `config.txt` file.
Note that running rboAnalyzer with `--download_rfam` will overwrite this manually installed file.

### BLAST database
The rboAnalyzer needs to get relevant 5' and 3' regions of subject sequence of HSPs, for this we use the BLAST database used in the BLAST search.
For each analysis you need to provide the nucleotide BLAST database containing whole sequences (complete genomes, etc.) for the sequence ids present in the BLAST output.

#### BLAST on NCBI web

If you used the BLAST using the NCBI web service against one of preformatted databased, you can download the whole database or use a `genomes_from_blast` script for download.

1. downloading whole database
The latest databases are provided here [NCBI LATEST](ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST).
Note that databases included in the BLAST database releases are not the latest ones.
This code snippet can be used to obtain and update the database:
    ```shell
    # cd to directory to which you want to download the database

    # for the "nt" database
    wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/nt*

    # for other databases provided by NCBI (insert database name without square brackets)
    wget -N ftp://ftp.ncbi.nih.gov/blast/db/cloud/LATEST/[database name]*
    ```

2. downloading only relevant sequences
  If you do not wish to download whole blastdb you may use prepared script
 `genomes_from_blast`, which downloads only the needed sequences
  (those in the BLAST output) and build the blastdb from them.
  This command will download all needed genomes and create blast database for you (if `makeblastdb` command is available).
    ```shell
    # The `YOUR_EMAIL` is needed so the NCBI would contact you in case of missuse of their resources.

    genomes_from_blast -email YOUR_EMAIL -blast_in BLAST_OUT_FILE -out FASTA_FILE_OUT
    ```

#### Custom BLAST database
If custom database was used for the BLAST search you need to ensure multiple things for the rboAnalyzer to find the sequences correctly.
1. custom database is nucleotide and it was created with `-parse_seqids` (this makes sequences retrievable by their ids).
2. provide regular expression capturing the sequence ids. By default the rboAnalyzer captures the Accession.Version as documented [here](https://www.ncbi.nlm.nih.gov/Sequin/acc.html).

### Installing UNAFold
Prediction methods using suboptimal structures need UNAFold software to work.
 It is available here <http://unafold.rna.albany.edu/>.
 Follow installation instructions. The pipeline uses the `hybrid-ss-min` program.
 Either add it to PATH or add path to directory containing the executable to `config.txt` file.

### Shell autocomplete (optional)
The pipeline is equipped with argument completion for bash shell.
To enable this feature you need to register the script (more info [here](https://pypi.org/project/argcomplete/)).

To get the autocomplete working run:
```shell
register-python-argcomplete rboAnalyzer >> ~/.bashrc
```

## Basic Usage
### Help
```shell
rboAnalyzer -h
```

### Usage
```shell
rboAnalyzer -in BLAST_OUTPUT.txt -db USED_DATABASE_PATH -q BLAST_QUERY.fasta
```

## Example
Examples are provided in example directory.

To try examples you will need to:

  1. Install the rboAnalyzer and download Rfam [how to](#rfam_download).

  2. Obtain a copy of `example` directory
   [here](https://github.com/cas-bioinf/rna_blast_analyze/tree/master/example).

  3. cd to `example` directory.

#### Example 1:

Analyzing subset of NCBI blast HITs for [6S RNA](https://doi.org/10.1038%2F229147a0).

1) Now you need to obtain a copy of the BLAST database with all
 accessions which are in the BLAST output.
 As above, you can either get the NCBI nt database or download only the sequences (genomes) from the BLAST output.
 Here we describe the variant with downloading only the necessary sequences.
 This is done by using the `genomes_from_blast` by calling:
    ```shell
    genomes_from_blast -email YOUR_EMAIL -in 6S_short.xml -out genomes.fasta
    ```
    The `YOUR_EMAIL` should be valid email on which NCBI staff could contact you
    if they need to. It is not saved nor logged by the tool.
    The 6S.fasta file is not needed if the BLAST database was created
    successfully and you can delete it.

    The BLAST database with name genomes.fasta.bdb was created for you if everything was successful.
    You will need it in the next step.

2) Now you can run the pipeline itself:
    ```shell
    rboAnalyzer -in 6S_short.xml -q 6S_query.fasta -db genomes.fasta.bdb -html 6S_out.html
    ```
3) The output is single html file. You can scroll through analyzed HSPs, show the genomic loci of the HSP and select data to export.

#### Example 2:
Analyzing possible remote homologs for [MS1 RNA](https://doi.org/10.1093%2Fnar%2Fgku793).

The BLAST was run with database where Streptomycetaceae and Mycobacteriaceae families where excluded.
As the MS1 RNA is primarily known from Mycobacteriaceae family we can expect
incomplete HITs and many false positives.

Also you can notice, that the BLAST output is now in text format, the rboAnalyzer accepts BLAST output in plain text and xml.

We get the BLAST database as in previous example. (If the `genomes_from_blast`
  script was used, you can now run the same command with `MS1_BLAST_output.txt`
  as input and only the sequences not in the database will be downloaded).

With the database (assume `genomes.fasta.bdb`) we can run the rboAnalyzer.
We can expect that the BLAST output contain many false positive HSPs,
so we choose prediction methods which do not rely on information in the BLAST output itself.
These are:
- rnafold
- rfam-Rc
- rfam-centroid
- Turbo-fast with "max_seqs_in_prediction" parameter set to 2 (`--turbo_fast_preset` flag)
- centroid-fast with "max_seqs_in_prediction" parameter set to 1 (`--centroid_fast_preset` flag)
- rfam-sub (If `UnaFold` is installed.)

In this case we will use the rnafold, rfam_rnafoldc and Turbofold_fast with preset.
```shell
rboAnalyzer -in MS1_BLAST_output -q MS1_query.fasta -db genomes.fasta.bdb -html MS1_out.html --prediction_method rnafold rfam-Rc Turbo-fast --turbo_fast_preset
```

## Solving issues:
- __One or more records not found__
  Reason: the blastdbcmd was not able to find sequence(s) with respective id(s) in provided database.
  This is due to inconsistency between the sequence accessions and the BLAST database.
  The inconsistency may rise from:
    1. __sequence is not in the database__
      Solution: Provide correct blast database (update current or create new with `genomes_from_blast`).
    2. __capturing regexp does not capture the accession number__
      Solution: Provide capturing regular expression (python 3 syntax) for capturing the sequence id from the fasta header (it must match the id to the BLAST database used)
    3. __the BLAST database was created without the `-parse_seqids` flag__
      Solution: Create new database from the sequences used to create new one, this time with `-parse_seqids` flag.

  Another option is to call pipeline with `--skip_missing` flag.
  This will skip the missing sequences.

  Note that no HSP for the missing sequence will be included in pipeline output
  and some prediction methods may be influenced by the missing sequence.

- __The `genomes_from_blast` failed__
  The `genomes_from_blast` script has build in handling of failed downloads,
  but by default it tries only 10 times. If you are on unstable connection
  you might get better results by setting the `--retry` to some larger number.
  Also check if NCBI entrez services are functional.

<!--
3. WARNING - refold.pl could not be located (not in PATH)
  Issued when refold.pl executable could not be found.
  If everything was installed with conda to same environment, this will be resolved automatically.
  If not, you must add refold.pl to PATH or add path to refold.pl to config.txt (see docs/config_how_to.md)

4. WARNING - The TurboFold is installed but the DATAPATH environment variable is not set nor present in config.txt
  Issued when TurboFold's DATAPATH env variable is not set.
  If everything was installed with conda to same environment, this will be resolved automatically.
  If not, you must add path to `RNAstructure/data_tables` to config.txt (see docs/config_how_to.md).
-->

## Notes
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

## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](docs/ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund - Project "ELIXIR-CZ: Budování kapacit" (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](docs/logolink_OP_VVV_hor_barva_eng.jpg)

This work was also supported by Czech Science Foundation GA15-00885S.k