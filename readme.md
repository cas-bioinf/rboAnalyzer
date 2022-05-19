# rboAnalyzer
A tool for analyzing BLAST search output for RNA sequences.

## Short description
rboAnalyzer is a tool complementing the BLAST algorithm when searching for a query
 sequence that is RNA with a secondary structure (which does not have to be known).

The high-scoring pairs (HSPs) in BLAST output are often incomplete
 (ie. the alignment in HSP does not cover the whole query sequence).
This is a major drawback when trying to characterize the potential ncRNA
 indicated by the HSP.

Therefore, rboAnalyzer tries to find full-length RNA sequences
 from the incomplete HSPs from the BLAST output and predict their secondary structures
 with one or more methods.
Score for similarity (as proxy to homology) between the estimated full-length sequence
 and query sequence is also computed.
The BLAST output is combined with computed data and presented in form of an interactive HTML page.

To achieve this, rboAnalyzer takes as input:
- the query sequence
- the BLAST output
- the BLAST database containing sequences within the output

The paper for this tool is available [here](https://doi.org/10.3389/fgene.2020.00675).

You can also try out the [rboanalyzer webserver](http://rboanalyzer.elixir-czech.cz) implementing interactive HSP analysis with selected secondary structure prediction methods.

## Installation
This tool was tested on 64-bit linux (Ubuntu 14, 18 and Centos 7). 

### Install via Conda
 The most convenient way to install this pipeline is to use Conda package manager. 
 
 If you don't have Conda for `Python3`, install it. We recommend the [miniconda3](https://conda.io/en/latest/miniconda.html).

The `rboanalyzer` package is available from `schwarz.marek` channel. You also need to include the `bioconda` and `conda-forge` channels when installing the `rboanalyzer`.

 
__Installation to new virtual environment__

The `rboAnalyzer` and its dependencies will be available only in the shell session for which the virtual environment was activated. You will need to `activate` the virtual environment for each session.

Following commands will install the rboAnalyzer into new conda environment named "rbo". If you don't wish to install rboAnalyzer to new environment, omit commands `1)` and `2)`.

```shell
# 1) update conda
conda update conda

# 2) create virtual environment
#conda create -n YOUR_VIRTUAL_ENV_NAME
conda create -n rbo

# 3) activate it
# conda activate YOUR_VIRTUAL_ENV_NAME
conda activate rbo

# 4) run installation
# ommit "-c conda-forge" and/or "-c bioconda" if those channels are in your .condarc file
conda install -c conda-forge -c bioconda -c schwarz.marek rboanalyzer
```

### Install from source

 __Prerequisites__
- python >= 3.4, <3.8 [link](https://www.python.org/downloads/)
  Verify that you have latest compatible `pip3`. If you try to use EOL python, please see the `setup.py` script on compatible `pip3` and `setuptools` versions and use them. 

- ncbi-blast+ >= 2.8.1, <2.10 [link](http://ftp.ncbi.nih.gov/blast/executables/blast+/2.9.0/)
  (The pipeline can use blast from version 2.6.0, however this version is not compatible with blast dbv5)
- locarna >= 1.9.2, <2 [link](https://github.com/s-will/LocARNA/releases/tag/v1.9.2.2)
- infernal >= 1.1, <1.2 [link](http://eddylab.org/infernal/)
- clustalo >= 1.2.4, <2 [link](http://www.clustal.org/omega/)
- muscle >= 3.8.31, <4 [link](https://www.drive5.com/muscle/downloads.htm)

For prediction:
- viennarna >=2.3.5, <3 (with refold.pl in PATH) [link](https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_3_x/ViennaRNA-2.3.5.tar.gz)
  Don't forget to add `refold.pl` to your `PATH`. The `refold.pl` script is located in the `ViennaRNA-[version]/src/Utils/`.
- centroid_homfold >= 0.0.15, <0.1 [link](https://github.com/satoken/centroid-rna-package/releases/tag/v0.0.15)
- RNAstructure >= 6.0, <7 (TurboFold - Text (Command Line) Interfaces ) [link](https://rna.urmc.rochester.edu/RNAstructure.html)
  Don't forget to set the `DATAPATH` environment variable [link](http://rna.urmc.rochester.edu/Text/Thermodynamics.html).

Optional (some prediction methods are not available without):
- UNAFold >= 3.8, <4 [link](http://www.unafold.org/)

Download this repository (or release), unpack it if needed. Go to directory with the source code for the rboAnalyzer and run

```shell
python3 setup.py install --user
```
Note that the `--user` switch puts the executables in `$HOME/.local/bin` and you may need to add it to `PATH`. 

The rboAnalyzer executable should be created.
To test it, restart terminal (close and open new) and run

```shell
rboAnalyzer --version
```
which should return the version number.

## Preparation

<a name="rfamdownload" id="rfamdownload"></a>

### Obtain Rfam database
For correct function the rboAnalyzer needs a copy of Rfam database.

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

### Installing UNAFold (optional)
Prediction methods using suboptimal structures need UNAFold software to work.
 It is available here [http://www.unafold.org](http://www.unafold.org).
 Follow installation instructions. The pipeline uses the `hybrid-ss-min` program.
 Either add it to PATH or add path to directory containing the executable to `config.txt` file.

### Shell autocomplete (optional)
The rboAnalyzer is equipped with argument completion for bash shell.
To enable this feature the `argcomplete` package needs to be installed (version >1.6, <2), and script registered with your shell.
It is installed by default in conda package. To install it manually run
```shell
pip3 install --user "argcomplete>1.6, <2a0"
```
Note that the `--user` switch puts the executables in `$HOME/.local/bin` and you may need to add it to `PATH`.

To register the script to your shell (more info [here](https://pypi.org/project/argcomplete/)).
```shell
register-python-argcomplete rboAnalyzer >> ~/.bashrc
```


### BLAST database
The rboAnalyzer needs to get relevant 5' and 3' regions of subject sequence of HSPs, for this we use the BLAST database used in the BLAST search.
For each analysis you need to provide the nucleotide BLAST database containing whole sequences (complete genomes, etc.) for the sequence ids present in the BLAST output.

The procedure on how to get the BLAST databases for the examples is described in the Example section. For the general information, please see the BLAST databases [section](#blastdatabase).

## Basic Usage
### Help
```shell
rboAnalyzer -h
```
<details>
<summary>help output</summary>

```
usage: rboAnalyzer [-h] -in PATH -q PATH (-db path | --entrez ENTREZ)
                   [--db_type {blastdb,fasta,gb,server}]
                   [--b_type {guess,xml,plain}] [--blast_regexp BLAST_REGEXP]
                   [--mode {simple,locarna,meta}] [--turbo_fast_preset]
                   [--centroid_fast_preset] [--html PATH] [--threads N]
                   [--csv CSV] [--json PATH] [--cm_file CM_file | --use_rfam]
                   [--download_rfam] [--version] [--config_file PATH]
                   [-pm [prediction_method_name [prediction_method_name ...]]]
                   [--pm_param_file PATH] [--logfile logfile]
                   [--subseq_window_locarna SUBSEQ_WINDOW_LOCARNA]
                   [--locarna_anchor_length LOCARNA_ANCHOR_LENGTH]
                   [--filter_by_eval FILTER_BY_EVAL | --filter_by_bitscore FILTER_BY_BITSCORE]
                   [-v] [--skip_missing] [--show_HSP]

rboAnalyzer - tool for analyzing BLAST search output for RNA sequences

optional arguments:
  -h, --help            show this help message and exit
  -db path, --blast_db path
                        Provide path to blast database, that is the complete
                        path with blast db name without any extension (*.nin,
                        nsd, nog, nsi, nhr, nsq, nal).
  --entrez ENTREZ       EMAIL - Indicate that you want to use NCBI Entrez
                        service to download required regions of sequences at
                        runtime. To comply with NCBI service rules you are
                        required to provide valid email address at which the
                        NCBI staff could contact you if they need to.

INPUT:
  -in PATH, --blast_in PATH
                        BLAST output file with hits to analyze.
  -q PATH, --blast_query PATH
                        The Blast query fasta file.

OUTPUT:
  --html PATH           Output html file with secondary structure pictures and
                        other useful stuff.
  --csv CSV             Output in csv table, infered sequence and structure
                        present.
  --json PATH           Dump all stored data to JSON (developer only - it is
                        possible to convert to all other output formats).

PARAMETERS:
  --mode {simple,locarna,meta}
                        Choose mode of hit elongation: simple (extend by
                        unaligned parts of query) locarna (run locarna
                        algorithm - uses secondary structure for better
                        alignment) meta (uses both methods and chooses the
                        alignment which has better RSEARCH score).
  --turbo_fast_preset   Act's as parameter preset for Turbo-fast setting the
                        max_seqs_in_prediction to 2. This means that only the
                        query sequence is used as homologous sequence. It is
                        useful if analyzing very distant BLAST HITs.
  --centroid_fast_preset
                        Parameter preset for centroid-fast. Set's the
                        max_seqs_in_prediction to 1. This means that only the
                        query sequence is used as homologous sequence for
                        prediction. It is useful if analyzing very distant
                        BLAST HITs.
  --config_file PATH    Provide config file if tools and data are in non-
                        default paths.
  -pm [prediction_method_name [prediction_method_name ...]], --prediction_method [prediction_method_name [prediction_method_name ...]]
                        Prediction method to use. Multiple prediction methods
                        are allowed. Possible values: C-A-U-r-Rc centroid
                        rnafold M-A-sub M-A-r-Rc fq-sub Turbo-fast TurboFold
                        centroid-fast C-A-sub C-A-r-Rc M-A-U-r-Rc rfam-sub
                        rfam-Rc rfam-centroid
  --pm_param_file PATH  Path to file with parameters for prediction methods in
                        JSON. Prediction methods not declared within provided
                        file are used with default values. File is in json
                        format. Default values (also example how to provide
                        parameters) are stored in '[install location]/rna_blas
                        t_analyze/BR_core/prediction_parameters.json'
  --subseq_window_locarna SUBSEQ_WINDOW_LOCARNA
                        N of nucleotides to add to expected start/end of
                        sequence before realignement. The unaligned
                        nucleotides are not included in reported sequence.
  --locarna_anchor_length LOCARNA_ANCHOR_LENGTH
                        Minimal number of adjacent matching bases in BLAST hit
                        to create an anchor for Locarna.

MISC:
  --db_type {blastdb,fasta,gb,server}
                        Type of a database provided. If 'fasta' or 'gb' then
                        --blast_db must be directory containing files with
                        names in accession.version format. Example '/home/my-
                        best-db/' with files like 'CP000001.1'.
  --b_type {guess,xml,plain}
  --blast_regexp BLAST_REGEXP
                        Provide python valid regular expression which capture
                        the index key to blastdb (usualy the accession.version
                        number).
  --threads N           Number of threads to use (default = N of logical cores
                        detected).
  --cm_file CM_file     Provided covariance model will be used for homology
                        inference instead of RSEARCH model.
  --use_rfam            Search in rfam database for covariance model to infer
                        homology with instead of RSEARCH model.
  --download_rfam       Retrieve RFAM covariance models database. Will
                        download only if new version avalible.
  --version             show program's version number and exit
  --logfile logfile     Path to where logfile should be written.
  --filter_by_eval FILTER_BY_EVAL
                        Filter the input blast by E-value. Only hits following
                        the rule will be kept. Example ">10e-10" will keep
                        only hits with eval greater then 10e-10. Interval can
                        be specified with "," e.g. ">10e-100, <10e-1". The
                        homologous sequences used with certain prediction
                        methods are taken from all hits (regardless of the
                        filtering).
  --filter_by_bitscore FILTER_BY_BITSCORE
                        Filter the input blast by bit score. Only hits
                        following the rule will be kept. Example "<20" will
                        keep only hits with bit score less then 20. Interval
                        can be specified with "," e.g. ">30, <45". The
                        homologous sequences used with certain prediction
                        methods are taken from all hits (regardless of the
                        filtering).
  -v, --verbose         output verbosity -> most detailed -vv (lot of output)
  --skip_missing        If given, the missing records in given blast database
                        will be skipped. This may alter the results of
                        bit_score computation (for homology prediction) and
                        secondary structure prediction for several methods.
  --show_HSP            Show HSP marker in NCBI sequence viewer.
```
</details>

### Usage
```shell
rboAnalyzer -in BLAST_OUTPUT.xml -db USED_DATABASE_PATH -q BLAST_QUERY.fasta
```
Note that only __one__ query sequence in the query FASTA file is expected.

## Example
Examples are provided in example directory.

To try examples you will need to:

  1. Install the rboAnalyzer and download Rfam. [How to?](#rfamdownload).

  2. Obtain a copy of `example` directory. If you've cloned or downloaded the program you should already have it.
   Otherwise it is [here](https://github.com/cas-bioinf/rna_blast_analyze/tree/master/example).

  3. cd to `example` directory.

#### Example 1:

Analyzing subset of NCBI blast HITs for [6S RNA](https://doi.org/10.1038%2F229147a0).

1) Now you need to obtain a copy of the BLAST database with all
 accessions which are in the BLAST output.
 As above, you can either get the NCBI nt database or download only the sequences (genomes) from the BLAST output.
 Here we describe the variant with downloading only the necessary sequences.
 This is done by using the `genomes_from_blast` by calling:
    ```shell
    genomes_from_blast -e YOUR_EMAIL_ADDRESS -in 6S_super_short.xml -o genomes.bdb
    ```
    The parameter `-e` `YOUR_EMAIL_ADDRESS` should be your valid email address on which NCBI staff could contact you
    if they need to. It is not logged by the tool.

    The parameter `-in` is path to file (in this example `6S_super_short.xml`) containing the BLAST output.

    The parameter `-o` is output file path. In this command the BLAST database with name `genomes.bdb` was created for you if everything was successful.
    You will need it in the next step.

    The intermediate file `genomes.bdb.fasta` was also created and contains all sequences added to the BLAST database.
    When another BLAST output is analyzed (and sequences are needed) then only those sequences not present in the intermediate FASTA file are downloaded (assuming same BLAST db name is used).

2) Now you can run the pipeline itself:
    ```shell
    rboAnalyzer -in 6S_super_short.xml -q 6S_query.fasta -db genomes.bdb --html 6S_out.html
    ```
3) The output is single html file. You can scroll through analyzed HSPs, show the genomic loci of the HSP and select data to export.

#### Example 2:
Analyzing possible remote homologs for [MS1 RNA](https://doi.org/10.1093%2Fnar%2Fgku793). This take about 10 minutes on average pc.

```shell
# update the genome database with new sequences
genomes_from_blast -e YOUR_EMAIL_ADDRESS -in MS1_BLAST_output.txt -o genomes.bdb

# run the rboAnalyzer
rboAnalyzer -in MS1_BLAST_output.txt -q MS1_query.fasta -db genomes.bdb --html MS1_out.html --prediction_method rnafold rfam-Rc Turbo-fast --turbo_fast_preset
```

The BLAST was run with database where Streptomycetaceae and Mycobacteriaceae families where excluded.
As the MS1 RNA is primarily known from Mycobacteriaceae family we can expect
  incomplete HITs and many false positives.

Also you can notice, that the BLAST output was in text format, the rboAnalyzer accepts BLAST output in plain text and xml.

The BLAST database is obtained with similar command as in previous example.
Since the output BLAST database file is the same as before the `genomes_from_blast`
 will check the intermediate fasta file (`genomes.bdb.fasta`) and will
 download only sequences which are not present.

We can expect that the BLAST output contain many false positive HSPs,
so we selected prediction methods from those, which do not rely on information in the BLAST output itself.
These are:
- rnafold
- rfam-Rc
- rfam-centroid
- Turbo-fast with "max_seqs_in_prediction" parameter set to 2 (`--turbo_fast_preset` flag)
- centroid-fast with "max_seqs_in_prediction" parameter set to 1 (`--centroid_fast_preset` flag)
- rfam-sub (If `UnaFold` is installed.)


### Solving issues:
- __BLAST txt output from WEB not recognized__    
  Reason: the NCBI changed the `txt` output for the WEB BLAST when they switched to new design. Our txt parser is compatible with commandline `txt` output (`-outfmt 0`) which was also the `txt` output for the WEB BLAST.   
  Solution: Download the `xml` output or switch to the "Traditional result page" and download `txt` there.
- __One or more records not found__   
  Reason: the blastdbcmd was not able to find sequence(s) with respective id(s) in provided database.
  This is due to inconsistency between the sequence accessions and the BLAST database.
  The inconsistency may rise from:
    1. __sequence is not in the database__   
      Solution: Provide correct BLAST database (update current or create new with `genomes_from_blast`).
    2. __capturing regexp does not capture the accession number__   
      Solution: Provide capturing regular expression (python 3 syntax) for capturing the sequence id from the fasta header (it must match the id to the BLAST database used)
    3. __the BLAST database was created without the `-parse_seqids` flag__   
      Solution: Create new database from the sequences used to create new one, this time with `-parse_seqids` flag.
    4. __inconsistent accession mapping__   
      Detailed cause: In certain NCBI's BLAST databases there are some sequences with accession numbers for which ENTREZ 
       return sequences with different accession numbers. (We observed this for several sequences with accession numbers starting with `GPS_`.) 
       The `genomes_from_blast` script depends on matching accession numbers and can't handle the inconsistency.   
      Solution: We recommend obtaining the database from NCBI which was used to generate the BLAST output. 
       Other option is to manually add the sequence in question to the blast database.  

  Another option is to call pipeline with `--skip_missing` flag.
  This will skip the missing sequences.

  Note that no HSP for the missing sequence will be included in pipeline output
  and some prediction methods may be influenced by the missing sequence.

- __The `genomes_from_blast` failed__   
  The `genomes_from_blast` script has build in handling of failed downloads,
  but by default it tries only 10 times. If you are on unstable connection
  you might get better results by setting the `--retry` to some larger number.
  Also check if NCBI ENTREZ services are functional.

- __ValueError: could not convert string to float__    
  (raised with traceback containing `infer_homology.py ... dtypes/cast.py)`    
  Please check, that there is only one query sequence in the input FASTA file. 

<a name="blastdatabase" id="blastdatabase"></a>

### BLAST databases

#### BLAST on NCBI web

If you used the BLAST using the NCBI web service against one of preformatted databased, you can download the whole database or use a `genomes_from_blast` script to download only the sequences in your blast output.

1. downloading whole database (~50GB)   
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
  This command will download all needed genomes and create BLAST database for you.
    ```shell
    # The `YOUR_EMAIL_ADDRESS` is needed so the NCBI would contact you in case of misuse of their resources.

    genomes_from_blast -e YOUR_EMAIL_ADDRESS -in BLAST_OUT_FILE -o BLAST_DATABASE_OUTFILE_NAME
    ```

#### Custom BLAST database
If custom database was used for the BLAST search you need to ensure multiple things for the rboAnalyzer to find the sequences correctly.
1. custom database is nucleotide and it was created with `-parse_seqids` (this makes sequences retrievable by their ids).
2. provide regular expression capturing the sequence ids. By default the rboAnalyzer captures the Accession.Version as documented [here](https://www.ncbi.nlm.nih.gov/Sequin/acc.html).

## Learn more
Description of the pipeline logic is avialable in the [publication](https://doi.org/10.3389/fgene.2020.00675).

Details can also be found [here](docs/help.md).

This readme and documentation in `docs` are valid for latest release 0.1.4 and dev 0.1.5a1.

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
 Methods in Molecular Biology™, vol 453. Humana Press. <https://doi.org/10.1007/978-1-60327-429-6_1>, [website](http://www.unafold.org)
- Biopython: Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics.
 Bioinformatics 2009 Jun 1; 25(11) 1422-3. <http://dx.doi.org/10.1093/bioinformatics/btp163>, [website](https://biopython.org/)
- NumPy [website](https://numpy.org)
- Pandas [website](https://pandas.pydata.org)
- Jinja2 [website](https://jinja.palletsprojects.com)
- Matplotlib [website](https://matplotlib.org)

## Citation
- If you find this software useful please cite this [paper](https://doi.org/10.3389/fgene.2020.00675).

  Schwarz, M., Vohradský, J., Modrák, M. and Pánek, J., 2020. rboAnalyzer: A Software to Improve Characterization of Non-coding RNAs From Sequence Database Search Output. Frontiers in genetics, 11, p.675.

  The rboAnalyzer version used for the paper is `0.1.4`.

- If you find helpfull the webserver implementation, please cite this [paper](https://doi.org/10.1093/bioinformatics/btab073).
                         
  Marek Schwarz, Jiří Vohradský, Josef Pánek, rboAnalyzer webserver: web service for non-coding RNA characterization from NCBI BLAST output, Bioinformatics, Volume 37, Issue 17, 1 September 2021, Pages 2755–2756, https://doi.org/10.1093/bioinformatics/btab073
                         
## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](docs/ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund - Project "ELIXIR-CZ: Budování kapacit" (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](docs/logolink_OP_VVV_hor_barva_eng.jpg)

This work was also supported by Czech Science Foundation GA15-00885S.k
