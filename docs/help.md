# rboAnalyzer
## Introduction
rboAnalyzer is a pipeline meant as a complementary step when BLAST algorithm
 was used to search for query that is non-coding RNA (ncRNA) with secondary structure (does not have to be known).
As the BLAST is general sequence alignment algorithm, it's results (output)
 is missing some very useful features in context of ncRNAs.

Also, the output apart from well scoring hits, usually contains sequence
 fragments (subject sequence only partially aligned to query).
Such hit may bring valuable information by capturing distant homology, or
it may be nothing.

With our rboAnalyzer we add information to such BLAST search to help researcher
 decide which hits are real ncRNAs and what their secondary structure might be.

### Searching sequence database with BLAST

Since BLAST scoring parameters can significantly influence the number and length of obtained HSPs we suggest to test out multiple BLAST parameter groups to find out which provides relevant results for your use-case. Analysis with rboAnalyzer is helpful by providing potential secondary structure and genomic context information to the HSPs of interest, especially for HSPs with partial alignment or alignment with many mismatches and/or gaps. Additionally it is possible to explore and combine multiple BLAST parameter settings in addition to the default ones, as no single ones are appropriate for all situations (ref. [1](#references)).

Here we list scoring parameters that where used for ncRNA searches (ref. [2-4](#references)) and recommended for cross-species exploration and RNA search (ref. [5](#references)):

| Parameter        | [ref. 2](#references) | [ref. 3](#references) | [ref. 4](#references) | [ref. 5](#references) |
|------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| MatchReward      | 5                     | 5                     | 5                     | 1                     |
| MismatchPenalty  | -4                    | -4                    | -4                    | -1                    |
| GapOpeningCost   | 10                    | 8                     | 10                    | 1                     |
| GapExtensionCost | 10                    | 6                     | 6                     | 2                     |
| WordSize         | 7                     | 4                     | 7                     | 9                     |

## Getting help
Commandline reference can be viewed [here](../readme.md#help) (click on the dropdown) or accessed by
```shell
rboAnalyzer -h
```

In case of any questions, comments, bugreports or suggestions, create an issue or write to `marek.schwarz AT biomed.cas.cz`.


## Functionality overview
<img src="RBA_pipeline_overview.svg" width="700px" />

The rboAnalyzer has 3 stages:
1) Estimation of full-length RNA sequence from HSP (extension).
2) Estimation of homology of estimated full-length sequences to query sequence.
3) Prediction of secondary structures.

Each of these stages has dedicated section.

### Methods for estimation of full-length sequences from HSPs
The pipeline has 3 methods for estimating the full-length sequences from BLAST HSPs.

1) __simple__

    This means that location of estimated full-length sequence is computed from
     unaligned parts of query sequence on 5' and 3' ends of the HSP.

2) __locarna__

    In this method, the loci containing hit with flanking regions at the subject sequence is realigned to
     the query sequence with Locarna algorithm. The sequence
     aligned to the query is considered to be the estimated full-length sequence.

3) __meta__

    Here the two aforementioned methods are combined and estimated full-length sequences
     are scored with covariance model. The better scoring sequence is chosen.

#### ad simple)
<img src="figure_blast_simple.svg" width="700px" />

In the __simple__ mode we compute the location of the extended subject sequence according 
to the unaligned parts of the query sequence, i.e. those that were not aligned in HSP 
and flank the HSP alignment at both 5’ and 3’ ends. In this toy example we have the Plus/Plus 
BLAST HSP with a section of the query sequence between nucleotides 10 and 21 aligned 
to a section of the subject sequence between nucleotides 1000 and 1009. 
The task is to extend the partial HSP subject sequence between nucleotides 1000 and 1009 
to the length of the query sequence.

Suppose that the query sequence (the red bar in the figure of the example) is 50 bases long. 
Then the length of the unaligned part of the query at 5’ end is 9 nucleotides 
(subtract nucleotide positions 10 - 1) and the length of the unaligned part of 
the query at 3’ end is 29 nucleotides (subtract positions 50 - 21). 
The positions of the extended subject sequence at the whole subject sequence is computed by 
adding/subtracting the lengths of the unaligned parts of the query sequence to/from 3’/5’ 
ends of the partial HSP subject sequence, respectively. 
Then, 5’ and 3’ ends of the extended HSP sequence lie at nucleotides 991 (1000 – 9) 
and 1038 (1009 + 29), respectively. Because there are 2 gaps in HPS subject sequence, 
the resulting extended subject sequence will be 2 nucleotides shorter than the query sequence. 

Theoretically, the 2 extra nucleotides could be added at 3’ end of the extended subject 
sequence to make it as long as the query sequence. 
But this might not be biologically relevant as the gap could occur naturally instead 
of being caused by the alignment.


#### ad locarna)
<img src="figure_blast_locarna.svg" width="700px" />

With __locarna__ mode we first extract so called _supersequence_, which is
 region on subject sequence as with __simple__,
 additionally padded on 5' and 3' ends by extra sequence from the subject sequence.
This _supersequence_ is then realigned with Locarna algorithm to obtain the estimated full-length sequence.

The Locarna algorithm utilises possible pairings in it's computations,
 thus it is better suited to align RNAs then BLAST algorithm.
The Locarna is by default called with `struct-local=0`,
 `sequ-local=0` and `free-endgaps=++++` parameters.
Additionally, the information about matching nucleotides from BLAST HSPs
is used to construct so called anchor for the Locarna algorithm.
The anchor defines columns of alignment which are considered aligned.
As the anchor we consider consecutive series of matches of length at
 least `L` in BLAST alignment.
The default value of `L` is 7.
This way the alignment is anchored and the Locarna algorithm can align
 query to the _supersequence_. With the `free-endgaps=++++` option,
 the algorithm does not put penalty to unaligned ends of _supersequence_.
The estimated full-length sequence is the continuous part of _supersequence_ aligned to the query sequence
 (i.e. the subject sequence between the bases, inclusive, on subject sequence matching to the 5' terminal and 3' terminal bases).

#### add meta)
This approach combines the __simple__ and __locarna__. It computes both and 
 for each HSPs it chooses the estimated full-length sequence with higher score to covariance model.

### Estimation of homology
Here we compute score for relation between the estimated full-length sequence and query sequence.
The computation is based on aligning covariance model (CM) to each estimated full-length
 sequence with `cmalign` program from the Infernal package.

We've implemented 3 options on how to provide covariance model:

1) build with RSEARCH (default)

    By default, we build the covariance model from the query sequence (secondary structure predicted by RNAfold) and RIBOSUM matrix.
    The RIBOSUM is RIBOSUM65 by default and it can be changed in [alternative](config_how_to.md) `config.txt` file.

2) supply your own model (the `--cm_file` option)

    If the covariance model is known, it can be provided with `--cm_file` option.
    Only one model per file is allowed.
    
    Note that if you provide the covariance model, it will also be used in all methods for prediction of secondary structures using covariance models (those starting with `rfam`).

3) infer from Rfam (the `--use_rfam` option)

    The Rfam database is searched with query sequence for the best matching
    model (`cmscan`).

### Prediction of secondary structures
The rboAnalyzer can use multiple approaches (prediction methods) to predict secondary structures.
The prediction methods can be (roughly) divided to following groups:

- Predict structure independently of other estimated full-length sequences
    The advantage for these methods is robustness to possible improper parameter choice.
    - [rnafold](prediction_methods.md#rnafold)
    - [fq-sub](prediction_methods.md#fq-sub)
    - [rfam-Rc](prediction_methods.md#rfam-Rc)
    - [rfam-centroid](prediction_methods.md#rfam-centroid)
    - [rfam-sub](prediction_methods.md#rfam-sub)

- Use of selected estimated full-length sequences as reference
    - [centroid](prediction_methods.md#centroid)
    - [TurboFold](prediction_methods.md#TurboFold)
    - [Turbo-fast](prediction_methods.md#Turbo-fast)
    - [centroid-fast](prediction_methods.md#centroid-fast)

- Use of selected estimated full-length sequences to build consensus secondary structure
    - [C-A-r-Rc](prediction_methods.md#C-A-r-Rc)
    - [M-A-r-Rc](prediction_methods.md#M-A-r-Rc)
    - [C-A-U-r-Rc](prediction_methods.md#C-A-U-r-Rc)
    - [M-A-U-r-Rc](prediction_methods.md#M-A-U-r-Rc)
    - [C-A-sub](prediction_methods.md#C-A-sub)
    - [M-A-sub](prediction_methods.md#M-A-sub)

## Output

### Output formats
The rboAnalyzer is able to produce several output formats, most handy
  being the `.html`.
- html
    Stand-alone web page containing estimated full-length sequences and predicted secondary structures.
    If internet connection is available, it can be used to view respective
    genome loci for each BLAST HSP using NCBI SeqViewer.
- json
    Json-readable rboAnalyzer output (contains all data).
- csv
    Output table in comma separated values. Contains all important information
    including original HSP data, estimated full-length sequence location,
    sequence and predicted secondary structure(s).

## HTML output
In the head section there is report on basic input data and name of Rfam
 covariance model with best score to the provided query sequence.

The html output is organized around BLAST output.
Each BLAST HSP gets it's separate section with five parts:

  1) the text representation of BLAST HSP

  2) rboAnalyzer report with estimated full-length sequence indices and RSEARCH bit score

  3) the estimated full-length sequence itself

  4) one or multiple predicted secondary structures

  5) NCBI Sequence viewer (optional - by default only load button is shown)

The `html` outputs offers sorting, selecting sequences and structures and
  their export to fasta format or fasta-like format with predicted secondary structures in dot-bracket notation.
If internet connection is available, the NCBI genome browser can be used
  to explore synteny and known features of current genome.

The header for each BLAST HSP contains Accession.Version number (based on provided regular expression).
The header is also color-coded on color scale from green to red based on RSEARCH score.
This allows rapid identification of interesting or suspicious HSPs differing from others.

### Example output ideal case
<img src="html_ba_1.png" width="695px" />

The black arrows points to the HSP header (color indicating homology) and control buttons respectively.

### Example output with  loaded NCBI sequence viewer
<img src="html_ba_2.png" width="824px" />

### Notes
#### Control buttons
At the bottom of the view there are general control buttons which allow
  selecting and deselecting of sequences and structures, sorting by E-value
  and bulk initialization of NCBI Sequence Viewer.
- Select/Unselect all Seqs. (will select/unselect (check checkbox) all estimated full-length sequences)
- Select/Unselect all Structs (will select/unselect (check checkbox) all predicted structures)
- Export sel. Structs (will trigger download of selected predicted secondary structures in fasta-like format)
- Export sel. Seqs (will trigger download of selected estimated full-length sequences in fasta-like format)
- Sort Eval desc/asc (will sort BLAST HSPs according to E-value Ascending or Descending)
- View all Regions (will trigger loading of NCBI viewer for all (not yet loaded) HSPs)

#### Report structure

1) Inputs: query input file and BLAST input file
  Best matching mode from Rfam

2) Estimated full-length sequences, predicted secondary structures and other data

3) Command and parameters
  - executed commandline string
  - date and time of run
  - parameters

#### The fasta-like format containing secondary structures
```
>uid:N|ACCESSION.VERSIONdirection-method_name START-END (genome location)
SEQUENCE
SECONDARY_STRUCTURE

# - where the N is serial number of BLAST HSP
# - direcion can be "fw" for plus strand and "rc" for minus strand
# - genome-location is then location of found sequence on original genome
#   in START-END format where START is always lower index then END (direction is defined by "direction")
# the sequence is always 5' to 3' direction

>uid:104|CP006976.1fw-rfam-Rc 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((.(((..(((((((..((((.((((((.(....((((......(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))...))))....).)))))).))))...))))))).))))))))))))))......
>uid:104|CP006976.1fw-rnafold 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((((((.((((((......)))......................(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))(((((((((....)))))))))......))).)))..))))))))))))......
```
#### The NCBI sequence Viewer
The NCBI sequence viewer works only if internet connection is available.
It may take some time to load (especially with large genomes) and when the report
 contains many BLAST hits it may require more substantial amount of RAM.
The data for the sequence viewer are not saved across browser sessions.

## References

1. Velandia-Huerto, C. A., Gittenberger, A. A., Brown, F. D., Stadler, P. F., & Bermúdez-Santana, C. I. (2016). Automated detection of ncRNAs in the draft genome sequence of a colonial tunicate: the carpet sea squirt Didemnum vexillum. BMC genomics, 17(1), 1-15.   
2. Freyhult, E. K., Bollback, J. P., & Gardner, P. P. (2007). Exploring genomic dark matter: a critical assessment of the performance of homology search methods on noncoding RNA. Genome research, 17(1), 117-125.   
3. Roshan, U., Chikkagoudar, S., & Livesay, D. R. (2008). Searching for evolutionary distant RNA homologs within genomic sequences using partition function posterior probabilities. BMC bioinformatics, 9(1), 1-9.   
4. Mount, S., & Nguyen, M.-C. (2006, December 14). BLASTN parameters for noncoding queries. Retrieved May 16, 2022, from http://stevemount.outfoxing.com/Posting0004.html  
5. Korf, I., Yandell, M., & Bedell, J. (2003). Blast (1st ed.). O’Reilly Media.”  


## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund - Project "ELIXIR-CZ: Budování kapacit" (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](logolink_OP_VVV_hor_barva_eng.jpg)
