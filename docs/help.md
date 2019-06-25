# rboAnalyzer
## Introduction
rboAnalyzer is pipeline meant as a complementary step when BLAST algorithm
 was used to search for query that is non-coding RNA (ncRNA) with secondary structure (does not have to be known).
As the BLAST is general sequence alignment algorithm, it's results (output)
 is missing some very useful features in context of ncRNAs.

Also the output apart from well scoring hits, usually contains sequence
 fragments (subject sequence only partially aligned to query).
Such hit may bring valuable information by capturing distant homology or
it may be nothing.

With our rboAnalyzer we add information to such BLAST search to help researcher
 decide which hits are real ncRNAs and what their secondary structure might be.

## Functionality overview
<img src="RBA_pipeline_overview.svg" width="700px" />

The pipeline has 3 stages:
1) Hit completion.
2) Inference of homology of completed hit to query sequence.
3) Secondary structure prediction.

Each of these stages has dedicated section.

### Hit completion
The pipeline has 3 modes for completing BLAST hits to full length of the query.

1) __simple__

    This means that location of completed sequence is computed from
     flanking unaligned parts of query sequence.

2) __locarna__

    In this mode, the loci containing hit with flanking regions at the subject sequence is realigned to
     the query sequence with Locarna algorithm. The sequence
     aligned to the query is considered to be completed sequence.

3) __joined__

    Here the two aforementioned modes are combined and extended sequences
     are scored with covariance model. The better scoring sequence
     is chosen.

#### Extension - simple
<img src="figure_blast_simple.svg" width="700px" />


In this mode we compute the completed seuence location by taking length of
 unaligned parts of query sequence in HSP (can be at start, end or both) and add or
 subtract it respectively from hit start/end index.
In the toy example we have the Plus/Plus BLAST HSP with query sequence
 aligned from 10 to 21 to subject sequence 1000 to 1009.
If query is 50 bases long, then length of unaligned query at start is 9
 and length of unaligned query at end is 29.
Then start of the extended sequence in subject is 1000 - 9 = 991 and end is
 1009 + 29 = 1038.


#### Extension - locarna
<img src="figure_blast_locarna.svg" width="700px" />

With __locarna__ mode we first extract so called _supersequence_. Which is
 region on subject sequence extended by unaligned query region with
 additional regions. This _supersequence_ is then realigned with Locarna
 algorithm to obtain extended sequence.
The Locarna algorithm utilises possible pairings in it's computations,
 thus it is better suited to align ncRNAs then BLAST.
The Locarna algorithm is by default called with `struct-local=0`,
 `sequ-local=0` and `free-endgaps=++++` parameters.
Additionaly the information about matching nucleotides from BLAST HSPs 
is used to construct so called anchor to the Locarna algorithm.
The anchor defines columns of alignment which are considered aligned.
As the anchor we consider consecutive series of matches of length at
 least `L` in BLAST alignment.
The default value of `L` is 7.
This way the alignment is anchored and the Locarna algorithm can align
 query to the _supersequence_. With the `free-endgaps=++++` option,
 the algorithm does not put penalty to unaligned ends of _supersequence_.

#### Extension - joined
This approach combines the __simple__ and __locarna__. It computes both and 
 for each HSPs it chooses completed sequence based on score to covariance model.

### Homology inference
Here we compute score for relation between completed sequence and query sequence.
The computation is based on aligning covariance model (CM) to each extended
 sequence with `cmalign` program from the Infernal package.

We've implemented 3 options on how to provide covariance model:

1) build with RSEARCH (default)

    By default we build the covariance model from the query sequence (secondary structure predicted by RNAfold) and RIBOSUM matrix.
    The RIBOSUM is RIBOSUM65 by default and it can be changed in [alternative](config_how_to.md) `config.txt` file.

2) supply your own model (the `--cm_file` option)

    If the covariance model is known, it can be provided with `--cm_file` option.
    Only one model per file is allowed (this is not checked by the pipeline).
    
    Note that provided covariance model will also be used in all prediction methods using covariance models (starting with `rfam`).

3) infer from RFAM (the `--use_rfam` option)

    The Rfam database is searched with query sequence for the best matching
    model (`cmscan`).

### Secondary structure prediction
The rboAnalyzer can use multiple approaches (prediction methods) to predict secondary structures.
The prediction methods can be (roughly) divided to following groups:

- Predict structure independently on other completed sequences
    The advantage for these methods is robustness to possible improper parameter choice.
    - [rnafold](prediction_methods.md#rnafold)
    - [subopt_fold_query](prediction_methods.md#subopt_fold_query)
    - [rfam_rnafoldc](prediction_methods.md#rfam_rnafoldc)
    - [rfam_centroid_homfold](prediction_methods.md#rfam_centroid_homfold)
    - [rfam_subopt](prediction_methods.md#rfam_subopt)

- Use _trusted_ completed sequences as reference
    - [centroid_homfold](prediction_methods.md#centroid_homfold)
    - [TurboFold](prediction_methods.md#TurboFold)
    - [TurboFold_fast](prediction_methods.md#TurboFold_fast)

- Use _trusted_ completed sequences to build consensus secondary structure
    - [clustalo_alifold_refold_rnafoldc](prediction_methods.md#clustalo_alifold_refold_rnafoldc)
    - [muscle_alifold_refold_rnafoldc](prediction_methods.md#muscle_alifold_refold_rnafoldc)
    - [clustalo_alifold_unpaired_conserved_refold_rnafoldc](prediction_methods.md#clustalo_alifold_unpaired_conserved_refold_rnafoldc)
    - [muscle_alifold_unpaired_conserved_refold_rnafoldc](prediction_methods.md#muscle_alifold_unpaired_conserved_refold_rnafoldc)
    - [subopt_fold_clustal_alifold](prediction_methods.md#subopt_fold_clustal_alifold)
    - [subopt_fold_muscle_alifold](prediction_methods.md#subopt_fold_muscle_alifold)

## Output

### Output formats
The pipeline is able to produce several output formats, most handy being
  being the `.html`.
- html
    Stand-alone web page containing sequences and predicted secondary structures.
    If internet connection is avalible, it can be used to view respective
    genome loci for each BLAST HSP using NCBI SeqViewer.
- json
    Complete json-readable pipeline output (all other output can be produced from this one.).
- csv
    Output table in comma separated values. Contains all important information
    including original HSP data, extended sequence location,
    sequence and secondary structure(s) itself.

## HTML output
The html output is organized around BLAST output. Each blast hit gets
it's separate section with five parts:

  1) the text representation of BLAST hit

  2) pipeline report with extended sequence indices and RSEARCH bit score

  3) the completed sequence itself (checkbox or avalible for direct copy)

  4) one or multiple predicted secondary structures, each with its own checkbox for export

  5) NCBI Sequence viewer (optional - by default only load button is shown)

The `html` otuput offers sorting, selecting sequences and structures and
  their export to fasta format or fasta-like format with dot-bracket secondary structures.
  Also if internet connection is avalible NCBI genome browser is avalible
  to explore synteny and known features of current genome.

### Example output ideal case
<img src="html_ba_1.png" width="695px" />

The black arrows points to the hit header (color indicating homology) and control buttons respectively.

### Example output with  loaded NCBI sequence viewer
<img src="html_ba_2.png" width="824px" />

### Notes
#### Control buttons
At the bottom of the view there are general control buttons which allow
  selecting and deselecting of sequences and structures, sorting by E-value
  and bulk NCBI Sequence Viewer loader.
- Select/Unselect all Seqs. (will select/unselect (check checkbox) all extended sequences)
- Select/Unselect all Structs (will select/unselect (check checkbox) all predicted structures)
- Export sel. Structs (will trigger download of selected structures in fasta-like format)
- Export sel. Seqs (will trigger donwload of selected sequences in fasta-like format)
- Sort Eval desc/asc (will sort BLAST hits according to E-value Ascending or Descending)
- View all Regions (will trigger loading of NCBI viewer for all (not yet loaded) hits)

#### Hit Header
The header for each blast hit contains Accession.Version number (based on provided regexp). The header is also color-coded on color scale from green to red based on
  RSEARCH score. This allows rapid identification of interesting or
  suspicious hits differing from others.

#### Report structures

1) Pipeline name and BLAST input file

2) Extended sequences report

3) Command and parameters
  - executed commandline string
  - date and time of run
  - paramaters

#### The fasta-like format containing secondary structures
```
>uid:N|accessionDIRECTION-METHOD_NAME genome location START-END
SEQUENCE
SECONDARY_STRUCTURE

# where the N is serial number of BLAST hit
#  and DIRECTION can be "fw" for plus strand and "rc" for minus strand
#  genome-location is then location of found sequence on original genome
#  in START-END format where START is always less then END (direction is defined by DIRECTION)
# the sequence is always 5' to 3' direction

>uid:104|CP006976.1fw-rfam_rnafoldc 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((.(((..(((((((..((((.((((((.(....((((......(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))...))))....).)))))).))))...))))))).))))))))))))))......
>uid:104|CP006976.1fw-rnafold 2175683-2175865
GAUUACCUGAGGUGUUUGCCAGUGGGUUAUGUCCCUGAGCCGAUACUUUUAUUUUAUGAAUCGGUUUCUAAUUGUUGGUGUGCAUGCUUAGCUUGACUAAGAAGCCUAAAAAUAGUUAUAACUGAUUCCCUUGAACCGUUGGGUUCAAGGACUGAGACUUGCAGCAGCAUCUCGGGUUCUUCC
....(((((((((((((((.((((((......)))......................(((((((((..((((((((..((.((.(.(((((.....)))))).)).))..)))))))).)))))))))(((((((((....)))))))))......))).)))..))))))))))))......
```
#### The NCBI sequence Viewer
The NCBI sequence viewer works only if Internet connection is available.
It may take some time to load (especially with large genomes) and when the report
 contains many BLAST hits it may required more substantial amount of RAM.
 The data for the sequence viewer are not saved across sessions (after you close the web page), and must be reloaded re-lunch.

## Funding

This work was supported by ELIXIR CZ research infrastructure project (MEYS Grant No: LM2015047) including access to computing and storage facilities.

![elixir logo](ELIXIR_CZECHREPUBLIC_white_background_small.png)

This work was supported from European Regional Development Fund-Project ELIXIR-CZ (No. CZ.02.1.01/0.0/0.0/16_013/0001777).

![msmt logo](logolink_OP_VVV_hor_barva_eng.jpg)