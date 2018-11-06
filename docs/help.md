# RNA BLAST ANALYZE pipeline
## Introduction
This pipeline is meant as a complementary step when using BLAST algorithm
 to search non-coding RNA (ncRNA) with secondary structure.
As the BLAST is general sequence alignment algorithm, it's results (output)
 is missing some very useful features in context of ncRNAs.

Also the output apart from well scoring hits, usually contains sequence
 fragments (subject sequence only partially aligned to query).
Such hit may bring valuable information by capturing distant homology or
it could be just ...

With our pipeline we add information to such BLAST search to help researcher
 decide which hits are real ncRNAs and what their secondary structure might be.

## Functionality overview
<img src="RBA_pipeline_overview.svg" width="700px" />

The pipeline has 3 stages:
1) Hit extension.
2) Inference of homology of extended hit to query sequence.
3) Secondary structure prediction.

Each of these stages has dedicated section.

### Hit extension
The pipeline has 3 modes for extending BLAST hits.

1) __simple__

    This means that location of extended sequence is computed from
     unaligned parts of query sequence.

2) __locarna__

    In this mode, the hit loci at the subject sequence is realigned to
     the query sequence with Locarna algorithm. The subject sequence
     aligned to the query is considered as extended sequence.

3) __joined__

    Here the two aforementioned mode are combined and extended sequences
     are scored with RSEARCH covariance model. The better scoring sequence
     is chosen.

#### Extension - simple
<img src="figure_blast_simple.svg" width="700px" />


In this mode we compute the extended hit location by taking length of
 unaligned query sequence (can be at start, end or both) and add or
 subtract it respectively from hit start/end index.
In the toy example we have the Plus/Plus BLAST hit with query sequence
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
The Locarna algorithm utilises secondary structure in it's computations,
 thus it is more suited to align ncRNAs.
The Locarna algorithm is by default called with `struct-local=0`,
 `sequ-local=0` and `free-endgaps=++++` parameters.
In addition the information from BLAST hits is utilized in form of
 providing so called anchor to the Locarna algorithm.
The anchor defines columns of alignment which are considered aligned.
As the anchor we consider consecutive series of matches of length at
 least `L` in BLAST alignment.
The default value of `L` is 7.
This way the alignment is anchored and the Locarna algorithm can align
 query to the _supersequence_. With the `free-endgaps=++++` option,
 the algorithm does not put penalty to unaligned ends of _supersequence_.

#### Extension - joined
This approach combines the __simple__ and __locarna__.

### Homology inference
Here we compute score for relation between extended sequence and query sequence.
The computation is based on aligning covariance model (CM) to each extended
 sequence with `cmalign` program from the Infernal package.

We implement 3 options on how to get covariance model:

1) build with RSEARCH (default)

    By default we build the covariance model from the query sequence and RIBOSUM65 matrix.
    The RIBOSUM file used can be changed by [alternative](config_how_to.md) `config.txt` file.

2) supply your own model (the `--cm_file` option)

    If the covariance model is known it can be provided with `--cm_file` option.
    Only one model per file is allowed (this is not checked by the pipeline).

3) infer from Rfam (the `--use_rfam` option)

    The Rfam database is searched with query sequence for the best matching
    model (`cmscan`).

### Secondary structure prediction
In the secondary structure prediction the pipeline can use multiple approaches (prediction methods).
The prediction methods can be (roughly) divided to following groups:

- Predict structure independently on other extended sequences
    The advantage here is robustness to possible improper parameter choice.
    - [rnafold](prediction_methods.md#rnafold)
    - [subopt_fold_query](prediction_methods.md#subopt_fold_query)
    - [rfam_rnafoldc](prediction_methods.md#rfam_rnafoldc)
    - [rfam_subopt](prediction_methods.md#rfam_subopt)

- Use _trusted_ extended sequences as reference
    - [centroid_homfold](prediction_methods.md#centroid_homfold)
    - [TurboFold](prediction_methods.md#TurboFold)
    - [TurboFold_fast](prediction_methods.md#TurboFold_fast)

- Use _trusted_ extended sequences to build consensus secondary structure
    - [alifold_refold_rnafold_c](prediction_methods.md#alifold_refold_rnafold_c)
    - [muscle_alifold_refold_rnafold_c](prediction_methods.md#muscle_alifold_refold_rnafold_c)
    - [tcoffee_rcoffee_alifold_refold_rnafoldc](prediction_methods.md#tcoffe_rcoffee_alifold_refold_rnafoldc)
    - [alifold_unpaired_conserved_refold](prediction_methods.md#alifold_unpaired_conserved_refold)
    - [muscle_alifold_unpaired_conserved_refold](prediction_methods.md#muscle_alifold_unpaired_conserved_refold)
    - [tcoffee_rcoffee_alifold_conserved_ss_rnafoldc](prediction_methods.md#tcoffee_rcoffee_alifold_conserved_ss_rnafoldc)
    - [subopt_fold_clustal_alifold](prediction_methods.md#subopt_fold_clustal_alifold)
    - [subopt_fold_muscle_alifold](prediction_methods.md#subopt_fold_muscle_alifold)

## Output

### Output formats
The pipeline is able to produce several output formats, most handy being
  being the `.html`.
- html
    Stand-alone web page containing sequences and predicted secondary structures.
    If internet connection is avalible, it can be used to view respective
    genome loci for each blast hit using NCBI SeqViewer.
- json
    Complete json-readable pipeline output (all other output can be produced from this one.).
- csv
    Output table in comma separated values. Contains all important information
    including original hit specification, extended sequence location,
    sequence and secondary structure itself.

## HTML output
The html output is organized around BLAST output. Each blast hit gets
it's separate section with five parts:

  1) the text representation of BLAST hit

  2) pipeline report with extended sequence indices and RSEARCH bit score

  3) the extended sequence itself (checkbox or avalible for direct copy)

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
