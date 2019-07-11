# Prediction methods in more detail

Structure of this document:

First parameters and their meaning is described.

Second, prediction methods are listed. Each prediction method is described
  (what and how the structures are predicted).
For each method, there is section called "parameters",
  which serves as an example for parameter definition in json format.
Parameters are specified for each method independently on other prediction methods.

## Parameters and their meaning

### Third party tools commandline argument specification
The pipeline offer a way to interact with most parts of tools used for
  secondary structure prediction.

These tools are used for specific tasks and the parameters controlling
  such behaviour are always set and are not exposed for user control.
  (io flags, mode flags, etc.)


#### RNAfold: "RNAFOLD PARAMETERS"
is string with commandline arguments for `RNAfold`
  (see RNAFold [documentation](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)).
  It must be specified with double quotes.

Default: No parameters specified.

#### alifold: "ALIFOLD PARAMETERS"
is string with commandline arguments for `RNAalifold`
  (see RNAalifold [documentation](https://www.tbi.univie.ac.at/RNA/RNAalifold.1.html)).
It must be specified with double quotes.

Default: No parameters specified.

#### cmscan: "CMSCAN PARAMETERS"
is string with commandline arguments for `cmscan`
  (see cmscan [documentation](http://eddylab.org/infernal/)).
  It must be specified with double quotes.

Default: No parameters specified.

#### cmalign: "CMALIGN PARAMETERS"
is string with commandline arguments for `cmalign`
  (see cmalign [documentation](http://eddylab.org/infernal/)).
  It must be specified with double quotes.

Default: No parameters specified.

#### mfold: \[P, W, X\]
specifies parameters for `hybrid-ss-min` program where P, W, and X are integers
  see UNAFold [documenation](http://unafold.rna.albany.edu/?q=unafold-man-pages/hybrid-ss-min).
  It must be specified without double quotes.

The hybrid-ss-min is additionally used with `--suffix=DAT --NA=RNA --noisolate` parameters.

#### clustalo: "CLUSTALO PARAMETERS"
is string of commandline arguments for `clustalo`
  (see clustal omega [documentation](www.clustal.org/omega))

Default: No parameters specified.

#### clustalo_profile: "CLUSTALO_PROFILE"
is string of commandline argument for `clustalo`. Allows specification
  of different `clustalo` parameters for the profile alignment stage.
  (see clustal omega [documentation](www.clustal.org/omega))

Default: `clustalo` called with `--profile1` only profile align compatible option can be specified.

#### centroid_homfold: "CENTROID_HOMFOLD PARAMETERS"
is string of commandline argument for `centroid_homfold`
  (see centroid rna package [documentation](https://github.com/satoken/centroid-rna-package/)).

By default `-g -1` is used. If `-g` or `-t` is specified by user then the `-g -1` is not used.

#### muscle: "MUSCLE PARAMETERS"
is a string argument for `muscle` aligner. For more details, see muscle
  [documentation](https://www.drive5.com/muscle/manual/).

By default `-seqtype rna -quiet -clwstrict` are used and cannot be changed.

### Selection of estimated full-length sequences to provide reference for secondary structure prediction (reference sequences)
Parameters used for selection of sequences similar to the query sequence.
The selection parameters influence how much each sequence can differ from covariance model,
  how much the sequences can be similar to each other and what is maximum
  accepted length difference between the estimated full-length sequence and the query sequence.

Selection of estimated full-length sequences:
All estimated full-length sequences are filtered based on cm bit-score,
  then the remaining sequences are filtered according to their length to
  be within specified length difference from query sequence.
  After that, the remaining sequences are filtered based on sequence similarity with each other
  (from pairwise sequence identity), so only those sequences with similarity
  less then defined similarity threshold are retained.

#### cmscore_percent: CMSCORE_PERCENT
Defines inclusion threshold in percent of bit-score value obtained by
  aligning query sequence to covariance model.
It serves as filter for not sufficiently related sequences.
The inclusion threshold is computed by `CMSCORE_PERCENT`&nbsp;*&nbsp;`QUERY_BITSCORE`&nbsp;/&nbsp;`100`,
  only estimated full-length sequences with CM alignment bit-score higher then inclusion
  threshold are considered for further computation.
`CMSCORE_PERCENT` ranges from 0 to 100, allowed values are integers.

The higher the threshold, the more conservative setting for trusted sequences
  (i.e. more similarity to cm is required).

#### query_max_len_diff: MAX_LEN_DIFF_to_QUERY
Defines a maximum length difference between the query sequence and the estimated full-length sequence.
This serves complementary to cmscore_percent to allow setting low cmscore_percent,
while preventing very short or very long estimated full-length sequences (for any reason)
 to be part of selected sequences set.

`MAX_LEN_DIFF_to_QUERY` ranges from 0 to 1, allowed valuse are floating point (with decimal dot).

The higher the threshold, the more difference is allowed.

#### pred_sim_threshold: PRED_SIM_THRESHOLD
Defines exclusion threshold in percent of sequence similarity.
This serves as a protection from populating set of selected sequences with
  too many too similar sequences (as it may skew alignment and other prediction methods).
`PRED_SIM_THRESHOLD` ranges from 0 to 100, allowed values are integers.

The query is included as first sequence and all sequences which are similar to it
  down to defined threshold are excluded. Then next sequence from the remaining sequences is
  analyzed in same manner.

The higher the similarity the more similar sequences are accepted.
  (i.e. `PRED_SIM_THRESHOLD = 100` means that only exact duplicate sequences are removed)

#### repred_unpaired_tr: REPRED_UNPAIRED_TR
Serves for prediction methods where the conserved unpaired positions in alignment are taken as constraints.
Defines how much MSA column must be conserved at to denote the position
  as single-strand constraint for `RNAfold`.
This parameter is always used together with conseq_conserved.

#### conseq_conserved: CONSEQ_CONSERVED
Serves for prediction method where the conserved unpaired positions in alignment are taken as constraints.
Defines how many bases in a row must be conserved to denote the position
  as singlestrand constraint for `RNAfold`.
This parameter is always used together with repred_unpaired_tr.

#### max_seqs_in_prediction: MAX_SEQS_IN_PREDICTION
Used with TurboFold to set required (and maximum) number of sequences in
  prediction. Allowed values are integers >= 2. Beware of setting this too high.
  Then the prediction is __very__ memory and time expensive.


## rnafold
Predict secondary structure for each estimated full-length sequence with RNAFold.

parameters:
```
    {"rnafold" : {
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## rfam-Rc
Obtain the related covariance model (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), align estimated full-length sequences to it (`cmalign`),
  then take conserved bp as constrains for for folding with `RNAfold -C`.

parameters:
```
  {"rfam-Rc": {
      "RNAfold": "RNAFOLD PARAMETERS",
      "cmscan": "CMSCAN PARAMETERS",
      "cmalign": "CMALIGN PARAMETERS"
    }
  }
```

Non optional parameters:
- `RNAFold`: `-C` (constraints)
- `cmscan`: `-g` (global aligment of models to query sequence)
- `cmalign`: `--notrunc` (disable truncated hits detection - speeds up computation).

## rfam-sub
Obtain the related CM (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), extract reference secondary structure,
  run `hybrid-ss-min` (UNAFold) for suboptimal structures and select the most similar
  one to the one from CM.

parameters:
```
    {"rfam-sub" : {
        "mfold": [P, W, X]
        }
    }

```

## fq-sub
Predict structure of query with RNAFold and take it as reference structure,
  then predict suboptimal structures with `hybrid-ss-min` (UNAFold) and select
  structure most similar to the predicted structure of query (by `RNAdistance` score).

parameters:
```
    {"fq-sub" : {
        "RNAfold": "RNAFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## C-A-sub
clustalo - RNAalifold - hybrid-ss-min
Select reference estimated full-length sequences, align them with `clustalo` (Clustal Omega),
  predict consensus structure with `RNAalifold`.
For each estimated full-length sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure (by `RNAdistance` score).

parameters:
```
    {"C-A-sub" : {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## M-A-sub
clustalo - RNAalifold - hybrid-ss-min
Select reference estimated full-length sequences, align them with `muscle`,
  predict consensus structure with `RNAalifold`.
For each estimated full-length sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure (by `RNAdistance` score).

parameters:
```
    {"M-A-sub" : {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## C-A-r-Rc
clustalo - RNAalifold - refold.pl - rnafold-C
Select reference estimated full-length sequences, compute alignment with `clustalo` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences
  and add the consensus structure to all estimated full-length sequences from profile alignment,
  finally run `refold.pl` and `RNAFold` for the result.

parameters:
```
    {"C-A-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## C-A-U-r-Rc
clustalo - RNAalifold - unpaired conserved - refold.pl - rnafold-C
Select reference estimated full-length sequences, compute alignment with `clustalo` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences,
  and add the consensus structure to all estimated full-length sequences from profile alignment.
Then select conserved parts of alignment where the consensus secondary
  structure annotation is single-strand (i.e. no base-pairs) and use them as constrains for `RNAFold -c`.

parameters:
```
    {"C-A-U-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```



## M-A-r-Rc
muscle - RNAalifold - refold.pl - rnafold -C
Select reference estimated full-length sequences, compute alignment with `muscle` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences
  and add the consensus structure to all estimated full-length sequences from profile alignment,
  finally run `refold.pl` and `RNAFold` for the result.

parameters:
```
    {"M-A-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## M-A-U-r-Rc
muscle - RNAalifold - unpaired conserved - refold.pl - rnafold-C
Select reference estimated full-length sequences, compute alignment with `muscle` (MSA1)
  and then compute consensus secondary structure with `RNAalifold` from the MSA1.
Then compute profile alignment of MSA1 with all estimated full-length sequences,
  and add the consensus structure to all estimated full-length sequences from profile alignment.
Then select conserved parts of alignment where the consensus secondary
  structure annotation is single-strand (i.e. no base-pairs) and use them as constrains for `RNAFold -c`.

parameters:
```
    {"M-A-U-r-Rc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## centroid
centroid_homfold
Select reference estimated full-length sequences and pass them to `centroid_homfold` as
  homologous sequences, predict secondary structure for all estimated full-length sequences.
The `centroid_homfold` is by default called with `-g -1` parameter.
With this setting the `centroid_homfold` predicts multiple structures
 (see [cetroid_homfold documentation](https://github.com/satoken/centroid-rna-package/)),
 the structure with best score is chosen.

parameters:
```
    {"centroid": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS"
        }
    }
```

## centroid-fast
With centroid-fast the reference estimated full-length sequences are selected as follows:
 we take non-redundant estimated full-length sequences without ambiguous
 bases with cm_score > 0 within allowed length-to-query-length difference.
Then we take up to `N` reference sequences as homologous sequences for `centroid_homfold`.
 The sequences are added in order of the original BLAST HSPs, starting with query.
For this method there is also parameter preset avalible which sets `N` to 1
 meaning that only query will be used as homologous sequence in prediction.

parameters:
```
    {"centroid-fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS",
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        }
    }
```

## rfam-centroid
Use covariance model (CM) to generate set of reference sequences for `centroid_homfold`.
By default the highest scoring CM in Rfam is found by `cmscan`
 (or provided one is used `--cm_file` option).
Then requested number of sequences (`N_SEQS`) is generated from the covariance model with `cmemit`.
The sequences are used as homologous sequences for `centroid_homfold`.
The generated sequences can differ between runs. If repeatable behaviour
 is desired the `cmemit` can be seeded (see it's options).

parameters:
```
    {"rfam-centroid": {
        "n_seqs": N_SEQS,
        "cmscan": "CMSCAN PARAMETERS",
        "cmemit": "CMEMIT PARAMETERS",
        }
    }
```

## Turbo-fast
With Turbo-fast the reference estimated full-length sequences are selected as follows:
we take non-redundant estimated full-length sequences without ambiguous
 bases with cm_score > 0 within allowed length-to-query-length difference.
For each estimated full-length sequence, we make non-redundant group of sequences
 consisting of the estimated full-length sequence for which we want to predict secondary structure
 and up to `N`-1 reference sequences.
The sequences are added in order of the original BLAST HSPs, starting with query.
That means that if `N` is 2 and query does not contain ambiguous bases,
 then each estimated full-length sequence secondary structure is computed with query sequence as a reference.
This setting is also available as commandline argument with `--turbo_fast_preset` flag and
 will override the `"max_seqs_in_prediction"` value from prediction_parameters file.

The `N` can be defined in parameters as `"max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION`.

parameters:
```
    {"Turbo-fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        "TurboFold": "TurboFold PARAMETERS"
        }
    }
```

## TurboFold
Select reference estimated full-length sequences without ambiguous bases.
Then continue as with Turbo-fast.

parameters:
```
    {"TurboFold": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        "TurboFold": "TurboFold PARAMETERS"
        }
    }
```
