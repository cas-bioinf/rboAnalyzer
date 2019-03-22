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

#### rcoffee: "RCOFFEE PARAMETERS"
is string argument for `t_coffee -mode rcoffee` program. All supplied
  arguments must be compatible with `-mode rcoffee` (a dedicated mode for
  aligning RNAs). More info in `t_coffee` [documentation](https://tcoffee.readthedocs.io/en/latest/).

By default if multiple threads are allowed, `t_coffee` is called with `-n_core=N-threads`.

### Reference set of trusted extended sequences generation parameters
Parameters used for selection of sequences related to the query sequence.
The selection parameters influence how much each sequence can differ from covariance model,
  how much the sequences can be similar to each other and what is maximum
  accepted length difference of extended and query sequence.

Standard trusted extended sequences processing works in this way:
All extended sequences are filtered based on cm bit-score, then selected
  sequences are filtered according to their length to be within specified
  length difference from query sequence. Remaining sequences are filtered
  based on sequence similarity (from pairwise sequence identity), so only
  those sequences with similarity less then defined similarity threshold
  are retained.

#### cmscore_percent: CMSCORE_PERCENT
Defines inclusion threshold in percent of bit-score value obtained by
  aligning query sequence to covariance model.
It serves as filter for not sufficiently related sequences.
The inclusion threshold is computed by `CMSCORE_PERCENT`&nbsp;*&nbsp;`QUERY_BITSCORE`&nbsp;/&nbsp;`100`,
  only extended sequences with CM alignment bit-score higher then inclusion
  threshold are considered for further computation.
`CMSCORE_PERCENT` ranges from 0 to 100, allowed values are integers.

The higher the threshold, the more conservative setting for trusted sequences
  (i.e. more similarity to cm is required).

#### query_max_len_diff: MAX_LEN_DIFF_to_QUERY
Defines a maximum length difference from query sequence for extended sequence.
This serves complementary to cmscore_percent to allow setting low cmscore_percent,
while preventing truncated hits (for any reason) to be part of trusted sequence set.

`MAX_LEN_DIFF_to_QUERY` ranges from 0 to 1, allowed valuse are floating point (with decimal dot).

The higher the threshold, the more difference is allowed.

#### pred_sim_threshold: PRED_SIM_THRESHOLD
Defines exclusion threshold in percent of sequence similarity.
This serves as a protection from populating set of trusted sequences with
  too many too similar sequences (as it may skew alignment and other prediction methods).
`PRED_SIM_THRESHOLD` ranges from 0 to 100, allowed values are integers.

The query is included as first sequence and all sequences which are similar to it
  down to defined threshold are excluded. Then next remaining sequences are
  analyzed in same manner.

The higher the similarity the more similar sequences are accepted.
  (i.e. `PRED_SIM_THRESHOLD = 100` means that only exact duplicate sequences are removed)

#### repred_unpaired_tr: REPRED_UNPAIRED_TR
Serves for prediction method where as constraints are taken the unpaired position.
Defines how much MSA column must be conserved at to denote the position
  as singlestrand constraint for `RNAfold`.
This parameter is always used together with conseq_conserved.

#### conseq_conserved: CONSEQ_CONSERVED
Serves for prediction method where as constraints are taken the unpaired position.
Defines how many bases in a row must be conserved to denote the position
  as singlestrand constraint for `RNAfold`.
This parameter is always used together with repred_unpaired_tr.

#### max_seqs_in_prediction: MAX_SEQS_IN_PREDICTION
Used with TurboFold to set required (and maximum) number of sequences in
  prediction. Allowed values are integers >= 2. Beware of setting this too high.
  Then the prediction is __very__ memory and time expensive.


## rnafold
Predict secondary structure for each extended sequence with RNAFold.

parameters:
```
    {"rnafold" : {
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## rfam_rnafoldc
Obtain the related CM (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), align extended sequences to it (`cmalign`),
  then take conserved bp as constrains for for folding with `RNAfold`.

parameters:
```
  {"rfam_rnafoldc": {
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

## rfam_subopt
Obtain the related CM (either find highest scoring one in Rfam by `cmscan`
  or use the one provided with `--cm_file`), extract reference structure,
  run `hybrid-ss-min` (UNAFold) for suboptimals and select the most similar
  one to the one in CM.

parameters:
```
    {"rfam_subopt" : {
        "mfold": [P, W, X]
        }
    }

```

## subopt_fold_query
Predict structure of query with RNAFold and take it as reference,
  then predict suboptimal structures with `hybrid-ss-min` (UNAFold) and select
  structure most similar to predicted structure of query (by `RNAdistance` score).

parameters:
```
    {"subopt_fold_query" : {
        "RNAfold": "RNAFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## subopt_fold_clustal_alifold
Select trusted extended sequences, align them with `clustalo` (Clustal Omega),
  predict consensus structure with `RNAalifold`.
For each extended sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure (by `RNAdistance` score).

parameters:
```
    {"subopt_fold_clustal_alifold" : {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## subopt_fold_muscle_alifold
Select trusted extended sequences, align them with `muscle`,
  predict consensus structure with `RNAalifold`.
For each extended sequence compute suboptimal structures with `hybrid-ss-min`
  and select structure from predicted suboptimal structures the one most
  similar to the consensus structure.

parameters:
```
    {"subopt_fold_muscle_alifold" : {
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "mfold": [P, W, X]
        }
    }
```

## clustalo_alifold_refold_rnafoldc
Select trusted sequences subset, compute alignment with `clustalo` for
  the trusted sequences, compute consensus structure with `RNAalifold`,
  then compute profile alignment of aligned trusted sequences with all extended sequences,
  then add the consensus structure to all sequences from profile alignment,
  and finally run `refold.pl` and `RNAFold` for the result alignment.

parameters:
```
    {"clustalo_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## clustalo_alifold_unpaired_conserved_refold_rnafoldc
Select trusted sequences subset, compute alignment with `clustalo` in those
  sequences, compute consensus structure with `RNAalifold`,
  compute profile alignment of selected sequences with all sequences,
  match the reference structure to all sequence from profile alignment,
  select conserved parts of alignment that have singlestrand consensus
  structure annotation and use them as constrains for `RNAFold -c`.

parameters:
```
    {"clustalo_alifold_unpaired_conserved_refold_rnafoldc": {
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



## muscle_alifold_refold_rnafoldc
Select trusted sequences subset, compute alignment with `muscle` for
  the trusted sequences, compute consensus structure with `RNAalifold`,
  then compute profile alignment of aligned trusted sequences with all extended sequences,
  then add the consensus structure to all sequences from profile alignment,
  and finally run `refold.pl` and `RNAFold -C`.

parameters:
```
    {"muscle_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "muscle": "MUSCLE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## muscle_alifold_unpaired_conserved_refold_rnafoldc
Select trusted sequences subset, compute alignment with `muscle` for these
  sequences, compute consensus structure with `RNAalifold`,
  compute profile alignment of selected sequences with all sequences,
  match the reference structure to all sequence from profile alignment,
  select conserved parts of alignment that have singlestrand consensus
  structure annotation and use them as constrains for `RNAFold -C`

parameters:
```
    {"muscle_alifold_refold_rnafoldc": {
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

## centroid_homfold
Select trusted sequences. Then pass the selected sequences to `centroid_homfold` as
  homologous sequence set and predict all sequences. The `centroid_homfold`
  is by default called with `-g -1` parameter. With this setting the
  `centroid_homfold` predict multiple structures - the structure with
  best score is chosen.

parameters:
```
    {"centroid_homfold": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS"
        }
    }
```

## centroid_homfold_fast
With centroid_homfold_fast we take non-redundant sequences without ambiguous
 bases with cm_score > 0 within allowed length query difference as trusted sequences.
 We take up to `N` trusted sequences as homologous sequences for centroid_homfold.
 The sequences are added in order of the original BLAST hits, starting with query.
 For this method there is also parameter preset avalible which sets `N` to 1
 meaning that only query will be used as homologous sequence in prediction.

parameters:
```
    {"centroid_homfold_fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "centroid_homfold": "CENTROID_HOMFOLD PARAMETERS",
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        }
    }
```

## rfam_centroid_homfold
Use covariance model (CM) to generate set of sequences considered homologous.
By default the highest scoring CM in Rfam is found by `cmscan`
 (or provided one is used `--cm_file` option). Then requested number of
 sequences is generated under the covariance model with `cmemit` which
 are used as homologous sequences.
 The generated sequences can differ between runs. If repeatable behaviour
 is desired the `cmemit` can be seeded (see it's options).

parameters:
```
    {"rfam_centroid_homfold": {
        "n_seqs": N_SEQS,
        "cmscan": "CMSCAN PARAMETERS",
        "cmemit": "CMEMIT PARAMETERS",
        }
    }
```

## TurboFold_fast
With TurboFold_fast we take non-redundant sequences without ambiguous
 bases with cm_score > 0 within allowed length query difference as trusted sequences.
 For each predicted sequence, we make non-redundant group of sequences
 consisting of the predicted sequence and up to `N` trusted sequences.
 The sequences are added in order of the original BLAST hits, starting with query.
 That means that if `N` is 2 and query does not contain ambiguous bases,
 then each extended sequence secondary structure is computed with query sequence as a reference.
 This setting is also available as commandline argument with `--turbofold_fast_preset` flag and
 will override the `"max_seqs_in_prediction"` value from prediction_parameters file.

 The `N` can be defined in parameters as `"max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION`.

parameters:
```
    {"TurboFold_fast": {
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "max_seqs_in_prediction": MAX_SEQS_IN_PREDICTION,
        "TurboFold": "TurboFold PARAMETERS"
        }
    }
```

## TurboFold
Select trusted sequences (unambiguous). Then continue as with Turbofold_fast.

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

## rcoffee_alifold_refold_rnafoldc
Select trusted sequences subset, compute alignment with `t_coffee -mode rcoffee` for
  the trusted sequences, compute consensus structure with `RNAalifold`,
  then compute profile alignment of aligned trusted sequences with all extended sequences,
  then add the consensus structure to all sequences from profile alignment,
  and finally run `refold.pl` and `RNAFold` for the result alignment.

parameters:
```
    {"rcoffee_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "rcoffee": "RCOFFEE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## rcoffee_alifold_unpaired_conserved_refold_rnafoldc
Select trusted sequences subset, compute alignment with `t_coffee -mode rcoffee` for
  the trusted sequences, compute consensus structure with `RNAalifold`,
  then compute profile alignment of aligned trusted sequences with all extended sequences,
  then match the computed reference structure to all sequence from profile alignment,
  select the conserved parts of alignment that have singlestrand consensus
  structure annotation and use them as constrains for `RNAFold -c`.

parameters:
```
    {"rcoffee_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "rcoffee": "RCOFFEE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## dh_rcoffee_alifold_refold_rnafoldc
Select trusted sequences, then compute alignment with `t_coffee -mode rcoffee` for these
  sequences, compute consensus structure with `RNAalifold`. Then run profile
  alignment separately for trusted and for remaining sequences.
  Then for each part add the consensus structure to sequences from profile alignment,
  and finally run `refold.pl` and `RNAFold -C`.

parameters:
```
    {"dh_rcoffee_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "clustalo": "CLUSTALO PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## dh_rcoffee_alifold_unpaired_conserved_rnafoldc
Select trusted sequences, then compute alignment with `t_coffee -mode rcoffee` for these
  sequences, compute consensus structure with `RNAalifold`. Then run profile
  alignment separately for trusted and for remaining sequences.
  Then for each part select conserved parts of alignment that have
  singlestrand consensus structure annotation and use them as constrains
  for `RNAFold -C`.

parameters:
```
    {"dh_rcoffee_alifold_unpaired_conserved_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "repred_unpaired_tr": REPRED_UNPAIRED_TR,
        "conseq_conserved": CONSEQ_CONSERVED,
        "rcoffee": "RCOFFEE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```

## dh_clustal_alifold_refold_rnafoldc
Select trusted sequences, then compute alignment with `clustalo` for these
  sequences, compute consensus structure with `RNAalifold`. Then run profile
  alignment separately for trusted and for remaining sequences.
  Then for each part add the consensus structure to sequences from profile alignment,
  and finally run `refold.pl` and `RNAFold -C`.

parameters:
```
    {"dh_clustal_alifold_refold_rnafoldc": {
        "cmscore_percent" : CMSCORE_PERCENT,
        "pred_sim_threshold": PRED_SIM_THRESHOLD,
        "query_max_len_diff": MAX_LEN_DIFF_to_QUERY,
        "rcoffee": "RCOFFEE PARAMETERS",
        "alifold": "ALIFOLD PARAMETERS",
        "RNAfold": "RNAFOLD PARAMETERS"
        }
    }
```


## dh_clustal_alifold_unpaired_conserved_rnafoldc
Select trusted sequences, then compute alignment with `clustalo` for these
  sequences, compute consensus structure with `RNAalifold`. Then run profile
  alignment separately for trusted and for remaining sequences.
  Then for each part select conserved parts of alignment that have
  singlestrand consensus structure annotation and use them as constrains
  for `RNAFold -C`

parameters:
```
    {"dh_clustal_alifold_unpaired_conserved_rnafoldc": {
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