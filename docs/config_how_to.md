# The configuration file
 The configuration file contains pipeline wide settings such are paths to tools
 and databases used.

 It should be used if non-default installation path is used.

 ## Configuration
 Configuration is done using the file `config.txt` placed in
 `INSTALL_LOCATION/rna_blast_analyze/BR_CORE/config.txt` or in custom location
 and running the pipeline with the `--config_file PATH_to_custom_LOCATION` argument.

 ## Specification

There are 3 sections

- TOOL_PATHS

  defines path to executable parent directory

- DATA

  defines paths to databases

- SSH

  used for development installation of RapidShapes inside virtual machine and
  made accessible through SSH (with sshpass).

Each section is specified by its name in square brackets.

## Example
 ```
[TOOL_PATHS]
clustal = /usr/bin/
[DATA]
rfam_dir = /home/user/rfamdir/
 ```

## Available settings
 Setting name with executable(s) which should be accessible from provided location (in brackets).

### TOOL_PATHS
 - refold (`refold.pl`)
 - infernal (`cmalign`, `cmbuild`, `cmfetch`, `cmscan`)
 - muscle (`muscle`)
 - clustal (`clustalo`)
 - locarna (`locarna`)
 - viennarna_bin (`RNAfold`, `RNAalifold`, `RNAplot`)
 - mfold (`hybrid-ss-min`)
 - centroid (`centroid_homfold`)
 - turbofold (`TurboFold`)
 - rnashapes
 - rapidshapes

### DATA
 - rsearch_ribosum

  specify full path to RSEARCH matrix file (default: `RIBOSUM65.mat`).

 - rfam_dir

  custom path where Rfam database should be stored.

 - rfam_url

  custom url from where the `rfam.cm.gz` (database dump) should be downloaded.
  This url is also checked when update is requested (wget timestamp is checked).
  (default: `INSTALL_LOCATION/rna_blast_analyze/3rd_party_source/rfamdb/`)

 - html_template_dir

  custom location for html jinja2 template (default: `INSTALL_LOCATION/rna_blast_analyze/BR_core/output/`)

 - rnastructure_datapath

  datapath for RNAstructure (see installation notes for RNAstructure)
