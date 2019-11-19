blast_minimal_version = [2, 6, 0]
blast_maximal_version = [3, 0, 0]
locarna_minimal_version = [1, 9, 2]
locarna_maximal_version = [1, 999, 999]
infernal_minimal_version = [1, 1, 2]
vrna_minimal_version = [2, 3, 5]
clustalo_minimal_version = [1, 2, 4]
muscle_minimal_version = [3, 8, 31]
centroid_homfold_minimal_version = [0, 0, 15]
turbofold_minimal_version = [6, 0]
mfold_minimal_version = [3, 6]


# ===== All possible prediction methods =====
#  list only methods which are not required by default
method_required_tools = {
    'locarna': {'locarna', 'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle', 'refold.pl', 'rfam'},
    'simple': {'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle', 'rfam'},
    'meta': {'locarna', 'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle', 'refold.pl', 'rfam'},
    'rfam-Rc': {'rfam'},
    'rfam-sub': {'mfold', 'rfam'},
    'rnafold': set(),
    'fq-sub': {'mfold', 'RNAdistance'},
    'C-A-sub': {'mfold', 'RNAalifold', 'RNAdistance'},
    'M-A-sub': {'mfold', 'RNAalifold', 'RNAdistance'},
    'C-A-r-Rc': {'RNAalifold', 'refold.pl'},
    'M-A-r-Rc': {'RNAalifold', 'refold.pl'},
    'C-A-U-r-Rc': {'RNAalifold', 'refold.pl'},
    'M-A-U-r-Rc': {'RNAalifold', 'refold.pl'},
    'centroid': {'centroid_homfold'},
    'Turbo-fast': {'turbofold'},
    'TurboFold': {'turbofold'},
    'centroid-fast': {'centroid_homfold'},
    'rfam-centroid': {'centroid_homfold', 'rfam'}
}

prediction_methods = set(method_required_tools.keys()) - {'simple', 'locarna', 'meta'}

pred_params = {
    'cmscore_percent',
    'cmscore_tr',
    'pred_sim_threshold',
    'query_max_len_diff',
    'repred_unpaired_tr',
    'conseq_conserved',
    'mfold',
    'muscle',
    'clustalo',
    'alifold',
    'turbofold_mode',
    'cmalign',
    'RNAfold',
    'centroid_homfold',
    'max_seqs_in_prediction',
    'n_seqs',
}

allowed_params = {
    'rfam-Rc': {
        'RNAfold', 'cmscan', 'cmalign'
    },
    'rfam-sub': {
        'mfold'
    },
    'rnafold': {
        'RNAfold'
    },
    'fq-sub': {
        'mfold', 'RNAfold'
    },
    'C-A-sub': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "clustalo", "alifold", "mfold"
    },
    'M-A-sub': {
        "muscle", "alifold", "mfold"
    },
    'C-A-r-Rc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "clustalo", "alifold", "RNAfold"
    },
    'M-A-r-Rc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "muscle", "alifold", "RNAfold"
    },
    'C-A-U-r-Rc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "clustalo", "alifold", "RNAfold"
    },
    'M-A-U-r-Rc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "muscle", "alifold", "RNAfold"
    },
    'centroid': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "centroid_homfold"
    },
    'Turbo-fast': {
        "query_max_len_diff", "max_seqs_in_prediction", "TurboFold"
    },
    'TurboFold': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "max_seqs_in_prediction", "TurboFold"
    },
    'centroid-fast': {
        "query_max_len_diff", "centroid_homfold", "max_seqs_in_prediction"
    },
    'rfam-centroid': {
        "n_seqs", "cmscan", "cmemit", "centroid_homfold"
    }
}

commandline_params = {"alifold", "RNAfold", "cmscan", "cmemit", "TurboFold", "clustalo", "muscle"}

allowed_pm_params_server = {k: allowed_params[k] - commandline_params for k in allowed_params.keys()}