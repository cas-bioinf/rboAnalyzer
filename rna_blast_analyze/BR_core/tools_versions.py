blast_minimal_version = [2, 6, 0]
locarna_minimal_version = [1, 9, 2]
infernal_minimal_version = [1, 1, 2]
vrna_minimal_version = [2, 3, 5]
clustalo_minimal_version = [1, 2, 4]
muscle_minimal_version = [3, 8, 31]
tcoffee_rcoffee_minimal_version = [10, 0]
centroid_homfold_minimal_version = [0, 0, 15]
turbofold_minimal_version = [6, 0]
mfold_minimal_version = [3, 6]

# ===== All possible prediction methods =====
#  list only methods which are not required by default
method_required_tools = {
    'locarna': {'locarna', 'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle', 'refold.pl'},
    'simple': {'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle'},
    'joined': {'locarna', 'infernal', 'blastdbcmd', 'clustalo', 'RNAfold', 'RNAplot', 'RNAdistance', 'muscle', 'refold.pl'},
    'rfam_rnafoldc': {'rfam'},
    'rfam_subopt': {'mfold', 'rfam'},
    'rnafold': set(),
    'subopt_fold_query': {'mfold', 'RNAdistance'},
    'subopt_fold_clustal_alifold': {'mfold', 'RNAalifold', 'RNAdistance'},
    'subopt_fold_muscle_alifold': {'mfold', 'RNAalifold', 'RNAdistance'},
    'clustalo_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'muscle_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'clustalo_alifold_unpaired_conserved_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'muscle_alifold_unpaired_conserved_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'dh_rcoffee_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl', 'rcoffee'},
    'dh_rcoffee_alifold_unpaired_conserved_rnafoldc': {'RNAalifold', 'refold.pl', 'rcoffee'},
    'dh_clustal_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'dh_clustal_alifold_unpaired_conserved_rnafoldc': {'RNAalifold', 'refold.pl'},
    'centroid_homfold': {'centroid_homfold'},
    'TurboFold_fast': {'turbofold'},
    'TurboFold': {'turbofold'},
    'rcoffee_alifold_refold_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl'},
    'rcoffee_alifold_unpaired_conserved_refold_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl'},
    'centroid_homfold_fast': {'centroid_homfold'},
    'rfam_centroid_homfold': {'centroid_homfold', 'rfam'}
}

prediction_methods = set(method_required_tools.keys()) - {'simple', 'locarna', 'joined'}

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
    'rcoffee',
    'rcoffee_profile',
    'cmalign',
    'rcoffee',
    'RNAfold',
    'centroid_homfold',
    'max_seqs_in_prediction',
    'n_seqs',
}

allowed_params = {
    'rfam_rnafoldc': {
        'RNAfold', 'cmscan', 'cmalign'
    },
    'rfam_subopt': {
        'mfold'
    },
    'rnafold': {
        'RNAfold'
    },
    'subopt_fold_query': {
        'mfold', 'RNAfold'
    },
    'subopt_fold_clustal_alifold': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "clustalo", "alifold", "mfold"
    },
    'subopt_fold_muscle_alifold': {
        "muscle", "alifold", "mfold"
    },
    'clustalo_alifold_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "clustalo", "alifold", "RNAfold"
    },
    'muscle_alifold_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "muscle", "alifold", "RNAfold"
    },
    'clustalo_alifold_unpaired_conserved_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "clustalo", "alifold", "RNAfold"
    },
    'muscle_alifold_unpaired_conserved_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "muscle", "alifold", "RNAfold"
    },
    'dh_rcoffee_alifold_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "rcoffee", "alifold", "RNAfold"
    },
    'dh_rcoffee_alifold_unpaired_conserved_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "rcoffee", "alifold", "RNAfold", "repred_unpaired_tr", "conseq_conserved"
    },
    'dh_clustal_alifold_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "clustalo", "alifold", "RNAfold"
    },
    'dh_clustal_alifold_unpaired_conserved_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "clustalo", "alifold", "RNAfold"
    },
    'centroid_homfold': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "centroid_homfold"
    },
    'TurboFold_fast': {
        "query_max_len_diff", "max_seqs_in_prediction", "TurboFold"
    },
    'TurboFold': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "max_seqs_in_prediction", "TurboFold"
    },
    'rcoffee_alifold_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "rcoffee", "alifold", "RNAfold"
    },
    'rcoffee_alifold_unpaired_conserved_refold_rnafoldc': {
        "cmscore_percent", "pred_sim_threshold", "query_max_len_diff", "repred_unpaired_tr", "conseq_conserved", "rcoffee", "alifold", "RNAfold"
    },
    'centroid_homfold_fast': {
        "query_max_len_diff", "centroid_homfold", "max_seqs_in_prediction"
    },
    'rfam_centroid_homfold': {
        "n_seqs", "cmscan", "cmemit", "centroid_homfold"
    }
}

commandline_params = {"rcoffee", "alifold", "RNAfold", "cmscan", "cmemit", "TurboFold", "clustalo", "muscle"}

allowed_pm_params_server = {k: allowed_params[k] - commandline_params for k in allowed_params.keys()}