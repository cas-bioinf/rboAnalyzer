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
pred_method_required_tools = {
    # ===== All possible prediction methods =====
    #  list only methods which are not required by default
    'rfam_rnafoldc': {'rfam'},
    'rfam_subopt': {'mfold'},
    'rfam_rapidshapes': {'rfam', 'rapidshapes'},
    'clustalo_alifold_rapidshapes': {'RNAalifold', 'rapidshapes'},
    'muscle_alifold_rapidshapes': {'RNAalifold', 'rapidshapes'},
    'rcoffee_alifold_rapidshapes': {'rcoffee', 'RNAalifold', 'rapidshapes'},
    'alifold_refold': {'RNAalifold', 'refold.pl'},
    'muscle_alifold_refold': {'RNAalifold', 'refold.pl'},
    'rnafold': set(),
    'subopt_fold_query': {'mfold', 'RNAdistance'},
    'subopt_fold_clustal_alifold': {'mfold', 'RNAalifold', 'RNAdistance'},
    'subopt_fold_muscle_alifold': {'mfold', 'RNAalifold', 'RNAdistance'},
    'alifold_refold_rnafold_c': {'RNAalifold', 'refold.pl'},
    'muscle_alifold_refold_rnafold_c': {'RNAalifold', 'refold.pl'},
    'alifold_unpaired_conserved_refold': {'RNAalifold', 'refold.pl'},
    'muscle_alifold_unpaired_conserved_refold': {'RNAalifold', 'refold.pl'},
    'dh_tcoffee_alifold_refold': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_tcoffee_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_tcoffee_alifold_conserved_ss_rnafoldc': {'RNAalifold', 'refold.pl', 'tcoffee'},
    'dh_clustal_alifold_refold': {'RNAalifold', 'refold.pl'},
    'dh_clustal_alifold_refold_rnafoldc': {'RNAalifold', 'refold.pl'},
    'dh_clustal_alifold_conserved_ss_rnafoldc': {'RNAalifold', 'refold.pl'},
    'pairwise_centroid_homfold': {'centroid_homfold'},
    'TurboFold_fast': {'turbofold'},
    'TurboFold': {'turbofold'},
    'tcoffee_rcoffee_alifold_refold': {'rcoffee', 'RNAalifold', 'refold.pl'},
    'tcoffee_rcoffee_alifold_refold_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl'},
    'tcoffee_rcoffee_alifold_conserved_ss_rnafoldc': {'rcoffee', 'RNAalifold', 'refold.pl'},
}
pred_params = {
    'cmscore_percent',
    'cmscore_tr',
    'pred_sim_threshold',
    'query_max_len_diff',
    'repred_unpaired',
    'conseq_conserved',
    'mfold',
    'muscle',
    'clustalo',
    'alifold',
    'turbofold_mode',
    'tcoffee_rcoffee',
    'tcoffee_profile',
    'cmalign',
    'rcoffee',
    'shape_level',
    'rapidshapes',
    'RNAfold',
    'centroid_homfold',
    'max_seqs_in_prediction',
}