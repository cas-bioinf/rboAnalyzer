pm2test = [
    'TurboFold',
    'alifold_refold',
    'alifold_refold_rnafold_c',
    'alifold_unpaired_conserved_refold',
    # 'clustalo_alifold_rapidshapes',
    'dh_clustal_alifold_conserved_ss_rnafoldc',
    'dh_clustal_alifold_refold',
    'dh_clustal_alifold_refold_rnafoldc',
    'dh_tcoffee_alifold_conserved_ss_rnafoldc',
    'dh_tcoffee_alifold_refold',
    'dh_tcoffee_alifold_refold_rnafoldc',
    # 'muscle_alifold_rapidshapes',
    'muscle_alifold_refold',
    'muscle_alifold_refold_rnafold_c',
    'muscle_alifold_unpaired_conserved_refold',
    'pairwise_centroid_homfold',
    # 'rcoffee_alifold_rapidshapes',
    # 'rfam_rapidshapes',
    'rfam_rnafoldc',
    'rfam_subopt',
    'rnafold',
    'subopt_fold_clustal_alifold',
    'subopt_fold_muscle_alifold',
    'subopt_fold_query',
    'tcoffee_rcoffee_alifold_conserved_ss_rnafoldc',
    'tcoffee_rcoffee_alifold_refold',
    'tcoffee_rcoffee_alifold_refold_rnafoldc',
]

for i in pm2test:
    print(
        """    def test_{}(self):
            self.run('{}')
        """.format(i, i)
    )