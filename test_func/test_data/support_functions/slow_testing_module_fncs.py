pm2test = [
    'TurboFold',
    # 'alifold_refold',
    'clustalo_alifold_refold_rnafoldc',
    'clustalo_alifold_unpaired_conserved_refold_rnafoldc',
    # 'clustalo_alifold_rapidshapes',
    'dh_clustal_alifold_unpaired_conserved_rnafoldc',
    # 'dh_clustal_alifold_refold',
    'dh_clustal_alifold_refold_rnafoldc',
    'dh_rcoffee_alifold_unpaired_conserved_rnafoldc',
    # 'dh_tcoffee_alifold_refold',
    'dh_rcoffee_alifold_refold_rnafoldc',
    # 'muscle_alifold_rapidshapes',
    # 'muscle_alifold_refold',
    'muscle_alifold_refold_rnafoldc',
    'muscle_alifold_unpaired_conserved_refold_rnafoldc',
    'centroid_homfold',
    # 'rcoffee_alifold_rapidshapes',
    # 'rfam_rapidshapes',
    'rfam_rnafoldc',
    'rfam_subopt',
    'rnafold',
    'subopt_fold_clustal_alifold',
    'subopt_fold_muscle_alifold',
    'subopt_fold_query',
    'rcoffee_alifold_unpaired_conserved_refold_rnafoldc',
    # 'tcoffee_rcoffee_alifold_refold',
    'rcoffee_alifold_refold_rnafoldc',
]

for i in pm2test:
    print(
        """    def test_{}(self):
            self.run('{}')
        """.format(i, i)
    )