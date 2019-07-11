pm2test = [
    'TurboFold',
    'C-A-r-Rc',
    'C-A-U-r-Rc',
    'M-A-r-Rc',
    'M-A-U-r-Rc',
    'centroid',
    'rfam-Rc',
    'rfam-sub',
    'rnafold',
    'C-A-sub',
    'M-A-sub',
    'fq-sub',
    'centroid-fast',
    'rfam-centroid',
]

for i in pm2test:
    print(
        """    def test_{}(self):
            self.run('{}')
        """.format(i, i)
    )