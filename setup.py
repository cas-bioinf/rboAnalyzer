from setuptools import setup

with open('rna_blast_analyze/VERSION', 'r') as o:
    version = o.read().strip()

package_data = {
    'rna_blast_analyze': [
        'BR_core/config.txt',
        'BR_core/prediction_parameters.json',
        'BR_core/output/*',
        '3rd_party_source/RSEARCH_matrices/*',
        'docs/*',
        'VERSION'
    ]
}

setup(
    name='rboAnalyzer',
    version=version,
    description='Analyze BLAST output for RNA query',
    author='Marek Schwarz',
    author_email='marek.schwarz@biomed.cas.cz',
    entry_points={
        'console_scripts': [
            'rboAnalyzer = rna_blast_analyze.BA:main',
            'genomes_from_blast = rna_blast_analyze.download_blast_genomes:main'
        ],
    },
    install_requires=[
        'matplotlib<3',
        'numpy>=1.14',
        'Jinja2>=2, <3',
        'pandas>=0.22',
        'biopython>=1.72, <1.74',
        'argcomplete',
    ],
    packages=['rna_blast_analyze', 'rna_blast_analyze.BR_core', 'rna_blast_analyze.BR_core.output'],
    package_data=package_data,
    include_package_data=True,
)
