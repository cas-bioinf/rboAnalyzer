from setuptools import setup, find_packages

# package_data
# - directly related to code
# - 3rd_party_data as data_files not working to expectations
package_data = {
    'rna_blast_analyze': [
        'BR_core/config.txt',
        'BR_core/prediction_parameters.json',
        'BR_core/output/*.html',
        'BR_core/output/style.css',
        '3rd_party_source/RSEARCH_matrices/*',
    ]
}

setup(
    name='rna_blast_analyze',
    version='0.0.4',
    description='Analyze BLAST output for RNA query',
    author='Marek Schwarz',
    author_email='marek.schwarz@biomed.cas.cz',
    entry_points={
        'console_scripts': [
            'rna_blast_analyze = rna_blast_analyze.BA:main',
            'genomes_from_blast = rna_blast_analyze.download_blast_genomes:main'
        ],
    },
    install_requires=[
        'matplotlib',
        'dill',
        'numpy',
        'Jinja2>=2, <3',
        'pandas>=0.22',
        'PyPDF2',
        'biopython',
        'argcomplete',
    ],
    packages=find_packages(),
    package_data=package_data,
    include_package_data=True,
)
