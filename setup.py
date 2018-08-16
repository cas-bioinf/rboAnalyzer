from setuptools import setup, find_packages

# package_data
# - directly related to code
# - 3rd_party_data as data_files not working to expectations
package_data = {
    'rna_blast_analyze': [
        'BR_core/config.txt',
        'BR_core/prediction_parameters.json',
        'BR_core/blast_bio/LICENSE',
        'BR_core/output/*.html',
        'BR_core/output/style.css',
        '3rd_party_source/RSEARCH_matrices/*',
        '3rd_party_source/refold/*',
    ]
}

setup(
    name='rna_blast_analyze',
    version='0.0.1',
    description='Analyze BLAST output for RNA query',
    author='Marek Schwarz',
    author_email='marek.schwarz@biomed.cas.cz',
    entry_points={
        'console_scripts': [
            'rna_blast_analyze = rna_blast_analyze.BA:main',
        ],
    },
    install_requires=[
        'matplotlib',
        'setuptools',
        'dill',
        'scipy',
        'numpy',
        'Jinja2',
        'pandas',
        'reportlab',
        'seaborn',
        'PyPDF2',
        'biopython',
    ],
    packages=find_packages(),
    package_data=package_data,
    include_package_data=True,
)
