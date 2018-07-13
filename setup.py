from setuptools import setup, find_packages
# from distutils.core import setup

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
    packages=find_packages()
)
