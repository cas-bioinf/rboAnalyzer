import os
from subprocess import call

fwd = os.path.dirname(os.path.dirname(__file__))
base = os.path.dirname(os.path.dirname(fwd))
blast_in = os.path.join(fwd, 'RF00001_short.blastout')
blast_query = os.path.join(fwd, 'RF00001.fasta')
blast_db = os.path.join(fwd, 'blastdb', 'RF00001-art.blastdb')

json_output = os.path.join(fwd, 'RF00001_output.json')

os.chdir(base)

try:
    os.remove(json_output)
except:
    pass

r = call(
    [
        'python3', '-m',
        'rna_blast_analyze.BA',
        '--blast_in', blast_in,
        '--blast_query', blast_query,
        '--blast_db', blast_db,
        '--mode', 'simple',
        '--blast_regexp', r'(?<=\|)[A-Z0-9]*\.?\d*$',
        '--b_type', 'plain',
        '--json', json_output,
        '--prediction_method', 'rnafold'
    ]
)

if r:
    raise ChildProcessError('Computation failed')
print('Computation succeeded.')