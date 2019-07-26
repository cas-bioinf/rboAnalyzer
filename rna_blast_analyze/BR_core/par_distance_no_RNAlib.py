import argparse
from multiprocessing import Pool
from subprocess import check_output

from rna_blast_analyze.BR_core.config import CONFIG


def f_parser():
    """
    Input parser
    :return: args structure
    """
    parser = argparse.ArgumentParser(description='BIF analysis')
    parser.add_argument('-s', help='file with structure')
    parser.add_argument('-c', help='file with structure to compare the rest with')
    args = parser.parse_args()

    return args


def vypocet(oi):
    """
    runs RNAdistance without need for RNAlib
    """
    try:
        ret = check_output(
            [
                '{}RNAdistance'.format(CONFIG.viennarna_path)
            ],
            input='{}\n{}'.format(oi[0], oi[1]).encode()
        )

        dist = int(ret.decode().strip().split(':')[1])

        return dist
    except ChildProcessError as e:
        print('RNAdistance failed. Returning large number ({})'.format(10**10))
        return 10**10


def the_main(fp):
    """input is list of dictionaries with fields
    name, seq, structure, consensus
    where consensus is the structure to which you want to compare the other one"""
    # dd = []
    # for i in fp:
    #     d = vypocet(i)
    #     dd.append(d)

    pool = Pool()
    if round(len(fp)/pool._processes) == 0:
        distances = pool.map(vypocet, fp)
    else:
        distances = pool.map(vypocet, fp, round(len(fp)/pool._processes))
    pool.close()
    return distances


def compute_distances(fp, threads=1):
    if threads == 1:
        dist = []
        for pair in fp:
            dist.append(vypocet(pair))
        return dist
    else:
        with Pool(processes=threads) as pool:
            return pool.map(vypocet, fp)


def two_files_input(fasta_structures, fasta_reference):
    fid = open(fasta_reference, 'r')
    lin = fid.read().splitlines()
    c_structure = lin[2].rsplit(' ')[0]
    fid.close()

    # read the file and parse args
    fp = []

    f = open(fasta_structures, 'r')
    c = 0
    names = []
    while True:
        c += 1
        k = f.readline()
        if k == "":
            break
        if k == '\n':
            continue
        k = k.rstrip('\n')

        if '>' == k[0]:
            names.append(k)
            seq = f.readline().rstrip('\n')
            structure = f.readline().rstrip('\n').rsplit(' ')[0]
            oi = (c_structure, structure)
        else:
            names.append(str(c))
            seq = k
            structure = f.readline().rstrip('\n').rsplit(' ')[0]
            oi = (c_structure, structure)
        fp.append(oi)
    f.close()

    distances = the_main(fp)
    return distances


if __name__ == '__main__':
    """
    begin computation and print structures
    now with multiple structures to fold at once
    """
    args = f_parser()

    distances = two_files_input(args.s, args.c)

    for j in range(len(distances)):
        print(distances[j])
