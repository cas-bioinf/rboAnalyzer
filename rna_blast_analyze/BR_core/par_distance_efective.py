import RNA
import argparse
from multiprocessing import Pool


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
    accepts list or tuple with two structures in dot bracket notation e.g. "(((...)).).."
    """

    c_exp = RNA.expand_Full(oi[0])    # consensus
    c_tree = RNA.make_tree(c_exp)

    s_exp = RNA.expand_Full(oi[1])    # structure
    s_tree = RNA.make_tree(s_exp)

    dist = int(RNA.tree_edit_distance(c_tree, s_tree))

    RNA.free_tree(c_tree)
    RNA.free_tree(s_tree)
    return dist


def the_main(fp):
    """input is list of dictionaries with fields
    name, seq, structure, consensus
    where consensus is the structure to which you want to compare the other one"""
    # dd = []
    # for i in fp:
    #     d = vypocet(i)
    #     dd.append(d)

    with Pool() as pool:
        distances = pool.map(vypocet, fp)
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


def RNAfold(seq):
    return RNA.fold(seq)


def compute_RNAfold(sequences, threads=1):
    if threads == 1:
        return [RNAfold(seq) for seq in sequences]
    else:
        with Pool(processes=threads) as pool:
            return pool.map(RNAfold, sequences)


if __name__ == '__main__':
    """
    begin computation and print structures
    now with multiple structures to fold at once
    """
    args = f_parser()

    distances = two_files_input(args.s, args.c)

    for j in range(len(distances)):
        print(distances[j])
