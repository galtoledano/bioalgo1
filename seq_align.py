import argparse
import numpy as np
from itertools import groupby
GAP = "-"

CONVERT_BASE_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}


def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def read_score(file):
    """reading the score table and process it to pandas matrix"""
    # score = pd.read_csv(file, sep='\t', index_col=0)
    return np.genfromtxt(file, delimiter='\t', skip_header=1, usecols=range(1, 6))


def score(score_matrix, i, j, ind_i, overlap):
    """
    gets two letters, returning it's score
    :param score_matrix:  the matrix with all scores
    :param i: the first letter
    :param j: the second letter
    :return: the score of i, j.
    """
    if overlap and (ind_i == 0):
        return 0
    return score_matrix[i][j]


def align(col, row, seq2, seq1, score_matrix, globaling, overlap):
    """
    Calculates the optimal sequence alignment according to the selected alignment type
    :param col: the numbers of letters at the first sequence
    :param row: the numbers of letters at the second sequence
    :param seq2: the second sequence
    :param seq1: the first sequence
    :param score_matrix: the matrix with all scores
    :param globaling: Decides the type of alignment. If true- global, if false- local
    :param overlap: Decides the type of alignment. If true- overlaping
    :return:
    """
    #  init the matrixes
    values = np.zeros((row, col), dtype=int)
    pointers = np.zeros((row, col), dtype=int)
    if globaling:
        init(col, pointers, row, score_matrix, seq1, seq2, values, overlap)
    build_matrix(col, pointers, row, score_matrix, seq1, seq2, values, globaling, overlap)
    # print(values)
    eli1, eli2, s = traceback(col, pointers,values, row, seq1, seq2, globaling, overlap)
    printall(eli1, eli2, pointers, values)
    return s


def build_matrix(col, pointers, row, score_matrix, seq1, seq2, values, globaling, overlap):
    """
    building the matrix according the formula
    :param col: the numbers of letters at the first sequence
    :param pointers: the pointer's matrix
    :param row: the numbers of letters at the second sequence
    :param score_matrix: the matrix with all scores
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param values: the value's matrix
    :param globaling: Decides the type of alignment. If true- global, if false- local
    :param overlap: Decides the type of alignment. If true- overlap.
    :return:
    """
    for i in range(1, row):
        print(i)
        for j in range(1, col):
            d = values[i - 1][j - 1] + score(score_matrix, CONVERT_BASE_TO_INT[seq1[i - 1]],
                                             CONVERT_BASE_TO_INT[seq2[j - 1]], i, overlap)  #diagonal
            h = values[i][j - 1] + score(score_matrix, CONVERT_BASE_TO_INT[GAP], CONVERT_BASE_TO_INT[seq2[j - 1]],
                                         i, overlap)  #horizontal
            v = values[i - 1][j] + score(score_matrix, CONVERT_BASE_TO_INT[seq1[i - 1]], CONVERT_BASE_TO_INT[GAP],
                                         i, overlap)  #vertical
            arr = np.array([v, h, d])
            if not globaling:
                arr = np.append(arr, 0)  #at local alignment, also can be 0
            pointers[i][j] = np.argmax(arr) + 1
            values[i][j] = np.max(arr)




def init(col, pointers, row, score_matrix, seq1, seq2, values, overlap):
    """
    Initializes two arrays. the first array used to hold formula results values. The second array is to save pointers
    for trace back the sequence.
    :param col: the numbers of letters at the first sequence
    :param pointers: the pointer's matrix
    :param row: the numbers of letters at the second sequence
    :param score_matrix: the matrix with all scores
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param values: the value's matrix
    :return:
    """
    for i in range(row):
        pointers[i][0] = 1
        values[i][0] = score(score_matrix, CONVERT_BASE_TO_INT[seq1[i - 1]], CONVERT_BASE_TO_INT[GAP], i, overlap) * i
    for j in range(col):
        values[0][j] = score(score_matrix, CONVERT_BASE_TO_INT[GAP], CONVERT_BASE_TO_INT[seq2[j - 1]], 0, overlap) * j
        pointers[0][j] = 2
    pointers[0][0] = 0



def printall(eli1, eli2, pointers, values):
    """
    prints the final alignments
    :param eli1:  the first alignment
    :param eli2: the second alignment
    :return:
    """
    eli1 = eli1[::-1]
    eli2 = eli2[::-1]
    # print("val : ")
    # print(values)
    # print("pointers: ")
    # print(pointers)
    # print("eli1 :")
    count = 0
    while len(eli1) > 50:
        print(''.join(eli1[:50]))
        eli1 = eli1[50:]
        print(''.join(eli2[:50]))
        eli2 = eli2[50:]
        print("")
    if len(eli1) >= 1:
        print(''.join(eli1))
        print(''.join(eli2))
        print("")


def traceback(col, pointers,values, row, seq1, seq2, glob, overlap):
    """
    trace back according to the pointers matrix and getting the best scored sequence

    :param col: the numbers of letters at the first sequence
    :param pointers: the pointer's matrix
    :param values: the value's matrix
    :param row: the numbers of letters at the second sequence
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param glob: Decides the type of alignment. If true- global, if false- local
    :return: the two final alignments
    """
    align1 = []
    align2 = []
    i = row - 1
    j = col - 1
    if not glob:
        # option to have an sub sequence
        vals = np.array(values)
        res = np.where(vals == np.amax(vals))
        res = np.array(res)
        i = res[0][0]
        j = res[1][0]
        print(pointers)
    if overlap:
        vals = np.array(values)
        res = np.where(vals == np.amax(vals[i]))
        res = np.array(res)
        j = res[1][0]
    p = pointers[i][j]
    s = values[i][j]
    # print("pointers:")
    # print(pointers)
    if overlap and (j < col - 1):
        align2 = [GAP] * (col - j - 1)
        temp = list(seq2)
        align1 = temp[j:col-1][::-1]
    while p != 0:
        if not glob and not overlap and p == 4 :
            break
        if int(p) == 1:
            align1.append(GAP)
            align2.append(seq1[i - 1])
            i -= 1
            p = pointers[i][j]
        elif int(p) == 2:
            align1.append(seq2[j - 1])
            align2.append(GAP)
            j -= 1
            p = pointers[i][j]
        elif int(p) == 3:
            align1.append(seq2[j - 1])
            align2.append(seq1[i - 1])
            i -= 1
            j -= 1
            p = pointers[i][j]
    return align1, align2, s



def main():
    parser = argparse.ArgumentParser()
    seq_a = parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    seq_b = parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    # align(4, 5, "AGC", "AAAC", score, True, False)  # to remember, len +1
    # align(8, 4, "AAAGCCG", "CCG", score, True, True)  # to remember, len +1
    # align(5, 6, "AAGA", "TTAAG", score, False, False)  # to remember, len +1

    command_args = parser.parse_args()
    seq_a = np.array(list(fastaread(command_args.seq_a).__next__()[1]))
    seq_b = np.array(list(fastaread(command_args.seq_b).__next__()[1]))
    score = read_score(command_args.score) #todo
    if command_args.align_type == 'global':
        s = align(len(seq_a) + 1, len(seq_b) + 1, seq_a, seq_b, score, True, False)
    elif command_args.align_type == 'local':
        s = align(len(seq_a) +1, len(seq_b) + 1, seq_a, seq_b, score, False, False)
    elif command_args.align_type == 'overlap':
        s = align(len(seq_a) +1, len(seq_b) + 1, seq_a, seq_b, score, True, True)
    # print the best alignment and score
    print(command_args.align_type + " : " + str(s))


if __name__ == '__main__':
    main()
