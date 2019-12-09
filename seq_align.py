#########################################importiong#########################################
import argparse
import numpy as np
from itertools import groupby

##################################### consts ################################################

GAP = "-"

CONVERT_BASE_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}

CONVERT_INT_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-'}

####################################### functions ############################################


def convert(seq1, seq2, d):
    """
    converting by dictionary dna to numbers and  numbers to dna
    :param seq1: the first alignment
    :param seq2: the second alignment
    :param d: the dictionary
    :return: the converted alignments
    """
    keys, choices = list(zip(*d.items()))
    seq_a = np.array(keys)[:, None, None] == seq1
    seq_b = np.array(keys)[:, None, None] == seq2
    seq_1 = np.select(seq_a, choices)[0]
    seq_2 = np.select(seq_b, choices)[0]
    return seq_1, seq_2


def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def read_score(file):
    """reading the score table and process it to pandas matrix"""
    return np.genfromtxt(file, delimiter='\t', skip_header=1, usecols=range(1, 6))


def align(col, row, seq2, seq1, score_matrix, glob, overlap):
    """
    Calculates the optimal sequence alignment according to the selected alignment type
    :param col: the numbers of letters at the first sequence
    :param row: the numbers of letters at the second sequence
    :param seq2: the second sequence
    :param seq1: the first sequence
    :param score_matrix: the matrix with all scores
    :param glob: Decides the type of alignment. If true- global, if false- local
    :param overlap: Decides the type of alignment. If true- overlaping
    :return: the final score
    """
    #  initialize the values and pointers matrixs
    values = np.zeros((row, col), dtype=int)
    pointers = np.zeros((row, col), dtype=int)
    #  converting the dna bases to numbers
    seq1, seq2 = convert(seq1, seq2, CONVERT_BASE_TO_INT)
    if glob:
        # initialize the matrixs if the type is not global.
        pointers, values = init(col, pointers, row, score_matrix, seq1, seq2, values, overlap)
    #  filling the values matrix
    build_matrix(col, pointers, row, score_matrix, seq1, seq2, values, glob, overlap)
    #  converting the number to  dna bases
    seq1, seq2 = convert(seq1, seq2, CONVERT_INT_TO_BASE)
    #  trace back the best alignment
    eli1, eli2, s = traceback(col, pointers, values, row, seq1, seq2, glob, overlap)
    #  prints the alignments
    print_all(eli1, eli2)
    return s


def build_matrix(col, pointers, row, score_matrix, seq1, seq2, values, glob, overlap):
    """
    filling the matrix with the best-fit values.
    :param col: the numbers of letters at the first sequence
    :param pointers: the pointer's matrix
    :param row: the numbers of letters at the second sequence
    :param score_matrix: the scores for each match
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param values: the value's matrix
    :param glob: Decides the type of alignment. If true- global or overlap, if false- local
    :param overlap: Decides the type of alignment. If true-overlap, if false- global or local
    """
    for i in range(1, row):
        #  splicing the latest full row
        upper_line = values[i-1]
        #  vertical values
        v = np.add(upper_line[1:], score_matrix[seq2, CONVERT_BASE_TO_INT[GAP]])
        #  horizontal values
        h = np.zeros(col)
        if overlap and (i == row - 1):
            #  promise no "cost"
            h[0] = values[i][0]
        else:
            h[0] = values[i][0] + score_matrix[CONVERT_BASE_TO_INT[GAP], seq1[i - 1]]
        #  diagonal values
        d = np.add(upper_line[: col-1], score_matrix[seq1[i-1], seq2])

        # choosing the best of three options
        for j in range(1, col):
            arr = [v[j-1], h[j-1], d[j-1]]
            if not glob:
                #  adding option on 0
                arr.append(0)
            values[i][j] = max(arr)
            if overlap and (i == row-1):
                h[j] = values[i][j]
            else:
                h[j] = values[i][j] + score_matrix[CONVERT_BASE_TO_INT[GAP], seq1[i-1]]
            pointers[i][j] = arr.index(max(arr)) + 1


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
    :param overlap: Decides the type of alignment. If true-overlap, if false- global or local
    """
    #  initialize the first column if it's not overlap match.
    if not overlap:
        for i in range(1, row):
            values[i][0] = score_matrix[CONVERT_BASE_TO_INT[GAP], seq1[i-1]] + values[i-1][0]
    for j in range(1, col):
        values[0][j] = score_matrix[CONVERT_BASE_TO_INT[GAP], seq2[j-1]] + values[0][j-1]
    #  initialize the first pointers row and column
    pointers[:, 1:] = 2
    pointers[1:, :] = 1
    pointers[0][0] = 0
    return pointers, values


def print_all(eli1, eli2):
    """
    prints the final alignments
    :param eli1:  the first alignment
    :param eli2: the second alignment
    """
    eli1 = eli1[::-1]
    eli2 = eli2[::-1]
    while len(eli1) > 50:
        print(''.join(eli1[:50]))
        print(''.join(eli2[:50]))
        print("")
        eli1 = eli1[50:]
        eli2 = eli2[50:]
    if len(eli1) >= 1:
        print(''.join(eli1))
        print(''.join(eli2))
        print("")


def traceback(col, pointers, values, row, seq1, seq2, glob, overlap):
    """
    trace back according to the pointers matrix and getting the best scored sequence

    :param col: the numbers of letters at the first sequence
    :param pointers: the pointer's matrix
    :param values: the value's matrix
    :param row: the numbers of letters at the second sequence
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param glob: Decides the type of alignment. If true- global or overlap, if false- local
    :param overlap: Decides the type of alignment. If true-overlap, if false- global or local
    :return: the two final alignments and score
    """
    #  init the aligns and indexes
    align1, align2, i, j, p, s = init_trace(col, glob, overlap, pointers, row, seq2, values)
    #  the trace back
    tracing(align1, align2, glob, i, j, overlap, p, pointers, seq1, seq2)
    return align1, align2, s


def tracing(align1, align2, glob, i, j, overlap, p, pointers, seq1, seq2):
    """
    do the tracing
    :param align1: first alignment
    :param align2: second alignment
    :param glob: flag if the type is global or not
    :param i: the i index
    :param j: the j index
    :param overlap: flag if the type is overlap  or not
    :param p: the position at the pointers matrix
    :param pointers: the pointers matrix
    :param seq1: the first seq
    :param seq2: the second seq
    """
    while p != 0:
        if not glob and not overlap and p == 4:
            break
        if int(p) == 1:
            align1.append(GAP)
            align2.append(seq1[i - 1])
            i -= 1
            p = pointers[i][j]
        elif int(p) == 2:
            align1.append(seq2[j - 1][0])
            align2.append(GAP)
            j -= 1
            p = pointers[i][j]
        elif int(p) == 3:
            align1.append(seq2[j - 1][0])
            align2.append(seq1[i - 1])
            i -= 1
            j -= 1
            p = pointers[i][j]


def init_trace(col, glob, overlap, pointers, row, seq2, values):
    """
    initilize the trace back lists and indexes
    :param col: number od columns
    :param glob: flag if the type is global or not
    :param overlap: flag if the type is overlap or not
    :param pointers: the pointers matrix
    :param row: numbers of rows
    :param seq2: the second seq
    :param values: the values matrix
    :return: both initilized aligns,  i, j - the first traceback index p - first pointer at the pointers matrix,
    s - the final score
    """
    align1 = []
    align2 = []
    i = row - 1
    j = col - 1
    if not glob:
        # option to have an sub sequence
        res = np.where(values == np.amax(values))
        res = np.array(res)
        i = res[0][0]
        j = res[1][0]
    if overlap:
        x = max(values[i])
        res = np.where(values[i] == x)
        res = np.array(res)
        j = res[0][0]
    p = pointers[i][j]
    s = values[i][j]
    #  init the gaps before the first match
    if overlap and (j < col - 1):
        align2 = [GAP] * (col - j - 1)
        align1 = list(seq2[j:col - 1][::-1])
    return align1, align2, i, j, p, s


def main():
    """
    the main function. parsing the progrem's arguments and running the match-algorithm to them.
    :print: the function prints the best match alignments, the matching type and the score
    """
    #  parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()
    seq_a = np.array(list(fastaread(command_args.seq_a).__next__()[1]))
    seq_b = np.array(list(fastaread(command_args.seq_b).__next__()[1]))
    score = np.array((read_score(command_args.score)))
    len_a = len(seq_a) + 1
    len_b = len(seq_b) + 1

    #  running the match- algorithm.
    if command_args.align_type == 'global':
        s = align(len_a, len_b, seq_a, seq_b, score, True, False)
    elif command_args.align_type == 'local':
        s = align(len_a, len_b, seq_a, seq_b, score, False, False)
    elif command_args.align_type == 'overlap':
        s = align(len_b, len_a, seq_b, seq_a, score, True, True)

    # print the best alignment and score
    print(command_args.align_type + " : " + str(s))


if __name__ == '__main__':
    main()

