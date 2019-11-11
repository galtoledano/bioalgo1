import argparse
import numpy as np
import pandas as pd
from itertools import groupby

GAP = "-"


def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def read_score(file):
    score = pd.read_csv(file, sep='\t', index_col=0)
    return score


def score(score_matrix, i, j):
    return score_matrix[i][j]


def global_align(col, row, seq2, seq1, score_matrix):
    values = np.zeros((row, col), dtype=int)
    pointers = np.zeros((row, col), dtype=int)
    init(col, pointers, row, score_matrix, seq1, seq2, values)
    build_matrix(col, pointers, row, score_matrix, seq1, seq2, values)
    eli1, eli2 = traceback(col, pointers, row, seq1, seq2)

    printall(eli1, eli2, pointers, values)


def build_matrix(col, pointers, row, score_matrix, seq1, seq2, values):
    for i in range(1, row):
        for j in range(1, col):
            d = values[i - 1][j - 1] + score(score_matrix, seq1[i - 1], seq2[j - 1])
            h = values[i - 1][j] + score(score_matrix, GAP, seq2[j - 1])
            v = values[i][j - 1] + score(score_matrix, seq1[i - 1], GAP)
            arr = np.array([v, h, d])
            pointers[i][j] = np.argmax(arr) + 1
            values[i][j] = np.max(arr)


def init(col, pointers, row, score_matrix, seq1, seq2, values):
    for i in range(row):
        values[i][0] = score(score_matrix, seq1[i - 1], GAP) * i
        pointers[i][0] = 1
    for j in range(col):
        values[0][j] = score(score_matrix, GAP, seq2[j - 1]) * j
        pointers[0][j] = 2
    pointers[0][0] = 0


def printall(eli1, eli2, pointers, values):
    print("val : ")
    print(values)
    print("pointers: ")
    print(pointers)
    print("eli1 :")
    for i in range(len(eli1)):
        print(eli1[len(eli1) - i - 1], end=" ")
    print("eli2 :")
    for i in range(len(eli2)):
        print(eli2[len(eli2) - i - 1], end=" ")


def traceback(col, pointers, row, seq1, seq2):
    eli1 = []
    eli2 = []
    i = row - 1
    j = col - 1
    p = pointers[i][j]
    while p != 0:
        if p == 1:
            eli1.append(GAP)
            eli2.append(seq2[j - 1])
            j -= 1
            p = pointers[i][j]
        elif p == 2:
            eli1.append(seq1[i - 1])
            eli2.append(GAP)
            i -= 1
            p = pointers[i][j]
        elif p == 3:
            eli1.append(seq1[i - 1])
            eli2.append(seq2[j - 1])
            i -= 1
            j -= 1
            p = pointers[i][j]
    return eli1, eli2


def main():
    score = read_score("score_matrix.tsv")
    global_align(4, 5, "AGC", "AAAC", score)
# parser = argparse.ArgumentParser()
    # parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    # parser.add_argument('seq_b', help='Path to second FASTA file')
    # parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    # parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    # command_args = parser.parse_args()
    # if command_args.align_type == 'global':
    #     raise NotImplementedError
    # elif command_args.align_type == 'local':
    #     raise NotImplementedError
    # elif command_args.align_type == 'overlap':
    #     raise NotImplementedError
    # print the best alignment and score


if __name__ == '__main__':
    main()
