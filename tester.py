import sys
import time
import ntpath

import seq_align

SCORE_MATRIX_TSV = "cbio_ex1_tests/score_matrix.tsv"

FILES = ["cbio_ex1_tests/fastas/f1.fasta", "cbio_ex1_tests/fastas/f2.fasta",
         "cbio_ex1_tests/fastas/f3.fasta", "cbio_ex1_tests/fastas/f4.fasta",
         "cbio_ex1_tests/fastas/f5.fasta", "cbio_ex1_tests/fastas/f6.fasta",
         "cbio_ex1_tests/fastas/f7.fasta", "cbio_ex1_tests/fastas/n1.fasta",
         "cbio_ex1_tests/fastas/n2.fasta"]

ALIGN_TYPE = ["global", "local"]
NUMBER_OF_RUNS = 1

if __name__ == '__main__':
    s = sys.argv
    time_counter = []
    score_counter = []
    score_f = open("socre.txt", "w")
    score_f.write("atype	score_mat	f1	f2	score \n")
    for align_type in ALIGN_TYPE:
        visited = set()
        for f1 in FILES:
            for f2 in FILES:
                if (f2, f1) in visited:
                    continue
                this_s = s.copy() + [f1, f2, "--align_type", align_type, "--score",
                                     SCORE_MATRIX_TSV]
                sys.argv = this_s
                start = time.perf_counter()
                for i in range(NUMBER_OF_RUNS):
                    score = seq_align.main()
                end = time.perf_counter()
                time_counter.append(str((end - start) / NUMBER_OF_RUNS))
                score_f.write(
                    "{0} {1} {2} {3} {4} \n".format(align_type, ntpath.basename(SCORE_MATRIX_TSV),
                                                    ntpath.basename(f1), ntpath.basename(f2),
                                                    score))
                visited.add((f1, f2))
    with open("timer.txt", "w") as timer:
        timer.write("\n".join(time_counter))

    score_f.close()
