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

    for i in range(row - 1):
        for j in range(col - 1):
    #     prior_row = values[row - 1, :]
    #     no_gap = np.add(prior_row[:-1], score_matrix[seq1[row - 1], seq2])
    #     seq2_aligned_with_gap = np.add(prior_row[1:], score_matrix[GAP, seq2])
    #     seq1_aligned_with_gap = np.zeros(col)
    #     seq1_aligned_with_gap[0] = values[row][0] + score_matrix[seq1[row - 1], GAP]
    #     for j in range(col):
    #         options = [no_gap[col - 1], seq1_aligned_with_gap[col - 1], seq2_aligned_with_gap[col - 1]]
    #         values[row][col] = max(options)
    #         seq1_aligned_with_gap[col] = values[row][col] + score_matrix[seq1[row - 1], GAP]
    #         pointers[row][col] = options.index(max(options))
            match_score = values[i][j] + score_matrix[seq1[i], seq2[j]]
            # if type == "overlap" and i == len(seq_a) - 1:
            #     a_gap_score = M[i + 1][j]
            # else:
            a_gap_score = values[i][j + 1] + score_matrix[seq1[i], CONVERT_BASE_TO_INT[GAP]]
            b_gap_score = values[i + 1][j] + score_matrix[CONVERT_BASE_TO_INT[GAP], seq2[j]]
            if match_score >= a_gap_score and match_score >= b_gap_score:
                values[i + 1][j + 1] = match_score
                pointers[i + 1][j + 1] = 3
            elif a_gap_score >= b_gap_score:
                values[i + 1][j + 1] = a_gap_score
                pointers[i + 1][j + 1] = 2
            else:
                values[i + 1][j + 1] = b_gap_score
                pointers[i + 1][j + 1] = 1
            # if type == "local":
                # if M[i + 1][j + 1] < 0:
                #     M[i + 1][j + 1] = 0
                #     pointer[i + 1][j + 1] = [0, 0]

    # for i in range(1, row):
    #     # d_row = values[i-1]
    #     # v_row = values[i-1]
    #     for j in range(1, col):
    #         if overlap and i == 0:
    #             d = values[i - 1][j - 1]
    #             h = values[i][j - 1]
    #             v = values[i - 1][j]
    #         else:
    #             d = values[i - 1][j - 1] + score_matrix[CONVERT_BASE_TO_INT[seq1[i-1]]][CONVERT_BASE_TO_INT[seq2[j-1]]]
    #             # score(score_matrix, CONVERT_BASE_TO_INT[seq1[i - 1]],
    #             #                                  CONVERT_BASE_TO_INT[seq2[j - 1]], i, overlap)  #diagonal
    #             h = values[i][j - 1] + score_matrix[CONVERT_BASE_TO_INT[GAP]][CONVERT_BASE_TO_INT[seq2[j-1]]]
    #                 # score(score_matrix, CONVERT_BASE_TO_INT[GAP], CONVERT_BASE_TO_INT[seq2[j - 1]],
    #                 #                          i, overlap)  #horizontal
    #             v = values[i - 1][j] + score_matrix[CONVERT_BASE_TO_INT[seq1[i-1]]][CONVERT_BASE_TO_INT[GAP]]
    #                 # score(score_matrix, CONVERT_BASE_TO_INT[seq1[i - 1]], CONVERT_BASE_TO_INT[GAP],
    #                 #                          i, overlap)  #vertical
    #             # arr = np.array([v, h, d])
    #         arr = [v, h, d]
    #         if not globaling:
    #             arr.append(0)
    #         values[i][j] = max(arr)
    #         pointers[i][j] = arr.index(max(arr)) + 1
