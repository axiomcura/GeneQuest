# ------------------------------
# alignment.py
# Erik Serrano
# erik.serrano@cuanschutz.edu
#
# Module containing functions for local sequence
# alignment
# ------------------------------
import numpy as np
import pandas as pd


def generate_scoring_matrix(n_rows, n_cols):
    """Generates a zero-matrix

    Parameters
    ----------
    n_rows : int
        number of rows
    col_size : int
        number of columns
    """
    matrix = np.zeros((n_rows, n_cols))
    return matrix


def trace_back(query, contig, core_matrix, top_score, score_loc):
    """

    Parameters
    ----------
    query : _type_
        _description_
    contig : _type_
        _description_
    core_matrix : _type_
        _description_
    top_score : _type_
        _description_
    score_loc : _type_
        _description_
    """
    # search index position in the matrix that contains the highest score

    # trace back by searching for the next score

    # repeat until reading "0" as the highest score

    # return alignment and postional data
    pass


def match_scoring(nuc1, nuc2, match_score=10, mismatch_score=-5):
    """Calculates the score when matching two nucleotide

    Parameters
    ----------
    nuc1 : str
        nucleotide 1
    nuc2 : str
        nucleotide 2
    match_score : int, optional
        matching score, by defauly 10
    mismatch_score : int, optional
        mismatch score, by default -5

    Returns
    -------
    int
        score
    """
    if nuc1 == nuc2:
        return match_score
    else:
        return mismatch_score


def local_align(query, contig, gap=-5, match=10, mismatch=-5):
    """Do a local alignment between x and y"""

    # create a zero-filled matrix
    score_matrix = generate_scoring_matrix(len(contig) + 1, len(query) + 1)

    # tracking alignment score and position
    best_score = 0
    best_score_loc = (0, 0)

    # populating matrix
    for i in range(1, len(contig) + 1):
        for j in range(1, len(query) + 1):

            # getting best score per row
            # -- we are maximizing the score
            score_matrix[i][j] = max(
                score_matrix[i][j - 1] + gap,
                score_matrix[i - 1][j] + gap,
                score_matrix[i - 1][j - 1] + match_scoring(contig[i - 1], query[j - 1]),
                0,
            )

            # tracking largest score
            if score_matrix[i][j] >= best_score:
                best_score = score_matrix[i][j]
                best_score_loc = (i, j)

    # convert into pandas dataframe
    col_idx = ["*"] + [n for n in query]
    indx = ["*"] + [n for n in contig]
    alignment_df = pd.DataFrame(score_matrix, columns=col_idx, index=indx)
    alignment_df.to_csv("alignment.csv")

    # obtain the sequence and positional data
    aligned_seq = trace_back(query, contig, score_matrix, best_score, best_score_loc)

    # return the opt score and the best location
    return best_score, best_score_loc, score_matrix
