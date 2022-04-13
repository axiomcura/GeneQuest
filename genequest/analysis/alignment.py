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


def trace_back(score_matrix: np.ndarray) -> list:
    """

    Parameters
    ----------
    query : str
        query sequence
    contig : str
        contig seq
    score_matrix : np.ndarray


    Returns
    -------
    list
        list containing contig_start, contig_end, query_start,
        query_end, alignment score
    """

    # initialize score and position tacker
    score = None
    track_back = []
    score_list = []

    # search index position in the matrix that contains the highest score
    start_x_pos, start_y_pos = np.unravel_index(
        score_matrix.argmax(), shape=score_matrix.shape
    )
    top_score = score_matrix[start_x_pos, start_y_pos]
    score = top_score
    track_back.append((start_x_pos, start_y_pos))
    score_list.append(top_score)

    # trace back by searching for the next score repeat until score == 0
    x = start_x_pos
    y = start_y_pos
    while score > 0:
        tracker = {}

        # getting values from scoring matrix
        left_pos = score_matrix[x - 1, y]
        up_pos = score_matrix[x, y - 1]
        diag_pos = score_matrix[x - 1, y - 1]

        # priority conditions
        # -- Top: want to prioritize diagonal position
        # -- second prior: want to prioritize which positions has the second highest score
        if diag_pos == left_pos and diag_pos == up_pos:
            left_pos = left_pos - 0.1
            up_pos = up_pos - 0.1
        elif diag_pos == left_pos:
            up_pos = left_pos - 0.1
        elif diag_pos == up_pos:
            up_pos = up_pos - 0.1
        elif up_pos == left_pos:
            up_x, up_y = (x, y - 1)
            up_max_score = max(
                score_matrix[up_x, up_y - 1],
                score_matrix[up_x - 1, up_y],
                score_matrix[up_x - 1, up_y - 1],
            )

            left_x, left_y = (x - 1, y)
            left_max_score = max(
                score_matrix[left_x, left_y - 1],
                score_matrix[left_x - 1, left_y],
                score_matrix[left_x - 1, left_y - 1],
            )

            if up_max_score > left_max_score:
                left_pos = left_pos - 0.1
            else:
                up_pos = up_pos - 0.1

        # tracking values
        tracker[f"{x-1}, {y}"] = left_pos
        tracker[f"{x}, {y-1}"] = up_pos
        tracker[f"{x-1}, {y-1}"] = diag_pos

        # finding the max
        max_score_key = max(tracker, key=tracker.get)
        max_score = tracker[max_score_key]
        score = max_score
        score_list.append(max_score)

        # setting new start position
        pos = tuple(int(i) for i in max_score_key.split(","))
        x = pos[0]
        y = pos[1]
        track_back.append((x, y))

    return (track_back, score_list)


def parse_traceback_scores(traceback_data):
    """parses traceback scores and returns a list containing alignment position and
    scores

    Parameters
    ----------
    traceback_data : tuple
        positions, scores obtained from the traceback function

    Returns
    --------
    List
       [contig_beg, contig_end, query_beg, query_end, aln score]
    """
    positions = traceback_data[0]
    scores = traceback_data[1]

    contig_indx = sorted([contig_pos for _, contig_pos in positions[:-1]])
    query_indx = sorted([query_pos for query_pos, _ in positions[:-1]])

    contig_beg, contig_end = min(contig_indx), max(contig_indx)
    query_beg, query_end = min(query_indx), max(query_indx)
    aln_score = sum(scores)

    return [contig_beg, contig_end, query_beg, query_end, aln_score]


def match_scoring(nuc1, nuc2, match_score, mismatch_score):
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


def score_alignment(
    query: str, contig: str, gap=-5, match=10, mismatch=-4
) -> np.ndarray:
    """Do a local alignment between x and y"""

    # create a zero-filled matrix
    score_matrix = generate_scoring_matrix(len(contig) + 1, len(query) + 1)

    # populating matrix
    for i in range(1, len(contig) + 1):
        for j in range(1, len(query) + 1):

            # getting best score per row
            # -- we are maximizing the score
            score_matrix[i][j] = max(
                score_matrix[i][j - 1] + gap,
                score_matrix[i - 1][j] + gap,
                score_matrix[i - 1][j - 1]
                + match_scoring(
                    contig[i - 1],
                    query[j - 1],
                    match_score=match,
                    mismatch_score=mismatch,
                ),
                0,
            )

    # convert into pandas dataframe
    return score_matrix


def covert_alignment_to_pandas(
    query: str, contig: str, score_matrix: str
) -> pd.DataFrame:
    """Converts local alignment scoring matrix into pandas

    Parameters
    ----------
    query : str
        query sequence
    contig : str
        contig sequence
    score_matrix : np.ndarray
        scoring matrix

    Returns
    -------
    pd.DataFrame
        local alignment score into a pandas dataframe
    """
    col_idx = ["*"] + [n for n in query]
    indx = ["*"] + [n for n in contig]
    alignment_df = pd.DataFrame(score_matrix, columns=col_idx, index=indx)
    return alignment_df


def run_local_alignment(contig, query, match=10, mismatch=-4, gap=-5):
    pass
