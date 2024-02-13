# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    testNW = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat',-10, -1)
    seq1Seq2Alignment = testNW.align(seq1, seq2)

    #Create the correct alignment matrices
    inf = np.inf
    epsilon = 0.00000001
    matchMatrix = np.array([[0, -inf, -inf, -inf], 
                             [-inf, 5, -12, -14], 
                             [-inf, -13, 4, -8],
                             [-inf, -13, -1, 5],
                             [-inf, -15, -6, 4]])
    gapAMatrix = np.array([[-11, -12, -13, -14], 
                             [-inf, -inf, -6, -7], 
                             [-inf, -inf, -24, -7],
                             [-inf, -inf, -24, -12],
                             [-inf, -inf, -26, -17]])
    gapBMatrix = np.array([[-11, -inf, -inf, -inf], 
                             [-12, -inf, -inf, -inf], 
                             [-13, -6, -23, -25],
                             [-14, -7, -7, -19],
                             [-15, -8, -8, -6]])
    
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):

            assert testNW._align_matrix[i][j] >= matchMatrix[i][j] - epsilon and testNW._align_matrix[i][j] <= matchMatrix[i][j] + epsilon

            assert testNW._gapA_matrix[i][j] >= gapAMatrix[i][j] - epsilon and testNW._gapA_matrix[i][j] <= gapAMatrix[i][j] + epsilon

            assert testNW._gapB_matrix[i][j] >= gapBMatrix[i][j] - epsilon and testNW._gapB_matrix[i][j] <= gapBMatrix[i][j] + epsilon
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    testNW = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat',-10, -1)
    seq3Seq4Alignment = testNW.align(seq3, seq4)
    epsilon = 0.00000001
    assert seq3Seq4Alignment[0] >= 17 - epsilon and seq3Seq4Alignment[0] <= 17 + epsilon
    assert seq3Seq4Alignment[1] == 'MAVHQLIRRP'
    assert seq3Seq4Alignment[2] == 'M---QLIRHP'




