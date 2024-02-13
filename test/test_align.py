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
    pass
    

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




