# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        len_A = len(seqA)
        len_B = len(seqB)

        inf = np.inf

        #intilize empty alignment, gap_A, and gap_B matrix
        self._align_matrix = np.empty((len_A + 1, len_B + 1))
        self._gapA_matrix  = np.empty((len_A + 1, len_B + 1))
        self._gapB_matrix = np.empty((len_A + 1, len_B + 1))

        #intilize empty backtrace matrices
        self._back = np.empty((len_A + 1, len_B + 1), dtype='U4')
        self._back_A = np.empty((len_A + 1, len_B + 1), dtype='U4')
        self._back_B = np.empty((len_A + 1, len_B + 1), dtype='U4')

        #M(0,0) = 0, I_A(0,0) = gap_open + gap_extend, I_B(0,0) = gap_open + gap_extend
        self._align_matrix[0][0] = 0
        self._gapA_matrix[0][0] = self.gap_open + self.gap_extend 
        self._gapB_matrix[0][0] = self.gap_open + self.gap_extend 

        #M(i,0) = M(0,j) = -inf i and j not equal to 0
        #gapA(i,0) = -inf, i not 0
        #gapA(0,j) = I_A(0,j-1) + gap_extend, j not 0
        #gapB(0,j) = -inf, j not 0
        #gapB(i,0) = I_B(i-1,0) + gap_extend, i not 0
        for i in range(1, len_A + 1):
            self._align_matrix[i][0] = -inf

            self._gapA_matrix[i][0] = -inf

            self._gapB_matrix[i][0] = self._gapB_matrix[i-1][0] + self.gap_extend 
            self._back_B[i][0] = 'gapB'

        for j in range(1, len_B + 1):
            self._align_matrix[0][j] = - inf

            self._gapA_matrix[0][j] = self._gapA_matrix[0][j-1] + self.gap_extend
            self._back_A[0][j] = 'gapA'

            self._gapB_matrix[0][j] = -inf
            
        # TODO: Implement global alignment here
        for i in range(1, len_A + 1):
            for j in range(1, len_B + 1):

                #calculate next alignment matrix entry
                s = self.sub_dict[(self._seqA[i-1], self._seqB[j-1])] #get match/mismatch value from sub matrix
                back_M = self._align_matrix[i-1][j-1] + s   #match/mismatch seqA and seqB position
                back_gapA = self._gapA_matrix[i-1][j-1] + s #insertion in seqA
                back_gapB = self._gapB_matrix[i-1][j-1] + s #insertion in seqB
                #Choose the max value of a match,a gap in A, or a gap in B to add to align matrix
                ##Keep track of where this max value comes from
                if back_M > back_gapA and back_M > back_gapB:
                    self._align_matrix[i][j] = back_M
                    self._back[i][j] = 'M'
                elif back_gapA > back_M and back_gapA > back_gapB:
                    self._align_matrix[i][j] = back_gapA
                    self._back[i][j] = 'gapA'
                else:
                    self._align_matrix[i][j] = back_gapB
                    self._back[i][j] = 'gapB'

                #calculate next gapA matrix entry
                back_M = self._align_matrix[i][j-1] + self.gap_open  + self.gap_extend     #open gap in seqA
                back_gapA = self._gapA_matrix[i][j-1] + self.gap_extend  #expand gap in seqA
                #Choose the max value of opening a gap in A or extending a gap in A to add to gapA matrix
                ##Keep track of where this max value comes from
                if back_M > back_gapA:
                    self._gapA_matrix[i][j] = back_M
                    self._back_A[i][j] = 'M'
                else:
                    self._gapA_matrix[i][j] = back_gapA
                    self._back_A[i][j] = 'gapA'

                #calculate next gapB matrix entry
                back_M = self._align_matrix[i-1][j] + self.gap_open + self.gap_extend     #open gap in seqB
                back_gapB = self._gapB_matrix[i-1][j] + self.gap_extend  #xpand gap in seqB
                #Choose the max value of opening a gap in A or extending a gap in A to add to gapA matrix
                ##Keep track of where this max value comes from
                if back_M > back_gapB:
                    self._gapB_matrix[i][j] = back_M
                    self._back_B[i][j] = 'M'
                else:
                    self._gapB_matrix[i][j] = back_gapB
                    self._back_B[i][j] = 'gapB'

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        #Need the length of seqA and B to start at the last mattrix entry
        len_A = len(self._seqA)
        len_B = len(self._seqB)

        #Find the start

        #get the values of the last entry for align, gapA, and gapB matrix
        align_mat_last = self._align_matrix[len_A][len_B]
        gapA_mat_last = self._gapA_matrix[len_A][len_B]
        gapB_mat_last = self._gapB_matrix[len_A][len_B]

        next_step = '' #Will hold which matrix the backtracing should go to next

        #Start at the largest of the last entry values, and set this value as the alignment score
        if align_mat_last > gapA_mat_last and align_mat_last > gapB_mat_last:
            next_step = 'M' #set the first step to match
            self.alignment_score = align_mat_last
        elif gapA_mat_last > align_mat_last and gapA_mat_last > gapB_mat_last:
            next_step = 'gapA' #set the first step to gap in A
            self.alignment_score = gapA_mat_last
        else:
            next_step = 'gapB' #set the first step to gap in B
            self._gapB_matrix_last

        #iterators to go back through the matrices
        i = len_A
        j = len_B

        #Will hold the rev pf the aligned sequences
        seqA_align_r = ''
        seqB_align_r = ''

        #go ack through the matirices until we hit the first entry
        while i > 0 or j > 0:
            #If the next step is match, match the two seqs, get the next step from _back, and update both i and j
            if next_step == 'M':
                next_step = self._back[i][j]
                seqA_align_r = seqA_align_r + self._seqA[i-1]
                seqB_align_r = seqB_align_r + self._seqB[j-1]

                i = i - 1
                j = j - 1

            #If the first step is gap in A, insert a gap in A, get the next step from _back_A, and update j
            elif next_step == 'gapA':
                next_step=self._back_A[i][j]
                seqA_align_r = seqA_align_r + '-'
                seqB_align_r = seqA_align_r + self._seqB[j-1]

                j = j - 1


            #If the next step is gap in B, insert a gap in B, get the next step from _back_B, and update i
            elif next_step == 'gapB':
                next_step=self._back_B[i][j]
                seqA_align_r = seqA_align_r + self._seqA[i-1]
                seqB_align_r = seqB_align_r + '-'

                i = i - 1

            #If the next step is not match, gap in A, or gap in B, something has gone wrong and raise an error
            else:
                raise(ValueError('Invalid step encountered in backtracing, quitting'))

        #set the aligned seqs
        self.seqA_align = seqA_align_r[::-1]
        self.seqB_align = seqB_align_r[::-1]
            

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
