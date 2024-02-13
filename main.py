# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")


    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    #Do all the alignments and store the alignment score for each in a dictionary
    BRD2_alignments = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat',-10, -1)
    seqs = [gg_seq, mm_seq, br_seq, tt_seq]
    species = ['Gallus gallus', 'Mus musculus', 'Balaeniceps rex', 'Tursiops truncatus']

    alignment_scores = {}

    for i in range(len(seqs)):
        hs_align_score = BRD2_alignments.align(hs_seq, seqs[i])[0]
        alignment_scores[species[i]] = hs_align_score

    #Sort the dictonary by values
    sorted_alingments = sorted(alignment_scores.items(), reverse=True)

    #print species and their score from most similar to least similar
    for alignment in sorted_alingments:
        print(alignment)



if __name__ == "__main__":
    main()
