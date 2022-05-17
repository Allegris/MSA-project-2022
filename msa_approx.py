import sys
import numpy as np
import time
import project2_linear as pa #pairwise alignment
import msa_common
import msa_sp_score_3k as msa #score for approx
import fasta_and_phylip as fp # Helper functions for reading/writing/parsing fasta and phylip files



'''
Computes a multiple sequence alignment (MSA) of the input nodes by:

- Computing the center string

- Uses the corresponding star tree as a "guide tree" when doing an approximation algorithm
  (Gusfield's approximation algorithm)

Returns the MSA
'''
def MSA_approx(S, center, sub_matrix, gap_cost):
    M = []
    S.remove(center)
    for s in S:
        A_matrix = pa.fill_table(center, s, sub_matrix, gap_cost)
        # optimal alignment
        A = pa.construct_alignment(A_matrix, center, s, sub_matrix, gap_cost)
        if(s != S[0]):
            M = msa_common.extend_M(M, A, row_idx = 0) # center string is row 0
        else:
            M = A
    return M


'''
Finds the center string of multiple sequences.
The center string is the sequences with lowest alignments score to the other sequences.
'''
def find_center_string(S, sub_matrix, gap_cost):
    # Contains pairwise distances from s to s
	score_matrix = np.full((len(S), len(S)), None)
    # Distances from s to s is 0
	for i in range(len(S)):
		score_matrix[i, i] = 0
    # Iterate through possible centers, S[i]
	for i in range(len(S)):
        # Score for S[i]
		sum_scores = 0
        # Iterate through all other strings, S[j]
		for j in range(len(S)):
            # If we have NOT already computed the distance from S[i] to S[j], do this
			if(score_matrix[i, j] == None):
				score = pa.fill_table(S[i], S[j], sub_matrix, gap_cost)[len(S[i]), len(S[j])]
                # Distance from S[i] to S[j] is equal to the distance from S[j] to S[i]
				score_matrix[i, j] = score
				score_matrix[j, i] = score
			sum_scores += score_matrix[i, j]
    # Calculate total scores for S[i]'s
	total_scores = [sum(score) for score in score_matrix]
    # The min score index
	best_score_index = np.argmin(total_scores)
    # Return string with min score (= center)
	return S[best_score_index]


##########################################################################
# Ordered Gusfield algorithm:
##########################################################################

'''
For sorted Gusfield variant:
Returns a sorted by score list of (idx, score) so the first entry is the center string,
the second entry should be the first to align with the center string,
the third should be the second to align with the center string, etc.
'''
def ordered_gusfield(S, sub_matrix, gap_cost):
    # Contains pairwise distances from s to s
	score_matrix = np.full((len(S), len(S)), None)
    # Distances from s to s is 0
	for i in range(len(S)):
		score_matrix[i, i] = 0
    # Iterate through possible centers, S[i]
	for i in range(len(S)):
        # Score for S[i]
		sum_scores = 0
        # Iterate through all other strings, S[j]
		for j in range(len(S)):
            # If we have NOT already computed the distance from S[i] to S[j], do this
			if(score_matrix[i, j] == None):
				score = pa.fill_table(S[i], S[j], sub_matrix, gap_cost)[len(S[i]), len(S[j])]
                # Distance from S[i] to S[j] is equal to the distance from S[j] to S[i]
				score_matrix[i, j] = score
				score_matrix[j, i] = score
			sum_scores += score_matrix[i, j]
    # Calculate total scores for S[i]'s
	total_scores = [sum(score) for score in score_matrix]
    # The min score index
	best_score_index = np.argmin(total_scores) # index of center string
	best_row = score_matrix[best_score_index,:] # distances to center string

	index_and_score = sorted(enumerate(best_row), key = lambda i: i[1])
	return index_and_score

def MSA_approx_ordered_gusfield(S, index_and_score, sub_matrix, gap_cost):
	M = []
	center_idx = index_and_score[0][0]
	center = S[center_idx]
	for i in range(1, len(index_and_score)):
		idx = index_and_score[i][0]
		s = S[idx]
		A_matrix = pa.fill_table(center, s, sub_matrix, gap_cost)
        # optimal alignment
		A = pa.construct_alignment(A_matrix, center, s, sub_matrix, gap_cost)
		if M:
			M = msa_common.extend_M(M, A, row_idx = 0) # center string is row 0
		else:
			M = A
	return M

##########################################################################
# Code to run
##########################################################################
'''
# Run from command line:
# python msa_approx.py sub_m.txt 5 brca.fasta
# Remenber to edit Storm's script, msa_sp_score_3k.py, to use the same sub_matrix and gap_cost

# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
S = fp.read_fasta_file(sys.argv[3])

# Get letters specified in substitution matrix file
letters = fp.parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all((c in letters for c in s) for s in S)):

	# ***** For normal gusfield only *****
	# Find center string and create MSA
	#center = find_center_string(S, sub_matrix, gap_cost)
	#seq_list = MSA_approx(S, center, sub_matrix, gap_cost)

	# ***** For ordered gusfield only *****
	index_and_score = ordered_gusfield(S, sub_matrix, gap_cost)
	seq_list = 	MSA_approx_ordered_gusfield(S, index_and_score, sub_matrix, gap_cost)

	# Print the SP score of the MSA, using Storm's script
	# Note that the command line specified sub_matrix and gap_cost should correspond to those in Storm's script
	fp.write_to_fasta_file("alignment.fasta", seq_list)
	#print(msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")
'''






