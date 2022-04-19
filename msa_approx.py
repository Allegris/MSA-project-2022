import sys
import numpy as np
import time
import seaborn as sns
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
#Fills out the M matrix with alignments found from backtracking
def MSA_approx(S, center, sub_matrix, gap_cost):
    M = []
    S.remove(center)
    for s in S:
        A_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, center, s)
        # optimal alignment
        A = pa.backtrack_nonrec(A_matrix, center, s, sub_matrix, gap_cost)
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
				score = pa.calculate_alignment_matrix(sub_matrix, gap_cost, S[i], S[j])[len(S[i]), len(S[j])]
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

	# Find center string and create MSA
	center = find_center_string(S, sub_matrix, gap_cost)
	seq_list = MSA_approx(S, center, sub_matrix, gap_cost)

	# Print the SP score of the MSA, using Storm's script
	# Note that the command line specified sub_matrix and gap_cost should correspond to those in Storm's script
	fp.write_to_fasta_file("alignment.fasta", seq_list)
	#print(msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")

'''

##########################################################################
# Measure running time
##########################################################################


'''
lst_time=[]
lst_length=[]
#print(len(str_A))
for i in range(1,len(str_A)//10):
    s=10*i
    start = time.time()
    #t3, p3 = calculate_alignment_matrix(sub_matrix, gap_cost, str_A[:s], str_B[:s])
    center = find_center_string(S)
    seq_list = MSA_approx(S, center)
    end = time.time()
    lst_time.append((end-start)/s**2)
    lst_length.append(s)
print(lst_time)

ax = sns.scatterplot(x = lst_length, y = lst_time)
ax.set(xlabel = "lenght of seq", ylabel = "Time (sec/n^2)")
figure = ax.get_figure()
figure.savefig("time_of_alg_linear_n2.png")
'''







