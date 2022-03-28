import sys
import time
import mst_msa_approx as mst_algo
import msa_approx as gusfield
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp


# Returns: (n, k), score, time
def evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, mst):
	letters = fp.parse_phylip(sys.argv[1], True) # Get letters specified in substitution matrix file
	# Check if sequences only contain allowed letters
	if(all((c in letters for c in s) for s in S)):
		start = time.time() # Start timer
		if mst:
			# Construct the MSA
		    seqs = mst_algo.MST_MSA_approx(S_idx, S, sub_matrix, gap_cost)
		else:
			# Find center string and create MSA
			center = gusfield.find_center_string(S, sub_matrix, gap_cost)
			seqs = gusfield.MSA_approx(S, center, sub_matrix, gap_cost)
		end = time.time() # End timer
		run_time = end - start
		fp.write_to_fasta_file("alignment.fasta", seqs)
		score = sp_score_msa.compute_sp_score("alignment.fasta")
	else:
		   print("Error: A letter in a sequence is not specified in the substitution matrix.")
	n = len(max(S, key = len)) # Longest string in S, i.e. "n"
	k = len(S) # Number of seqs
	return (n, k), score, run_time


##########################################################################
# Code to run
##########################################################################

# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
S = fp.read_fasta_file(sys.argv[3])
mst = bool(sys.argv[3])
# Assign indices to the strings in node_strings
S_idx = list(range(len(S)))

print(evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, mst))









