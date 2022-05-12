import sys
import prim
import project2_linear as pa # Pairwise alignment / SP score
from project2_linear import calculate_alignment_matrix
from project2_linear import backtrack_nonrec
from prim import MST_prim
#import msa_common as msa
from msa_common import extend_M
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp # Helper functions for reading/writing/parsing fasta and phylip files


##########################################################################
# Compute MSA using MST
##########################################################################


'''
Computes a multiple sequence alignment (MSA) of the input nodes by:

- Computing a minimum spanning tree (MST) of the input nodes, using Prim's algorithm

- Uses this MST as a "guide tree" when doing an approximation algorithm much like Gusfield's 2-approximation algorithm,
  but where the matrix M is extended with the pairwise alignments in the order that they are added to the MST

- Sorts the MSA st. the row order corresponds to the order in node_strings

Returns the MSA
'''
def MST_MSA_approx(seq_indices, seqs, sub_matrix, gap_cost):
	M = []
	# Keep track of seq index vs row index in M
	seq_idx_to_row = [None] * len(seq_indices)
	# Iterate MST edges
	MST = MST_prim(seq_indices, seqs, sub_matrix, gap_cost)
	for i in range(len(MST)):
		# Edge (pi[u], u)
		parent = MST[i][0]
		node = MST[i][1]
		# Fill out dyn. prog. table for pairwise alignment
		table = calculate_alignment_matrix(sub_matrix, gap_cost, seqs[parent], seqs[node])
		# Construct opt. align. from table
		A = backtrack_nonrec(table, seqs[parent], seqs[node], sub_matrix, gap_cost)
		if i == 0:
			seq_idx_to_row[parent] = 0
			seq_idx_to_row[node] = 1
			M = A
		else:
			seq_idx_to_row[node] = len(M)
			M = extend_M(M, A, seq_idx_to_row[parent])
	return M
	# Used for testing: does not change score of M
	# Sort M, s.t. the row order correspond
	# to the order of input seqs
	'''
	sorted_M = [None] * len(M)
	for i in range(len(str_idx_to_row)):
		row_idx = str_idx_to_row[i]
		sorted_M[i] = M[row_idx]
	return sorted_M
	'''



# Function for testing
# Input is two rows from alignment
# Removes all gap columns from alignment and returns the SP-score of the resulting alignment
def induced_pair_score(seq_1, seq_2):
	res_1, res_2 = "", ""
	for i in range(len(seq_1)):
		if seq_1[i] == "-" and seq_2[i] == "-":
			continue
		else:
			res_1 += seq_1[i]
			res_2 += seq_2[i]
	fp.write_to_fasta_file("alignment.fasta", [res_1, res_2])
	# Compute the score of the induced alignment
	pair_score = sp_score_msa.compute_sp_score("alignment.fasta")
	return pair_score


##########################################################################
# Code to run
##########################################################################
'''
# Run from command line:
# python mst_msa_approx.py sub_m.txt 5 brca.fasta
# Remenber to edit Storm's script, msa_sp_score_3k.py, to use the same sub_matrix and gap_cost

# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
node_strings = fp.read_fasta_file(sys.argv[3])

# Assign indices to the strings in node_strings
nodes = list(range(len(node_strings)))
#print("Sequences:", node_strings, "\nSequence indices:", nodes)

# Get letters specified in substitution matrix file
letters = fp.parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all((c in letters for c in s) for s in node_strings)):
    # Construct the MSA
    seqs = MST_MSA_approx(nodes, node_strings, sub_matrix, gap_cost)
	# Write MSA to file "alignment.fasta"
    fp.write_to_fasta_file("alignment.fasta", seqs)
	# Print the SP score of the MSA, using Storm's script
	# Note that the command line specified sub_matrix and gap_cost should correspond to those in Storm's script
    print(sp_score_msa.compute_sp_score("alignment.fasta"))



    # CORRECTNESS TEST: Are scores the same for induced pair alignments and pair alignments of pairs in MST?
    tests_true = True
    mst = prim.MST_prim(nodes, node_strings, sub_matrix, gap_cost)
    for pair in mst:
        str_1 = node_strings[pair[0]]
        str_2 = node_strings[pair[1]]
        pair_score = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_1, str_2)[len(str_1), len(str_2)]
        induced_score = induced_pair_score(seqs[pair[0]], seqs[pair[1]])
        if(induced_score == pair_score):
            #print("Scores ARE the same for", pair[0], "and", pair[1])
            tests_true = True
        else:
            tests_true = False
            #print("*Scores are NOT the same for", pair[0], "and", pair[1], "....")
            #print("Score OPT:", pair_score)
            #print("Score induced:", induced_score)
    print("***DID ALL TESTS SUCCEED:", tests_true)


else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")
'''

'''
# Example to run from this file (not a very good example, rather run from commandline)

#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))

seqs = MSA_approx(nodes, node_strings, sub_matrix, gap_cost)
#print(seqs)

print_alignment_to_file(seqs)
print(sp_score_msa.compute_sp_score("alignment.fasta"))
'''










