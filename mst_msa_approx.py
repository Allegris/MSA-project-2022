import sys
import prim
import project2_linear as pa # Pairwise alignment / SP score
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp # Helper functions for reading/writing/parsing fasta and phylip files


##########################################################################
# Compute MSA using MST
##########################################################################


'''
Computes a multiple sequence alignment (MSA) of the input nodes by:

- Computing a minimum spanning tree (MST) of the input nodes, using Prim's algorithm
  (use_center_string is whether the "center string" should be used as start node in Prim's algorithm)

- Uses this MST as a "guide tree" when doing an approximation algorithm much like Gusfield's 2-approximation algorithm,
  but where the matrix M is extended with the pairwise alignments in the order that they are added to the MST

- Sorts the MSA st. the row order corresponds to the order in node_strings

Returns the MSA
'''
def MST_MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string):
    MST_pairs_to_align = prim.MST_prim(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
    M = []
	# Contains, for string index i, the row index in M that corresponds to this string
	# This is later used for sorting the rows in M s.t. the string with index 0 is first, then index 1, etc.
    str_idx_to_row = [None] * len(nodes)

    for i in range(len(MST_pairs_to_align)):
		# The pair of strings to align in this iteration
        pair = MST_pairs_to_align[i]
        str_1 = node_strings[pair[0]]
        str_2 = node_strings[pair[1]]
		# Make alignment matrix for the pair of strings
        pair_align_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_1, str_2) # Score of alignment A
		# Create an optimal alignment, by backtracking the matrix
        pair_opt_align = pa.backtrack_nonrec(pair_align_matrix, str_1, str_2, sub_matrix, gap_cost)
        if i == 0:
            M = pair_opt_align
			# The first pair are the first two rows in M
            str_idx_to_row[pair[0]] = 0
            str_idx_to_row[pair[1]] = 1
        else:
			# We now add string with index pair[1] to the M matrix as row index len(M) - record this info
            str_idx_to_row[pair[1]] = len(M)
			# Extend matrix M with the new pairwise optimal alignment
            M = extend_M(M, pair_opt_align, pair, str_idx_to_row)
	# Sort M so that the row order correspond to the order in node_strings
    sorted_M = [None] * len(M)
    for i in range(len(str_idx_to_row)):
        row_idx = str_idx_to_row[i]
        sorted_M[i] = M[row_idx]
    return sorted_M


'''
Extends the current M matrix with a new optimal pairwise alignment

Returns the new M matrix (which contains one more row than the input M)
'''
# Extend the matrix M with a new optimal alignment, pair_align
def extend_M(M, pair_align, pair_idx, str_idx_to_row):
	# Contains the new M string
	new_M_str = ""
	# Let i point to a column in M
	# (just pick M[0] as range, all M strings should be of same length)
	i = 0
	# Let j point to a column in pair_opt_align
	# (just pick pair_opt_align[0] as range, both strings should be of same length)
	j = 0
	while i < len(M[0]) and j < len(pair_align[0]):
		# This is the symbol in M that we are interested in (the row corresponding to string with index pair[0])
		M_symbol = M[str_idx_to_row[pair_idx[0]]][i]
		# Upper and lower symbol in the pairwise alignment column
		upper_symbol = pair_align[0][j]
		lower_symbol = pair_align[1][j]

		# Now we have four cases of how the two columns look.
		# Case 1: M_symbol is "-" and upper_symbol is "-"
		if M_symbol == "-" and upper_symbol == "-":
			# We add lower_symbol to new_M_str and increment i and j
			new_M_str += lower_symbol
			i += 1
			j += 1
		# Case 2: M_symbol is "-", but upper_symbol is not "-"
		elif M_symbol == "-" and upper_symbol != "-":
			# We add "-" to new_M_str and increment i
			new_M_str += "-"
			i += 1
		# Case 3: M_symbol is not "-", but upper_symbol is "-"
		elif M_symbol != "-" and upper_symbol == "-":
			# We create new column with lower_symbol in last row and only "-"'s above it, and increment j
			# We also increment i, since we add a new column, and we don't want to look at this column in next iteration
			M = [s[:i] + "-" + s[i:] for s in M]
			new_M_str += lower_symbol
			i += 1 # This is here, because we insert a new column to the left of col i - and we still want to look at col i and not this new column
			j += 1
		# Case 4: M_symbol is not "-" and upper_symbol is not "-"
		elif M_symbol != "-" and upper_symbol != "-" and M_symbol == upper_symbol:
			new_M_str += lower_symbol
			i += 1
			j += 1
	# If we broke while loop, because we ran out of M-string,
	# add "-"'s to M-strings, except for the upper string in pair_align and new_M_str,
	# for these, just add the rest of the pair_align
	if i >= len(M[0]):
		new_M_str += pair_align[1][j:]
		row_of_upper = str_idx_to_row[pair_idx[0]]
		M[row_of_upper] += pair_align[0][j:]
		M = [M[i] + "-" * len(pair_align[0][j:]) if i != row_of_upper else M[i] for i in range(len(M))]
	# If we broke while loop, because we ran out of pair_align string,
	# add "-"'s to new_M_str
	elif j >= len(pair_align[0]):
		# Diff in length of M strings and new_M_str
		diff_len = len(M[0]) - len(new_M_str)
		new_M_str += "-" * diff_len
	# new_M_str is done - append it to M and return M
	M.append(new_M_str)
	return M


##########################################################################
# Code to run
##########################################################################


# Run from command line:
# python mst_msa_approx.py sub_m.txt 5 brca.fasta True
# (use False if center string should not be used as start node in MST, but rather just the first string in node_strings)
# Remenber to edit Storm's script, msa_sp_score_3k.py, to use the same sub_matrix and gap_cost

# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
node_strings = fp.read_fasta_file(sys.argv[3])
use_center_string = bool(sys.argv[4]) # Should center string be used as start node in MST?
# TODO: It seems that finding the center string yields the same SP scores as when not, so maybe delete this option later...

# Assign indices to the strings in node_strings
nodes = list(range(len(node_strings)))

# Get letters specified in substitution matrix file
letters = fp.parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all((c in letters for c in s) for s in node_strings)):
    # Construct the MSA
    seqs = MST_MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
	# Write MSA to file "alignment.fasta"
    fp.write_to_fasta_file(seqs)
	# Print the SP score of the MSA, using Storm's script
	# Note that the command line specified sub_matrix and gap_cost should correspond to those in Storm's script
    print(sp_score_msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")



'''
# Example to run from this file (not a very good example, rather run from commandline)

#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))

# Should we use the center string as start node in MST
use_center_string = True

seqs = MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
#print(seqs)

print_alignment_to_file(seqs)
print(sp_score_msa.compute_sp_score("alignment.fasta"))
'''










