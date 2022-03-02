import prim
import project2_linear as pa #pairwise alignment / SP score


# Fills out the M matrix with alignments found from backtracking
def MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string):
    MST_pairs_to_align = prim.MST_prim(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
    M = []
	# Contains, for string index i, the row index in M that corresponds to this string
	# This is later used for sorting the rows in M s.t. the string with index 0 is first, then index 1, etc.
    str_idx_to_row = [None] * len(nodes)

    for i in range(len(MST_pairs_to_align)):
		# The pair of strings to align in this iteration
        pair = MST_pairs_to_align[i] # TODO: maybe refactor this var (delete)
        str_1 = node_strings[pair[0]]
        str_2 = node_strings[pair[1]]
		# Make alignment matrix for the pair of strings
        pair_align_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_1, str_2) # Score of alignment A
		# Create an optimal alignment, by backtracking the matrix
        pair_opt_align = pa.backtrack_nonrec(pair_align_matrix, str_1, str_2, sub_matrix, gap_cost, "", "", len(str_1), len(str_2)) # TODO: refactor the last 2 params in other file
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
    return M


# Extend the matrix M with a new optimal alignment, pair_align
def extend_M(M, pair_align, pair_idx, str_idx_to_row):
	# Contains the new M string
	new_M_str = ""
	# Let i point to a column in M
	# (just pick M[0] as range, all M strings should be of same length)
	# Let j point to a column in pair_opt_align
	# (just pick pair_opt_align[0] as range, both strings should be of same length)
	i = j = 0
	while i < len(M[0]) and j < len(pair_align[0]):  # TODO: Should this be OR??? What to do if we run out of string in one of them?
		# Now we have four cases of how the two columns look.
		# This is the symbol in M that we are interested in (the row corresponding to string with index pair[0])
		M_symbol = M[str_idx_to_row[pair_idx[0]]][i]
		# Upper and lower symbol in the pairwise alignment column
		upper_symbol = pair_align[0][j]
		lower_symbol = pair_align[1][j]
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
			M = [s + "-" for s in M]
			new_M_str += lower_symbol
			i += 1 # TODO: Should this be here?? (se above comment)
			j += 1
		# Case 4: M_symbol is not "-" and upper_symbol is not "-"
		# TODO: This case is the same as case 1: we can combine them to test if M_symbol == upper_symbol
		elif M_symbol != "-" and upper_symbol != "-":
			new_M_str = lower_symbol
			i += 1
			j += 1
	# If we have run through all of pair_align, but now M, then add the corresponding number of "-"'s to new_M_str
	diff_len = len(M[0]) - len(new_M_str)
	new_M_str += "-" * diff_len
	# new_M_str is done - append it to M and return M
	M.append(new_M_str)
	return M




##### Code to run #####

#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))

# Should we use the center string as start node in MST
use_center_string = True

print(MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string))










