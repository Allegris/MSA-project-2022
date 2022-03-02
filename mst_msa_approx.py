import prim
import project2_linear as pa #pairwise alignment / SP score


#Fills out the M matrix with alignments found from backtracking
def OLD_multiple_align(S, center):
    M = []
    S.remove(center)
    for s in S:
        A_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, center, s)
        # optimal alignment
        A = pa.backtrack_nonrec(A_matrix, center, s, sub_matrix, gap_cost, "", "", len(center), len(s))
        if(s != S[0]):
            M = OLD_extend_M_with_A(M, A)
        else:
            M = A
    return M

#Extend a matrix with two or more alignments with one extra aligment.
def OLD_extend_M_with_A(M, A):
    # List for holding new M-strings
    new_M = ["" for i in range(len(M))]
    # String holding the last M-string for the new M-list
    new_M_str = ""
    # Iterate through columns in A (opt alignment)
    for j in range(len(A[0])):
        # Upper and lower symbol in the j'th column
        upper_symbol = A[0][j]
        lower_symbol = A[1][j]
        # Case: insertion
        if upper_symbol == "-":
            # We want to add a column to M consisting of only "-"'s except
            # the last row, which should contain lower_symbol:
            # Add "-" to new M-strings (except for last M-string)
            new_M = [old_str + "-" for old_str in new_M]
            # Add lower_symbol to last M-string
            new_M_str += lower_symbol
        # Case: (mis)match or deletion
        else:
            # Find first occurrance of upper_symbol in first M-string
            # (we have cut off the part the string that we have already added to the new M-strings)
            pos = M[0].find(upper_symbol)
            # Add everything up until this symbol to all the new M-strings (except the last)
            prefixes_M = [m_str[:pos+1] for m_str in M]
            new_M = [new_M[i] + prefixes_M[i] for i in range(len(new_M))]
            # To the last string, add correspondingly many gaps (minus 1) and lower_symbol
            new_M_str += "-"*pos + lower_symbol
            # Update M to only contain suffix that we have not yet "transferred" to new_M yet
            M = [m_str[pos+1:] for m_str in M]
    # Check if there is any suffix of the M-strings remaining and add these to new_M (and corresponding gaps to new_M_str)
    suffix_length_M = len(M[0])
    if suffix_length_M > 0:
        new_M = [new_M[i] + M[i] for i in range(len(new_M))]
        new_M_str += "-"*suffix_length_M
    # Add the last M-string to the updated M-strings in new_M
    new_M.append(new_M_str)
    # Return entire new M
    return new_M

# Fills out the M matrix with alignments found from backtracking
def multiple_align(nodes, node_strings, sub_matrix, gap_cost, use_center_string):
    MST_pairs_to_align = prim.MST_prim(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
    M = []
	# Contains, for string index i, the row index in M that corresponds to this string
	# This is later used for sorting the rows in M s.t. the string with index 0 is first, then index 1, etc.
	str_idx_to_row = None * nodes

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
            M = extend_M(M, pair_opt_align, pair)
    return M

# Extend the matrix M with a new optimal alignment, pair_align
def extend_M(M, pair_align, pair_idx, str_idx_to_row):
	# Contains the new M string
	new_M_str = ""
	# Let i point to a column in M
	# (just pick M[0] as range, all M strings should be of same length)
	while i < len(M[0]):
		# Let j point to a column in pair_opt_align
		# (just pick pair_opt_align[0] as range, both strings should be of same length)
		while j < len(pair_opt_align[0]):
			# Now we have four cases of how the two columns look.
			# This is the symbol in M that we are interested in (the row corresponding to string with index pair[0])
			M_symbol = M[str_idx_to_row[pair[0]]][i]
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
				M = [s + "-" for s in M]
				new_M_str += lower_symbol
				j += 1
			elif M_symbol != "-" and upper_symbol != "-":



# Extend the M matrix with a new alignment, A
def OLD_extend_M(M, A):
    # List for holding new M-strings
    new_M = ["" for i in range(len(M))]
    # String holding the last M-string for the new M-list
    new_M_str = ""
    # Iterate through columns in A (opt alignment)
    for j in range(len(A[0])):
        # Upper and lower symbol in the j'th column
        upper_symbol = A[0][j]
        lower_symbol = A[1][j]
        # Case: insertion
        if upper_symbol == "-":
            # We want to add a column to M consisting of only "-"'s except
            # the last row, which should contain lower_symbol:
            # Add "-" to new M-strings (except for last M-string)
            new_M = [old_str + "-" for old_str in new_M]
            # Add lower_symbol to last M-string
            new_M_str += lower_symbol
        # Case: (mis)match or deletion
        else:
            # Find first occurrance of upper_symbol in first M-string
            # (we have cut off the part the string that we have already added to the new M-strings)
            pos = M[0].find(upper_symbol)
            # Add everything up until this symbol to all the new M-strings (except the last)
            prefixes_M = [m_str[:pos+1] for m_str in M]
            new_M = [new_M[i] + prefixes_M[i] for i in range(len(new_M))]
            # To the last string, add correspondingly many gaps (minus 1) and lower_symbol
            new_M_str += "-"*pos + lower_symbol
            # Update M to only contain suffix that we have not yet "transferred" to new_M yet
            M = [m_str[pos+1:] for m_str in M]
    # Check if there is any suffix of the M-strings remaining and add these to new_M (and corresponding gaps to new_M_str)
    suffix_length_M = len(M[0])
    if suffix_length_M > 0:
        new_M = [new_M[i] + M[i] for i in range(len(new_M))]
        new_M_str += "-"*suffix_length_M
    # Add the last M-string to the updated M-strings in new_M
    new_M.append(new_M_str)
    # Return entire new M
    return new_M




##### Code to run #####

#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))












