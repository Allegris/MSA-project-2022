
'''
Extends the current M matrix with a new optimal pairwise alignment, A

row_idx is the row index in M corresponding to the string in the upper row of A

Returns the new M matrix (which contains one more row than the input M)
'''
def extend_M(M, A, row_idx):
	# Contains the new M string
	new_M_str = ""
	# Let i point to a column in M
	# (just pick M[0] as range, all M strings should be of same length)
	i = 0
	# Let j point to a column in pair_align
	# (just pick pair_opt_align[0] as range, both strings should be of same length)
	j = 0
	while i < len(M[0]) and j < len(A[0]):
		# This is the symbol in M that we are interested in (the row corresponding to string with index pair[0])
		M_symbol = M[row_idx][i]
		# Upper and lower symbol in the pairwise alignment column
		upper_symbol = A[0][j]
		lower_symbol = A[1][j]

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
	# add "-"'s to M-strings, except for the upper string in pair_align and new_M_str;
	# for these, just add the rest of the pair_align
	if i >= len(M[0]):
		new_M_str += A[1][j:]
		M[row_idx] += A[0][j:]
		M = [M[i] + "-" * len(A[0][j:]) if i != row_idx else M[i] for i in range(len(M))]
	# If we broke while loop, because we ran out of pair_align string,
	# add "-"'s to new_M_str
	elif j >= len(A[0]):
		# Diff in length of M strings and new_M_str
		diff_len = len(M[0]) - len(new_M_str)
		new_M_str += "-" * diff_len
	# new_M_str is done - append it to M and return M
	M.append(new_M_str)
	return M