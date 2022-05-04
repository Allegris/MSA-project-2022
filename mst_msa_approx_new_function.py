def MST_MSA_approx(seq_indices, seqs, sub_matrix, gap_cost):
    MST = prim.MST_prim(seq_indices, seqs, sub_matrix, gap_cost)
    M = []
	# Contains, for seq with index i,
	# the row index in M that corresponds to this seq
	# This is later used for sorting the rows in M,
	# s.t. the seq with index 0 is first, then index 1, etc.
    seq_idx_to_row = [None] * len(seq_indices)

    for i in range(len(MST)):
		# The pair of seqs to align in this iteration
        pair = MST[i]
        seq_1 = seqs[pair[0]]
        seq_2 = seqs[pair[1]]
		# Create (optimal) pairwise alignment matrix for the pair of seqs
        pair_align_matrix = pa.calculate_alignment_matrix
							(sub_matrix, gap_cost, seq_1, seq_2) # Alignment A
		# Get score of the optimal alignment
        pair_opt_align = pa.backtrack_nonrec
						(pair_align_matrix, seq_1, seq_2, sub_matrix, gap_cost)
        if i == 0:
			# The first pair are the first two rows in M
            M = pair_opt_align
            seq_idx_to_row[pair[0]] = 0
            seq_idx_to_row[pair[1]] = 1
        else:
			# We now add seq with index pair[1] to the M matrix
			# as row with index len(M) - record this info
            seq_idx_to_row[pair[1]] = len(M)
			# Row index in M of the upper seq in pairwise alignment
            row_idx = seq_idx_to_row[pair[0]]
			# Extend matrix M with the pairwise alignment
            M = msa_common.extend_M(M, pair_opt_align, row_idx)
	# Sort M, s.t. the row order correspond
	# to the order in input seqs
    sorted_M = [None] * len(M)
    for i in range(len(seq_idx_to_row)):
        row_idx = seq_idx_to_row[i]
        sorted_M[i] = M[row_idx]
    return sorted_M
