import prim
import project2_linear as pa #pairwise alignment / SP score


#Fills out the M matrix with alignments found from backtracking
def multiple_align(S, center):
    M = []
    S.remove(center)
    for s in S:
        A_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, center, s)
        # optimal alignment
        A = pa.backtrack_nonrec(A_matrix, center, s, sub_matrix, gap_cost, "", "", len(center), len(s))
        if(s != S[0]):
            M = extend_M_with_A(M, A)
        else:
            M = A
    return M

#Extend a matrix with two or more alignments with one extra aligment.
def extend_M_with_A(M, A):
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
print(prim.MST_prim(nodes, node_strings, sub_matrix, gap_cost, True))











