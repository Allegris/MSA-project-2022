from Bio import SeqIO
import sys
import prim
import project2_linear as pa #pairwise alignment / SP score
import msa_sp_score_3k as sp_score_msa # Storm's script


'''
Reads a fasta file and changes all nucleic symbols not in [A, C, G, T] to A

Returns a list containing the strings/sequences in the fasta file
'''
def read_fasta_file(filename):
    rec_list = []
    nucleic_list = ["U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]
    for record in SeqIO.parse(filename, "fasta"):
        corrected_seq = str(record.seq)
        for symbol in nucleic_list:
            corrected_seq = corrected_seq.replace(symbol, "A")
        rec_list.append(corrected_seq)
    return rec_list


'''
Reads a file of this format:
4
A  10  2  5  2
C  2  10  2  5
G  5  2  10  2
T  2  5  2  10

representing a substitution matrix and returns a dictionary corresponding to the substitutionmatrix

Returning a dictionary of this format (if getAlphabet = False):
{"A": {"A": 10, "C": 2, "G": 5, "T": 2},
 "C": {"A": 2, "C": 10, "G": 2, "T": 5},
 "G": {"A": 5, "C": 2, "G": 10, "T": 2},
 "T": {"A": 2, "C": 5, "G": 2, "T": 10}}

If getAlphabet = True, we instead return a list of the alphabet letters:
['A', 'C', 'G', 'T']
'''
def parse_phylip(filename, getAlphabet = False):
    f= open(filename, "r")
    f1 = f.readlines()
    f2 = list()
    for x in f1:
        f2.append(x.split())
    alph_size = int(f2[0][0])

    letters = list()
    for i in range(1, alph_size+1):
        letters.insert(i, f2[i][0])

    sub_matrix = dict()
    for i in range(len(letters)):
        inner_dict = dict()
        for j in range(len(letters)):
            inner_dict[letters[j]] = int(f2[i+1][j+1])
        sub_matrix[letters[i]] = inner_dict
    if(getAlphabet):
        return letters
    else:
        return sub_matrix


'''
Computes a multiple sequence alignment (MSA) of the input nodes by:

- Computing a minimum spanning tree (MST) of the input nodes, using Prim's algorithm
  (use_center_string is whether the "center string" should be used as start node in Prim's algorithm)

- Uses this MST as a "guide tree" when doing an approximation algorithm much like Gusfield's 2-approximation algorithm,
  but where the matrix M is extended with the pairwise alignments in the order that they are added to the MST

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
		# TODO: SORT M!!!
    return M


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
	# Let j point to a column in pair_opt_align
	# (just pick pair_opt_align[0] as range, both strings should be of same length)
	i = j = 0
	while i < len(M[0]) and j < len(pair_align[0]):
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
			M = [s[:i] + "-" + s[i:] for s in M]
			new_M_str += lower_symbol
			i += 1 # This is here, because we insert a new column to the left of col i - and we still want to look at col i and not this new column
			j += 1
		# Case 4: M_symbol is not "-" and upper_symbol is not "-"
		elif M_symbol != "-" and upper_symbol != "-":
			new_M_str += lower_symbol
			i += 1
			j += 1
	# If we have run through all of pair_align, but now M, then add the corresponding number of "-"'s to new_M_str
	diff_len = len(M[0]) - len(new_M_str)
	new_M_str += "-" * diff_len
	# new_M_str is done - append it to M and return M
	M.append(new_M_str)
	return M


'''
Writes the strings/sequences in seq_list to the file "alignment.fasta" with headers seq1, seq2 etc.
'''
def print_alignment_to_file(seq_list):
    x = open("alignment.fasta", "w")
    for i in range(len(seq_list)):
        x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
    x.close()


##### Code to run #####

# Run from command line:
# python mst_msa_approx.py sub_m.txt 5 brca.fasta True
# (use False if center string should not be used as start node in MST, but rather just the first string in node_strings)
# Remenber to edit Storm's script, msa_sp_score_3k.py, to use the same sub_matrix and gap_cost

# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
node_strings = read_fasta_file(sys.argv[3])
use_center_string = bool(sys.argv[4]) # Should center string be used as start node in MST?

# Assign indices to the strings in node_strings
nodes = list(range(len(node_strings)))

# Get letters specified in substitution matrix file
letters = parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all((c in letters for c in s) for s in node_strings)):
    # Construct the MSA
    seqs = MST_MSA_approx(nodes, node_strings, sub_matrix, gap_cost, use_center_string)
	# Write MSA to file "alignment.fasta"
    print_alignment_to_file(seqs)
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










