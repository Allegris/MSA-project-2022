import numpy as np
from Bio import SeqIO
import project2_linear as pa #pairwise alignment
import sys 
import msa_sp_score_3k as msa #score for approx


# Read fasta files
def read_fasta_file(filename):
    rec_list = []
    nucleic_list = ["U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]
    for record in SeqIO.parse(filename, "fasta"):
        corrected_seq = str(record.seq)
        for symbol in nucleic_list:
            corrected_seq = corrected_seq.replace(symbol, "A")
        rec_list.append(corrected_seq)
    return rec_list

"""
HELPER FUNCTION
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
"""
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

#Finds the center string of multiple sequences. The center string is the sequences with lowest alignments score to the other sequences. 
def find_center_string_fast(S):
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

#Writes a fasta file with the aligned sequences 
def print_alignment_to_file(seq_list):
    x = open("alignment.fasta", "w")
    for i in range(len(seq_list)):    
        x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
    x.close()



# Code we run from command line
# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
S = read_fasta_file(sys.argv[3])


# Get letters specified in substitution matrix file
letters = parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters 
if(all((c in letters for c in s) for s in S)):
    # Calculate alignment matrix and print optimal cost and write fasta file
    center = find_center_string_fast(S)
    seq_list = multiple_align(S, center)
    print_alignment_to_file(seq_list)
    print(msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")

