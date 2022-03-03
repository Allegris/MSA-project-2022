import numpy as np
from Bio import SeqIO
import sys
import msa_sp_score_3k as sp_score_msa # Storm's script


#Calculate cost of an optimal alignment for string str_A and str_B with substitution matrix sm and gap cost gc
def calculate_alignment_matrix(sub_m, gap_cost, strA, strB, strC):
    # Global vars
    global str_A
    global str_B
    global str_C
    global sm
    global gc
    global T

    # Set global vars
    str_A = strA
    str_B = strB
    str_C = strC
    sm = sub_m
    gc = gap_cost

    T = np.full((len(str_A) + 1, len(str_B) + 1, len(str_C) + 1), None)
    #iterate through rows
    for i in range(0, len(str_A) + 1):
        #iterate through columns
        for j in range(0, len(str_B) + 1):
            for k in range(0, len(str_C) + 1):
                T[i, j, k] = calc_cost_nonrec(i, j, k)
    #print(T[len(str_A),len(str_B),len(str_C)])
    return T

#Non reccursive calculation of cost
def calc_cost_nonrec(i, j, k):
    if(T[i,j,k] is None):
        v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")
        #get diagonal value
        if(i > 0 and j > 0 and k > 0):
            v1 = T[i-1, j-1, k-1] + sm[str_A[i-1]][str_B[j-1]] + sm[str_A[i-1]][str_C[k-1]] + sm[str_B[j-1]][str_C[k-1]]
        #get above value
        if(i > 0 and j > 0 and k>=0):
            v2 = T[i-1, j-1, k] + sm[str_A[i-1]][str_B[j-1]] + 2*gc
        #get left value
        if(i > 0 and j >= 0 and k>0):
            v3 = T[i-1, j, k-1] + sm[str_A[i-1]][str_C[k-1]] + 2*gc
        if(i >= 0 and j > 0 and k>0):
            v4 = T[i, j-1, k-1] + sm[str_B[j-1]][str_C[k-1]] + 2*gc
        if(i > 0 and j >= 0 and k >= 0):
            v5 = T[i-1, j, k] + 2*gc
        if(i >= 0 and j > 0 and k >= 0):
            v6 = T[i, j-1, k] + 2*gc
        if(i >= 0 and j >= 0 and k > 0):
            v7 = T[i, j,k-1] + 2*gc
        #Left top corner
        if(i == 0 and j == 0 and k ==0):
            v0 = 0
        min_val =  min(v0, v1, v2, v3, v4, v5, v6, v7)
        return min_val
    else:
        return T[i,j,k]


#Non recursive backtracking
def backtrack_nonrec(T, str_A, str_B, str_C):
    res_str_A = ""
    res_str_B = ""
    res_str_C = ""
    i = len(str_A)
    j = len(str_B)
    k = len(str_C)
    while(i >= 0 and j >= 0 and k>= 0):
        cell = T[i, j, k]
        #diagonal cell - substitution
        if (i > 0 and j > 0 and k > 0 and cell == T[i-1, j-1, k-1] + sm[str_A[i-1]][str_B[j-1]] + sm[str_A[i-1]][str_C[k-1]] + sm[str_B[j-1]][str_C[k-1]]):
            res_str_A += str_A[i-1]
            res_str_B += str_B[j-1]
            res_str_C += str_C[k-1]
            i -= 1
            j -= 1
            k -= 1
        #upper cell - insertion
        elif (i > 0 and j > 0 and k>= 0 and cell == T[i-1, j-1, k] + sm[str_A[i-1]][str_B[j-1]] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += str_B[j-1]
            res_str_C += "-"
            i -= 1
            j -= 1

        elif (i > 0 and j >= 0 and k> 0 and cell == T[i-1, j, k-1] + sm[str_A[i-1]][str_C[k-1]] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += "-"
            res_str_C += str_C[k-1]
            i -= 1
            k -= 1
        elif (i >= 0 and j > 0 and k> 0 and cell == T[i, j-1, k-1] + sm[str_B[j-1]][str_C[k-1]] + 2*gc):
            res_str_A += "-"
            res_str_B += str_B[j-1]
            res_str_C += str_C[k-1]
            j -= 1
            k -= 1

        elif (i > 0 and j >= 0 and k>= 0 and cell == T[i-1, j, k] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += "-"
            res_str_C += "-"
            i -= 1
        elif (i >= 0 and j > 0 and k>= 0 and cell == T[i, j-1, k] + 2*gc):
            res_str_A += "-"
            res_str_B += str_B[j-1]
            res_str_C += "-"
            j -= 1
        elif (i >= 0 and j >= 0 and k> 0 and cell == T[i, j, k-1] + 2*gc):
            res_str_A += "-"
            res_str_B += "-"
            res_str_C += str_C[k-1]
            k -= 1

        elif (i==0 and j==0 and k==0):
            return [res_str_A[::-1], res_str_B[::-1], res_str_C[::-1]]


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
Writes the strings/sequences in seq_list to the file "alignment.fasta" with headers seq1, seq2 etc.
'''
def write_to_fasta_file(seq_list):
    x = open("alignment.fasta", "w")
    for i in range(len(seq_list)):
        x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
    x.close()


'''
Parses a phylip file by:

Reading a file of this format:
4
A  10  2  5  2
C  2  10  2  5
G  5  2  10  2
T  2  5  2  10

representing a substitution matrix and returns a dictionary corresponding to the substitutionmatrix

Returns a dictionary of this format (if getAlphabet = False):
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



# Code we run from command line
# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
str_A = read_fasta_file(sys.argv[3])[0]
str_B = read_fasta_file(sys.argv[3])[1]
str_C = read_fasta_file(sys.argv[3])[2]

# Get letters specified in substitution matrix file
letters = parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all(c in letters for c in str_A) and all(c in letters for c in str_B) and all(c in letters for c in str_C)):
    # Calculate alignment matrix and print optimal cost
    t = calculate_alignment_matrix(sub_matrix, gap_cost, str_A, str_B, str_C)
    # If we want to backtrack, write optimal alignment in file alignment.fasta
    if len(sys.argv)==5 and sys.argv[4]=="True":
        b = backtrack_nonrec(t, str_A, str_B, str_C)
        write_to_fasta_file(b)
        print(sp_score_msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")


