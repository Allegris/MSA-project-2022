import sys
import numpy as np
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp # Helper functions for reading/writing/parsing fasta and phylip files


##########################################################################
# Compute an optimal MSA for 3 sequences (exact method)
##########################################################################


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





# Code we run from command line
# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
str_A = fp.read_fasta_file(sys.argv[3])[0]
str_B = fp.read_fasta_file(sys.argv[3])[1]
str_C = fp.read_fasta_file(sys.argv[3])[2]

# Get letters specified in substitution matrix file
letters = fp.parse_phylip(sys.argv[1], True)

# Check if sequences only contain allowed letters
if(all(c in letters for c in str_A) and all(c in letters for c in str_B) and all(c in letters for c in str_C)):
    # Calculate alignment matrix and print optimal cost
    t = calculate_alignment_matrix(sub_matrix, gap_cost, str_A, str_B, str_C)
    # If we want to backtrack, write optimal alignment in file alignment.fasta
    if len(sys.argv)==5 and sys.argv[4]=="True":
        b = backtrack_nonrec(t, str_A, str_B, str_C)
        fp.write_to_fasta_file(b)
        print(sp_score_msa.compute_sp_score("alignment.fasta"))
else:
    print("Error: A letter in a sequence is not specified in the substitution matrix.")


