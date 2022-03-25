import sys
import project2_linear as pa # Pairwise alignment / SP score
import fasta_and_phylip as fp

str_1 = "TTTAACTGGGTTCGAAAAATAGCAATTCCTCGTGGAAGGCGGGCTACCAAGCGTCGCTGATTCTATCCGATGTCTGCTAGCGCCCGCAGTC"
str_2 = "TCGAGGACAATGGTGGGTCACATGCAATGGCTACCGTATCGTGTTTCCACCATGCGTACGTCGGTAGCTGGCCAGTTTCGGCTTCCCTATGCTATCGGC"
str_3 = "TTTGATCCAACTCTGTGTCAGCGAATGACAGCTACGGGCATCCTACTCGGTTGAGTCCAATTTCATGGGTTCCTTCTGCCGCAGGGCTC"

sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])

def score(str_1, str_2):
	A = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_1, str_2) # Alignment A
	return A[len(str_1), len(str_2)]

print(score(str_1, str_2))
print(score(str_1, str_3))
print(score(str_2, str_3))
