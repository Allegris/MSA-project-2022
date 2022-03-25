from Bio import SeqIO


##########################################################################
# Helper functions for reading and writing to files
##########################################################################


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
def write_to_fasta_file(filename, seq_list):
    x = open(filename, "w") # eg. "alignment.fasta"
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


