import numpy as np
from Bio import SeqIO

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


