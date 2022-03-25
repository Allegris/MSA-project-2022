from numpy.random import choice

alpha = ["A", "C", "G", "T"]
probs = [0.25] * len(alpha)



def common_ancestor(n):
	return "".join(choice(alpha, n, probs))

def descendants(anc, alpha, sub_rate, indel_rate, k):
	res = []
	for i in range(k):
		res.append(descendant(anc, alpha, sub_rate, indel_rate))
	return res

def descendant(anc, alpha, sub_rate, indel_rate):
	res = ""
	for char in anc:
		#print("****CHAR", char)
		# Insert symbol before char?
		insert = binary_choice(indel_rate)
		if insert:
			res += random_symbol(alpha)
			#print("INS before", char, ":", res[-1])
		# Delete char?
		delete = binary_choice(indel_rate)
		if delete:
			#print("DEL", char)
			continue
		else:
			# Substitute char?
			sub = binary_choice(sub_rate)
			if sub:
				res += random_symbol(alpha)
				#print("SUB", char, "for", res[-1])
			else:
				#print("Adding char", char)
				res += char
	return res


def random_symbol(alpha):
	return "".join(choice(alpha, 1))

def binary_choice(rate):
	return choice([True, False], 1, [rate, 1 - rate])

#Writes a fasta file with the aligned sequences
def print_seqs_to_file(path, seq_list, idx):
	title = path + str(idx) + ".fasta"
	x = open(title, "w")
	for i in range(len(seq_list)):
		x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
	x.close()



##########################################################################
# Code to run
##########################################################################


n = 100
k = 3
anc = common_ancestor(n)
sub_rate = 0.75 # 0.75 corresponds to random strings
indel_rate = sub_rate / 5
path = "simulated_data\\test_seqs_" # start of path to put sequences in
idx = 1 # used in file name

'''
# Small test: will a division of the sequences in two clusters influence results
desc1 = descendants(anc, alpha, sub_rate, indel_rate, k)
desc2 = descendants(anc, alpha, sub_rate*5, indel_rate*5, k)
print_seqs_to_file(path, desc1 + desc2, idx)
'''

for i in range(10):
	d = descendants(anc, alpha, sub_rate, indel_rate, k)
	print_seqs_to_file(path, seq_list = d, idx = i)



'''
print("anc:", anc)
for i in range(len(desc)):
	print("desc:", i, ":", desc[i])
'''

