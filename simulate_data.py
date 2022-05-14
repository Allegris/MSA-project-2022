from numpy.random import choice
import numpy

alpha = ["A", "C", "G", "T"]
probs = [0.25] * len(alpha)



def common_ancestor(n):
	return "".join(choice(alpha, n, p=probs))

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
	return choice([True, False], 1, p=[rate, 1 - rate])

#Writes a fasta file with the aligned sequences
def print_seqs_to_file(seq_list, n, k):
	path = "simulated_data\\n" + str(n) + "\\"
	title = path + "test_seq_" + "n" + str(n) + "_k" + str(k) + ".fasta"
	x = open(title, "w")
	for i in range(len(seq_list)):
		x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
	x.close()



##########################################################################
# Code to run
##########################################################################


ns = [10, 50, 100, 150]
sub_rate = 0.05 # 1 corresponds to random strings
indel_rate = sub_rate / 5

'''
# Small test: will a division of the sequences in two clusters influence results
desc1 = descendants(anc, alpha, sub_rate, indel_rate, k)
desc2 = descendants(anc, alpha, sub_rate*5, indel_rate*5, k)
print_seqs_to_file(path, desc1 + desc2, idx)
'''

'''
OLD CODE
anc1 = common_ancestor(n)
anc2 = common_ancestor(n)
#anc2 = descendant(anc1, alpha, sub_rate = 1, indel_rate = 0)
#anc2 = common_ancestor(n)
#anc2 = descendant(anc1, alpha, sub_rate = 1, indel_rate = 0.25) #"A" * n
for i in range(1, 11):
	d1 = descendants(anc1, alpha, sub_rate, indel_rate, i)
	d2 = descendants(anc2, alpha, sub_rate, indel_rate, i)
	print_seqs_to_file(path, d1 + d2 + [anc2], n, k = i*2)
'''


for n in ns:
	anc1 = common_ancestor(n)
	anc2 = common_ancestor(n)
	for k in range(1, 11):
		d1 = descendants(anc1, alpha, sub_rate, indel_rate, k)
		d2 = descendants(anc2, alpha, sub_rate, indel_rate, k)
		print_seqs_to_file(d1 + d2, n, k*2)
'''


for n in ns:
	anc = common_ancestor(n)
	for k in range(2, 21):
		d = descendants(anc, alpha, sub_rate, indel_rate, k)
		print_seqs_to_file(d, n, k)



'''






