import numpy as np
import itertools as it


def func_numpy(str_ls):
	strings_len = [len(s)+1 for s in str_ls]

	dimensions = [*range(len(str_ls))]
	print("DIMENSION: ", dimensions)

	T = np.full(strings_len, None)
	for index, values in np.ndenumerate(T):
		T[index] = new_func(list(index))

	combs = []
	for i in range(1, len(strings_len)+1):
		# Decrement i of the index values in every possible way
		for comb in it.combinations(dimensions, i):
			combs.append(comb)

	#print(combs)





def new_func(index):
	return 4


strings = ["aaa", "ac", "atgg"]
func_numpy(strings)


my_list = [0,0,1,1,0]
print(list(it.permutations(my_list, 2)))


'''
def calc_cost_nonrec_new(index):
	if(T[index] is None):
		# We want to find all adjacent cells that have already been filled out
		prev_nb_cells = []
		for i in range(len(index)):
			# Decrement i of the index values in every possible way
			for numbers_to_modify in it.combinations(index, repeat=i)
			# Do not exceed matrix bounds
			'''