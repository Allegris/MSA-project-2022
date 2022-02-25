import sys
import numpy as np
import project2_linear as pa #pairwise alignment




# Find the string with minimum SP-score with a string in MST
def min_SP_score_index(keys, MST, non_MST):
		min_score, min_index = sys.maxsize, -1

		for idx in range(len(non_MST)):
			if keys[idx] < min_score and MST[idx] == False:
				min_score = keys[idx]
				min_index = idx
		return min_index


def printMST(MST, score_matrix):
	print ("Edge \tWeight")
	for i in range(1, len(MST)):
		print (MST[i], "-", i, "\t", score_matrix[i][MST[i]])

def MST_prim(nodes):

	# Construct score matrix and center string
	score_matrix = construct_score_matrix(nodes)
	start_node = find_center_string_index(score_matrix)

	# Contains all the nodes currently in MST
	MST = [start_node]
	# Contains all the nodes currently NOT in MST
	non_MST = nodes
	# Remove start node from non_MST
	if start_node in non_MST:
		non_MST.remove(start_node)

	# Contains a key for each node (minimum weight on edge to MST)
	keys = [sys.maxsize] * len(nodes)

	# Set key for start_node to 0, so we pick this as first node
	keys[start_node] = 0

	# While we still have nodes in non_MST, we add the minimizing node to MST
	while len(non_MST) > 0:
		# Get the min node
		min_node = min_SP_score_index(keys, MST, non_MST)
		# Add it to MST and remove from non_MST
		MST.append(min_node)
		non_MST.remove(min_node)
		# Update keys, now that a new node is in MST
		for i in range(len(nodes)):
			if i not in MST:
				score = score_matrix[min_node][i]
				if score > 0 and score < keys[i]:
					keys[i] = score

	printMST(MST, score_matrix)


#Finds the center string of multiple sequences. The center string is the sequences with lowest alignments score to the other sequences.
def construct_score_matrix(S):
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
    return score_matrix


def find_center_string_index(score_matrix):
	 total_scores = [sum(score) for score in score_matrix]
	 best_score_index = np.argmin(total_scores)
	 return best_score_index




global submatrix, gap_cost
#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

nodes = ["AACG", "AAAA", "CCCC", "GGGG"]
MST_prim(nodes)





# May not need this:
'''

print("Score matrix:\n", m)
print("Center string index: ", s1)
print("Center string: ", nodes[s1])

str_idx = {}

for i in range(len(nodes)):
	str_idx[i] = nodes[i]

#print("str_idx: ", str_idx)

'''




