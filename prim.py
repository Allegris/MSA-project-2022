import sys
import numpy as np
import project2_linear as pa #pairwise alignment / SP score



# Find the node in non_MST with smallest SP score to a node in MST
def min_SP_score_index(keys, MST, non_MST):
	min_score, min_index = sys.maxsize, -1

	# Run through all nodes in non_MST and pick the one with the smallest SP score to a node in MST
	for idx in non_MST:
		if keys[idx] < min_score and idx not in MST and idx in non_MST:
			min_score = keys[idx]
			min_index = idx
	return min_index


# PRIM's algorithm for computing a minimum spanning tree for a set of strings, nodes
# (nodes are the string indices, the actual strings can be accessed as node_strings[i] for i in nodes)
def MST_prim(nodes, use_center_string):
	# Construct score matrix
	score_matrix = construct_score_matrix(nodes)
	# If we want to use the center string as the start node, do that
	# Else just choose the first node in the list nodes as the start node
	if use_center_string:
		start_node = find_center_string_index(score_matrix)
	else:
		start_node = 0
	# Contains all the nodes currently in MST
	MST = []
	# Contains all the nodes currently NOT in MST
	non_MST = nodes.copy()
	# Contains a key for each node (minimum SP score to node in MST, ie. minimum weight on edge to MST)
	keys = [sys.maxsize] * len(nodes)
	# Set key for start_node to 0, so we pick this as first node
	keys[start_node] = 0

	# While we still have nodes in non_MST, we add the minimizing node to MST
	while len(non_MST) > 0:
		# Get the min node
		min_node = min_SP_score_index(keys, MST, non_MST)
		# Add it to MST and remove from non_MST
		MST.append(min_node)
		non_MST.remove(min_node) # min_node should be present in non_MST, as we have just picked it from there
		# Update keys, now that a new node is in MST
		for i in range(len(nodes)):
			# For each node in non_MST
			if i not in MST and i in non_MST:
				# Get the SP score of min_node and node i
				score = score_matrix[min_node][i]
				# If this score is smaller than the current key, update the key
				if score > 0 and score < keys[i]:
					keys[i] = score
	return MST


# Construct a score matrix for the strings in list S
# Uses pairwise alignment scores (SP scores)
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
                str_A = node_strings[S[i]]
                str_B = node_strings[S[j]]
                score = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_A, str_B)[len(str_A), len(str_B)]
                # Distance from S[i] to S[j] is equal to the distance from S[j] to S[i]
                score_matrix[i, j] = score
                score_matrix[j, i] = score
            sum_scores += score_matrix[i, j]
    return score_matrix

# Finds the center string of multiple sequences.
# The center string is the sequences with lowest SP score with the other sequences.
# Input is the score matrix for the strings
def find_center_string_index(score_matrix):
	 total_scores = [sum(score) for score in score_matrix]
	 best_score_index = np.argmin(total_scores)
	 return best_score_index



##### Code to run #####

global submatrix, gap_cost
#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

global node_strings
node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))
print(MST_prim(nodes, True))




