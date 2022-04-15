import sys
import numpy as np
import project2_linear as pa #pairwise alignment / SP score


# Constructs a MST from a set of nodes and cost matrix (as sub_matrix and gap_cost)
# node_strings is a list of the actual strings
# nodes is a corresponding numbering of the strings: [0, 1, 2,...]
def MST_prim(nodes, node_strings, sub_matrix, gap_cost):
	V = nodes.copy() # The set of nodes
	# Construct score matrix
	cost = construct_score_matrix(nodes, node_strings, sub_matrix, gap_cost)
	mst = [] # Will contain the MST edges (i.e. the resulting MST)
	key = [sys.maxsize] * len(V) # The value of the lightest weight connecting a node to the current MST
	parent = [None] * len(V) # The parent node that a given node was connected to the MST using

	key[0] = 0 # Set the key of node 0 to 0, ensuring that we pick this as the first node

	# While not all nodes are in MST
	while len(mst) < len(V):
		# Find the node that can be added to the MST giving
		# the lowest possible weight addition
		u = find_min_node(V, mst, key)
		mst.append(u)

		# Update keys and parents of nodes not in MST, since MST has changed
		for v in V:
			if v not in mst and cost[u][v] < key[v]:
				key[v] = cost[u][v]
				parent[v] = u
	# Return MST as edge pairs (u, v)
	# The edges appear in the order that they are added to the MST
	mst = [(parent[v], v) for v in mst[1:]] # Remove the first entry, since this is "None node to the root", i.e. (None, 0)
	return mst


# Returns the node in the set *V - mst* with the lowest key value
def find_min_node(V, mst, key):
	min_val = sys.maxsize
	min_node = None
	for v in V:
		if v not in mst and key[v] < min_val:
			min_node = v
			min_val = key[v]
	return min_node


# Construct a score matrix for the string indices in list S (the actual strings are in node_strings)
# Uses pairwise alignment scores (SP scores)
def construct_score_matrix(S, node_strings, sub_matrix, gap_cost):
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


##### Code to run #####
'''
#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))
print(MST_prim(nodes, node_strings, sub_matrix, gap_cost))
'''



