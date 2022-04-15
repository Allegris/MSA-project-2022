import sys
import numpy as np
import project2_linear as pa #pairwise alignment / SP score


# PRIM's algorithm for computing a minimum spanning tree for a set of strings, nodes
# (nodes are the string indices, the actual strings can be accessed as node_strings[i] for i in nodes)
def old_MST_prim(nodes, node_strings, sub_matrix, gap_cost):
	# Construct score matrix
	cost = construct_score_matrix(nodes, node_strings, sub_matrix, gap_cost)
	mst = [] # Will contain the MST
	in_MST = [0]
	not_in_MST = nodes.copy()
	not_in_MST.remove(0)
	# We want all nodes in the MST - when we have this, we are done
	while len(not_in_MST) > 0:
		min_cost = sys.maxsize # Min edge to MST
		pair = (-1, -1)
		# Check all possible edges (pair of nodes)
		for u in in_MST:
			for v in not_in_MST:
				# If the edge has lower cost than min_cost
				if cost[u][v] < min_cost:
					min_cost = cost[u][v]
					pair = (u, v)
		if pair != (-1, -1):
			# Add the found edge to MST
			mst.append(pair)
			in_MST.append(pair[1])
			not_in_MST.remove(pair[1])
	return mst


def find_min_node(V, mst, key):
	min_val = sys.maxsize
	min_node = None
	for v in V:
		if v not in mst and key[v] < min_val:
			min_node = v
			min_val = key[v]
	return min_node


def MST_prim(nodes, node_strings, sub_matrix, gap_cost):
	V = nodes.copy()
	# Construct score matrix
	cost = construct_score_matrix(nodes, node_strings, sub_matrix, gap_cost)
	mst = [] # Will contain the MST edges (i.e. the resulting MST)
	key = [sys.maxsize] * len(V)
	parent = [None] * len(V)

	key[0] = 0

	#for _ in range(len(V)):
	while len(mst) < len(V):
		u = find_min_node(V, mst, key)
		mst.append(u)

		for v in V:
			if v not in mst and cost[u][v] < key[v]:
				key[v] = cost[u][v]
				parent[v] = u

	mst = [(parent[v], v) for v in mst[1:]] #remove the first entry, since this is "None node to the root", i.e. (-1, 0)
	return mst





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

#sub_matrix = {"A": {"A": 5, "C": 5, "G": 5, "T": 5}, "C": {"A": 5, "C": 5, "G": 5, "T": 5}, "G": {"A": 5, "C": 5, "G": 5, "T": 5}, "T": {"A": 5, "C": 5, "G": 5, "T": 5}}
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, "C": {"A": 2, "C": 10, "G": 2, "T": 5}, "G": {"A": 5, "C": 2, "G": 10, "T": 2}, "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = 5

node_strings = ["AACG", "AAAA", "CCCC", "GGGG"]
nodes = list(range(len(node_strings)))
print(MST_prim(nodes, node_strings, sub_matrix, gap_cost))




