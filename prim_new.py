import sys
import numpy as np
import project2_linear as pa #pairwise alignment / SP score


def MST_prim

def allowed_edge(u, v, in_MST):
	# If nodes are the same, return false
	if u == v:
		return False
	# If neither node is in MST, return false
	if in_MST[u] == False and in_MST[v] == False:
		return False
	# If both nodes are in MST, return false
	elif in_MST[u] == True and in_MST[v] == True:
		return False
	# Else, return true
	return True

def MST_prim_complete(nodes, node_strings, sub_matrix, gap_cost):
	# Construct score matrix
	cost = construct_score_matrix(nodes, node_strings, sub_matrix, gap_cost)
	mst = [] # Will contain the MST
	V = len(nodes) # Number of vertices in input graph
	in_MST = [False] * V # Is a node of given index already in the MST
	in_MST[0] = True # Let node 0 be the root of the MST
	MST_edge_no = 0 # Number of edges currently in MST
	# We want V-1 edges in the MST - when we have this, we are done
	while MST_edge_no < V - 1:
		min_cost = sys.maxsize # Min edge to MST
		u, v = -1, -1
		# Check all possible edges (pair of nodes)
		for i in range(V):
			for j in range(V):
				# If the edge has lower cost than min_cost
				if cost[i][j] < min_cost:
					# Check if one node is in MST and the other is not
					if allowed_edge(i, j, in_MST):
						min_cost = cost[i][j] # Update min_cost
						u, v = i, j
		if u != -1 and v != -1:
			# Add the found edge to MST
			mst.append((u, v))
			in_MST[u], in_MST[v] = True, True
			MST_edge_no += 1 # We have added one more edge to MST
	return mst












# Find the node in non_MST with smallest SP score to a node in MST
def min_SP_score_index(keys, MST, non_MST):
	min_score, min_index = sys.maxsize, -1

	# Run through all nodes in non_MST and pick the one with the smallest SP score to a node in MST
	for idx in non_MST:
		if keys[idx][0] < min_score and idx not in MST and idx in non_MST:
			min_score = keys[idx][0]
			min_index = idx
	return min_index


# PRIM's algorithm for computing a minimum spanning tree for a set of strings, nodes
# (nodes are the string indices, the actual strings can be accessed as node_strings[i] for i in nodes)
def MST_prim(nodes, node_strings, sub_matrix, gap_cost, use_center_string):
	# Construct score matrix
	score_matrix = construct_score_matrix(nodes, node_strings, sub_matrix, gap_cost)
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
	# Contains a key for each node (minimum SP score to node in MST, ie. minimum weight on edge to MST) and the node in MST
	keys = [(sys.maxsize, None)] * len(nodes)
	# Set key for start_node to 0, so we pick this as first node
	keys[start_node] = (0, None)

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
				# If this score is smaller than the current key value, update the key value
				# and the node at the other end of the edge that we add to MST
				if score > 0 and score < keys[i][0]:
					keys[i] = (score, min_node)
	# Now we have the order that the nodes were added to the MST (in list, MST)
	# And we have the final keys and their corresponding nodes
	# Combining these, we can get the pairs in the order that we want to align them in MST MSA
	pairs_to_align = [(MST[0], MST[1])] # First align start node and the next node added to MST
	# For each other node, n, in MST (in the inserting order of MST), add the pair consisting of (m, n)
	# where the edge (m, n) were the edge that connected node n to the MST in the algorithm
	for i in range(2, len(MST)):
		new_MST_node = MST[i]
		min_connected_MST_node = keys[new_MST_node][1]
		pairs_to_align.append((min_connected_MST_node, new_MST_node))
	return pairs_to_align


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
print(MST_prim(nodes, node_strings, sub_matrix, gap_cost, True))
print(MST_prim_complete(nodes, node_strings, sub_matrix, gap_cost))




