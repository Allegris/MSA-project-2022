import sys
import os
import time
import matplotlib.pyplot as plt
import mst_msa_approx as mst_algo
import msa_approx as gusfield
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp


# Returns: (n, k), score, time
def evaluate_MSA_algo(folder, sub_matrix_filename, gap_cost, n):

	d = sorted(os.listdir(folder))

	k_list = []
	score_diffs = []
	time_diffs = []

	g_score_list = []
	m_score_list = []

	g_time_list = []
	m_time_list = []

	for filename in d:
		print(filename)
		S = fp.read_fasta_file(folder + "\\" + filename)
		k = len(S) # Number of seqs
		k_list.append(k)
		sub_matrix = fp.parse_phylip(sub_matrix_filename)
		letters = fp.parse_phylip(sub_matrix_filename, True) # Get letters specified in substitution matrix file
		# Check if sequences only contain allowed letters
		if(all((c in letters for c in s) for s in S)):
			##### MSA #####
			S_idx = list(range(len(S)))

			# Run algo 10 times and use average running time
			m_times = []
			for i in range(10):
				m_start = time.time() # Start timer
				seqs = mst_algo.MST_MSA_approx(S_idx, S, sub_matrix, gap_cost)
				m_end = time.time() # Stop timer
				m_times.append(m_end - m_start)

			# Average running time
			m_time = sum(m_times)/(len(m_times))
			m_time_list.append(m_time)

			# Score
			fp.write_to_fasta_file("alignment_" + str(k) +".fasta", seqs)
			m_score = sp_score_msa.compute_sp_score("alignment_" + str(k) +".fasta")
			m_score_list.append(m_score)


			##### Gusfield #####
			g_times = []
			for i in range(10):
				g_start = time.time() # Start timer
				center = gusfield.find_center_string(S, sub_matrix, gap_cost)
				seqs = gusfield.MSA_approx(S, center, sub_matrix, gap_cost)
				g_end = time.time() # Stop timer
				g_times.append(g_end - g_start)

			# Average running time
			g_time = sum(g_times)/(len(g_times))
			g_time_list.append(g_time)

			# Score
			fp.write_to_fasta_file("alignment.fasta", seqs)
			g_score = sp_score_msa.compute_sp_score("alignment.fasta")
			g_score_list.append(g_score)


			# Differences between MST and Gusfield
			'''
			diff_score = ((g_score - m_score)/g_score)*100
			score_diffs.append(diff_score)
			diff_time = ((g_time - m_time)/g_time)*100
			time_diffs.append(diff_time)
			'''
		else:
			   print("Error: A letter in a sequence is not specified in the substitution matrix.")


	# Time plot
	print("Time for n = ", n)
	plt.scatter(k_list, g_time_list, color = "red") # Gusfield times
	plt.scatter(k_list, m_time_list, color="blue") # MST times
	#plt.xticks(range(1,20))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec)", fontsize = 13)
	plt.savefig('res_time_random_n100_kmax20.png')
	plt.show()

	# Clear plot
	plt.clf()

	# Score plot
	print("Score for n = ", n)
	plt.scatter(k_list, g_score_list, color = "red") # Gusfield scores
	plt.scatter(k_list, m_score_list, color="blue") # MST scores
	#plt.xticks(range(1,20))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("SP-score", fontsize = 13)
	plt.savefig('res_score_random_n100_kmax20.png')
	plt.show()

	return score_diffs, time_diffs #(n, k), diff_score, diff_time



##########################################################################
# Code to run
##########################################################################
'''
# Get sub matrix, gap cost, and sequences from command line variables
sub_matrix = fp.parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
S = fp.read_fasta_file(sys.argv[3])
mst = bool(sys.argv[3])
# Assign indices to the strings in node_strings
S_idx = list(range(len(S)))
'''
folder = "simulated_data"
sub_matrix_filename = "sub_m.txt"
gap_cost = 5
n = 1000

print(evaluate_MSA_algo(folder, sub_matrix_filename, gap_cost, n))

#mst_res = evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, True)
#gusfield_res = evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, False)

#print(evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, mst))









