import sys
import os
import time
import numpy as np
import prim
import matplotlib.pyplot as plt
import mst_msa_approx as mst_algo
import msa_approx as gusfield
import msa_sp_score_3k as sp_score_msa # Storm's script
import fasta_and_phylip as fp


# Returns: (n, k), score, time
def evaluate_MSA_algo(folder, sub_matrix_filename, gap_cost, n):

	d = sorted(os.listdir(folder))

	k_list = []

	g_score_list = []
	m_score_list = []
	og_score_list = []

	p_time_list = [] # prims algo only
	g_time_list = []
	m_time_list = []
	og_time_list = []

	p_exptime_list = [] # prims algo only
	g_exptime_list = []
	m_exptime_list = []
	og_exptime_list = []

	for filename in d:
		print(filename)
		S = fp.read_fasta_file(folder + "\\" + filename)
		k = len(S) # Number of seqs
		k_list.append(k)
		sub_matrix = fp.parse_phylip(sub_matrix_filename)
		letters = fp.parse_phylip(sub_matrix_filename, True) # Get letters specified in substitution matrix file
		# Check if sequences only contain allowed letters
		if(all((c in letters for c in s) for s in S)):

			S_idx = list(range(len(S)))


			##### PRIM ALGO ONLY #####
			p_times = []
			for i in range(5):
				p_start = time.time() # Start timer
				prim.MST_prim(S_idx, S, sub_matrix, gap_cost)
				p_end = time.time() # Stop timer
				p_times.append(p_end - p_start)
			# Average running time
			p_time = sum(p_times)/(len(p_times))
			# Min time, not currently used
			#p_time = min(p_times)
			p_time_list.append(p_time)

			# Expected runningtime: time / k^2
			p_exptime_list.append(p_time/(k**2))


			'''
			##### MSA #####
			# Run algo 5 times and use average running time
			m_times = []
			for i in range(5):
				m_start = time.time() # Start timer
				seqs = mst_algo.MST_MSA_approx(S_idx, S, sub_matrix, gap_cost)
				m_end = time.time() # Stop timer
				m_times.append(m_end - m_start)

			# Average running time
			m_time = sum(m_times)/(len(m_times))
			m_time_list.append(m_time)

			# Expected runningtime: time / k^2
			m_exptime_list.append(m_time/((k**2)*(n**2)))

			# Score
			fp.write_to_fasta_file("alignment_" + str(k) +".fasta", seqs)
			m_score = sp_score_msa.compute_sp_score("alignment_" + str(k) +".fasta")
			m_score_list.append(m_score)
			'''

			'''
			#CORRECTNESS TEST: Are scores the same for induced pair alignments and pair alignments of pairs in MST?
			tests_true = True
			mst = prim.MST_prim(S_idx, S, sub_matrix, gap_cost)
			for pair in mst:
			    str_1 = S[pair[0]]
			    str_2 = S[pair[1]]
			    pair_score = pa.calculate_alignment_matrix(sub_matrix, gap_cost, str_1, str_2)[len(str_1), len(str_2)]
			    induced_score = induced_pair_score(seqs[pair[0]], seqs[pair[1]])
			    if(induced_score == pair_score):
			        #print("Scores ARE the same for", pair[0], "and", pair[1])
			        tests_true = True
			    else:
			        tests_true = False
			        #print("*Scores are NOT the same for", pair[0], "and", pair[1], "....")
			        #print("Score OPT:", pair_score)
			        #print("Score induced:", induced_score)
			print("***DID ALL TESTS SUCCEED:", tests_true)
			'''

			'''
			##### Ordered Gusfield #####
			og_times = []
			for _ in range(10):
				og_start = time.time() # Start timer
				index_and_score = gusfield.ordered_gusfield(S.copy(), sub_matrix, gap_cost)
				seqs = 	gusfield.MSA_approx_ordered_gusfield(S.copy(), index_and_score, sub_matrix, gap_cost)
				og_end = time.time() # Stop timer
				og_times.append(og_end - og_start)

			# Average running time
			og_time = sum(og_times)/(len(og_times))
			og_time_list.append(og_time)

			# Expected runningtime: time / k^2
			og_exptime_list.append(og_time/((k**2)*(n**2)))

			# Score
			fp.write_to_fasta_file("alignment.fasta", seqs)
			og_score = sp_score_msa.compute_sp_score("alignment.fasta")
			og_score_list.append(og_score)
			'''

			'''
			##### Gusfield #####
			g_times = []
			for _ in range(10):
				g_start = time.time() # Start timer
				center = gusfield.find_center_string(S.copy(), sub_matrix, gap_cost)
				seqs = gusfield.MSA_approx(S.copy(), center, sub_matrix, gap_cost)
				g_end = time.time() # Stop timer
				g_times.append(g_end - g_start)

			# Average running time
			g_time = sum(g_times)/(len(g_times))
			g_time_list.append(g_time)

			# Expected runningtime: time / k^2
			g_exptime_list.append(g_time/((k**2)*(n**2)))

			# Score
			fp.write_to_fasta_file("alignment.fasta", seqs)
			g_score = sp_score_msa.compute_sp_score("alignment.fasta")
			g_score_list.append(g_score)
			'''

		else:
			   print("Error: A letter in a sequence is not specified in the substitution matrix.")


	# Prim time plot
	print("Prim time for n = ", n)
	plt.scatter(k_list, p_time_list, color = "green") # Prim times
	plt.xticks(range(5,41,5)) # plt.xticks(range(1, 21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec)", fontsize = 13)
	plt.savefig('res_prim_time_random_n' + str(n) + '_kmax40.png')
	plt.show()
	plt.clf() # Clear plot

	# Prim expected time plot
	print("Prim exp time for n = ", n)
	plt.scatter(k_list, p_exptime_list, color = "green") # Prim times
	plt.xticks(range(5,41,5)) # plt.xticks(range(1, 21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec) / k^2", fontsize = 13)
	plt.savefig('res_prim_exptime_random_n' + str(n) + '_kmax40.png')
	plt.show()
	plt.clf() # Clear plot

	'''
	# Time plot
	print("Time for n = ", n)
	plt.scatter(k_list, g_time_list, color = "red", alpha = 0.5) # Gusfield times
	plt.scatter(k_list, m_time_list, color="blue", alpha = 0.5) # MST times
	plt.xticks(range(1,21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec)", fontsize = 13)
	plt.savefig('res_time_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()
	plt.clf() # Clear plot

	# Expected time plot
	print("Time/(k^2*n^2) for n = ", n)
	plt.scatter(k_list, g_exptime_list, color = "red", alpha = 0.5) # Gusfield times
	plt.scatter(k_list, m_exptime_list, color="blue", alpha = 0.5) # MST times
	plt.xticks(range(1,21))
	#plt.ylim(-0.005, 0.005)
	plt.ylim(0, 10*(10**(-6)))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec) / k^2*n^2", fontsize = 13)
	plt.savefig('res_exptime_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()
	plt.clf() # Clear plot

	# Score plot
	print("Score for n = ", n)
	plt.scatter(k_list, g_score_list, color = "red", alpha = 0.5) # Gusfield scores
	plt.scatter(k_list, m_score_list, color="blue", alpha = 0.5) # MST scores
	plt.xticks(range(1,21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("SP-score", fontsize = 13)
	plt.savefig('res_score_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()


	# Time plot
	print("2Time for n = ", n)
	plt.scatter(k_list, g_time_list, color = "red", alpha = 0.5) # Gusfield times
	plt.scatter(k_list, m_time_list, color="blue", alpha = 0.5) # MST times
	plt.scatter(k_list, og_time_list, color="green", alpha = 0.5) # ordered Gusfield times
	plt.xticks(range(1,21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec)", fontsize = 13)
	plt.savefig('res2_time_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()
	plt.clf() # Clear plot

	# Expected time plot
	print("2Time/(k^2*n^2) for n = ", n)
	plt.scatter(k_list, g_exptime_list, color = "red", alpha = 0.5) # Gusfield times
	plt.scatter(k_list, m_exptime_list, color="blue", alpha = 0.5) # MST times
	plt.scatter(k_list, og_exptime_list, color="green", alpha = 0.5) # ordered Gusfield times
	plt.xticks(range(1,21))
	#plt.ylim(-0.005, 0.005)
	plt.ylim(0, 10*(10**(-6)))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("Time (sec) / k^2*n^2", fontsize = 13)
	plt.savefig('res2_exptime_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()
	plt.clf() # Clear plot

	# Score plot
	print("2Score for n = ", n)
	plt.scatter(k_list, g_score_list, color = "red", alpha = 0.5) # Gusfield scores
	plt.scatter(k_list, m_score_list, color="blue", alpha = 0.5) # MST scores
	plt.scatter(k_list, og_score_list, color="green", alpha = 0.5) # ordered Gusfield times
	plt.xticks(range(1,21))
	plt.title("n = " + str(n))
	plt.xlabel("Number of sequences, k", fontsize = 13)
	plt.ylabel("SP-score", fontsize = 13)
	plt.savefig('res2_score_low2anc_n' + str(n) + '_kmax20.png')
	plt.show()

	# Differences between MST and Gusfield
	# List of S_m / S_g
	m_g_frac_score = [m_score_list[i]/g_score_list[i] for i in range(len(g_score_list))  if g_score_list[i] != 0.0]
	# List of T_m / T_g
	m_g_frac_time = [m_time_list[i]/g_time_list[i] for i in range(len(g_time_list)) if g_time_list[i] != 0.0]

	score_frac_avg = np.average(m_g_frac_score)
	score_frac_var = np.var(m_g_frac_score)
	time_frac_avg = np.average(m_g_frac_time)
	time_frac_var = np.var(m_g_frac_time)

	x = open('res_score_low2anc_n' + str(n) + '_kmax20.txt', "w")
	x.write("S_m/S_g average: " + str(score_frac_avg) + "\nS_m/S_g variance: " + str(score_frac_var) +
		 "\nT_m/T_g average: " + str(time_frac_avg) + "\nT_m/T_g variance: " + str(time_frac_var))
	x.close()


	return score_frac_avg, score_frac_var, time_frac_avg, time_frac_var
	'''
	return "DONE"

# Function for testing
# Input is two rows from alignment
# Removes all gap columns from alignment and returns the SP-score of the resulting alignment
def induced_pair_score(seq_1, seq_2):
	res_1, res_2 = "", ""
	for i in range(len(seq_1)):
		if seq_1[i] == "-" and seq_2[i] == "-":
			continue
		else:
			res_1 += seq_1[i]
			res_2 += seq_2[i]
	fp.write_to_fasta_file("alignment.fasta", [res_1, res_2])
	# Compute the score of the induced alignment
	pair_score = sp_score_msa.compute_sp_score("alignment.fasta")
	return pair_score



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


#folder = "simulated_data"
sub_matrix_filename = "sub_m.txt"
gap_cost = 5
#n = 10

ns = [10, 50, 100, 150]#[150, 200, 250, 300]

for n in ns:
	print("**************N:", n)
	folder = "simulated_data\\n" + str(n)
	print(evaluate_MSA_algo(folder, sub_matrix_filename, gap_cost, n))




#mst_res = evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, True)
#gusfield_res = evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, False)

#print(evaluate_MSA_algo(S_idx, S, sub_matrix, gap_cost, mst))









