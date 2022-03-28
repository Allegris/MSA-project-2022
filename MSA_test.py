import os
import subprocess


def time_msa():

	folder = "simulated_data"
	d = sorted(os.listdir(folder))

	time_diffs = []
	score_diffs = []

	for file in d:
		print(file)
		res_mst = subprocess.check_output("python msa_eval.py sub_m.txt 5 " + folder + "\\" + file + " True", shell=True, text=True)
		res_gusfield = subprocess.check_output("python msa_eval.py sub_m.txt 5 " + folder + "\\" + file + " False", shell=True, text=True)

		print("MST:", res_mst)
		print("Gus:", res_gusfield)

		# ((n, k), diff)
		#score_diffs.append(tuple(res_mst[0]), int(res_mst[1]) - int(res_gusfield[1]))
		#time_diffs.append(tuple(res_mst[0]), int(res_mst[2]) - int(res_gusfield[2]))
	return score_diffs, time_diffs



##########################################################################
# Code to run
##########################################################################

# Run from command line:
# python msa_approx.py sub_m.txt 5 brca.fasta
# Remenber to edit Storm's script, msa_sp_score_3k.py, to use the same sub_matrix and gap_cost

time_msa()

'''
diffs = time_msa()

print("Score diffs:", diffs[0])
print("Time diffs:", diffs[1])
'''







