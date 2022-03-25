import os
from Bio import SeqIO
import subprocess

d = sorted(os.listdir("./testseqs"))


'''

for file in d:
	print(file)
	result_mst = subprocess.check_output("python mst_msa_approx.py sub_m.txt 5 ./testseqs/" + file, shell=True, text=True)
	print(result_mst.strip())
'''

# Test for diff between MSA and MST_MSA


approx_score = []
mst_approx_score = []

for file in d:
	print(file)
	result_gusfield = subprocess.check_output("python msa_approx.py sub_m.txt 5 ./testseqs/" + file + " True", shell=True, text=True)
	approx_score.append(int(result_gusfield.strip()))
	print("Gusfield: ", approx_score)
	#print(result_gusfield.strip())

	result_mst = subprocess.check_output("python mst_msa_approx.py sub_m.txt 5 ./testseqs/" + file, shell=True, text=True)
	mst_approx_score.append(int(result_mst.strip()))
	print("MST: ", mst_approx_score)
	#print(result_mst.strip())

print(approx_score)
print(mst_approx_score)

for i in range(len(approx_score)):
	print("Diff", approx_score[i] - mst_approx_score[i])
