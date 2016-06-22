import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sys import argv

def main():
	print 
	print "MR Plinker (v2.1) - Collider Bias"
	print

	fname = argv[1]
	num_sims = int(argv[2])
	print "Running plink on %d simulations for %s:" % (num_sims, fname)
	print "=================================================="
	print
	runTests(fname, num_sims)


# plinkSim - once mr_predictor has been run, runs a speficic plink analysis on a given simulation, and returns an array of the results.
# @inputs: string testname which is the name used in the Mr_predictor output files, int simNum, which is the number associated with the simulation,
#          and bools cov_SBP and T2D_adj for whether you want plink to use SBP as a covariate or exclude T2D cases, respectively
# @outputs: an array which is the assoc.linear output from plink.
def plinkSim(testname, sim, cov_SBP, T2D_adj):
	#runs Plink on the .map, .pheno, and .ped files saved with testname.  Returns an array of SNPs with Beta values read from the test.assoc file;

	sim = "_" + str(sim)
	cmd = "plink --noweb --ped "+testname+sim+".ped --map "+testname+".map --pheno "+testname+sim+".pheno --pheno-name BMI --linear "
	skip = 1;
	outname = testname

	if(cov_SBP):
		cmd = cmd + "--covar " + testname+sim+ ".pheno --covar-name SBP "
		outname = outname + "_cov_SBP"
		skip = 0

	if(T2D_adj == -1):
		subprocess.call(["awk '{if ($9=='1') print $1,$2;}' " + testname + sim + ".pheno > healthy.txt"], shell=True)
		cmd = cmd + "--remove healthy.txt "
		outname = outname + "_T2D_only"
	
	if(T2D_adj == 1):
		subprocess.call(["awk '{if ($9=='2') print $1,$2;}' " + testname + sim + ".pheno > T2D_cases.txt"], shell=True)
		cmd = cmd + "--remove T2D_cases.txt "
		outname = outname + "_exT2D"

	outname = outname + sim

	FNULL = open(os.devnull, 'w')
	subprocess.call([cmd+ "--out " + outname], stdout=FNULL, stderr=subprocess.STDOUT, shell=True)

	subprocess.call(["rm " + outname + ".log"], shell = True)

	outname = outname + ".assoc.linear"

	if(cov_SBP):
		subprocess.call(["awk '{if ($5==\"ADD\") print}' " + outname + " >" + outname+ ".temp"], shell=True)
		subprocess.call(["rm " + outname], shell = True)
		outname += ".temp"

	if(T2D_adj == -1):
		subprocess.call(["rm healthy.txt"], shell = True)

	if(T2D_adj == 1):
		subprocess.call(["rm T2D_cases.txt"], shell = True)

	result = readArray(outname, skip)

	subprocess.call(["rm " + outname], shell = True)

	print result
	print

	return result

# plinkAllSims - once mr_predictor has been run, runs a speficic plink analysis on all simulations, and returns an array of the results.
# @inputs: string testname which is the name used in the Mr_predictor output files, int num_sims, which is the number of simulations to run Plink for,
#          and bools cov_SBP and T2D_adj for whether you want plink to use SBP as a covariate or exclude T2D cases, respectively
# @outputs: array of the beta values predicted by plink, organized by SNP (cols) and simulation (rows)
def plinkAllSims(testname, num_sims, cov_SBP, T2D_adj):
	#initialize the array with the results from the first sim and a columnn of rsids.
	printTestType(num_sims, cov_SBP, T2D_adj)
	result = plinkSim(testname, 1, cov_SBP, T2D_adj)
	
	i = 2;
	while (i < num_sims+1):
		sim_Result = plinkSim(testname, i, cov_SBP, T2D_adj)
		result = np.concatenate((result, sim_Result[:,1:]), axis = 1)

		i+=1

	return result

# runTests - once mr_predictor has been run, runs four different plink analyses on all simulations:
#			1. Standard linear association with no covariates or exclusions
#			2. Lin assoc. with SBP as covariate
#			3. Lin assoc. with T2D cases excluded
#			4. Lin assoc. with SBP as covariate and T2D cases excluded
# @inputs: string testname which is the name used in the Mr_predictor output files, int num_sims, which is the number of simulations mr_predictor generated
# @outputs: void. Writes an array to file for each SNP which contains the Beta values predicted in each simulation for each of the four analyses.
def runTests(testname, num_sims):

	norm = plinkAllSims(testname, num_sims, False, 0)
	cSBP = plinkAllSims(testname, num_sims, True, 0)
	xT2D = plinkAllSims(testname, num_sims, False, 1)
	cSBP_xT2D = plinkAllSims(testname, num_sims, True, 1)
	T2D_only = plinkAllSims(testname, num_sims, False, -1)
	cSBP_T2D_only = plinkAllSims(testname, num_sims, True, -1)

	num_SNPs = norm.shape[0]

	for i in range(0, num_SNPs):
		SNP = cSBP[i, 0]
		result = norm[i:i+1, :]
		result = np.concatenate ((result, cSBP[i:i+1,:]), axis = 0)
		result = np.concatenate ((result, xT2D[i:i+1,:]), axis = 0)
		result = np.concatenate ((result, cSBP_xT2D[i:i+1,:]), axis = 0)
		result = np.concatenate ((result, T2D_only[i:i+1,:]), axis = 0)
		result = np.concatenate ((result, cSBP_T2D_only[i:i+1,:]), axis = 0)
		result = np.delete(result, 0, axis = 1)	
		
		print "Saving results for " + SNP + " to " + testname + "_" + SNP + "_phenoBMI"+ ".csv"
		print
		rowNames = ["norm", "cSBP", "xT2D", "cSBP_xT2D", "T2D_adj", "cSBP_T2D_adj"]
		df = pd.DataFrame(result, index = rowNames)
		df.to_csv(testname + "_" + SNP + "_phenoBMI.csv" , index = True, sep = " ")
		
		print "Saving boxplot for " + SNP + " to " + testname + "_" + SNP + "_phenoBMI" +  ".png"
		print
		savePlot(testname + "_" + SNP + "_phenoBMI", SNP)

	return

# savePlot - reads in the simulation betas from a .csv and makes a boxplot of it
# @inputs: fname, the name of the .csv file (without the file extension)
# @outputs: none.  Saves a plot with filename fname under file extension ".png"
def savePlot(fname,SNP):
		arr = np.genfromtxt(fname + ".csv", dtype = float, skip_header = 1)
		arr = np.delete(arr, 0, axis = 1)
		arr = arr.T
		
		plt.clf()
		
		plt.title(SNP)
		plt.xlabel('Method')
		plt.ylabel('Association with BMI')

		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		bp = ax.boxplot(arr, patch_artist = True)

		for box in bp['boxes']:
			box.set(color = '#7570b3', linewidth = 2)
			box.set(facecolor = '#1b9e77')

		for whisker in bp['whiskers']:
			whisker.set(color = '#7570b3')

		for cap in bp['caps']:
			cap.set(color = '#7570b3')

		ax.set_xticklabels(['Unadjusted', 'SBP Covariate', 'Exclude T2D', 'CovSBP & ExT2D', 'T2D only', 'CovSBP & T2D only'])

		fig.savefig((fname + ".png"))

def printTestType(num_sims, cov_SBP, T2D_adj):
	border = "=================================================="
	print 
	print border
	print
	print "Running plink for " + str(num_sims) + " simulations."
	if (cov_SBP):
		print "--Using SBP as a covariate."
	if (T2D_adj == 1):
		print "--Excluding T2D cases."
	if (T2D_adj == -1):
		print "--T2D cases only."
	print 
	print border
	print

def readArray(fname, skip_header):
	arr = []
	with open(fname) as f:
		if skip_header:
			next(f)
		for line in f:
			row = line.split()
			arr.append([row[1], float(row[6])])

	return np.asarray(arr)

if __name__ == '__main__':
	main()