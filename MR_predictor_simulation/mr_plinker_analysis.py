import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sys import argv

# Main - calculates error in plink association statistics from a set of mr_plinker output files
# USAGE: 'mr_plinker_analysis.py BETAS_FILE TEST_NAME NUM_SIMS'
def main():
	types = {1:'unadjusted', 2:'conditioning on BMI', 3:'exluding diabetes cases', 4:'conditioning on BMI and exluding diabetes cases'}
	print 
	print "Error Analysis"
	print
	print "=================================================="
	print
	betas = getParams(argv[1])
	num_sims = int(argv[3])

	for i in range(4):
		print "Analyzing results for %d simulation analyses %s" % (num_sims, types[1])
		oneTest(betas, argv[2], num_sims, i+1)
		print "-------------------------------------------------"
		print


# getParams - reads in a mr_predictor scorefile and returns a list of (SNP,beta) pairs 
# @inputs: fname, the name of the scorefile
# @outputs: betas, a list of (SNP,beta) pairs
def getParams(fname):
	print "Reading in beta values... "
	betas = []
	with open(fname) as f:
		for line in f:
			row = line.split()
			betas.append((row[0],float(row[2])))
	print "Complete."
	return betas


# oneTest - calculates and plots error for a single test type
# @inputs: betas, a list of (SNP,beta) pairs; tname, the root name for test files; sim_num, the number of simulations run;
#  		   and ttype, an int (1 for unadjusted, 2 for covBMI, 3 for xT2D, 4 for covBMI_xT2D)
# @outputs: none, calls errHist to save a histogram
def oneTest(betas, tname, num_sims, ttype):
	#add a file extension to show what corrections were made in finding linear associations
	extensions = {1:'_unadj', 2:'_cBMI', 3:'_xT2D', 4:'_cBMI_xT2D'}

	data = getData(betas, tname, num_sims, ttype)
	outname = tname + extensions[ttype]
	errs = getErr(betas,data)
	print
	print "Saving results for to " + outname + ".txt"
	print
	np.savetxt(outname + ".txt", errs)
	errHist(errs, outname, ttype)


# getData - given a dictionary of SNP:beta values, a test name, and a test type, reads in the appropriate txt files 
#           and constructs an array with SNPs as rows and simulation results as columns
# @inputs: betas, a list of (SNP,beta) pairs; tname, the root name for test files; sim_num, the number of simulations run;
#  		   and ttype, an int (1 for unadjusted, 2 for covBMI, 3 for xT2D, 4 for covBMI_xT2D)
# @outputs: data, an array with SNPs as rows and simulation results as columns
def getData(betas, tname, num_sims, ttype):
	print "Reading in simulation results... "

	data = np.zeros((len(betas), num_sims))

	#read in data from file
	for i in range(len(betas)):
		with open(tname+"_"+betas[i][0]+".csv") as f:
			myfile = f.readlines()
		line = myfile[ttype].split()
		line.pop(0)
		data[i,:] = map(float, line)

	print "Complete."
	return data


# errHist - makes and saves a histogram of a 1D numpy array
# @inputs: arr, a 1D array of errors for many simulations; fname, the name the histogram should be saved under
# @outputs: none.  Saves a histogram of the array as a .png
def errHist(arr, fname, ttype):
	plt.clf()
	types = {1:'Unadjusted', 2:'Conditioning on BMI', 3:'Exluding diabetes cases', 4:'Conditioning on BMI and exluding diabetes cases'}

	plt.title(types[ttype])
	plt.xlabel('Mean Error')
	plt.ylabel('Count')

	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	bp = ax.hist(arr)

	print 'Saving a histogram of these errors to %s_hist.png.' % fname

	fig.savefig((fname + "_hist.png"))
	

# getErr - given an array, calculates meanErr by column and returns a 1D np array of error
# @inputs: data, an array with SNPs as rows and simulation results as columns, and betas, a list of (SNP,beta) pairs.
# @outputs: a 1D np array of errors by column
def getErr(betas, data):
	print "Calculating error over all SNPs for each simulation... ",
	errors = []
	params = np.asarray([item[1] for item in betas])
	for i in range(data.shape[1]):
		errors.append(meanErr(params,data[:,i]))
	print "Complete."
	return np.asarray(errors)


# meanErr - given a parameter value and a np array of simulated estimates, returns the mean error
# @inputs: param, a float; stats, a 1D np array of estimates
# @outputs: err, a float
def meanErr(params, stats):
	n = len(stats)

	#correcting some stats for sign
	swaps = [8,9,11]
	for i in range(stats.size):
		if i in swaps:
			stats[i] = -stats[i]
			
	return np.sum(stats-params)/n


# savePlot - reads in the simulation betas from a .txt and makes a boxplot of it
# @inputs: fname, the name of the .txt file (without the file extension)
# @outputs: none.  Saves a plot with filename fname under file extension ".png"
def savePlot(fname,SNP):
	arr = np.genfromtxt(fname + ".txt", dtype = float, skip_header = 1)
	arr = np.delete(arr, 0, axis = 1)
	arr = arr.T
	
	plt.clf()
	
	plt.title(SNP)
	plt.xlabel('Method')
	plt.ylabel('Association with SBP')

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

	ax.set_xticklabels(['Unadjusted', 'BMI Covariate', 'Exclude T2D', 'CovBMI & ExT2D'])

	fig.savefig((fname + ".png"))

if __name__ == '__main__':
	main()