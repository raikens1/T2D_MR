\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Causal and Bias estimation by Egger Regression 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Rachael Aikens and Ben Voight
Questions and Comments to Ben Voight at bvoight@upenn.edu
-----------------------------------------------

PIPELINE:

Breifly, in order to replicate the analyses shown in the paper, once should:

1) Generate n=1000 MR_predictor simulations of 150,000 individuals using the input files provided in the Input_files subdirectory

	- using the "A" type phenofile and idfile will generate simulations based on a model in which BMI and SBP both increase T2D risk
	- using the "B" type phenofile and idfile will generate simulations based on a model in which only BMI increases T2D risk

2) Run mr_plinker on all simulations.
	
	- optionally, run mr_plinker_collider_bias on the "A" model simulations to replicate our illustration of collider bias.  

3) Run mr_plinker_analysis to generate box plots and mean error data files across all simulations.

 -----------------------------------------------

 This folder includes:
 - README.txt
 - mr_plinker
 - mr_plinker_collider_bias.py
 - mr_plinker_analysis.py
 - Input_files (directory)

 -----------------------------------------------

mr_plinker.py: python script designed to run PLINK on a set of MR_predictor files to estimate SBP association and generate an output file of the results.

	mr_plinker will run the PLINK linear association tool on each dataset in a set of simulations in using four different adjustments:
		1) No adjustments
		2) Use BMI as a covariate
		3) Exclude type 2 diabetics
		4) Use BMI as a covariate and exclude type 2 diabetics

	USAGE: "mr_plinker_v2.1.py TEST_NAME NUMBER_SIMULATIONS"

	INPUT FILES: mr_plinker works on output files directly from MR_predictor.  They should have the naming schemes "testname_simNumber.pheno", "testname_simNumber.ped," and "testname.pheno."  This should follow simply from the MR_predictor output; for a test called "mytest", simply set "mytest" as the output file for MR_predictor, and MR_predictor will number simulation output for you.

	OUTPUT: mr_plinker will print an output .csv file for each SNP in the set of 13.  Each column of this file contains the simulation number followed by the SBP association for that SNP predicted by PLINK using the each of the four options above (in units of SD increase in SBP). mr_plinker will also generate a histogram of the results for each SNP.
mr_plinker_collider_bias.py: python script designed to run PLINK on a set of MR_predictor files to estimate BMI association and generate an output file of the results

mr_plinker_collider_bias.py will run the PLINK linear association tool on each dataset in a set of simulations in using six different adjustments:
		1) No adjustments
		2) Use SBP as a covariate
		3) Exclude type 2 diabetics
		4) Use SBP as a covariate and exclude type 2 diabetics
		5) Include type 2 diabetics only
		6) Use SBP as a covariate and include type 2 diabetics only

	Using the MR_predictor input files provided (Model A), this analysis will generate an example of collider bias.  None of the simulated SNPs are BMI related, yet conditioning on SBP (which is affected by both BMI and genetics), will cause SBP-associated SNPs to falsely associate with BMI.

	USAGE: mr_plinker_collider_bias.py TEST_NAME NUMBER_SIMULATIONS"

	INPUT FILES: mr_plinker_collider_bias works on output files directly from MR_predictor.  They should have the naming schemes "testname_simNumber.pheno", "testname_simNumber.ped," and "testname.pheno."  This should follow simply from the MR_predictor output; for a test called "mytest", simply set "mytest" as the output file for MR_predictor, and MR_predictor will number simulation output for you.

	OUTPUT: mr_plinker_collider_bias will print an output .csv file for each SNP in the set of 13.  Each column of this file contains the simulation number followed by the SBP association for that SNP predicted by PLINK using the each of the four options above (in units of SD increase in SBP).mr_plinker_collider_bias will also generate a histogram of the results for each SNP

mr_plinker_analysis.py: python script designed to calculate mean error over all SNPs for a set of MR_plinker output files.  Calculates mean error for each analysis type and saves txt files and a boxplot of the results.

	USAGE: "mr_plinker_v2.1.py TEST_NAME NUMBER_SIMULATIONS"

	INPUT: mr_plinker_analysis.py is designed to work directly from MR_plinker output files for the same testname

	OUTPUT: mr_plinker_analysis.py will generate a single boxplot the distribution of mean error for the four analysis types.  It will additionally save one text file for each analysis type which contains the mean errors over all simulations.

Input_files(directory): this directory contains MR_predictor input files used to run the simualtions which appear in this report.  An additional file in this directory, Input_File_Info.txt, contains additional information on these files.

-----------------------------------------------
