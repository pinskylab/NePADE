This repository contains data, scripts and figures for a project investigating how genetic diversity and effective population size have changed over time for summer flounder

## Data

### Main analyses

• **PADE_stock_assessment16.xlsx**: number of fish at age, estimated fishing mortality at age, natural mortality at age & proportion mature at age for 1982-2015

• **Ne280_1196loci_nomaf_MAFpop0.obs**: msfs for 2008 cohort

• **Ne280_1196loci_nomaf_MAFpop1.obs**: msfs for 1997 cohort

• **Ne280_1196loci_nomaf_MAFpop2.obs**: msfs for 1994 cohort

• **SNP.DP3g95nomaf.FIL.FIL.recode.140trimmed.280fish.firstsnp.genepop.gen**: genepop file containing full datset of 3821 SNPs for 280 larvae

• **Ne_PADE_1196loci_complete.gen**: genepop file containing subset of full dataset with no missing data, 1196 SNPs across 280 larvae, used for demographic modeling

### Sensitivity analyses: no minor allele count (--mac) filter

• **Ne_PADE_1084loci_complete.gen**: genepop file containing subset of full dataset with no missing data, 1084 SNPs across 284 larvae, used for sensitivity analysis of demographic modeling

• **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen**: genepop file containing full dataset of 3905 SNPs for 284 larvae

## Scripts

• **best_lhoods.R** script for determining maximum-likelihood of each demographic scenario, calculating AIC and plotting demographic modeling results

• **demo_modeling_allele_counts.R**: script to determine the allele counts of all alleles used in demographic modeling (main analyses & no --mac filter sensitivity analyses)

• **generation_length.R** script calculates generation length for females and males over time using the 2016 stock assessment data, age-length relationships from Penttila et al. (1989) and age-fecundity curves from Morse (1981)

• **PADE_stock_assessment16.R** script imports data from 2016 summer flounder stock assessment, calculates the number of potential breeders at each age and year and plots census and breeding abundance over time

• **Ne.R** script for figuring out which summer flounder to resequence

• **Ne_analyses.R** exploratory script for plotting PCAs and calculating diversity statistics

• **Ne_genotyping_results.R** script for calculating coverage statistics

### **demo_modeling** directory
This directory contains the input (including the observed multiSFS of 1196 SNPs across 280 larvae) and output files from demographic modeling using the fastsimcoal program. Additional scripts used for data preparation and the script to submit the job to a slurm scheduler on a shared computing cluster is also provided.

   • **Ne_boxplots.png** shows ML point estimates and CIs for the best model  
   • **Ne_observedSFS.R** plots observed SFS and expected SFSs from demographic modeling  
   • **all_demo_scenarios.pdf** figure of tested demographic scenarios  
   • **model6_lineplot_a_and_b.png** shows how point estimates and CIs of Ne from the fastsimcoal demographic modeling change over time for Model 6  
   • **obs_sfs_polyonly.png** observed SFS for each larval cohort  
   • **obs_sfs_polyonly_nomafnomac.png** observed SFS for each larval cohort resulting from no --mac filter    
   • **obs_sfs_polyonly_nosibs.png** observed SFS for each larval cohort with one individual from putative sibling pairs removed

``` fsc_models/ ```
Contains the .est and .tpl files necessary to run each of the six models for the main analyses in fastsimcoal. Also contains the observed MSFS whose name must match the .est and .tpl file names when running fastsimcoal. Example slurm scripts for model selection and CI estimation, as well as helper scripts for concatenating fsc .bestlhoods files when fsc is run many times. These files can serve as templates for running the sensitivity analyses for dataset resulting from no --mac filter.

   • **cat_bestlhoods.sh** an example of how to concatenate fsc results for downstream model selection   
   • **cat_cis.sh** an example of how to concatenate fsc results when fsc is run on many simulated SFSs  
   • **cat_max_summary.sh** an example of how to select the ML run for many simulated SFSs and concatenate for CI estimation  
   • **run_fsc_Ne_CI_singlethread.sh** an example slurm script for running fsc on many simulated SFSs for CI estimation    
   • **run_model6_singlethread.sh** an example of how to submit a fastsimcoal run to a shared computing cluster with each job on a single thread (this is faster than multithreading)  
   • **run_model6_multithread.sh** an example of how to submit a fastsimcoal run to a shared computing cluster (this can be slow) 

``` model_results/ ```
Contains the fastsimcoal results for each model and the CIs for the best model resulting from the main demographic modeling analyses.

``` nomac_models/ ```
Contains the fastsimcoal results for each model and the CIs for the best model resulting from sensitivity analyses using the dataset when no --mac filter was applied.

``` nosibs_models/ ```
Contains the fastsimcoal results for each model resulting from sensitivity analyses when putative sibs were removed from the dataset used in the main analyses.

``` power/ ```
Contains a script to read in and plot the results of simulations using pseudo-observed datasets. Also contains the model fits for each of the PODs simulated using the ML parameters of each model. These simulations were based off of the dataset used in the main analyses.

``` sim_sfs/ ```
Contains a figure comparing the observed SFS to those of the best demographic models, with SFSs comparisons broken out by larval fish cohort. An example slurm script is also provided for generating parametric bootstrapped SFSs from the ML parameters of a model. Based off of the dataset used in the main analyses.

## Manuscript Figures

• [Nc_and_Nb_overtime](https://github.com/pinskylab/NePADE/blob/master/Nc_and_Nb_overtime.png) plot of census and breeding abundance over time calculated from the 2016 stock assessment (Figure 1)

• [all_demo_scenarios.pdf](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/all_demo_scenarios.pdf) figure of tested demographic scenarios (Figure 2)

• [model6_lineplot_a_and_b.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model6_lineplot_a_and_b.png) shows how point estimates and CIs of Ne from the fastsimcoal demographic modeling change over time for Model 6 (Figure 3)

• [Ne_boxplots.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/Ne_boxplots.png) shows ML point estimates and CIs for the best model (Figure S1)

• [sfs.comparison.top3.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/sfs.comparison.top3.png) shows SFSs by cohort comparing observed vs simulated SFSs of three best-fit models (Figure S2)

• [/demo_modeling/power/model_confusion_matrix.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/model_confusion_matrix.png) confusion matrix showing the power in the dataset to differentiate between demographic models (Figure S4)
