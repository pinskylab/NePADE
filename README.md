This repository contains data, scripts and figures for the analyses in "Genetic decline and recovery of a demographically rebuilt fishery species, " which investigates how genetic diversity and effective population size have changed over time for summer flounder.

*These data, methods and scripts are provided in the interests of open science. If you have questions or find errors, please let me know.*

Citable as

Contact: Jennifer Hoey [(jahoey13@gmail.com)](mailto:jahoey13@gmail.com)

## Data

**Note**: FASTQ files for individuals involved in analyses can be downloaded from the NCBI Sequency Read Archive under BioProject: PRJNA750099. Downstream genetic data files, as well as stock assessment data used to calculate summer flounder abundance, can be found in the ``` data/ ``` directory as follows:

### Within ``` data/ ```:

``` bayescan_output/ ```
Contains output from Bayescan looking for temporal outlier loci in dataset used for demographic modeling (1068 loci across 279 summer flounder)

• **Ne279_1068loci.arp**: Arlequin input file for dataset of 1068 loci across 279 summer flounder

• **Ne279_1068loci_MAFpop0.obs**: msfs for 2008 cohort

• **Ne279_1068loci_MAFpop1.obs**: msfs for 1997 cohort

• **Ne279_1068loci_MAFpop2.obs**: msfs for 1994 cohort

• **Ne279_1068loci_MSFS.obs**: multi msfs used for demographic modeling in fastsimcoal

• **Ne_279PADE_1068loci_complete.gen**: genepop file containing subset of full dataset with no missing data, 1068 SNPs across 279 larvae, used for demographic modeling

• **Ne_279PADE_1068loci_complete_correctnot.txt**: same as **Ne_279PADE_1068loci_complete.gen**, but as a text file, used to generate .arp file for Arlequin

• **Ne_279PADE_3749loci_missingallowed.gen**: genepop file after appling HW filter based on **Ne_HWP_test.txt**

• **Ne_HWP_test.txt**: results of HW test in **02_additional_filters.R**

• **PADE_stock_assessment16.xlsx**: number of fish at age, estimated fishing mortality at age, natural mortality at age & proportion mature at age for 1982-2015

• **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.vcf**: vcf files containing all available SNPs on a contig

• **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen**: genepop file generated from SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.gen by retaining only first SNP on a contig, full datset is composed of 3905 SNPs for 284 larvae


## Scripts & Results
### The numbered scripts indicate the order of their usage for the demographic analysis in fastsimcoal2

• **00_Ne_genotyping_results.R** script for calculating coverage statistics

• **01_keepfirstsnponly.R**: Reads in full vcf and writes a new vcf containing only the first SNP on a contig

• **02_additional_filters.R**: Removes several fish based on high individual heterozygosity. Then removes loci not in HWE.

• **02a_Ne_diversity_analyses.R**: script for plotting PCAs and calculating diversity statistics

• **03_prepNe_observedSFS.R**: additional filtering steps necessary for demographic modeling

• **03a_bayescan_plot_R.R**: plots output from Bayescan analysis looking for temporal outlier loci

• **03b_plot_observedSFS.R**: plots the observed msfs for each larval cohort

• **04_best_lhoods.R** script for determining maximum-likelihood of each demographic scenario, calculating AIC and plotting demographic modeling results

• **04a_plot_SFScomparisions.R**: script to plot the observed msfs for each larval cohort against the top three demographic models

• **04b_plot_SFScomparisions_temporal_test_240fish.R**: script to plot the expected msfs for each larval cohort under the best-fit model with equal sampling across larval cohorts

• **generation_length.R** script calculates generation length for females and males over time using the 2016 stock assessment data, age-length relationships from Penttila et al. (1989) and age-fecundity curves from Morse (1981)

• **PADE_stock_assessment16.R** script imports data from 2016 summer flounder stock assessment, calculates the number of potential breeders at each age and year and plots census and breeding abundance over time

• **Ne.R** script for figuring out which summer flounder to resequence

### **demo_modeling** directory
This directory contains the input and output files from demographic modeling using the fastsimcoal2 program. Additional scripts used for data preparation and the script to submit the job to a slurm scheduler on a shared computing cluster are also provided.

   • **Ne_boxplots.png** shows ML point estimates and CIs for the best model  
   • **all_demo_scenarios.pdf** figure of tested demographic scenarios  
   • **model6_lineplot_a_and_b.png** shows how point estimates and CIs of Ne from the fastsimcoal2 demographic modeling change over time for Model 6  
   • **obs_sfs_polyonly.png** observed SFS for each larval cohort  

``` fsc_models/ ```
Contains the .est and .tpl files necessary to run each of the demographic models for the main analyses in fastsimcoal2. The observed MSFS whose name must match the .est and .tpl file names when running fastsimcoal2 is in the ``` data/ ``` directory. Example slurm scripts for model selection and CI estimation, as well as helper scripts for concatenating fsc .bestlhoods files when fsc is run many times.

   • **cat_bestlhoods.sh** an example of how to concatenate fsc results for downstream model selection   
   • **cat_cis.sh** an example of how to concatenate fsc results when fsc is run on many simulated SFSs (this was done for CI estimation)  
   • **cat_max_summary.sh** an example of how to select the ML run from many simulated SFSs and concatenate for CI estimation   
   • **run_fsc_Ne_CI_singlethread.sh** an example slurm script for running a fsc job array on many simulated SFSs for CI estimation    
   • **run_model_Ne279_1068loci.sh** example slurm script for submitting fsc jobs based on the 1068 loci across 279 dataset   
   • **run_sbatch_Ne.sh** helper script to submit a script to SLURM many times

``` model_results/ ```
Contains the fastsimcoal2 results for each model and the CIs for the best model resulting from the main demographic modeling analyses.

``` power/ ```
Contains a script to read in and plot the results of simulations using pseudo-observed datasets. Also contains the model fits for each of the PODs simulated using the ML parameters of each model. These simulations were based off of the dataset used in the main analyses.

``` sim_sfs/ ```
Contains a figure comparing the observed SFS to those of the best demographic models, with SFSs comparisons broken out by larval fish cohort. An example slurm script is also provided for generating parametric bootstrapped SFSs from the ML parameters of a model. Based off of the dataset used in the main analyses.

## Manuscript Figures

• [Nc_and_Nb_overtime.png](https://github.com/pinskylab/NePADE/blob/master/Nc_and_Nb_overtime.png) plot of census and breeding abundance over time calculated from the 2016 stock assessment (Figure 1)

• [all_demo_scenarios.pdf](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/all_demo_scenarios.pdf) figure of tested demographic scenarios (Figure 2)

• [model6_lineplot_a_and_b.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model6_lineplot_a_and_b.png) shows how point estimates and CIs of Ne from the fastsimcoal2 demographic modeling change over time for Model 6 (Figure 3)

• [pade_generation_length.png](https://https://github.com/pinskylab/NePADE/blob/master/pade_generation_length.png) shows calculated summer flounder generation length over time (Figure S1)

• [Ne_boxplots.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/Ne_boxplots.png) shows ML point estimates and CIs for the best model (Figure S2)

• [/demo_modeling/power/model_confusion_matrix.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/model_confusion_matrix.png) confusion matrix showing the power in the dataset to differentiate between demographic models (Figure S3)  

• [/demo_modeling/power/PODs_fit.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/PODs_fit.png) Ne estimates over time for each model described in Figure 2 and the inferred history of 10 PODs based on each model's ML parameter estimates (Figure S4)

• [/demo_modeling/power_even_sample_size/PODs_fit_even_sampling.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_fit_even_sampling.png) Ne estimates over time for the best-fit model and the inferred history of 50 PODs with equal sampling across larval cohorts (Figure S5)

• [sfs.comparison.top3.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/sfs.comparison.top3.png) shows SFSs by cohort comparing observed vs simulated SFSs of three best-fit models (Figure S6)

• [/demo_modeling/sim_sfs/exp_sfs_model6.png](NePADE/demo_modeling/sim_sfs/exp_sfs_model6.png) shows expected msfs by larval cohort based on the best-fit model (Figure S7)

• [/demo_modeling/sim_sfs/exp_sfs_model6_temporal_test_240fish_1068loci.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/exp_sfs_model6_temporal_test_240fish_1068loci.png) shows expected msfs by larval cohort based on the best-fit model with equal sampling sizes across cohorts (Figure S8)
