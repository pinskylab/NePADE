This repository contains data, scripts and figures for the analyses in "Genetic decline and recovery of a demographically rebuilt fishery species, " which investigates how genetic diversity and effective population size have changed over time for summer flounder.

*These data, methods and scripts are provided in the interests of open science. If you have questions or find errors, please let me know.*

Citable as

Contact: Jennifer Hoey [(jahoey13@gmail.com)](mailto:jahoey13@gmail.com)

## Data

**Note**: FASTQ files for individuals involved in analyses can be downloaded from the NCBI Sequency Read Archive under BioProject: PRJNA750099. Downstream genetic data files, as well as stock assessment data used to calculate summer flounder abundance, can be found in the ``` data/ ``` directory as follows:

### Within ``` data/ ```:

• [Ne279_1068loci.arp](https://github.com/pinskylab/NePADE/blob/master/data/Ne279_1068loci.arp): Arlequin input file for dataset of 1068 loci across 279 summer flounder

• [Ne279_1068loci_MAFpop0.obs](https://github.com/pinskylab/NePADE/blob/master/data/Ne279_1068loci_MAFpop0.obs): observed mSFS for 2008 cohort

• [Ne279_1068loci_MAFpop1.obs](https://github.com/pinskylab/NePADE/blob/master/data/Ne279_1068loci_MAFpop1.obs): observed mSFS for 1997 cohort

• [Ne279_1068loci_MAFpop2.obs](https://github.com/pinskylab/NePADE/blob/master/data/Ne279_1068loci_MAFpop2.obs): observed mSFS for 1994 cohort

• [Ne279_1068loci_MSFS.obs](https://github.com/pinskylab/NePADE/blob/master/data/Ne279_1068loci_MSFS.obs): observed multi mSFS used for demographic modeling in fastsimcoal2 (sometimes referred to as fsc in this README too)

• [Ne_279PADE_1068loci_complete.gen](https://github.com/pinskylab/NePADE/blob/master/data/Ne_279PADE_1068loci_complete.gen): genepop file containing subset of full dataset with no missing data, 1068 SNPs across 279 larvae, used for demographic modeling

• [Ne_279PADE_1068loci_complete_correctnot.txt](https://github.com/pinskylab/NePADE/blob/master/data/Ne_279PADE_1068loci_complete_correctnot.txt): same as [Ne_279PADE_1068loci_complete.gen](https://github.com/pinskylab/NePADE/blob/master/data/Ne_279PADE_1068loci_complete.gen), but as a text file, used to generate .arp file for Arlequin

• [Ne_279PADE_3749loci_missingallowed.gen](https://github.com/pinskylab/NePADE/blob/master/data/Ne_279PADE_3749loci_missingallowed.gen): genepop file after appling HW filter based on [Ne_HWP_test.txt](https://github.com/pinskylab/NePADE/blob/master/data/Ne_HWP_test.txt)

• [Ne_HWP_test.txt](https://github.com/pinskylab/NePADE/blob/master/data/Ne_HWP_test.txt): results of HW test in [02_additional_filters.R](https://github.com/pinskylab/NePADE/blob/master/02_additional_filters.R)

• [PADE_stock_assessment16.xlsx](https://github.com/pinskylab/NePADE/blob/master/data/PADE_stock_assessment16.xlsx): number of fish at age, estimated fishing mortality at age, natural mortality at age & proportion mature at age for 1982-2015

• **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.vcf**: vcf file containing all available SNPs on a contig (3905 SNPs for 284 larvae)

• **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.vcf**: vcf file containing the first SNP on a contig

• [SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen](https://github.com/pinskylab/NePADE/blob/master/data/SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen): genepop file generated from SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.vcf. This genepop file was then used to evaluate individuals with high heterozygosity and loci not in HWP.

#### ``` bayescan_output/ ``` <br>

Contains output from Bayescan looking for temporal outlier loci in dataset used for demographic modeling (1068 loci across 279 summer flounder). Bayescan was run on Rutgers School of Environmental & Biological Sciences Annotate server.

## Methods

[Ne_bioinformatics.txt](https://github.com/pinskylab/NePADE/blob/master/Ne_bioinformatics.txt) contains bioinformatic steps taken to go from sequencing reads to SNP genotypes, as well as steps taken to produce filtered SNP genotypes in **SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.vcf**

## Scripts & Results
### The numbered scripts indicate the order of their usage for the demographic analysis in fastsimcoal2.
#### These R scripts were run using R v.4.0.3 (at least) on a MacBook Pro running OS Catalina (at least).

• [00_Ne_genotyping_results.R](https://github.com/pinskylab/NePADE/blob/master/00_Ne_genotyping_results.R): script for calculating coverage statistics

• [01_keepfirstsnponly.R](https://github.com/pinskylab/NePADE/blob/master/01_keepfirstsnponly.R): reads in full vcf and writes a new vcf containing only the first SNP on a contig

• [02_additional_filters.R](https://github.com/pinskylab/NePADE/blob/master/02_additional_filters.R): removes several fish based on high individual heterozygosity. Then removes loci not in HWE.

• [02a_Ne_diversity_analyses.R](https://github.com/pinskylab/NePADE/blob/master/02a_Ne_diversity_analyses.R): script for plotting PCAs and calculating diversity statistics

• [03_prepNe_observedSFS.R](https://github.com/pinskylab/NePADE/blob/master/03_prepNe_observedSFS.R): additional filtering steps necessary for demographic modeling

• [03a_bayescan_plot_R.R](https://github.com/pinskylab/NePADE/blob/master/03a_bayescan_plot_R.R): plots output from Bayescan analysis looking for temporal outlier loci

• [03b_plot_observedSFS.R](https://github.com/pinskylab/NePADE/blob/master/03b_plot_observedSFS.R): plots the observed msfs for each larval cohort

• [04_best_lhoods.R](https://github.com/pinskylab/NePADE/blob/master/04_best_lhoods.R): script for determining maximum-likelihood of each demographic scenario, calculating AIC and plotting demographic modeling results

• [04a_plot_SFScomparisions.R](https://github.com/pinskylab/NePADE/blob/master/04a_plot_SFScomparisions.R): script to plot the observed msfs for each larval cohort against the top three demographic models

• [04b_plot_SFScomparisions_temporal_test_240fish.R](https://github.com/pinskylab/NePADE/blob/master/04b_plot_SFScomparisions_temporal_test_240fish.R): script to plot the expected msfs for each larval cohort under the best-fit model with equal sampling across larval cohorts

• [generation_length.R](https://github.com/pinskylab/NePADE/blob/master/generation_length.R): script calculates generation length for females and males over time using the 2016 stock assessment data, age-length relationships from Penttila et al. (1989) and age-fecundity curves from Morse (1981)

• [PADE_stock_assessment16.R](https://github.com/pinskylab/NePADE/blob/master/PADE_stock_assessment16.R): script imports data from 2016 summer flounder stock assessment, calculates the number of potential breeders at each age and year and plots census and breeding abundance over time

### Within ``` demo_modeling/ ```: 
This directory contains the input and output files from demographic modeling using the fastsimcoal2 program. Additional scripts to submit the job to a slurm scheduler on a shared computing cluster and for manipulating output files are also provided. All demographic modeling analyses were performed on Rutgers Amarel computing cluster.

   • [Ne_boxplots.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/Ne_boxplots.png): shows ML point estimates and CIs for the best model (Figure S2)    
   • [all_demo_scenarios.pdf](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/all_demo_scenarios.pdf): figure of tested demographic scenarios (Figure 2)  
   • [model6_lineplot_a_and_b.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model6_lineplot_a_and_b.png): shows how point estimates and CIs of Ne from the fastsimcoal2 demographic modeling change over time for Model 6 (Figure 3)  
   • [obs_sfs_polyonly.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/obs_sfs_polyonly.png): observed SFS for each larval cohort using dataset of 1068 loci across 279 fish  

#### ``` fsc_models/ ``` <br>
Each of the ``` modelX/``` directories, where X refers to a model described in Figure 2, contains the .est and .tpl files necessary to run that particular demographic model in fastsimcoal2. Models are based off of our observed dataset of 279 summer flounder larvae and 1068 loci. The observed mSFS whose name must match the .est and .tpl file names when running fastsimcoal2 is in the ``` data/ ``` directory. Example slurm scripts for model selection and CI estimation, as well as helper scripts for concatenating fastsimcoal2 .bestlhoods files when fastsimcoal2 is run many times. See below for details:

   • [cat_bestlhoods.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_bestlhoods): an example of how to concatenate fsc results for downstream model selection   
   • [cat_cis.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_cis): an example of how to concatenate fsc results when fsc is run on many simulated SFSs (this was done for CI estimation)  
   • [cat_max_summary.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/cat_max_summary): an example of how to select the ML run from many simulated SFSs and concatenate for CI estimation   
   • [run_fsc_Ne_CI_singlethread.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_fsc_Ne_CI_singlethread.sh): an example slurm script for running a fsc job array on many simulated SFSs for CI estimation    
   • [run_model_Ne279_1068loci.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_model_Ne279_1068loci.sh): example slurm script for submitting fsc jobs based on the 1068 loci across 279 larvae dataset   
   • [run_sbatch_Ne.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/fsc_models/run_sbatch_Ne.sh): helper script to submit a script to SLURM many times

#### ``` model_results/ ``` <br>
This directory contains the fastsimcoal2 demographic modeling results for each of the seven models (plus any associated sensitivity analyses) and the CIs for the best-fit model.

   • [model1.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model1.bestlhoods.summary.txt): contains Model 1 estimated parameters     
   • [model2.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model2.bestlhoods.summary.txt): contains Model 2 estimated parameters    
   • [model3.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model3.bestlhoods.summary.txt): contains Model 3 estimated parameters    
   • [model4.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model4.bestlhoods.summary.txt): contains Model 4 estimated parameters    
   • [model5.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model5.bestlhoods.summary.txt): contains Model 5 estimated parameters    
   • [model6.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model6.bestlhoods.summary.txt): contains Model 6 estimated parameters    
   • [model6.nonparametric.ci.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model6.nonparametric.ci.summary.txt): contains Model 6 estimated parameters based on 100 SFSs generated via nonparametric bootstrapping   
   • [model6_fixTLEN.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model6_fixTLEN.bestlhoods.summary.txt): contains Model 6 estimated parameters for when TLEN was fixed at three generations   
   • [model6_wider_priors.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model6_wider_priors.bestlhoods.summary.txt): contains Model 6 estimated parameters for when priors were allowed to be wider  
   • [model7.bestlhoods.summary.txt](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model_results/model7.bestlhoods.summary.txt): contains Model 7 estimated parameters     

#### ``` power/ ``` <br>
This directory represents a series of fastsimcoal2 simulations to determine the power within our data for distinguishing among the demographic hypotheses. Contains the parameter estimates resulting from fitting the seven demographic models to each of the pseudo-observed datasets (POD) simulated from the ML parameters of each model, a script to determine the best-fit model for each POD, and a scripts to visualize these results. These simulations were based off of our observed dataset of 279 summer flounder larvae and 1068 loci. Model fits for each POD are housed within a directory where the model number is idenfified with an X: ``` modXpods/ ```. Additional details about files in ``` power/ ``` are described below:

   • [PODs_fit.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/PODs_fit.png): shows Ne estimates over time for each model described in Figure 2 and the inferred history of 10 PODs based on each model's ML parameter estimates (Figure S4)  
   • [PODs_fits.R](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/PODs_fits.R): script plots the ML parameter estimates for each model over time, reads in the results from each ``` modXpods/ ``` directory, calculates the best-fit model for each POD, and plots these over time. The result is [PODs_fit.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/PODs_fit.png).   
   • [model_confusion_matrix.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/model_confusion_matrix.png): confusion matrix showing the number of times the inferred model corresponded with the known generating model (Figure S3)   
   • [power.R](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/power.R): script reads in simulation results and plots the number of times the inferred model corresponded with the known generating model. The result is [model_confusion_matrix.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power/model_confusion_matrix.png).
   
#### ``` power_even_sample_size/ ``` <br>
This directory represents another series of simulations designed to assess the power for inferring the correct demographic model when equal numbers of fish were sampled across larval cohorts. As a result, this directory contains results from fitting the seven demographic models to 50 PODs simulated from the best-fit model. 80 diploids were sampled from each cohort in this analysis. Each POD is represented by a unique file, where X represents the POD number: model6_maxL_X.summary.txt. Additional details about files in ``` power_even_sample_size/ ``` are described below:

   • [PODs_determine_AIC.R](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_determine_AIC.R): script reads in the results from model-fitting and selects the best-fit model for each POD using AIC   
   • [PODs_fit_even_sampling.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_fit_even_sampling.png): shows Ne estimates over time for the best-fit model and the inferred history of 50 pseudo-observed datasets simulated from the best-fit model with equal sample sizes across cohorts (Figure S5)     
   • [PODs_plot_best_fits.R](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_plot_best_fits.R): script is similar to [PODs_determine_AIC.R](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_determine_AIC.R), but also plots [PODs_fit_even_sampling.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/power_even_sample_size/PODs_fit_even_sampling.png).
   
#### ``` sim_sfs/ ``` <br>
This directory contains scripts for simulating and visualizing mSFSs under particular demographic models and/or conditions. Unless otherwise specified, SFSs are based off of the observed dataset of 279 summer flounder larvae and 1068 loci. See below for more details:

   • [exp_sfs_model6.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/exp_sfs_model6.png): shows the expected minor SFS by larval cohort based on the best-fit model (Figure S7)     
   • [exp_sfs_model6_temporal_test_240fish_1068loci.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/exp_sfs_model6_temporal_test_240fish_1068loci.png): shows the expected mSFS by larval cohort based on the best-fit model with equal sampling 80 diploids in each cohort  
   • [run_fsc_boot.sh](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/run_fsc_boot.sh): is an example SLURM script for generating parametric bootstrapped minor SFSs based on ML parameters for the best-fit model using fastsimcoal2   
   • [sfs.comparison.top3.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/sfs.comparison.top3.png): shows the mSFSs for the observed dataset and the three best-fit demographic scenarios by larval cohort  

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
