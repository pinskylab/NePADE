This repository contains data, scripts and figures for a project investigating how genetic diversity and effective population size has changed over time for summer flounder

## Scripts

• **best_lhoods.R** script for determining maximum-likelihood of each demographic scenario and calculating AIC

• **generation_length.R** script calculates generation length for females and males over time using the 2016 stock assessment data, age-length relationships from Penttila et al. (1989) and age-fecundity curves from Morse (1981)

• **PADE_stock_assessment16.R** script imports data from 2016 summer flounder stock assessment and calculates the number of potential breeders at each age and year

• **Ne.R** script for figuring out which summer flounder to resequence

• **Ne_analyses.R** exploratory script for plotting PCAs and calculating diversity statistics

• **Ne_genotyping_results.R** script for calculating coverage statistics

### **demo_modeling** directory
This directory contains the input and output files from demographic modeling using the fastsimcoal program. Additional scripts used for data preparation and the script to submit the job to a slurm scheduler on a shared computing cluster is also provided.

#### fsc_model directory 
Contains the .est and .tpl files necessary to run each model in fastsimcoal. Also contains the observed MSFS whose name must match the .est and .tpl file names when running fastsimcoal.

• **run_model4.sh** an example of how to submit a fastsimcoal run to a shared compusing cluster 

#### model_results directory
Contains the fastsimcoal results for each model and the CIs for the best model 

#### sim_sfs directory 
Contains a figure comparing the observed SFS to those of the best demographic models, with SFSs comparisons broken out by larval fish cohort

   • **Ne_boxplots.png** shows ML point estimates and CIs for the best model  
   • **Ne_observedSFS.R** prepares SNP data for SFS required by fastsimcoal  
   • **all_demo_scenarios.pdf** figure of tested demographic scenarios  
   • **model4_lineplot_a_and_b.png** shows how point estimates and CIs of Ne from the fastsimcoal demographic modeling change over time for Model 4  
   • **obs_sfs_polyonly** observed SFS for each larval cohort  

## Figures

• **Nc_and_Nb_overtime** plot of census and breeding abundance over time calculated from the 2016 stock assessment
