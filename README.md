This repository contains data, scripts and figures for a project investigating how genetic diversity and effective population size has changed over time for summer flounder

## Scripts

• **best_lhoods.R** script for determining maximum-likelihood of each demographic scenario, calculating AIC and plotting demographic modeling results

• **generation_length.R** script calculates generation length for females and males over time using the 2016 stock assessment data, age-length relationships from Penttila et al. (1989) and age-fecundity curves from Morse (1981)

• **PADE_stock_assessment16.R** script imports data from 2016 summer flounder stock assessment, calculates the number of potential breeders at each age and year and plots census and breeding abundance over time

• **Ne.R** script for figuring out which summer flounder to resequence

• **Ne_analyses.R** exploratory script for plotting PCAs and calculating diversity statistics

• **Ne_genotyping_results.R** script for calculating coverage statistics

### **demo_modeling** directory
This directory contains the input and output files from demographic modeling using the fastsimcoal program. Additional scripts used for data preparation and the script to submit the job to a slurm scheduler on a shared computing cluster is also provided.

   • **Ne_boxplots.png** shows ML point estimates and CIs for the best model  
   • **Ne_observedSFS.R** plots observed SFS and expected SFSs from demographic modeling  
   • **all_demo_scenarios.pdf** figure of tested demographic scenarios  
   • **model6_lineplot_a_and_b.png** shows how point estimates and CIs of Ne from the fastsimcoal demographic modeling change over time for Model 6  
   • **obs_sfs_polyonly** observed SFS for each larval cohort  
   • **obs_sfs_polyonly_nosibs** observed SFS for each larval cohort with one individual from putative sibling pairs removed

``` fsc_model/ ```
Contains the .est and .tpl files necessary to run each of the six models in fastsimcoal. Also contains the observed MSFS whose name must match the .est and .tpl file names when running fastsimcoal.

   • **run_model6_singlethread.sh** an example of how to submit a fastsimcoal run to a shared computing cluster with each job on a single thread (this is faster than multithreading)  
   • **run_model6_multithread.sh** an example of how to submit a fastsimcoal run to a shared computing cluster (this can be slow) 

``` model_results ```
Contains the fastsimcoal results for each model and the CIs for the best model 

``` sim_sfs ```
Contains a figure comparing the observed SFS to those of the best demographic models, with SFSs comparisons broken out by larval fish cohort. An example slurm script is also provived for generating parametric bootstrapped SFSs from the ML parameters of a model

## Manuscript Figures

• [Nc_and_Nb_overtime](https://github.com/pinskylab/NePADE/blob/master/Nc_and_Nb_overtime.png) plot of census and breeding abundance over time calculated from the 2016 stock assessment (Figure 1)

• [all_demo_scenarios.pdf](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/all_demo_scenarios.pdf) figure of tested demographic scenarios (Figure 2)

• [model6_lineplot_a_and_b.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/model6_lineplot_a_and_b.png) shows how point estimates and CIs of Ne from the fastsimcoal demographic modeling change over time for Model 6 (Figure 3)

• [Ne_boxplots.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/Ne_boxplots.png) shows ML point estimates and CIs for the best model (Figure S1)

• [sfs.comparison.top3.png](https://github.com/pinskylab/NePADE/blob/master/demo_modeling/sim_sfs/sfs.comparison.top3.png) shows SFSs by cohort comparing observed vs simulated SFSs of three best-fit models (Figure S2)
