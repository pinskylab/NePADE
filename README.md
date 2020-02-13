# NePADE
Scripts for investigating Ne over time in summer flounder

• **best_lhoods.R** script for determining maximum-likelihood of each demographic scenario and calculating AIC

• **PADE_stock_assessment16.R** script imports data from 2016 summer flounder stock assessment and calculates the number of potential breeders at each age and year

• **Ne.R** script for figuring out which summer flounder to resequence

#### **demo_modeling** directory
This directory contains figure tested demographic scenarios, scripts for data preparation, input files necessary to run fastsimcoal and the associated script to submit the job to a slurm scheduler on a shared computing cluster  

   • **from_amarel/growing_nomaf** fastsimcoal input files for exponential growth following a popuation bottleneck  
   • **from_amarel/nobot_nomaf_oneparam** fastsimcoal input files for constant population size  
   • **from_amarel/shrinking_nomaf** fastsimcoal input files for exponential decay following a popuation bottleneck  
   • **from_amarel/stable_instant_nomaf_Nlowerlimit100** fastsimcoal input files for instantaneous recovery following population bottleneck  
   • **Ne_observedSFS.R** prepares SNP data for SFS required by fastsimcoal  
   • **all_demo_scenarios.pdf** figure of demographic scenarios  
   • **stable_nomaf_bootstrap_cis.R** script for parametric bootstrapped CIs  
