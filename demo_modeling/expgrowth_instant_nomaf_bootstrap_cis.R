#### Plotting the results from the fsc26 bootstrapping ####

# Read in the parameters estimated from simulated SFS
cis <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/cis_sfs_summary.txt', header = TRUE) # Read in ML bootstrapped parameters following 50 runs for each simulated SFS

# Read in ML parameters from best fit model
expgrowth_instant.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.expgrowth_instant_nomaf.txt", header = TRUE)# Read in ML parameters from instantaneous recovery model
max(expgrowth_instant.best.nomaf$MaxEstLhood)
best <- expgrowth_instant.best.nomaf[which(expgrowth_instant.best.nomaf$MaxEstLhood == max(expgrowth_instant.best.nomaf$MaxEstLhood)),]

# Histograms of parameters estimated from simulated SFS, plus a line for the point estimates
hist(log10(cis$NPOP08), main = 'Population size in 2008', xlab = 'log10(NPOP08)')
abline(v = log10(best$NPOP08), col = 'tomato')
hist(log10(cis$NANC), main = 'Ancestral population size', xlab = 'log10(NANC)')
abline(v=log10(best$NANC), col = 'tomato')
hist(log10(cis$NBOT), main = 'Population size during bottleneck', xlab = 'log10(NBOT')
abline(v=log10(best$NBOT), col = 'tomato')
hist(cis$TBOT, main = 'Time after bottleneck', xlab = 'Generations (2 years/generation)')
abline(v=best$TBOT, col = 'tomato')
hist(cis$TLEN, main = 'Length of bottleneck', xlab = 'Generations (2 years/generation)', breaks = 5)
abline(v=best$TLEN, col = 'tomato')
hist(cis$RANC, main = 'NANC growth rate', xlab = 'Growth rate of NANC going back in time\n(Negative means positive growth)', xlim = c(-0.007,-0.0001))
abline(v=best$RANC, col = 'tomato')
hist(cis$MaxEstLhood)
hist(cis$MaxObsLhood)

# Boxplots of population size from simulated SFS & point estimates of the ML parameters
boxplot(log10(cis$NANC), log10(cis$NBOT), log10(cis$NPOP08), ylab = expression('log'[10]*'(N'[e]* 'estimate)'), names = c('NANC', 'NBOT', 'NPOP08'))
points(c(1,2,3), c(log10(best$NANC), log10(best$NBOT), log10(best$NPOP08)), col = 'tomato', pch = 19)

# Calculate 95% confidence intervals
# NPOP08
# mean(cis$NPOP08) + 1.960*(sd(cis$NPOP08)/sqrt(100))
# mean(cis$NPOP08) - 1.960*(sd(cis$NPOP08)/sqrt(100))
quantile(cis$NPOP08, c(0.025, 0.975))

# NANC
quantile(cis$NANC, c(0.025, 0.975))

# NBOT
quantile(cis$NBOT, c(0.025, 0.975))

# TBOT
quantile(cis$TBOT, c(0.025, 0.975))

# TLEN
quantile(cis$TLEN, c(0.025, 0.975))

# RANC
quantile(cis$RANC, c(0.025, 0.975))

