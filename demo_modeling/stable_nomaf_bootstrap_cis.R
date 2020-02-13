#### Plotting the results from the fsc26 bootstrapping ####

# Read in the parameters estimated from simulated SFS
cis <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/cis_summary.txt', header = TRUE)

# Read in ML parameters from best fit model
stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.stable_nomaf.txt", header = TRUE)
max(stable.best.nomaf$MaxEstLhood)
best <- stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]

# Histograms of parameters estimated from simulated SFS
hist(log10(cis$NPOP08), main = 'Population size in 2008', xlab = 'log10(NPOP08)')
hist(log10(cis$NANC), main = 'Ancestral population size', xlab = 'log10(NANC')
hist(log10(cis$NBOT), main = 'Population size during bottleneck', xlab = 'log10(NBOT')
hist(cis$TBOT, main = 'Time after bottleneck', xlab = 'Years (2 x # generations)')
hist(cis$TLEN, main = 'Length of bottleneck', xlab = 'Years (2 x # generations)', breaks = c(1:4))
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


