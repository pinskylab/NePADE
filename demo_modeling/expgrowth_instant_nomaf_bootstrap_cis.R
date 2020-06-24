#### Plotting the results from the fsc26 bootstrapping ####

# Read in the parameters estimated from simulated SFS (nonparametric)
cis <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/nonparametric_ci_summary.txt', header = TRUE) # Read in ML bootstrapped parameters following 10 runs for each simulated SFS

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
hist(log10(cis$NBOT/cis$NPOP08), xlab = 'log10(NBOT/NPOP08)', main = 'NBOT/NPOP08 ratio') # same as hist(log10(cis$NBOT)-log10(cis$NPOP08))
abline(v = log10(best$NBOT/best$NPOP08), col = 'tomato')
hist(log10(cis$NANC/cis$NBOT), xlab = 'log10(NANC/NBOT)', main = 'NANC/NBOT ratio')
abline(v = log10(best$NANC/best$NBOT), col = 'tomato')

# Boxplots of population size from simulated SFS & point estimates of the ML parameters
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_boxplots.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# boxplot(log10(cis$NANC), log10(cis$NBOT), log10(cis$NPOP08), ylab = expression('log'[10]*'(N'[e]* 'estimate)'), names = c('NANC', 'NBOT', 'NPOP08')) # haploid
# points(c(1,2,3), c(log10(best$NANC), log10(best$NBOT), log10(best$NPOP08)), col = 'tomato', pch = 19)
boxplot(log10(cis$NANC/2), log10(cis$NBOT/2), log10(cis$NPOP08/2), ylab = expression('log'[10]*'(N'[e]* ' estimate)'), names = c('NANC', 'NBOT', 'NPOP08')) # convert Ne to diploid by dividing haploid number by 2
points(c(1,2,3), c(log10(best$NANC/2), log10(best$NBOT/2), log10(best$NPOP08/2)), col = 'tomato', pch = 19)

dev.off()

# Calculate 95% confidence intervals
# NPOP08
# mean(cis$NPOP08) + 1.960*(sd(cis$NPOP08)/sqrt(100))
# mean(cis$NPOP08) - 1.960*(sd(cis$NPOP08)/sqrt(100))
quantile(cis$NPOP08/2, c(0.025, 0.975)) # diploid

# NANC
quantile(cis$NANC/2, c(0.025, 0.975)) # diploid

# NBOT
quantile(cis$NBOT/2, c(0.025, 0.975)) # diploid

# TBOT
quantile(cis$TBOT, c(0.025, 0.975))

# TLEN
quantile(cis$TLEN, c(0.025, 0.975))

# RANC
quantile(cis$RANC, c(0.025, 0.975))

# NPOP08/NBOT ratio
best$NPOP08/best$NBOT # point estimate, same regardless of haploid or diploid
quantile((cis$NPOP08/cis$NBOT), c(0.025, 0.975))

# NBOT/NANC ratio
best$NBOT/best$NANC # point estimate, same regardless of haploid or diploid
quantile((cis$NBOT/cis$NANC), c(0.025, 0.975))

