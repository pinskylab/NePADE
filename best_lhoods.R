#### Fsc26 has been run 50 times and the maximum likelihoods and associated parameters have been concatenated ####
# Using a SFS that summarizes 280 larvae across three cohorts, 1196 loci and no MAF filter
# Read in the estimated parameters for each model
mod1 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model1.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model2.bestlhoods.summary.txt", header = TRUE)
mod3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model3.bestlhoods.summary.txt", header = TRUE)
mod4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model4.bestlhoods.summary.txt", header = TRUE)
mod5 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model5.bestlhoods.summary.txt", header = TRUE)
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.bestlhoods.summary.txt", header = TRUE)

# Now find the ML run for each model
# Model 1, no bottleneck
mod1_ml <- max(mod1$MaxEstLhood)
mod1[which(mod1$MaxEstLhood == max(mod1$MaxEstLhood)),]
mod1[19,]$MaxEstLhood - mod1[19,]$MaxObsLhood
mod1_param_no <- 1 #should be hard-coded because complex parameters are sometimes reported. Number of parameters can be found in .est files

# Model 2, instantaneous bottleneck
mod2_ml <- max(mod2$MaxEstLhood)
mod2[which(mod2$MaxEstLhood == max(mod2$MaxEstLhood)),]
mod2[21,]$MaxEstLhood - mod2[21,]$MaxObsLhood # difference between MaxEst and MaxObs
mod2_param_no <- 5

# Model 3, exponential change after bottleneck
mod3_ml <- max(mod3$MaxEstLhood)
mod3[which(mod3$MaxEstLhood == max(mod3$MaxEstLhood)),]
mod3[34,]$MaxEstLhood - mod3[34,]$MaxObsLhood # difference between MaxEst and MaxObs
mod3_param_no <- 5

# Model 4, exponential growth, then bottleneck, then instantaneous recovery
mod4_ml <- max(mod4$MaxEstLhood)
mod4[which(mod4$MaxEstLhood == max(mod4$MaxEstLhood)),]
mod4[41,]$MaxEstLhood - mod4[41,]$MaxObsLhood # difference between MaxEst and MaxObs
mod4_param_no <- 6

# Model 5, two bottlenecks
mod5_ml <- max(mod5$MaxEstLhood)
mod5[which(mod5$MaxEstLhood == max(mod5$MaxEstLhood)),]
mod5[37,]$MaxEstLhood - mod5[37,]$MaxObsLhood # difference between MaxEst and MaxObs
mod5_param_no <- 9

# Model 6, exponential change in pop size before and after bottleneck
mod6_ml <- max(mod6$MaxEstLhood)
mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]
mod6[2,]$MaxEstLhood - mod6[2,]$MaxObsLhood
mod6_param_no <- 6

#### Plot ML for parameters over fsc iterations ####
plot(mod6$MaxEstLhood, ylab = 'Maximum likelihood', xlab = 'Iteration')
lines(mod6$MaxEstLhood)
points(2,-3703.332, col = 'tomato', pch = 19)

plot(mod6)

mod6$ML_diff <- mod6$MaxEstLhood - (-3703.332)
plot(mod6$ML_diff)
lines(mod6$ML_diff)
barplot(mod6$ML_diff, xlab = 'Iteration', ylab = "Difference from ML")

#### AIC calculations ####
# a test
# a <- c(-918.395, -687.045, -740.019, -782.598	) # MaxEstLhood
# b <- c(4, 6, 7, 7) # number of estimated parameters

# Data from my models
a <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml) # MaxEstLhood. These are log10 likelihoods
aa <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml)*2.303 # Convert from log10 to ln
b <- c(mod1_param_no, mod2_param_no, mod3_param_no, mod4_param_no, mod5_param_no, mod6_param_no) # number of estimated parameters

aic <- 2*b-2*aa

delta_aic <-round(aic,0) - min(round(aic,0))

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

# for (i in 1:length(aic)) {
#   w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))))
# }

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))), exp(-0.5*(aic[5]-max(aic))), exp(-0.5*(aic[6]-max(aic))))
}

#### Plot all fsc runs from nonparametric bootstrapping, plus the best fit model ####
# Read in ML and nonparametric bootstrapped parameters from exponential growth of NANC, then bottleneck & exponential recovery model (Model 6)
boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.nonparametric.ci.summary.txt', header = TRUE) # Read in ML nonparametric bootstrapped parameters following 10 runs for each simulated SFS
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.bestlhoods.summary.txt", header = TRUE) # Read in ML parameters from Model 6

# Find maximum from best model
# Exponential growth of NANC, then bottleneck & exponential recovery (Model 6)
max(mod6$MaxEstLhood)
mod6_ml_parameters <- mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]

#### Examine parameter point estimates and CIs ####
# Histograms of parameters estimated from simulated SFS, plus a line for the point estimates
hist(log10(boot$NPOP08), main = 'Population size in 2008', xlab = 'log10(NPOP08)')
abline(v = log10(mod6_ml_parameters$NPOP08), col = 'tomato')
hist(log10(boot$NPREBOT), main = 'Pre-bottleneck population size', xlab = 'log10(NPREBOT)')
abline(v=log10(mod6_ml_parameters$NPREBOT), col = 'tomato')
hist(log10(boot$NBOT), main = 'Population size during bottleneck', xlab = 'log10(NBOT)')
abline(v=log10(mod6_ml_parameters$NBOT), col = 'tomato')
hist(boot$TBOT, main = 'Time after bottleneck', xlab = 'Generations (2 years/generation)')
abline(v=mod6_ml_parameters$TBOT, col = 'tomato')
hist(boot$TLEN, main = 'Length of bottleneck', xlab = 'Generations (2 years/generation)', breaks = 4)
abline(v=mod6_ml_parameters$TLEN, col = 'tomato')
hist(boot$RANC, main = 'NANC growth rate', xlab = 'Growth rate of NANC going back in time\n(Negative means positive growth)')
abline(v=mod6_ml_parameters$RANC, col = 'tomato')
hist(log10(boot$NANC), main = 'Population size prior to pre-bottleneck growth', xlab = 'log10(NANC))')
abline(v=log10(mod6_ml_parameters$NANC), col = 'tomato')
hist(boot$MaxEstLhood)
hist(boot$MaxObsLhood)
hist(log10(boot$NBOT/boot$NPOP08), xlab = 'log10(NBOT/NPOP08)', main = 'NBOT/NPOP08 ratio') # same as hist(log10(boot$NBOT)-log10(boot$NPOP08))
abline(v = log10(mod6_ml_parameters$NBOT/mod6_ml_parameters$NPOP08), col = 'tomato')
hist(log10(boot$NPREBOT/boot$NBOT), xlab = 'log10(NPREBOT/NBOT)', main = 'NPREBOT/NBOT ratio')
abline(v = log10(mod6_ml_parameters$NPREBOT/mod6_ml_parameters$NBOT), col = 'tomato')

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
boxplot(log10(boot$NPREBOT/2), log10(boot$NBOT/2), log10(boot$NPOP08/2), range = 0, ylab = expression('log'[10]*'('*italic(N[e])* ' estimate)'), names = c('NPREBOT', 'NBOT', 'NPOP08')) # convert Ne to diploid by dividing haploid number by 2; italic Ne

points(c(1,2,3), c(log10(mod6_ml_parameters$NPREBOT/2), log10(mod6_ml_parameters$NBOT/2), log10(mod6_ml_parameters$NPOP08/2)), col = 'tomato', pch = 19)

dev.off()

# Calculate 95% confidence intervals
# NPOP08
# mean(cis$NPOP08) + 1.960*(sd(cis$NPOP08)/sqrt(100))
# mean(cis$NPOP08) - 1.960*(sd(cis$NPOP08)/sqrt(100))
quantile(boot$NPOP08/2, c(0.025, 0.975)) # diploid

# NPREBOT
quantile(boot$NPREBOT/2, c(0.025, 0.975)) # diploid

# NBOT
quantile(boot$NBOT/2, c(0.025, 0.975)) # diploid

# TBOT
quantile(boot$TBOT, c(0.025, 0.975))

# TLEN
quantile(boot$TLEN, c(0.025, 0.975))

# RANC
quantile(boot$RANC, c(0.025, 0.975))

# NANC
quantile(boot$NANC/2, c(0.025, 0.975)) # diploid

# R2008
quantile(r, c(0.025, 0.975))

# NPOP08/NBOT ratio
mod6_ml_parameters$NPOP08/mod6_ml_parameters$NBOT # point estimate, same regardless of haploid or diploid
quantile((boot$NPOP08/boot$NBOT), c(0.025, 0.975))

# NBOT/NPREBOT ratio
mod6_ml_parameters$NBOT/mod6_ml_parameters$NPREBOT # point estimate, same regardless of haploid or diploid
quantile((boot$NBOT/boot$NPREBOT), c(0.025, 0.975))

# NPOP08/NPREBOT ratio
mod6_ml_parameters$NPOP08/mod6_ml_parameters$NPREBOT # point estimate, same regardless of haploid or diploid
quantile((boot$NPOP08/boot$NPREBOT), c(0.025, 0.975))

#### Manipulate ML and bootstrapped data so that it can be plotted over time ####
# Model 6 & convert haploid numbers to diploid by dividing by 2; check RANC by calculating log(N_0/N_t)/(t_0/t_t)
# Calculate recent growth rate R2008 from population size for later
r2008 <- (log(mod6_ml_parameters$NBOT/mod6_ml_parameters$NPOP08)/(mod6_ml_parameters$TBOT)) # check that this is correct; growth rate going back in time

max <- data.frame(matrix(NA, nrow = 12, ncol = 4))
max[,1] <- c(0, mod6_ml_parameters$TBOT-8, mod6_ml_parameters$TBOT-6, mod6_ml_parameters$TBOT-4, mod6_ml_parameters$TBOT-3, mod6_ml_parameters$TBOT-1, mod6_ml_parameters$TBOT, (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+2), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+5), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
max[,2] <- 2008 - 2*max[,1]
max[,3] <- c(mod6_ml_parameters$NPOP08, (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-8))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-6))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-4))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-3))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-1))), mod6_ml_parameters$NBOT, mod6_ml_parameters$NBOT, mod6_ml_parameters$NPREBOT, (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 2)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 5)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 1000))) #haploid; when t is replaced by 9000 generations, NANC = 2877
max[,4] <- max[,3]/2 #diploid

# Nonparametric bootstrapped data
# First calculate recent r (R2008) for each SFS
r <- vector(length = 100)
for (i in 1:length(r)) {
  r[i] <- (log(boot$NBOT[i]/boot$NPOP08[i])/(boot$TBOT[i])) # check a few to make sure they're correct
}

# Fewer discrete points, curve is rougher, but avoids plotting beyond 2008
boot.max <- array(numeric(), c(10,4,100))
for (i in 1:nrow(boot)) {
  boot.max[,1,i] <- c(0, boot$TBOT[i]-3, boot$TBOT[i]-2, boot$TBOT[i]-1, boot$TBOT[i], (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]+2), (boot$TBOT[i]+boot$TLEN[i]+5), (boot$TBOT[i]+boot$TLEN[i]+1000)) #generations going back in time
  boot.max[,2,i] <- 2008 -2*boot.max[,1,i] #convert to years assuming summer flounder generation time is 2 years
  boot.max[,3,i] <- c(boot$NPOP08[i], (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-3))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-2))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-1))), boot$NBOT[i], boot$NBOT[i], boot$NPREBOT[i], (boot$NPREBOT[i]*exp(boot$RANC[i] * 2)), (boot$NPREBOT[i]*exp(boot$RANC[i] * 5)), (boot$NPREBOT[i]*exp(boot$RANC[i] * 1000))) #haploid
  boot.max[,4,i] <- boot.max[,3,i]/2 #diploid
}

# More discrete points for smoother exponential curve
# boot.max <- array(numeric(), c(14,4,100))
# for (i in 1:nrow(boot)) {
#   boot.max[,1,i] <- c(0, boot$TBOT[i]-12, boot$TBOT[i]-10, boot$TBOT[i]-8, boot$TBOT[i]-6, boot$TBOT[i]-4, boot$TBOT[i]-2, boot$TBOT[i]-1, boot$TBOT[i], (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]+2), (boot$TBOT[i]+boot$TLEN[i]+5), (boot$TBOT[i]+boot$TLEN[i]+1000)) #generations going back in time
#   boot.max[,2,i] <- 2008 -2*boot.max[,1,i] #convert to years assuming summer flounder generation time is 2 years
#   boot.max[,3,i] <- c(boot$NPOP08[i], (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-12))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-10))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-8))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-6))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-4))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-2))), (boot$NPOP08[i]*exp(r[i] * (boot$TBOT[i]-1))), boot$NBOT[i], boot$NBOT[i], boot$NPREBOT[i], (boot$NPREBOT[i]*exp(boot$RANC[i] * 2)), (boot$NPREBOT[i]*exp(boot$RANC[i] * 5)), (boot$NPREBOT[i]*exp(boot$RANC[i] * 1000))) #haploid
#   boot.max[,4,i] <- boot.max[,3,i]/2 #diploid
# }


# Set up plot and plot bootstrapped data
options(scipen = 5)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model6_lineplot_a_and_b.png",width=6, height=10, res=300, units="in")
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model6_lineplot.png",width=6, height=5, res=300, units="in")
par(mar=c(4, 6.5, 1.1, 1)+0.1, # panel margin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12,
    mfrow = c(2,1)
)

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid numbers
cols <- adjustcolor('gray70', alpha.f = 0.5)
plot(max$X2, max$X4, xlab = '', ylab = '', type = 'n', xlim = c(1980,2008), ylim = c(0,65000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)
mtext('(a)', 2,4.5, cex = 1.4, las = 1, at = 65000)

# for (l in 1:100) {
#   lines(boot.max[,,l], col = 'gray90')
# }

# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4, lwd = 1.8)

# dev.off()

# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model6_lineplot_deeptime.png",width=6, height=5, res=300, units="in")
# par(mar=c(4, 6, 1, 1), # panel margin size in "line number" units
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=12
# )

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid
plot(max$X2, max$X4, xlab = '', ylab = '', type = 'n', xlim = c(1500,2008), ylim = c(0,64000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)
mtext('(b)', 2,4.5, cex = 1.4, las = 1, at = 65000)


# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4, lwd = 1.8)

dev.off()


#### Using a SFS that summarizes 276 larvae across three cohorts, 1196 loci and no MAF filter ####
# Read in the estimated parameters for each model
mod1 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model1.bestlhoods.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model2.bestlhoods.txt", header = TRUE)
mod3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model3.bestlhoods.txt", header = TRUE)
mod4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model4.bestlhoods.txt", header = TRUE)
mod5 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model5.bestlhoods.txt", header = TRUE)
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nosibs_models/model6.bestlhoods.txt", header = TRUE)

# Now find the ML run for each model
mod1_ml <- max(mod1$MaxEstLhood)
mod2_ml <- max(mod2$MaxEstLhood)
mod3_ml <- max(mod3$MaxEstLhood)
mod4_ml <- max(mod4$MaxEstLhood)
mod5_ml <- max(mod5$MaxEstLhood)
mod6_ml <- max(mod6$MaxEstLhood)

# AIC calculation
a <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml) # MaxEstLhood. These are log10 likelihoods
aa <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml)*2.303 # Convert from log10 to ln
b <- c(1,5,5,6,9,6) # number of estimated parameters

aic <- 2*b-2*aa

delta_aic <-round(aic,0) - min(round(aic,0))

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))), exp(-0.5*(aic[5]-max(aic))), exp(-0.5*(aic[6]-max(aic))))
}

##########################################################################################################
#### Using a SFS that summarizes 284 larvae across three cohorts, 1084 loci and no MAF or MAC filters ####
# Read in the estimated parameters for each model
mod1 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model1.bestlhoods.txt", header = TRUE) # constant population size for comparison to all the bottlenecks, but only estimating NPOP08
mod2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model2.bestlhoods.txt", header = TRUE)
mod3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model3.bestlhoods.txt", header = TRUE)
mod4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model4.bestlhoods.txt", header = TRUE)
mod5 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model5.bestlhoods.txt", header = TRUE)
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model6.bestlhoods.txt", header = TRUE)

# Now find the ML run for each model
mod1_ml <- max(mod1$MaxEstLhood)
mod2_ml <- max(mod2$MaxEstLhood)
mod3_ml <- max(mod3$MaxEstLhood)
mod4_ml <- max(mod4$MaxEstLhood)
mod5_ml <- max(mod5$MaxEstLhood)
mod6_ml <- max(mod6$MaxEstLhood)

# AIC calculation
a <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml) # MaxEstLhood. These are log10 likelihoods
aa <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml, mod6_ml)*2.303 # Convert from log10 to ln
b <- c(1,5,5,6,9,6) # number of estimated parameters

aic <- 2*b-2*aa

delta_aic <-round(aic,0) - min(round(aic,0))

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))), exp(-0.5*(aic[5]-max(aic))), exp(-0.5*(aic[6]-max(aic))))
}

#### Plot all fsc runs from nonparametric bootstrapping, plus the best fit model ####
# Read in ML and nonparametric bootstrapped parameters from exponential growth of NANC, then bottleneck & exponential recovery model (Model 6)
boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model6.nonparametric.ci.summary.txt', header = TRUE) # Read in ML nonparametric bootstrapped parameters following 30 runs for each simulated SFS
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/nomac_models/model6.bestlhoods.txt", header = TRUE) # Read in ML parameters from Model 6

# Find maximum from best model
max(mod6$MaxEstLhood)
mod6_ml_parameters <- mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]

#### Examine parameter point estimates and CIs ####
# Histograms of parameters estimated from simulated SFS, plus a line for the point estimates
hist(log10(boot$NPOP08), main = 'Population size in 2008', xlab = 'log10(NPOP08)')
abline(v = log10(mod6_ml_parameters$NPOP08), col = 'tomato')
hist(log10(boot$NPREBOT), main = 'Pre-bottleneck population size', xlab = 'log10(NPREBOT)')
abline(v=log10(mod6_ml_parameters$NPREBOT), col = 'tomato')
hist(log10(boot$NBOT), main = 'Population size during bottleneck', xlab = 'log10(NBOT)')
abline(v=log10(mod6_ml_parameters$NBOT), col = 'tomato')
hist(boot$TBOT, main = 'Time after bottleneck', xlab = 'Generations (2 years/generation)')
abline(v=mod6_ml_parameters$TBOT, col = 'tomato')
hist(boot$TLEN, main = 'Length of bottleneck', xlab = 'Generations (2 years/generation)', breaks = 4)
abline(v=mod6_ml_parameters$TLEN, col = 'tomato')
hist(boot$RANC, main = 'NANC growth rate', xlab = 'Growth rate of NANC going back in time\n(Negative means positive growth)')
abline(v=mod6_ml_parameters$RANC, col = 'tomato')
hist(log10(boot$NANC), main = 'Population size prior to pre-bottleneck growth', xlab = 'log10(NANC))')
abline(v=log10(mod6_ml_parameters$NANC), col = 'tomato')
hist(log10(boot$NBOT/boot$NPOP08), xlab = 'log10(NBOT/NPOP08)', main = 'NBOT/NPOP08 ratio') # same as hist(log10(boot$NBOT)-log10(boot$NPOP08))
abline(v = log10(mod6_ml_parameters$NBOT/mod6_ml_parameters$NPOP08), col = 'tomato')
hist(log10(boot$NPREBOT/boot$NBOT), xlab = 'log10(NPREBOT/NBOT)', main = 'NPREBOT/NBOT ratio')
abline(v = log10(mod6_ml_parameters$NPREBOT/mod6_ml_parameters$NBOT), col = 'tomato')

# Calculate 95% confidence intervals
# NPOP08
quantile(boot$NPOP08/2, c(0.025, 0.975)) # diploid

# NPREBOT
quantile(boot$NPREBOT/2, c(0.025, 0.975)) # diploid

# NBOT
quantile(boot$NBOT/2, c(0.025, 0.975)) # diploid

# TBOT
quantile(boot$TBOT, c(0.025, 0.975))

# TLEN
quantile(boot$TLEN, c(0.025, 0.975))

# NANC
quantile(boot$NANC/2, c(0.025, 0.975)) # diploid

# NBOT/NPOP08
quantile((boot$NBOT/2)/(boot$NPOP08/2), c(0.025, 0.975))

# PREBOT/NBOT
quantile((boot$NPREBOT/2)/(boot$NBOT/2), c(0.025, 0.975))

# PREBOT/NPOP08
quantile((boot$NPREBOT/2)/(boot$NPOP08/2), c(0.025, 0.975))

