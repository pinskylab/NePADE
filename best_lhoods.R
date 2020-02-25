#### Fsc26 has been run 50 times and the maximum likelihoods and associated parameters have been concatinated ####
# Using a SFS where a MAF < 0.05 was used
# Instantaneous recovery following a bottleneck
stable.best <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing.res/best.lhood.summary.stableinstant.txt", header = TRUE)

plot(stable.best)
plot(stable.best$NBOT ~ stable.best$TLEN, xlab = 'Length of bottleneck', ylab = 'Ne bottleneck')

max(stable.best$MaxEstLhood)
stable.best[which(stable.best$MaxEstLhood == max(stable.best$MaxEstLhood)),]

#####################################################################################################################
# Malin is worried that the length and mangnitude of the bottleneck are correlated --> holding TLEN = 2
stable_fixed.best <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixedTLEN/best.lhood.summary.stableinstant_fixedtlen.txt", header = TRUE)
plot(stable_fixed.best)

max(stable_fixed.best$MaxEstLhood)
stable_fixed.best[which(stable_fixed.best$MaxEstLhood == max(stable_fixed.best$MaxEstLhood)),]

# Malin is worried that the length and mangnitude of the bottleneck are correlated --> holding TBOT = 7
stable_fixed2.best <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixedTBOT/best.lhood.summary.stable_fixedTBOT.txt", header = TRUE)

max(stable_fixed2.best$MaxEstLhood)
stable_fixed2.best[which(stable_fixed2.best$MaxEstLhood == max(stable_fixed2.best$MaxEstLhood)),]

#####################################################################################################################
# Instantaneous recovery following a bottleneck using the SFS with no MAF
stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing_nomaf2.res/best.lhood.summary.stable_nomaf.txt", header = TRUE)

plot(stable.best.nomaf)
plot(stable.best.nomaf$NBOT ~ stable.best.nomaf$TLEN, xlab = 'Length of bottleneck', ylab = 'Ne bottleneck')

max(stable.best.nomaf$MaxEstLhood)
stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]

#######################################################################################################################################################################################################
# Needed to play around with the parameters, but here are the best ML for instantaneous recovery, exp growth and exp decay following a bottleneck. The lower limit for N is 100 (50 diploid individuals)
stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.stable_nomaf.txt", header = TRUE)
growing.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/growing_nomaf/best.lhood.summary.growing_nomaf_Nlowerlimit100.txt", header = TRUE)
shrinking.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/shrinking_nomaf/best.lhood.summary.shrinking_nomaf.txt", header = TRUE)
# nobot.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/nobot_nomaf/best.lhood.summary.nobot_nomaf.txt", header = TRUE) # constant population size for comparision to all the bottlenecks
nobot.best.nomaf2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/nobot_nomaf_oneparam/best.lhood.summary.nobot_nomaf.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
tbot.fixed.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/instant_nomaf_fixedTBOT/best.lhood.summary.stable_nomaf_fixedTBOT.txt", header = TRUE)

# Stable
max(stable.best.nomaf$MaxEstLhood)
stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]
stable.best.nomaf[23,]$MaxEstLhood - stable.best.nomaf[23,]$MaxObsLhood # difference between MaxEst and MaxObs

# Growing
max(growing.best.nomaf$MaxEstLhood)
growing.best.nomaf[which(growing.best.nomaf$MaxEstLhood == max(growing.best.nomaf$MaxEstLhood)),]
growing.best.nomaf[16,]$MaxEstLhood - growing.best.nomaf[16,]$MaxObsLhood # difference between MaxEst and MaxObs

# Shrinking
max(shrinking.best.nomaf$MaxEstLhood)
shrinking.best.nomaf[which(shrinking.best.nomaf$MaxEstLhood == max(shrinking.best.nomaf$MaxEstLhood)),]
shrinking.best.nomaf[9,]$MaxEstLhood - shrinking.best.nomaf[9,]$MaxObsLhood # difference between MaxEst and MaxObs

# Constant
# max(nobot.best.nomaf$MaxEstLhood)
# nobot.best.nomaf[which(nobot.best.nomaf$MaxEstLhood == max(nobot.best.nomaf$MaxEstLhood)),]
# nobot.best.nomaf[28,]$MaxEstLhood - nobot.best.nomaf[28,]$MaxObsLhood # difference between MaxEst and MaxObs

# Constant, one parameter only
max(nobot.best.nomaf2$MaxEstLhood)
nobot.best.nomaf2[which(nobot.best.nomaf2$MaxEstLhood == max(nobot.best.nomaf2$MaxEstLhood)),]
nobot.best.nomaf2[5,]$MaxEstLhood - nobot.best.nomaf2[5,]$MaxObsLhood

# Fixed TBOT to 7
max(tbot.fixed.nomaf$MaxEstLhood)
tbot.fixed.nomaf[which(tbot.fixed.nomaf$MaxEstLhood == max(tbot.fixed.nomaf$MaxEstLhood)),]
tbot.fixed.nomaf[5,]$MaxEstLhood - tbot.fixed.nomaf[5,]$MaxObsLhood

#### AIC calculations ####
# a test
# a <- c(-918.395, -687.045, -740.019, -782.598	) # MaxEstLhood
# b <- c(4, 6, 7, 7) # number of estimated parameters

# Data from my models
a <- c(-4016.987, -4021.143, -4038.039, -4054.514) # MaxEstLhood
b <- c(5, 5, 5, 1) # number of estimated parameters

aic <- 2*b-2*a

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

# for (i in 1:length(aic)) {
#   w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))))
# }

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))))
}

#### Plot all fsc runs from parametric bootstrapping, plus the best fit model ####
boot <- read.table() # Read in ML bootstrapped parameters
stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.stable_nomaf.txt", header = TRUE)# Read in ML parameters from best model

# Find maximum from best model
max(stable.best.nomaf$MaxEstLhood)
instant <- stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]

# Manipulate data so that its plottable
# Best fit model
max <- data.frame(matrix(NA, nrow = 6, ncol = 2))
max[,1] <- c(2008, 2008-(2*instant$TBOT), 2008-(2*instant$TBOT), 2008-(2*instant$TBOT+instant$TLEN), 2008-(2*instant$TBOT+instant$TLEN), 1980)
max[,2] <- c(instant$NPOP08, instant$NPOP08, instant$NBOT, instant$NBOT, instant$NANC, instant$NANC)

# Parametric bootstrapped data
boot.max <- array(numeric(), c(6,2,100))
for (i in 1:nrow(boot)) {
  boot.max[,1,i] <- c(2008, 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]+boot$TLEN[i]), 2008-(2*boot$TBOT[i]+boot$TLEN[i]), 1980)
  boot.max[,2,i] <- c(boot$NPOP08[i], boot$NPOP08[i], boot$NBOT[i], boot$NBOT[i], boot$NANC[i], boot$NANC[i])
}

# Set up plot and plot bootstrapped data
plot(max$X1, max$X2, xlab = 'Year', ylab = 'Effective population size', type = 'n', ylim = c(0,100000))
for (l in 1:100) {
  lines(jitter(boot.max[,1,l], factor = 0.3), boot.max[,2,l], col = 'gray90')
}

# for (l in 1:100) {
#   lines(boot.max[,,l], col = 'gray90')
# }

# parameters from best fit model
lines(max$X1, max$X2)

