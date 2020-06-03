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
# growing.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/growing_nomaf/best.lhood.summary.growing_nomaf_Nlowerlimit100.txt", header = TRUE)
# shrinking.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/shrinking_nomaf/best.lhood.summary.shrinking_nomaf.txt", header = TRUE)
# nobot.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/nobot_nomaf/best.lhood.summary.nobot_nomaf.txt", header = TRUE) # constant population size for comparision to all the bottlenecks
nobot.best.nomaf2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/nobot_nomaf_oneparam/best.lhood.summary.nobot_nomaf.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
#tbot.fixed.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/instant_nomaf_fixedTBOT/best.lhood.summary.stable_nomaf_fixedTBOT.txt", header = TRUE)
expgrowth_instant <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.expgrowth_instant_nomaf.txt", header = TRUE)
two_bots <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/two_bots/best.lhood.summary.two_bots_nomaf.txt", header = TRUE)
exp_change <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/exp_change_nomaf/best.lhood.summary.exp_change_nomaf.txt", header = TRUE)

# Stable
stable_ml <- max(stable.best.nomaf$MaxEstLhood)
stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]
stable.best.nomaf[23,]$MaxEstLhood - stable.best.nomaf[23,]$MaxObsLhood # difference between MaxEst and MaxObs
stable_param_no <- ncol(stable.best.nomaf)-2

# Growing
# max(growing.best.nomaf$MaxEstLhood)
# growing.best.nomaf[which(growing.best.nomaf$MaxEstLhood == max(growing.best.nomaf$MaxEstLhood)),]
# growing.best.nomaf[16,]$MaxEstLhood - growing.best.nomaf[16,]$MaxObsLhood # difference between MaxEst and MaxObs

# Shrinking
# max(shrinking.best.nomaf$MaxEstLhood)
# shrinking.best.nomaf[which(shrinking.best.nomaf$MaxEstLhood == max(shrinking.best.nomaf$MaxEstLhood)),]
# shrinking.best.nomaf[9,]$MaxEstLhood - shrinking.best.nomaf[9,]$MaxObsLhood # difference between MaxEst and MaxObs

# Constant
# max(nobot.best.nomaf$MaxEstLhood)
# nobot.best.nomaf[which(nobot.best.nomaf$MaxEstLhood == max(nobot.best.nomaf$MaxEstLhood)),]
# nobot.best.nomaf[28,]$MaxEstLhood - nobot.best.nomaf[28,]$MaxObsLhood # difference between MaxEst and MaxObs

# Constant, one parameter only
nobot_ml <- max(nobot.best.nomaf2$MaxEstLhood)
nobot.best.nomaf2[which(nobot.best.nomaf2$MaxEstLhood == max(nobot.best.nomaf2$MaxEstLhood)),]
nobot.best.nomaf2[5,]$MaxEstLhood - nobot.best.nomaf2[5,]$MaxObsLhood
nobot_param_no <- ncol(nobot.best.nomaf2)-2

# Fixed TBOT to 7
# max(tbot.fixed.nomaf$MaxEstLhood)
# tbot.fixed.nomaf[which(tbot.fixed.nomaf$MaxEstLhood == max(tbot.fixed.nomaf$MaxEstLhood)),]
# tbot.fixed.nomaf[5,]$MaxEstLhood - tbot.fixed.nomaf[5,]$MaxObsLhood

# Exponential growth, then bottleneck, then instantaneous recovery
expgrowth_instant_ml <- max(expgrowth_instant$MaxEstLhood)
expgrowth_instant[which(expgrowth_instant$MaxEstLhood == max(expgrowth_instant$MaxEstLhood)),]
expgrowth_instant[43,]$MaxEstLhood - expgrowth_instant[43,]$MaxObsLhood # difference between MaxEst and MaxObs
expgrowth_instant_param_no <- ncol(expgrowth_instant)-2

# Two bottlenecks
two_bots_ml <- max(two_bots$MaxEstLhood)
two_bots[which(two_bots$MaxEstLhood == max(two_bots$MaxEstLhood)),]
two_bots[21,]$MaxEstLhood - two_bots[21,]$MaxObsLhood # difference between MaxEst and MaxObs
two_bots_param_no <- ncol(two_bots)-2

# Exponential change after bottleneck
exp_change_ml <- max(exp_change$MaxEstLhood)
exp_change[which(exp_change$MaxEstLhood == max(exp_change$MaxEstLhood)),]
exp_change[34,]$MaxEstLhood - exp_change[34,]$MaxObsLhood # difference between MaxEst and MaxObs
exp_change_param_no <- ncol(exp_change)-2

#### AIC calculations ####
# a test
# a <- c(-918.395, -687.045, -740.019, -782.598	) # MaxEstLhood
# b <- c(4, 6, 7, 7) # number of estimated parameters

# Data from my models
a <- c(nobot_ml, stable_ml, exp_change_ml, expgrowth_instant_ml, two_bots_ml) # MaxEstLhood
b <- c(nobot_param_no, stable_param_no, exp_change_param_no, expgrowth_instant_param_no, two_bots_param_no) # number of estimated parameters

aic <- 2*b-2*a

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

# for (i in 1:length(aic)) {
#   w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))))
# }

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))), exp(-0.5*(aic[5]-max(aic))))
}

#### Plot all fsc runs from parametric bootstrapping, plus the best fit model ####
# Read in ML and bootstrapped parameters from instananeous recovery model
# boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/max.summary.txt', header = TRUE) # Read in ML bootstrapped parameters following 50 runs for each simulated SFS
# stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.stable_nomaf.txt", header = TRUE)# Read in ML parameters from instantaneous recovery model

# Read in ML and bootstrapped parameters from exponential growth of NANC, then bottleneck & instananeous recovery model
boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/cis_sfs_summary.txt', header = TRUE) # Read in ML bootstrapped parameters following 50 runs for each simulated SFS
expgrowth_instant.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.expgrowth_instant_nomaf.txt", header = TRUE)# Read in ML parameters from instantaneous recovery model

# Find maximum from best models
# Instantaneous recovery
# max(stable.best.nomaf$MaxEstLhood)
# instant <- stable.best.nomaf[which(stable.best.nomaf$MaxEstLhood == max(stable.best.nomaf$MaxEstLhood)),]

# Exponential growth of NANC, then bottleneck & instantaneous recovery
max(expgrowth_instant.best.nomaf$MaxEstLhood)
expgrowth_instant <- expgrowth_instant.best.nomaf[which(expgrowth_instant.best.nomaf$MaxEstLhood == max(expgrowth_instant.best.nomaf$MaxEstLhood)),]

#### Manipulate data so that its plottable ####
# Instantaneous recovery model
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
options(scipen = 5)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/instant_recovery_lineplot.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Plots bootstrapped 50 runs used to estimate 95% CI
plot(max$X1, max$X2, xlab = 'Year', ylab = expression('N'[e]), type = 'n', ylim = c(0,100000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,1,l], factor = 0.3), boot.max[,2,l], col = 'gray90')
}

# for (l in 1:100) {
#   lines(boot.max[,,l], col = 'gray90')
# }

# Plots parameters from best fit model
lines(max$X1, max$X2)

dev.off()

# Exponenetial growth of NANC then bottleneck & instantaneous recovery model & convert haploid numbers to diploid by dividing by 2; check RANC by calculating log(N_0/N_t)/(t_0/t_t)
max <- data.frame(matrix(NA, nrow = 8, ncol = 4))
max[,1] <- c(0, expgrowth_instant$TBOT, expgrowth_instant$TBOT, (expgrowth_instant$TBOT+expgrowth_instant$TLEN), (expgrowth_instant$TBOT+expgrowth_instant$TLEN), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+2), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+5), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+1000)) #generations going back in time
# max[,2] <- c(2008, 2008-(2*expgrowth_instant$TBOT), 2008-(2*expgrowth_instant$TBOT), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN)), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN)), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-10, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-20, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-30) #convert to years assuming summer flounder generation time of 2 years
max[,2] <- 2008 - 2*max[,1]
max[,3] <- c(expgrowth_instant$NPOP08, expgrowth_instant$NPOP08, expgrowth_instant$NBOT, expgrowth_instant$NBOT, expgrowth_instant$NANC, expgrowth_instant$NANC/exp(expgrowth_instant$RANC * -2), expgrowth_instant$NANC/exp(expgrowth_instant$RANC * -5), expgrowth_instant$NANC/exp(expgrowth_instant$RANC * -1000)) #haploid
max[,4] <- max[,3]/2 #diploid

# Parametric bootstrapped data
boot.max <- array(numeric(), c(8,4,100))
for (i in 1:nrow(boot)) {
  boot.max[,1,i] <- c(0, boot$TBOT[i], boot$TBOT[i], (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]+2), (boot$TBOT[i]+boot$TLEN[i]+5), (boot$TBOT[i]+boot$TLEN[i]+1000)) #generations going back in time
  boot.max[,2,i] <- 2008 -2*boot.max[,1,i] #convert to years assuming summer flounder generation time is 2 years
  #boot.max[,2,i] <- c(2008, 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]), 2008-((2*boot$TBOT[i])+(2*boot$TLEN[i])), 2008-((2*boot$TBOT[i])+(2*boot$TLEN[i])), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-10, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-20, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-30)
  boot.max[,3,i] <- c(boot$NPOP08[i], boot$NPOP08[i], boot$NBOT[i], boot$NBOT[i], boot$NANC[i], boot$NANC[i]/exp(boot$RANC[i] * -2), boot$NANC[i]/exp(boot$RANC[i] * -5), boot$NANC[i]/exp(boot$RANC[i] * -1000)) #haploid
  boot.max[,4,i] <- boot.max[,3,i]/2 #diploid
}

# Set up plot and plot bootstrapped data
options(scipen = 5)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/expgrowth_instant_lineplot.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid numbers
cols <- adjustcolor('gray70', alpha.f = 0.5)
plot(max$X2, max$X4, xlab = 'Year', ylab = expression('N'[E]), type = 'n', xlim = c(1980,2008), ylim = c(0,70000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}

# for (l in 1:100) {
#   lines(boot.max[,,l], col = 'gray90')
# }

# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4)

dev.off()

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/expgrowth_instant_lineplot_deeptime.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid
plot(max$X2, max$X4, xlab = 'Year', ylab = expression('N'[E]), type = 'n', xlim = c(1500,2008), ylim = c(0,70000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}

# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4)

dev.off()
