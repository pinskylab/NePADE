#### Fsc26 has been run 50 times and the maximum likelihoods and associated parameters have been concatenated ####
# Using a SFS that summarizes 280 larvae across three cohorts, 1196 loci and no MAF filter
# Read in the estimated parameters for each model
mod1 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model1.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model2.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model3.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model4.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod5 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model5.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08

# Now find the ML run for each model
# Model 1, no bottleneck
mod1_ml <- max(mod1$MaxEstLhood)
mod1[which(mod1$MaxEstLhood == max(mod1$MaxEstLhood)),]
mod1[19,]$MaxEstLhood - mod1[19,]$MaxObsLhood
mod1_param_no <- ncol(mod1)-2

# Model 2, instantaneous bottleneck
mod2_ml <- max(mod2$MaxEstLhood)
mod2[which(mod2$MaxEstLhood == max(mod2$MaxEstLhood)),]
mod2[21,]$MaxEstLhood - mod2[21,]$MaxObsLhood # difference between MaxEst and MaxObs
mod2_param_no <- ncol(mod2)-2

# Model 3, exponential change after bottleneck
mod3_ml <- max(mod3$MaxEstLhood)
mod3[which(mod3$MaxEstLhood == max(mod3$MaxEstLhood)),]
mod3[8,]$MaxEstLhood - mod3[8,]$MaxObsLhood # difference between MaxEst and MaxObs
mod3_param_no <- ncol(mod3)-2

# Model 4, exponential growth, then bottleneck, then instantaneous recovery
mod4_ml <- max(mod4$MaxEstLhood)
mod4[which(mod4$MaxEstLhood == max(mod4$MaxEstLhood)),]
mod4[2,]$MaxEstLhood - mod4[2,]$MaxObsLhood # difference between MaxEst and MaxObs
mod4_param_no <- ncol(mod4)-2

# Model 5, two bottlenecks
mod5_ml <- max(mod5$MaxEstLhood)
mod5[which(mod5$MaxEstLhood == max(mod5$MaxEstLhood)),]
mod5[37,]$MaxEstLhood - mod5[37,]$MaxObsLhood # difference between MaxEst and MaxObs
mod5_param_no <- ncol(mod5)-2

#### Plot ML for parameters over fsc iterations ####
plot(mod4$MaxEstLhood, ylab = 'Maximum likelihood', xlab = 'Iteration')
lines(mod4$MaxEstLhood) # there's a local maximum when RANC is positive
points(2,-3713.266, col = 'tomato', pch = 19)

plot(mod4)

mod4$ML_diff <- mod4$MaxEstLhood - (-3713.266)
plot(mod4$ML_diff)
lines(mod4$ML_diff)
barplot(mod4$ML_diff, xlab = 'Iteration', ylab = "Difference from ML")

#### AIC calculations ####
# a test
# a <- c(-918.395, -687.045, -740.019, -782.598	) # MaxEstLhood
# b <- c(4, 6, 7, 7) # number of estimated parameters

# Data from my models
a <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml) # MaxEstLhood. These are log10 likelihoods
aa <- c(mod1_ml, mod2_ml, mod3_ml, mod4_ml, mod5_ml)*2.303 # Convert from log10 to ln
b <- c(mod1_param_no, mod2_param_no, mod3_param_no, mod4_param_no, mod5_param_no) # number of estimated parameters

aic <- 2*b-2*aa

delta_aic <-aic - min(aic)

# Akaike weight = provides relative weight of evidence for each model. Probability that model i is the best model for the observed data, given the candidate set of models
w <- vector()

# for (i in 1:length(aic)) {
#   w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))))
# }

for (i in 1:length(aic)) {
  w[i] <- (exp(-0.5*(aic[i]-max(aic))))/sum(exp(-0.5*(aic[1]-max(aic))), exp(-0.5*(aic[2]-max(aic))), exp(-0.5*(aic[3]-max(aic))), exp(-0.5*(aic[4]-max(aic))), exp(-0.5*(aic[5]-max(aic))))
}

#### Plot all fsc runs from nonparametric bootstrapping, plus the best fit model ####
# Read in ML and bootstrapped parameters from instantaneous recovery model
# boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/max.summary.txt', header = TRUE) # Read in ML bootstrapped parameters following 50 runs for each simulated SFS
# stable.best.nomaf <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/stable_instant_nomaf_Nlowerlimit100/best.lhood.summary.stable_nomaf.txt", header = TRUE)# Read in ML parameters from instantaneous recovery model

# Read in ML and nonparametric bootstrapped parameters from exponential growth of NANC, then bottleneck & instananeous recovery model
boot <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/from_amarel/expgrowth_then_stable_instant_nomaf_Nlowerlimit100/nonparametric_ci_summary.txt', header = TRUE) # Read in ML nonparametric bootstrapped parameters following 10 runs for each simulated SFS
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
# max <- data.frame(matrix(NA, nrow = 6, ncol = 2))
# max[,1] <- c(2008, 2008-(2*instant$TBOT), 2008-(2*instant$TBOT), 2008-(2*instant$TBOT+instant$TLEN), 2008-(2*instant$TBOT+instant$TLEN), 1980)
# max[,2] <- c(instant$NPOP08, instant$NPOP08, instant$NBOT, instant$NBOT, instant$NANC, instant$NANC)
# 
# # Parametric bootstrapped data
# boot.max <- array(numeric(), c(6,2,100))
# for (i in 1:nrow(boot)) {
#   boot.max[,1,i] <- c(2008, 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]+boot$TLEN[i]), 2008-(2*boot$TBOT[i]+boot$TLEN[i]), 1980)
#   boot.max[,2,i] <- c(boot$NPOP08[i], boot$NPOP08[i], boot$NBOT[i], boot$NBOT[i], boot$NANC[i], boot$NANC[i])
# }
# 
# # Set up plot and plot bootstrapped data
# options(scipen = 5)
# 
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/instant_recovery_lineplot.png",width=6, height=5, res=300, units="in")
# par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=12
# )
# 
# # Plots bootstrapped 50 runs used to estimate 95% CI
# plot(max$X1, max$X2, xlab = 'Year', ylab = expression('N'[e]), type = 'n', ylim = c(0,100000), las = 1)
# for (l in 1:100) {
#   lines(jitter(boot.max[,1,l], factor = 0.3), boot.max[,2,l], col = 'gray90')
# }
# 
# # for (l in 1:100) {
# #   lines(boot.max[,,l], col = 'gray90')
# # }
# 
# # Plots parameters from best fit model
# lines(max$X1, max$X2)
# 
# dev.off()

# Exponential growth of NANC then bottleneck & instantaneous recovery model & convert haploid numbers to diploid by dividing by 2; check RANC by calculating log(N_0/N_t)/(t_0/t_t)
max <- data.frame(matrix(NA, nrow = 8, ncol = 4))
max[,1] <- c(0, expgrowth_instant$TBOT, expgrowth_instant$TBOT, (expgrowth_instant$TBOT+expgrowth_instant$TLEN), (expgrowth_instant$TBOT+expgrowth_instant$TLEN), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+2), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+5), (expgrowth_instant$TBOT+expgrowth_instant$TLEN+1000)) #generations going back in time
# max[,2] <- c(2008, 2008-(2*expgrowth_instant$TBOT), 2008-(2*expgrowth_instant$TBOT), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN)), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN)), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-10, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-20, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-30) #convert to years assuming summer flounder generation time of 2 years
max[,2] <- 2008 - 2*max[,1]
max[,3] <- c(expgrowth_instant$NPOP08, expgrowth_instant$NPOP08, expgrowth_instant$NBOT, expgrowth_instant$NBOT, expgrowth_instant$NANC, (expgrowth_instant$NANC*exp(expgrowth_instant$RANC * 2)), (expgrowth_instant$NANC*exp(expgrowth_instant$RANC * 5)), (expgrowth_instant$NANC*exp(expgrowth_instant$RANC * 1000))) #haploid
max[,4] <- max[,3]/2 #diploid

# Parametric bootstrapped data
boot.max <- array(numeric(), c(8,4,100))
for (i in 1:nrow(boot)) {
  boot.max[,1,i] <- c(0, boot$TBOT[i], boot$TBOT[i], (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]), (boot$TBOT[i]+boot$TLEN[i]+2), (boot$TBOT[i]+boot$TLEN[i]+5), (boot$TBOT[i]+boot$TLEN[i]+1000)) #generations going back in time
  boot.max[,2,i] <- 2008 -2*boot.max[,1,i] #convert to years assuming summer flounder generation time is 2 years
  #boot.max[,2,i] <- c(2008, 2008-(2*boot$TBOT[i]), 2008-(2*boot$TBOT[i]), 2008-((2*boot$TBOT[i])+(2*boot$TLEN[i])), 2008-((2*boot$TBOT[i])+(2*boot$TLEN[i])), 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-10, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-20, 2008-((2*expgrowth_instant$TBOT)+(2*expgrowth_instant$TLEN))-30)
  boot.max[,3,i] <- c(boot$NPOP08[i], boot$NPOP08[i], boot$NBOT[i], boot$NBOT[i], boot$NANC[i], (boot$NANC[i]*exp(boot$RANC[i] * 2)), (boot$NANC[i]*exp(boot$RANC[i] * 5)), (boot$NANC[i]*exp(boot$RANC[i] * 1000))) #haploid
  boot.max[,4,i] <- boot.max[,3,i]/2 #diploid
}

# Set up plot and plot bootstrapped data
options(scipen = 5)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/expgrowth_instant_lineplot.png",width=6, height=5, res=300, units="in")
par(mar=c(4, 6, 1, 1), # panel margin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid numbers
cols <- adjustcolor('gray70', alpha.f = 0.5)
plot(max$X2, max$X4, xlab = '', ylab = '', type = 'n', xlim = c(1980,2008), ylim = c(0,70000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression('N'[E]), 2, 3.7, cex = 1.2)

# for (l in 1:100) {
#   lines(boot.max[,,l], col = 'gray90')
# }

# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4, lwd = 1.8)

dev.off()

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/expgrowth_instant_lineplot_deeptime.png",width=6, height=5, res=300, units="in")
par(mar=c(4, 6, 1, 1), # panel margin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Plots bootstrapped 50 runs used to estimate 95% CI. Specify haploid or diploid
plot(max$X2, max$X4, xlab = '', ylab = '', type = 'n', xlim = c(1500,2008), ylim = c(0,72000), las = 1)
for (l in 1:100) {
  lines(jitter(boot.max[,2,l], factor = 0.2), boot.max[,4,l], col = cols)
}
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression('N'[E]), 2, 3.7, cex = 1.2)

# Plots parameters from best fit model. Specify haploid or diploid
lines(max$X2, max$X4, lwd = 1.8)

dev.off()
