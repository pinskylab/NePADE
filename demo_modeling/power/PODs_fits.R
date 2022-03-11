#### This script plots the best-fit model for Models 1-7 as lines. The best fit models were then used to generate 10 PODs for each model, and Models 1-7 were fit to each POD. Then line plots of the best-fit model for each POD were plotted ####
# 1. Reads in the parameters for Models 1-7 resulting from 50 fastsimcoal runs
# 2. Determines best fit model and creates a data frame of points to plot as a line
# 3. Reads in the parameters resulting from fitting Models 1-7 to the PODS for each Model
# 4. Determines the best-fit model for each POD
# 5. Creates an array of points to plot each best-fit model to a POD as a line
# 6. Plots points from steps 2 and 5 to create a 7 panel figure, with one panel for each Model

# Read in fsc data to select ML parameters for each model and plot
# Read in the estimated parameters for each model
mod1 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model1.bestlhoods.summary.txt", header = TRUE) # constant population size for comparision to all the bottlenecks, but only estimating NPOP08
mod2 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model2.bestlhoods.summary.txt", header = TRUE)
mod3 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model3.bestlhoods.summary.txt", header = TRUE)
mod4 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model4.bestlhoods.summary.txt", header = TRUE)
mod5 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model5.bestlhoods.summary.txt", header = TRUE)
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.bestlhoods.summary.txt", header = TRUE)
mod7 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model7.bestlhoods.summary.txt", header = TRUE)

# Now find the ML run for each model
# Model 1, no bottleneck
mod1_ml <- max(mod1$MaxEstLhood)
mod1_ml_parameters <- mod1[which(mod1$MaxEstLhood == max(mod1$MaxEstLhood)),]

# Model 2, instantaneous bottleneck
mod2_ml <- max(mod2$MaxEstLhood)
mod2_ml_parameters <- mod2[which(mod2$MaxEstLhood == max(mod2$MaxEstLhood)),]

# Model 3, exponential change after bottleneck
mod3_ml <- max(mod3$MaxEstLhood)
mod3_ml_parameters <- mod3[which(mod3$MaxEstLhood == max(mod3$MaxEstLhood)),]

# Model 4, exponential change in population size, then bottleneck, then instantaneous recovery
mod4_ml <- max(mod4$MaxEstLhood)
mod4_ml_parameters <- mod4[which(mod4$MaxEstLhood == max(mod4$MaxEstLhood)),]

# Model 5, two bottlenecks
mod5_ml <- max(mod5$MaxEstLhood)
mod5_ml_parameters <- mod5[which(mod5$MaxEstLhood == max(mod5$MaxEstLhood)),]

# Model 6, exponential change in pop size before and after bottleneck
mod6_ml <- max(mod6$MaxEstLhood)
mod6_ml_parameters <- mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]

# Model 7, ancestral change in population size up to carrying capacity
mod7_ml <- max(mod7$MaxEstLhood)
mod7_ml_parameters <- mod7[which(mod7$MaxEstLhood == max(mod7$MaxEstLhood)),]

#### Determine best-fit model for each POD ####
b <- c(1, 5, 5, 6, 9, 6, 3) # number of parameters in each of the models going from model 1 to model 7


#### Plot ML parameters for each model as a line. Then plot the ML model for each of the 10 PODS ####
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/PODs_fit.png",width=8, height=11, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel margin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12,
    mfrow = c(4,2)
)

#### Model 1 ####
max.mod1 <- data.frame(matrix(NA, nrow = 13, ncol = 4))
max.mod1[,1] <- c(seq(0, 23, 2),100) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
max.mod1[,2] <- 2008 - 2*max.mod1[,1]
max.mod1[,3] <- rep(mod1_ml_parameters$NPOP08, 13)
max.mod1[,4] <- max.mod1[,3]/2 #diploid

# Model 1 PODs
mod1pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod1pods", pattern = "*summary.txt", full.names = TRUE)
mod1_files <- lapply(mod1pods_names, function(x)read.table(x, header=F, fill = T))
mod1_files_array <- array(as.numeric(unlist(mod1_files)), dim=c(7, 11, 10))

mod1_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod1_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod1_files_array_sub[1,,i] <- mod1_files_array[1,c(1,4:5,2:3),i]
  mod1_files_array_sub[2,,i] <- mod1_files_array[2,c(1:3,6:7),i]
  mod1_files_array_sub[3,,i] <- mod1_files_array[3,c(1:3,6:7),i]
  mod1_files_array_sub[4,,i] <- mod1_files_array[4,c(1:3,8:9),i]
  mod1_files_array_sub[5,,i] <- mod1_files_array[5,c(1:3,10:11),i]
  mod1_files_array_sub[6,,i] <- mod1_files_array[6,c(1:3,8:9),i]
  mod1_files_array_sub[7,,i] <- mod1_files_array[7,c(1,7:8,5:6),i]
}

# Model 1
mod1.aic <- 2*b-2*(mod1_files_array_sub[,4,] *2.303) # Convert from log10 to ln, first
mod1.ass <- apply(mod1.aic,2,which.min)

# Get ML parameter for each fit POD
mod1_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod1_pods_mlfits[i,] <- mod1_files_array[mod1.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod1_pods_mod1_fit <- mod1_pods_mlfits[which(mod1.ass == 1),-c(4:11)]
colnames(mod1_pods_mod1_fit) <- colnames(mod1_ml_parameters)
mod1_pods_mod2_fit <- mod1_pods_mlfits[which(mod1.ass == 2),-c(8:11)]
colnames(mod1_pods_mod2_fit) <- colnames(mod2_ml_parameters)
mod1_pods_mod3_fit <- mod1_pods_mlfits[which(mod1.ass == 3),-c(8:11)]
colnames(mod1_pods_mod3_fit) <- colnames(mod3_ml_parameters)
mod1_pods_mod5_fit <- mod1_pods_mlfits[which(mod1.ass == 5),]
colnames(mod1_pods_mod5_fit) <- colnames(mod5_ml_parameters)
mod1_pods_mod7_fit <-mod1_pods_mlfits[which(mod1.ass == 7),-c(7:11)]
colnames(mod1_pods_mod7_fit) <- colnames(mod7_ml_parameters)

# Make curves for fitted PODs
# Model 1
mod1pod.mod1fit.coord <- array(numeric(), c(12,4,nrow(mod1_pods_mod1_fit)))
for (i in 1:nrow(mod1_pods_mod1_fit)) {
  mod1pod.mod1fit.coord[,1,i] <- c(seq(0, 21, 2),2000) #generations going back in time
  mod1pod.mod1fit.coord[,2,i] <- 2008 - 2*mod1pod.mod1fit.coord[,1,i]
  mod1pod.mod1fit.coord[,3,i] <- rep(mod1_pods_mod1_fit$NPOP08[i], 12)
  mod1pod.mod1fit.coord[,4,i] <- mod1pod.mod1fit.coord[,3,i]/2 #diploid
}

# Model 2
mod1pod.mod2fit.coord <- array(numeric(), c(10,4,nrow(mod1_pods_mod2_fit)))
for (i in 1:nrow(mod1_pods_mod2_fit)) {
  mod1pod.mod2fit.coord[,1,i] <- c(mod1_pods_mod2_fit$TBOT[i]-mod1_pods_mod2_fit$TBOT[i], mod1_pods_mod2_fit$TBOT[i]-3, mod1_pods_mod2_fit$TBOT[i]-2, mod1_pods_mod2_fit$TBOT[i]-1, mod1_pods_mod2_fit$TBOT[i], mod1_pods_mod2_fit$TBOT[i], (mod1_pods_mod2_fit$TBOT[i]+mod1_pods_mod2_fit$TLEN[i]), (mod1_pods_mod2_fit$TBOT[i]+mod1_pods_mod2_fit$TLEN[i]), (mod1_pods_mod2_fit$TBOT[i]+mod1_pods_mod2_fit$TLEN[i]+2), (mod1_pods_mod2_fit$TBOT[i]+mod1_pods_mod2_fit$TLEN[i]+50)) #generations going back in time
  mod1pod.mod2fit.coord[,2,i] <- 2008 - 2*mod1pod.mod2fit.coord[,1,i]
  mod1pod.mod2fit.coord[,3,i] <- c(mod1_pods_mod2_fit$NPOP08[i], mod1_pods_mod2_fit$NPOP08[i], mod1_pods_mod2_fit$NPOP08[i], mod1_pods_mod2_fit$NPOP08[i], mod1_pods_mod2_fit$NPOP08[i], mod1_pods_mod2_fit$NBOT[i], mod1_pods_mod2_fit$NBOT[i], mod1_pods_mod2_fit$PREBOT[i], mod1_pods_mod2_fit$NPREBOT[i], mod1_pods_mod2_fit$NPREBOT[i], mod1_pods_mod2_fit$NPREBOT[i]) #haploid
  mod1pod.mod2fit.coord[,4,i] <- mod1pod.mod2fit.coord[,3,i]/2 #diploid
}

# Model 3
r2008 <- (log(mod1_pods_mod3_fit$NBOT/mod1_pods_mod3_fit$NPOP08)/(mod1_pods_mod3_fit$TBOT)) # check that this is correct; growth rate going back in time

mod1pod.mod3fit.coord <- array(numeric(), c(13,4,nrow(mod1_pods_mod3_fit)))
for (i in 1:nrow(mod1_pods_mod3_fit)) {
  mod1pod.mod3fit.coord[,1,i] <- c(mod1_pods_mod3_fit$TBOT[i]-mod1_pods_mod3_fit$TBOT[i], mod1_pods_mod3_fit$TBOT[i]-6, mod1_pods_mod3_fit$TBOT[i]-5, mod1_pods_mod3_fit$TBOT[i]-4, mod1_pods_mod3_fit$TBOT[i]-3, mod1_pods_mod3_fit$TBOT[i]-2, mod1_pods_mod3_fit$TBOT[i]-1, mod1_pods_mod3_fit$TBOT[i], (mod1_pods_mod3_fit$TBOT[i]+mod1_pods_mod3_fit$TLEN[i]), (mod1_pods_mod3_fit$TBOT[i]+mod1_pods_mod3_fit$TLEN[i]), (mod1_pods_mod3_fit$TBOT[i]+mod1_pods_mod3_fit$TLEN[i]+2), (mod1_pods_mod3_fit$TBOT[i]+mod1_pods_mod3_fit$TLEN[i]+5), (mod1_pods_mod3_fit$TBOT[i]+mod1_pods_mod3_fit$TLEN[i]+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
  mod1pod.mod3fit.coord[,2,i] <- 2008 - 2*mod1pod.mod3fit.coord[,1,i]
  mod1pod.mod3fit.coord[,3,i] <- c(mod1_pods_mod3_fit$NPOP08[i], (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-6))), (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-5))), (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-4))), (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-3))), (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-2))), (mod1_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod1_pods_mod3_fit$TBOT[i]-1))), mod1_pods_mod3_fit$NBOT[i], mod1_pods_mod3_fit$NBOT[i], mod1_pods_mod3_fit$NPREBOT[i], mod1_pods_mod3_fit$NPREBOT[i], mod1_pods_mod3_fit$NPREBOT[i], mod1_pods_mod3_fit$NPREBOT[i]) #haploid
  mod1pod.mod3fit.coord[,4,i] <- mod1pod.mod3fit.coord[,3,i]/2 #diploid
}

# Model 5
mod1pod.mod5fit.coord <- array(numeric(), c(10,4,nrow(mod1_pods_mod5_fit)))
for (i in 1:nrow(mod1_pods_mod5_fit)) {
  mod1pod.mod5fit.coord[,1,i] <- c((mod1_pods_mod5_fit$TBOTTWO[i]-mod1_pods_mod5_fit$TBOTTWO[i]), mod1_pods_mod5_fit$TBOTTWO[i], mod1_pods_mod5_fit$TBOTTWO[i], (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i] + mod1_pods_mod5_fit$TBOTONE[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i] + mod1_pods_mod5_fit$TBOTONE[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i] + mod1_pods_mod5_fit$TBOTONE[i] + mod1_pods_mod5_fit$TLENONE[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i] + mod1_pods_mod5_fit$TBOTONE[i] + mod1_pods_mod5_fit$TLENONE[i]), (mod1_pods_mod5_fit$TBOTTWO[i] + mod1_pods_mod5_fit$TLENTWO[i] + mod1_pods_mod5_fit$TBOTONE[i] + mod1_pods_mod5_fit$TLENONE[i] + 100)) #generations going back in time
  mod1pod.mod5fit.coord[,2,i] <- 2008 - 2*mod1pod.mod5fit.coord[,1,i]
  mod1pod.mod5fit.coord[,3,i] <- c(mod1_pods_mod5_fit$NPOP08[i], mod1_pods_mod5_fit$NPOP08[i], mod1_pods_mod5_fit$NBOTTWO[i], mod1_pods_mod5_fit$NBOTTWO[i], mod1_pods_mod5_fit$NPREBOT[i], mod1_pods_mod5_fit$NPREBOT[i], mod1_pods_mod5_fit$NBOTONE[i], mod1_pods_mod5_fit$NBOTONE[i], mod1_pods_mod5_fit$NANC[i], mod1_pods_mod5_fit$NANC[i]) #haploid
  mod1pod.mod5fit.coord[,4,i] <- mod1pod.mod5fit.coord[,3,i]/2 #diploid
}

# Model 7
mod1pod.mod7fit.coord <- array(numeric(), c(4,4,nrow(mod1_pods_mod7_fit)))
for (i in 1:nrow(mod1_pods_mod7_fit)) {
  mod1pod.mod7fit.coord[,1,i] <- c(mod1_pods_mod7_fit$TCAR[i]-mod1_pods_mod7_fit$TCAR[i], mod1_pods_mod7_fit$TCAR[i], (mod1_pods_mod7_fit$TCAR[i]+5), (mod1_pods_mod7_fit$TCAR[i]+1000)) #generations going back in time
  mod1pod.mod7fit.coord[,2,i] <- 2008 -2*mod1pod.mod7fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod1pod.mod7fit.coord[,3,i] <- c(mod1_pods_mod7_fit$NPOP08[i], mod1_pods_mod7_fit$NPOP08[i], (mod1_pods_mod7_fit$NANC[i]*exp(mod1_pods_mod7_fit$RANC[i] * 5)), (mod1_pods_mod7_fit$NANC[i]*exp(mod1_pods_mod7_fit$RANC[i] * 1000))) #haploid
  mod1pod.mod7fit.coord[,4,i] <- mod1pod.mod7fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod1$X2, max.mod1$X4, xlab = '', ylab = '', type = 'n', xlim = c(1965,2008), ylim = c(0,30000), las = 1)
lines(max.mod1$X2, max.mod1$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:5) {
  lines(jitter(mod1pod.mod1fit.coord[,2,l], factor = 0.2), mod1pod.mod1fit.coord[,4,l], col = cols) # plots 5 line: Model 1 fits to Model 1 PODs
}
for (l in 1:1) {
  lines(jitter(mod1pod.mod2fit.coord[-2,2,l], factor = 0.2), mod1pod.mod2fit.coord[-2,4,l], col = cols) # plots 2 line: Model 2 fits to Model 1 PODs; plot separately so that the first line stops at 2008
}
for (l in 2:2) {
  lines(jitter(mod1pod.mod2fit.coord[,2,l], factor = 0.2), mod1pod.mod2fit.coord[,4,l], col = cols) # plots 2 line: Model 2 fits to Model 1 PODs
}
for (l in 1:1) {
  lines(jitter(mod1pod.mod3fit.coord[-c(2:5),2,l], factor = 0.2), mod1pod.mod3fit.coord[-c(2:5),4,l], col = cols) # plots 1 line: Model 3 fits to Model 1 PODs; manually make line stop at 2008
}
for (l in 1:1) {
  lines(jitter(mod1pod.mod5fit.coord[,2,l], factor = 0.2), mod1pod.mod5fit.coord[,4,l], col = cols) # plots 1 line: Model 5 fits to Model 1 PODs
}
for (l in 1:1) {
  lines(jitter(mod1pod.mod7fit.coord[,2,l], factor = 0.2), mod1pod.mod7fit.coord[,4,l], col = cols) # plots 1 line: Model 7 fits to Model 1 PODs
}

#### Model 2 ####
max.mod2 <- data.frame(matrix(NA, nrow = 10, ncol = 4))
max.mod2[,1] <- c(mod2_ml_parameters$TBOT-7, mod2_ml_parameters$TBOT-5, mod2_ml_parameters$TBOT-3, mod2_ml_parameters$TBOT-1, mod2_ml_parameters$TBOT, mod2_ml_parameters$TBOT, (mod2_ml_parameters$TBOT+mod2_ml_parameters$TLEN), (mod2_ml_parameters$TBOT+mod2_ml_parameters$TLEN), (mod2_ml_parameters$TBOT+mod2_ml_parameters$TLEN+2), (mod2_ml_parameters$TBOT+mod2_ml_parameters$TLEN+50)) #generations going back in time
max.mod2[,2] <- 2008 - 2*max.mod2[,1]
max.mod2[,3] <- c(mod2_ml_parameters$NPOP08, mod2_ml_parameters$NPOP08, mod2_ml_parameters$NPOP08, mod2_ml_parameters$NPOP08, mod2_ml_parameters$NPOP08, mod2_ml_parameters$NBOT, mod2_ml_parameters$NBOT, mod2_ml_parameters$PREBOT, mod2_ml_parameters$NPREBOT, mod2_ml_parameters$NPREBOT, mod2_ml_parameters$NPREBOT) #haploid
max.mod2[,4] <- max.mod2[,3]/2 #diploid

# Model 2 POD fits
mod2pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod2pods", pattern = "*summary.txt", full.names = TRUE)
mod2_files <- lapply(mod2pods_names, function(x)read.table(x, header=F, fill = T))
mod2_files_array <- array(as.numeric(unlist(mod2_files)), dim=c(7, 11, 10))

mod2_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod2_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod2_files_array_sub[1,,i] <- mod2_files_array[1,c(1,4:5,2:3),i]
  mod2_files_array_sub[2,,i] <- mod2_files_array[2,c(1:3,6:7),i]
  mod2_files_array_sub[3,,i] <- mod2_files_array[3,c(1:3,6:7),i]
  mod2_files_array_sub[4,,i] <- mod2_files_array[4,c(1:3,8:9),i]
  mod2_files_array_sub[5,,i] <- mod2_files_array[5,c(1:3,10:11),i]
  mod2_files_array_sub[6,,i] <- mod2_files_array[6,c(1:3,8:9),i]
  mod2_files_array_sub[7,,i] <- mod2_files_array[7,c(1,7:8,5:6),i]
}

mod2.aic <- 2*b-2*(mod2_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod2.ass <- apply(mod2.aic,2,which.min)

# Get ML parameter for each fit POD
mod2_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod2_pods_mlfits[i,] <- mod2_files_array[mod2.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod2_pods_mod1_fit <- mod2_pods_mlfits[which(mod2.ass == 1),-c(4:11)]
colnames(mod2_pods_mod1_fit) <- colnames(mod1_ml_parameters)
mod2_pods_mod2_fit <- mod2_pods_mlfits[which(mod2.ass == 2),-c(8:11)]
colnames(mod2_pods_mod2_fit) <- colnames(mod2_ml_parameters)
mod2_pods_mod3_fit <- mod2_pods_mlfits[which(mod2.ass == 3),-c(8:11)]
colnames(mod2_pods_mod3_fit) <- colnames(mod3_ml_parameters)
mod2_pods_mod4_fit <-mod2_pods_mlfits[which(mod2.ass == 4),-c(10:11)]
colnames(mod2_pods_mod4_fit) <- colnames(mod4_ml_parameters)
mod2_pods_mod5_fit <- mod2_pods_mlfits[which(mod2.ass == 5),]
colnames(mod2_pods_mod5_fit) <- colnames(mod5_ml_parameters)
mod2_pods_mod6_fit <- mod2_pods_mlfits[which(mod2.ass == 6),-c(10:11)]
colnames(mod2_pods_mod6_fit) <- colnames(mod6_ml_parameters)

# Make curves for fitted PODs
# Model 1
mod2pod.mod1fit.coord <- array(numeric(), c(12,4,nrow(mod2_pods_mod1_fit)))
for (i in 1:nrow(mod2_pods_mod1_fit)) {
  mod2pod.mod1fit.coord[,1,i] <- c(seq(0, 21, 2),2000) #generations going back in time
  mod2pod.mod1fit.coord[,2,i] <- 2008 - 2*mod2pod.mod1fit.coord[,1,i]
  mod2pod.mod1fit.coord[,3,i] <- rep(mod2_pods_mod1_fit$NPOP08[i], 12)
  mod2pod.mod1fit.coord[,4,i] <- mod2pod.mod1fit.coord[,3,i]/2 #diploid
}

# Model 2
mod2pod.mod2fit.coord <- array(numeric(), c(10,4,nrow(mod2_pods_mod2_fit)))
for (i in 1:nrow(mod2_pods_mod2_fit)) {
  mod2pod.mod2fit.coord[,1,i] <- c(mod2_pods_mod2_fit$TBOT[i]-mod2_pods_mod2_fit$TBOT[i], mod2_pods_mod2_fit$TBOT[i]-3, mod2_pods_mod2_fit$TBOT[i]-2, mod2_pods_mod2_fit$TBOT[i]-1, mod2_pods_mod2_fit$TBOT[i], mod2_pods_mod2_fit$TBOT[i], (mod2_pods_mod2_fit$TBOT[i]+mod2_pods_mod2_fit$TLEN[i]), (mod2_pods_mod2_fit$TBOT[i]+mod2_pods_mod2_fit$TLEN[i]), (mod2_pods_mod2_fit$TBOT[i]+mod2_pods_mod2_fit$TLEN[i]+2), (mod2_pods_mod2_fit$TBOT[i]+mod2_pods_mod2_fit$TLEN[i]+50)) #generations going back in time
  mod2pod.mod2fit.coord[,2,i] <- 2008 - 2*mod2pod.mod2fit.coord[,1,i]
  mod2pod.mod2fit.coord[,3,i] <- c(mod2_pods_mod2_fit$NPOP08[i], mod2_pods_mod2_fit$NPOP08[i], mod2_pods_mod2_fit$NPOP08[i], mod2_pods_mod2_fit$NPOP08[i], mod2_pods_mod2_fit$NPOP08[i], mod2_pods_mod2_fit$NBOT[i], mod2_pods_mod2_fit$NBOT[i], mod2_pods_mod2_fit$PREBOT[i], mod2_pods_mod2_fit$NPREBOT[i], mod2_pods_mod2_fit$NPREBOT[i], mod2_pods_mod2_fit$NPREBOT[i]) #haploid
  mod2pod.mod2fit.coord[,4,i] <- mod2pod.mod2fit.coord[,3,i]/2 #diploid
}

# Model 3
r2008 <- (log(mod2_pods_mod3_fit$NBOT/mod2_pods_mod3_fit$NPOP08)/(mod2_pods_mod3_fit$TBOT)) # check that this is correct; growth rate going back in time

mod2pod.mod3fit.coord <- array(numeric(), c(13,4,nrow(mod2_pods_mod3_fit)))
for (i in 1:nrow(mod2_pods_mod3_fit)) {
  mod2pod.mod3fit.coord[,1,i] <- c(mod2_pods_mod3_fit$TBOT[i]-mod2_pods_mod3_fit$TBOT[i], mod2_pods_mod3_fit$TBOT[i]-6, mod2_pods_mod3_fit$TBOT[i]-5, mod2_pods_mod3_fit$TBOT[i]-4, mod2_pods_mod3_fit$TBOT[i]-3, mod2_pods_mod3_fit$TBOT[i]-2, mod2_pods_mod3_fit$TBOT[i]-1, mod2_pods_mod3_fit$TBOT[i], (mod2_pods_mod3_fit$TBOT[i]+mod2_pods_mod3_fit$TLEN[i]), (mod2_pods_mod3_fit$TBOT[i]+mod2_pods_mod3_fit$TLEN[i]), (mod2_pods_mod3_fit$TBOT[i]+mod2_pods_mod3_fit$TLEN[i]+2), (mod2_pods_mod3_fit$TBOT[i]+mod2_pods_mod3_fit$TLEN[i]+5), (mod2_pods_mod3_fit$TBOT[i]+mod2_pods_mod3_fit$TLEN[i]+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
  mod2pod.mod3fit.coord[,2,i] <- 2008 - 2*mod2pod.mod3fit.coord[,1,i]
  mod2pod.mod3fit.coord[,3,i] <- c(mod2_pods_mod3_fit$NPOP08[i], (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-6))), (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-5))), (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-4))), (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-3))), (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-2))), (mod2_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod2_pods_mod3_fit$TBOT[i]-1))), mod2_pods_mod3_fit$NBOT[i], mod2_pods_mod3_fit$NBOT[i], mod2_pods_mod3_fit$NPREBOT[i], mod2_pods_mod3_fit$NPREBOT[i], mod2_pods_mod3_fit$NPREBOT[i], mod2_pods_mod3_fit$NPREBOT[i]) #haploid
  mod2pod.mod3fit.coord[,4,i] <- mod2pod.mod3fit.coord[,3,i]/2 #diploid
}

# Model 4
mod2pod.mod4fit.coord <- array(numeric(), c(12,4,nrow(mod2_pods_mod4_fit)))
for (i in 1:nrow(mod2_pods_mod4_fit)) {
  mod2pod.mod4fit.coord[,1,i] <- c(mod2_pods_mod4_fit$TBOT[i]-mod2_pods_mod4_fit$TBOT[i], mod2_pods_mod4_fit$TBOT[i]-3, mod2_pods_mod4_fit$TBOT[i]-2, mod2_pods_mod4_fit$TBOT[i]-1, mod2_pods_mod4_fit$TBOT[i], mod2_pods_mod4_fit$TBOT[i], (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]), (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]), (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]+2), (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]+5), (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]+10), (mod2_pods_mod4_fit$TBOT[i]+mod2_pods_mod4_fit$TLEN[i]+1000)) #generations going back in time
  mod2pod.mod4fit.coord[,2,i] <- 2008 -2*mod2pod.mod4fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod2pod.mod4fit.coord[,3,i] <- c(mod2_pods_mod4_fit$NPOP08[i], mod2_pods_mod4_fit$NPOP08[i], mod2_pods_mod4_fit$NPOP08[i], mod2_pods_mod4_fit$NPOP08[i], mod2_pods_mod4_fit$NPOP08[i], mod2_pods_mod4_fit$NBOT[i], mod2_pods_mod4_fit$NBOT[i], mod2_pods_mod4_fit$NPREBOT[i], (mod2_pods_mod4_fit$NPREBOT[i]*exp(mod2_pods_mod4_fit$RANC[i] * 2)), (mod2_pods_mod4_fit$NPREBOT[i]*exp(mod2_pods_mod4_fit$RANC[i] * 5)), (mod2_pods_mod4_fit$NPREBOT[i]*exp(mod2_pods_mod4_fit$RANC[i] * 10)), (mod2_pods_mod4_fit$NPREBOT[i]*exp(mod2_pods_mod4_fit$RANC[i] * 1000))) #haploid
  mod2pod.mod4fit.coord[,4,i] <- mod2pod.mod4fit.coord[,3,i]/2 #diploid
}

# Model 5
mod2pod.mod5fit.coord <- array(numeric(), c(10,4,nrow(mod2_pods_mod5_fit)))
for (i in 1:nrow(mod2_pods_mod5_fit)) {
  mod2pod.mod5fit.coord[,1,i] <- c((mod2_pods_mod5_fit$TBOTTWO[i]-mod2_pods_mod5_fit$TBOTTWO[i]), mod2_pods_mod5_fit$TBOTTWO[i], mod2_pods_mod5_fit$TBOTTWO[i], (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i] + mod2_pods_mod5_fit$TBOTONE[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i] + mod2_pods_mod5_fit$TBOTONE[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i] + mod2_pods_mod5_fit$TBOTONE[i] + mod2_pods_mod5_fit$TLENONE[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i] + mod2_pods_mod5_fit$TBOTONE[i] + mod2_pods_mod5_fit$TLENONE[i]), (mod2_pods_mod5_fit$TBOTTWO[i] + mod2_pods_mod5_fit$TLENTWO[i] + mod2_pods_mod5_fit$TBOTONE[i] + mod2_pods_mod5_fit$TLENONE[i] + 100)) #generations going back in time
  mod2pod.mod5fit.coord[,2,i] <- 2008 - 2*mod2pod.mod5fit.coord[,1,i]
  mod2pod.mod5fit.coord[,3,i] <- c(mod2_pods_mod5_fit$NPOP08[i], mod2_pods_mod5_fit$NPOP08[i], mod2_pods_mod5_fit$NBOTTWO[i], mod2_pods_mod5_fit$NBOTTWO[i], mod2_pods_mod5_fit$NPREBOT[i], mod2_pods_mod5_fit$NPREBOT[i], mod2_pods_mod5_fit$NBOTONE[i], mod2_pods_mod5_fit$NBOTONE[i], mod2_pods_mod5_fit$NANC[i], mod2_pods_mod5_fit$NANC[i]) #haploid
  mod2pod.mod5fit.coord[,4,i] <- mod2pod.mod5fit.coord[,3,i]/2 #diploid
}

# Model 6
# First calculate recent r (R2008) for each SFS
r <- vector(length = nrow(mod2_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod2_pods_mod6_fit$NBOT[i]/mod2_pods_mod6_fit$NPOP08[i])/(mod2_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod2pod.mod6fit.coord <- array(numeric(), c(10,4,nrow(mod2_pods_mod6_fit)))
for (i in 1:nrow(mod2_pods_mod6_fit)) {
  mod2pod.mod6fit.coord[,1,i] <- c(0, mod2_pods_mod6_fit$TBOT[i]-3, mod2_pods_mod6_fit$TBOT[i]-2, mod2_pods_mod6_fit$TBOT[i]-1, mod2_pods_mod6_fit$TBOT[i], (mod2_pods_mod6_fit$TBOT[i]+mod2_pods_mod6_fit$TLEN[i]), (mod2_pods_mod6_fit$TBOT[i]+mod2_pods_mod6_fit$TLEN[i]), (mod2_pods_mod6_fit$TBOT[i]+mod2_pods_mod6_fit$TLEN[i]+2), (mod2_pods_mod6_fit$TBOT[i]+mod2_pods_mod6_fit$TLEN[i]+5), (mod2_pods_mod6_fit$TBOT[i]+mod2_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod2pod.mod6fit.coord[,2,i] <- 2008 -2*mod2pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod2pod.mod6fit.coord[,3,i] <- c(mod2_pods_mod6_fit$NPOP08[i], (mod2_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod2_pods_mod6_fit$TBOT[i]-3))), (mod2_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod2_pods_mod6_fit$TBOT[i]-2))), (mod2_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod2_pods_mod6_fit$TBOT[i]-1))), mod2_pods_mod6_fit$NBOT[i], mod2_pods_mod6_fit$NBOT[i], mod2_pods_mod6_fit$NPREBOT[i], (mod2_pods_mod6_fit$NPREBOT[i]*exp(mod2_pods_mod6_fit$RANC[i] * 2)), (mod2_pods_mod6_fit$NPREBOT[i]*exp(mod2_pods_mod6_fit$RANC[i] * 5)), (mod2_pods_mod6_fit$NPREBOT[i]*exp(mod2_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod2pod.mod6fit.coord[,4,i] <- mod2pod.mod6fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod2$X2, max.mod2$X4, xlab = '', ylab = '', type = 'n', xlim = c(1976,2008), ylim = c(0,30000), las = 1)
lines(max.mod2$X2, max.mod2$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:1) {
  lines(jitter(mod2pod.mod1fit.coord[,2,l], factor = 0.2), mod2pod.mod1fit.coord[,4,l], col = cols) # plots 1 line: Model 1 fits to Model 2 PODs
}
for (l in 1:1) {
  lines(jitter(mod2pod.mod2fit.coord[,2,l], factor = 0.2), mod2pod.mod2fit.coord[,4,l], col = cols) # plots 1 line: Model 2 fits to Model 2 PODs
}
for (l in 1:1) {
  lines(jitter(mod2pod.mod3fit.coord[-c(2:5),2,l], factor = 0.2), mod2pod.mod3fit.coord[-c(2:5),4,l], col = cols) # plots 2 line: Model 3 fits to Model 2 PODs; plot these separately to stop the line from extending past 2008
}
for (l in 2:2) {
  lines(jitter(mod2pod.mod3fit.coord[,2,l], factor = 0.2), mod2pod.mod3fit.coord[,4,l], col = cols) # plots 2 line: Model 3 fits to Model 2 PODs
}
for (l in 1:1) {
  lines(jitter(mod2pod.mod4fit.coord[-c(2:3),2,l], factor = 0.2), mod2pod.mod4fit.coord[-c(2:3),4,l], col = cols) # plots 1 line: Model 4 fits to Model 2 PODs; manually stop line from extending past 2008
}
for (l in 1:2) {
  lines(jitter(mod2pod.mod5fit.coord[,2,l], factor = 0.2), mod2pod.mod5fit.coord[,4,l], col = cols) # plots 2 line: Model 5 fits to Model 2 PODs
}
for (l in 1:3) {
  lines(jitter(mod2pod.mod6fit.coord[,2,l], factor = 0.2), mod2pod.mod6fit.coord[,4,l], col = cols) # plots 3 line: Model 6 fits to Model 2 PODs
}

#### Model 3 ####
r2008 <- (log(mod3_ml_parameters$NBOT/mod3_ml_parameters$NPOP08)/(mod3_ml_parameters$TBOT)) # check that this is correct; growth rate going back in time

max.mod3 <- data.frame(matrix(NA, nrow = 13, ncol = 4))
max.mod3[,1] <- c(0, mod3_ml_parameters$TBOT-10, mod3_ml_parameters$TBOT-8, mod3_ml_parameters$TBOT-6, mod3_ml_parameters$TBOT-4, mod3_ml_parameters$TBOT-3, mod3_ml_parameters$TBOT-1, mod3_ml_parameters$TBOT, (mod3_ml_parameters$TBOT+mod3_ml_parameters$TLEN), (mod3_ml_parameters$TBOT+mod3_ml_parameters$TLEN), (mod3_ml_parameters$TBOT+mod3_ml_parameters$TLEN+2), (mod3_ml_parameters$TBOT+mod3_ml_parameters$TLEN+5), (mod3_ml_parameters$TBOT+mod3_ml_parameters$TLEN+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
max.mod3[,2] <- 2008 - 2*max.mod3[,1]
max.mod3[,3] <- c(mod3_ml_parameters$NPOP08, (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-10))), (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-8))), (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-6))), (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-4))), (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-3))), (mod3_ml_parameters$NPOP08*exp(r2008 * (mod3_ml_parameters$TBOT-1))), mod3_ml_parameters$NBOT, mod3_ml_parameters$NBOT, mod3_ml_parameters$NPREBOT, mod3_ml_parameters$NPREBOT, mod3_ml_parameters$NPREBOT, mod3_ml_parameters$NPREBOT) #haploid
max.mod3[,4] <- max.mod3[,3]/2 #diploid

# Model 3 POD fits
mod3pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod3pods", pattern = "*summary.txt", full.names = TRUE)
mod3_files <- lapply(mod3pods_names, function(x)read.table(x, header=F, fill = T))
mod3_files_array <- array(as.numeric(unlist(mod3_files)), dim=c(7, 11, 10))

mod3_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod3_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod3_files_array_sub[1,,i] <- mod3_files_array[1,c(1,4:5,2:3),i]
  mod3_files_array_sub[2,,i] <- mod3_files_array[2,c(1:3,6:7),i]
  mod3_files_array_sub[3,,i] <- mod3_files_array[3,c(1:3,6:7),i]
  mod3_files_array_sub[4,,i] <- mod3_files_array[4,c(1:3,8:9),i]
  mod3_files_array_sub[5,,i] <- mod3_files_array[5,c(1:3,10:11),i]
  mod3_files_array_sub[6,,i] <- mod3_files_array[6,c(1:3,8:9),i]
  mod3_files_array_sub[7,,i] <- mod3_files_array[7,c(1,7:8,5:6),i]
}

mod3.aic <- 2*b-2*(mod3_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod3.ass <- apply(mod3.aic,2,which.min)

# Get ML parameter for each fit POD
mod3_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod3_pods_mlfits[i,] <- mod3_files_array[mod3.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod3_pods_mod1_fit <- mod3_pods_mlfits[which(mod3.ass == 1),-c(4:11)]
colnames(mod3_pods_mod1_fit) <- colnames(mod1_ml_parameters)
mod3_pods_mod2_fit <- mod3_pods_mlfits[which(mod3.ass == 2),-c(8:11)]
colnames(mod3_pods_mod2_fit) <- colnames(mod2_ml_parameters)
mod3_pods_mod3_fit <- mod3_pods_mlfits[which(mod3.ass == 3),-c(8:11)]
colnames(mod3_pods_mod3_fit) <- colnames(mod3_ml_parameters)
mod3_pods_mod5_fit <- mod3_pods_mlfits[which(mod3.ass == 5),]
colnames(mod3_pods_mod5_fit) <- colnames(mod5_ml_parameters)
mod3_pods_mod6_fit <- mod3_pods_mlfits[which(mod3.ass == 6),-c(10:11)]
colnames(mod3_pods_mod6_fit) <- colnames(mod6_ml_parameters)
mod3_pods_mod7_fit <-mod3_pods_mlfits[which(mod3.ass == 7),-c(7:11)]
colnames(mod3_pods_mod7_fit) <- colnames(mod7_ml_parameters)

# Make curves for fitted PODs
# Model 1
mod3pod.mod1fit.coord <- array(numeric(), c(12,4,nrow(mod3_pods_mod1_fit)))
for (i in 1:nrow(mod3_pods_mod1_fit)) {
  mod3pod.mod1fit.coord[,1,i] <- c(seq(0, 21, 2),2000) #generations going back in time
  mod3pod.mod1fit.coord[,2,i] <- 2008 - 2*mod3pod.mod1fit.coord[,1,i]
  mod3pod.mod1fit.coord[,3,i] <- rep(mod3_pods_mod1_fit$NPOP08[i], 12)
  mod3pod.mod1fit.coord[,4,i] <- mod3pod.mod1fit.coord[,3,i]/2 #diploid
}

# Model 2
mod3pod.mod2fit.coord <- array(numeric(), c(10,4,nrow(mod3_pods_mod2_fit)))
for (i in 1:nrow(mod3_pods_mod2_fit)) {
  mod3pod.mod2fit.coord[,1,i] <- c(mod3_pods_mod2_fit$TBOT[i]-mod3_pods_mod2_fit$TBOT[i], mod3_pods_mod2_fit$TBOT[i]-3, mod3_pods_mod2_fit$TBOT[i]-2, mod3_pods_mod2_fit$TBOT[i]-1, mod3_pods_mod2_fit$TBOT[i], mod3_pods_mod2_fit$TBOT[i], (mod3_pods_mod2_fit$TBOT[i]+mod3_pods_mod2_fit$TLEN[i]), (mod3_pods_mod2_fit$TBOT[i]+mod3_pods_mod2_fit$TLEN[i]), (mod3_pods_mod2_fit$TBOT[i]+mod3_pods_mod2_fit$TLEN[i]+2), (mod3_pods_mod2_fit$TBOT[i]+mod3_pods_mod2_fit$TLEN[i]+50)) #generations going back in time
  mod3pod.mod2fit.coord[,2,i] <- 2008 - 2*mod3pod.mod2fit.coord[,1,i]
  mod3pod.mod2fit.coord[,3,i] <- c(mod3_pods_mod2_fit$NPOP08[i], mod3_pods_mod2_fit$NPOP08[i], mod3_pods_mod2_fit$NPOP08[i], mod3_pods_mod2_fit$NPOP08[i], mod3_pods_mod2_fit$NPOP08[i], mod3_pods_mod2_fit$NBOT[i], mod3_pods_mod2_fit$NBOT[i], mod3_pods_mod2_fit$PREBOT[i], mod3_pods_mod2_fit$NPREBOT[i], mod3_pods_mod2_fit$NPREBOT[i], mod3_pods_mod2_fit$NPREBOT[i]) #haploid
  mod3pod.mod2fit.coord[,4,i] <- mod3pod.mod2fit.coord[,3,i]/2 #diploid
}

# Model 3
r2008 <- (log(mod3_pods_mod3_fit$NBOT/mod3_pods_mod3_fit$NPOP08)/(mod3_pods_mod3_fit$TBOT)) # check that this is correct; growth rate going back in time

mod3pod.mod3fit.coord <- array(numeric(), c(13,4,nrow(mod3_pods_mod3_fit)))
for (i in 1:nrow(mod3_pods_mod3_fit)) {
  mod3pod.mod3fit.coord[,1,i] <- c(mod3_pods_mod3_fit$TBOT[i]-mod3_pods_mod3_fit$TBOT[i], mod3_pods_mod3_fit$TBOT[i]-6, mod3_pods_mod3_fit$TBOT[i]-5, mod3_pods_mod3_fit$TBOT[i]-4, mod3_pods_mod3_fit$TBOT[i]-3, mod3_pods_mod3_fit$TBOT[i]-2, mod3_pods_mod3_fit$TBOT[i]-1, mod3_pods_mod3_fit$TBOT[i], (mod3_pods_mod3_fit$TBOT[i]+mod3_pods_mod3_fit$TLEN[i]), (mod3_pods_mod3_fit$TBOT[i]+mod3_pods_mod3_fit$TLEN[i]), (mod3_pods_mod3_fit$TBOT[i]+mod3_pods_mod3_fit$TLEN[i]+2), (mod3_pods_mod3_fit$TBOT[i]+mod3_pods_mod3_fit$TLEN[i]+5), (mod3_pods_mod3_fit$TBOT[i]+mod3_pods_mod3_fit$TLEN[i]+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
  mod3pod.mod3fit.coord[,2,i] <- 2008 - 2*mod3pod.mod3fit.coord[,1,i]
  mod3pod.mod3fit.coord[,3,i] <- c(mod3_pods_mod3_fit$NPOP08[i], (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-6))), (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-5))), (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-4))), (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-3))), (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-2))), (mod3_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod3_pods_mod3_fit$TBOT[i]-1))), mod3_pods_mod3_fit$NBOT[i], mod3_pods_mod3_fit$NBOT[i], mod3_pods_mod3_fit$NPREBOT[i], mod3_pods_mod3_fit$NPREBOT[i], mod3_pods_mod3_fit$NPREBOT[i], mod3_pods_mod3_fit$NPREBOT[i]) #haploid
  mod3pod.mod3fit.coord[,4,i] <- mod3pod.mod3fit.coord[,3,i]/2 #diploid
}

# Model 5
mod3pod.mod5fit.coord <- array(numeric(), c(10,4,nrow(mod3_pods_mod5_fit)))
for (i in 1:nrow(mod3_pods_mod5_fit)) {
  mod3pod.mod5fit.coord[,1,i] <- c((mod3_pods_mod5_fit$TBOTTWO[i]-mod3_pods_mod5_fit$TBOTTWO[i]), mod3_pods_mod5_fit$TBOTTWO[i], mod3_pods_mod5_fit$TBOTTWO[i], (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i] + mod3_pods_mod5_fit$TBOTONE[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i] + mod3_pods_mod5_fit$TBOTONE[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i] + mod3_pods_mod5_fit$TBOTONE[i] + mod3_pods_mod5_fit$TLENONE[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i] + mod3_pods_mod5_fit$TBOTONE[i] + mod3_pods_mod5_fit$TLENONE[i]), (mod3_pods_mod5_fit$TBOTTWO[i] + mod3_pods_mod5_fit$TLENTWO[i] + mod3_pods_mod5_fit$TBOTONE[i] + mod3_pods_mod5_fit$TLENONE[i] + 100)) #generations going back in time
  mod3pod.mod5fit.coord[,2,i] <- 2008 - 2*mod3pod.mod5fit.coord[,1,i]
  mod3pod.mod5fit.coord[,3,i] <- c(mod3_pods_mod5_fit$NPOP08[i], mod3_pods_mod5_fit$NPOP08[i], mod3_pods_mod5_fit$NBOTTWO[i], mod3_pods_mod5_fit$NBOTTWO[i], mod3_pods_mod5_fit$NPREBOT[i], mod3_pods_mod5_fit$NPREBOT[i], mod3_pods_mod5_fit$NBOTONE[i], mod3_pods_mod5_fit$NBOTONE[i], mod3_pods_mod5_fit$NANC[i], mod3_pods_mod5_fit$NANC[i]) #haploid
  mod3pod.mod5fit.coord[,4,i] <- mod3pod.mod5fit.coord[,3,i]/2 #diploid
}

# Model 6
# First calculate recent r (R2008) for each SFS
r <- vector(length = nrow(mod3_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod3_pods_mod6_fit$NBOT[i]/mod3_pods_mod6_fit$NPOP08[i])/(mod3_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod3pod.mod6fit.coord <- array(numeric(), c(10,4,nrow(mod3_pods_mod6_fit)))
for (i in 1:nrow(mod3_pods_mod6_fit)) {
  mod3pod.mod6fit.coord[,1,i] <- c(0, mod3_pods_mod6_fit$TBOT[i]-3, mod3_pods_mod6_fit$TBOT[i]-2, mod3_pods_mod6_fit$TBOT[i]-1, mod3_pods_mod6_fit$TBOT[i], (mod3_pods_mod6_fit$TBOT[i]+mod3_pods_mod6_fit$TLEN[i]), (mod3_pods_mod6_fit$TBOT[i]+mod3_pods_mod6_fit$TLEN[i]), (mod3_pods_mod6_fit$TBOT[i]+mod3_pods_mod6_fit$TLEN[i]+2), (mod3_pods_mod6_fit$TBOT[i]+mod3_pods_mod6_fit$TLEN[i]+5), (mod3_pods_mod6_fit$TBOT[i]+mod3_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod3pod.mod6fit.coord[,2,i] <- 2008 -2*mod3pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod3pod.mod6fit.coord[,3,i] <- c(mod3_pods_mod6_fit$NPOP08[i], (mod3_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod3_pods_mod6_fit$TBOT[i]-3))), (mod3_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod3_pods_mod6_fit$TBOT[i]-2))), (mod3_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod3_pods_mod6_fit$TBOT[i]-1))), mod3_pods_mod6_fit$NBOT[i], mod3_pods_mod6_fit$NBOT[i], mod3_pods_mod6_fit$NPREBOT[i], (mod3_pods_mod6_fit$NPREBOT[i]*exp(mod3_pods_mod6_fit$RANC[i] * 2)), (mod3_pods_mod6_fit$NPREBOT[i]*exp(mod3_pods_mod6_fit$RANC[i] * 5)), (mod3_pods_mod6_fit$NPREBOT[i]*exp(mod3_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod3pod.mod6fit.coord[,4,i] <- mod3pod.mod6fit.coord[,3,i]/2 #diploid
}

# Model 7
mod3pod.mod7fit.coord <- array(numeric(), c(4,4,nrow(mod3_pods_mod7_fit)))
for (i in 1:nrow(mod3_pods_mod7_fit)) {
  mod3pod.mod7fit.coord[,1,i] <- c(mod3_pods_mod7_fit$TCAR[i]-mod3_pods_mod7_fit$TCAR[i], mod3_pods_mod7_fit$TCAR[i], (mod3_pods_mod7_fit$TCAR[i]+5), (mod3_pods_mod7_fit$TCAR[i]+1000)) #generations going back in time
  mod3pod.mod7fit.coord[,2,i] <- 2008 -2*mod3pod.mod7fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod3pod.mod7fit.coord[,3,i] <- c(mod3_pods_mod7_fit$NPOP08[i], mod3_pods_mod7_fit$NPOP08[i], (mod3_pods_mod7_fit$NANC[i]*exp(mod3_pods_mod7_fit$RANC[i] * 5)), (mod3_pods_mod7_fit$NANC[i]*exp(mod3_pods_mod7_fit$RANC[i] * 1000))) #haploid
  mod3pod.mod7fit.coord[,4,i] <- mod3pod.mod7fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod3$X2, max.mod3$X4, xlab = '', ylab = '', type = 'n', xlim = c(1900,2008), ylim = c(0,35000), las = 1)
lines(max.mod3$X2, max.mod3$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:1) {
  lines(jitter(mod3pod.mod1fit.coord[,2,l], factor = 0.2), mod3pod.mod1fit.coord[,4,l], col = cols) # plots 1 line: Model 1 fits to Model 3 PODs
}
for (l in 1:2) {
  lines(jitter(mod3pod.mod2fit.coord[,2,l], factor = 0.2), mod3pod.mod2fit.coord[,4,l], col = cols) # plots 2 lines: Model 2 fits to Model 3 PODs
}
for (l in 1:4) {
  lines(jitter(mod3pod.mod3fit.coord[,2,l], factor = 0.2), mod3pod.mod3fit.coord[,4,l], col = cols) # plots 4 lines: Model 3 fits to Model 3 PODs
}
for (l in 1:1) {
  lines(jitter(mod3pod.mod5fit.coord[,2,l], factor = 0.2), mod3pod.mod5fit.coord[,4,l], col = cols) # plots 1 line: Model 5 fits to Model 3 PODs
}
for (l in 1:1) {
  lines(jitter(mod3pod.mod6fit.coord[,2,l], factor = 0.2), mod3pod.mod6fit.coord[,4,l], col = cols) # plots 1 line: Model 6 fits to Model 3 PODs
}
for (l in 1:1) {
  lines(jitter(mod3pod.mod7fit.coord[,2,l], factor = 0.2), mod3pod.mod7fit.coord[,4,l], col = cols) # plots 1 line: Model 7 fits to Model 3 PODs
}

#### Model 4 ####
max.mod4 <- data.frame(matrix(NA, nrow = 12, ncol = 4))
max.mod4[,1] <- c(0, mod4_ml_parameters$TBOT-6, mod4_ml_parameters$TBOT-4, mod4_ml_parameters$TBOT-3, mod4_ml_parameters$TBOT-1, mod4_ml_parameters$TBOT, mod4_ml_parameters$TBOT, (mod4_ml_parameters$TBOT+mod4_ml_parameters$TLEN), (mod4_ml_parameters$TBOT+mod4_ml_parameters$TLEN), (mod4_ml_parameters$TBOT+mod4_ml_parameters$TLEN+2), (mod4_ml_parameters$TBOT+mod4_ml_parameters$TLEN+5), (mod4_ml_parameters$TBOT+mod4_ml_parameters$TLEN+1000)) #generations going back in time
max.mod4[,2] <- 2008 - 2*max.mod4[,1]
max.mod4[,3] <- c(mod4_ml_parameters$NPOP08, mod4_ml_parameters$NPOP08, mod4_ml_parameters$NPOP08, mod4_ml_parameters$NPOP08, mod4_ml_parameters$NPOP08, mod4_ml_parameters$NPOP08, mod4_ml_parameters$NBOT, mod4_ml_parameters$NBOT, mod4_ml_parameters$NPREBOT, (mod4_ml_parameters$NPREBOT*exp(mod4_ml_parameters$RANC * 2)), (mod4_ml_parameters$NPREBOT*exp(mod4_ml_parameters$RANC * 5)), (mod4_ml_parameters$NPREBOT*exp(mod4_ml_parameters$RANC * 1000))) #haploid
max.mod4[,4] <- max.mod4[,3]/2 #diploid

# Model 4 POD fits
mod4pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod4pods", pattern = "*summary.txt", full.names = TRUE)
mod4_files <- lapply(mod4pods_names, function(x)read.table(x, header=F, fill = T))
mod4_files_array <- array(as.numeric(unlist(mod4_files)), dim=c(7, 11, 10))

mod4_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod4_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod4_files_array_sub[1,,i] <- mod4_files_array[1,c(1,4:5,2:3),i]
  mod4_files_array_sub[2,,i] <- mod4_files_array[2,c(1:3,6:7),i]
  mod4_files_array_sub[3,,i] <- mod4_files_array[3,c(1:3,6:7),i]
  mod4_files_array_sub[4,,i] <- mod4_files_array[4,c(1:3,8:9),i]
  mod4_files_array_sub[5,,i] <- mod4_files_array[5,c(1:3,10:11),i]
  mod4_files_array_sub[6,,i] <- mod4_files_array[6,c(1:3,8:9),i]
  mod4_files_array_sub[7,,i] <- mod4_files_array[7,c(1,7:8,5:6),i]
}

mod4.aic <- 2*b-2*(mod4_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod4.ass <- apply(mod4.aic,2,which.min)

# Get ML parameter for each fit POD
mod4_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod4_pods_mlfits[i,] <- mod4_files_array[mod4.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod4_pods_mod4_fit <- mod4_pods_mlfits[which(mod4.ass == 4),-c(10:11)]
colnames(mod4_pods_mod4_fit) <- colnames(mod4_ml_parameters)
mod4_pods_mod6_fit <- mod4_pods_mlfits[which(mod4.ass == 6),-c(10:11)]
colnames(mod4_pods_mod6_fit) <- colnames(mod6_ml_parameters)
mod4_pods_mod7_fit <-mod4_pods_mlfits[which(mod4.ass == 7),-c(7:11)]
colnames(mod4_pods_mod7_fit) <- colnames(mod7_ml_parameters)

# Make curves for fitted PODs
# Model 4
mod4pod.mod4fit.coord <- array(numeric(), c(12,4,nrow(mod4_pods_mod4_fit)))
for (i in 1:nrow(mod4_pods_mod4_fit)) {
  mod4pod.mod4fit.coord[,1,i] <- c(mod4_pods_mod4_fit$TBOT[i]-mod4_pods_mod4_fit$TBOT[i], mod4_pods_mod4_fit$TBOT[i]-3, mod4_pods_mod4_fit$TBOT[i]-2, mod4_pods_mod4_fit$TBOT[i]-1, mod4_pods_mod4_fit$TBOT[i], mod4_pods_mod4_fit$TBOT[i], (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]), (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]), (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]+2), (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]+5), (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]+10), (mod4_pods_mod4_fit$TBOT[i]+mod4_pods_mod4_fit$TLEN[i]+1000)) #generations going back in time
  mod4pod.mod4fit.coord[,2,i] <- 2008 -2*mod4pod.mod4fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod4pod.mod4fit.coord[,3,i] <- c(mod4_pods_mod4_fit$NPOP08[i], mod4_pods_mod4_fit$NPOP08[i], mod4_pods_mod4_fit$NPOP08[i], mod4_pods_mod4_fit$NPOP08[i], mod4_pods_mod4_fit$NPOP08[i], mod4_pods_mod4_fit$NBOT[i], mod4_pods_mod4_fit$NBOT[i], mod4_pods_mod4_fit$NPREBOT[i], (mod4_pods_mod4_fit$NPREBOT[i]*exp(mod4_pods_mod4_fit$RANC[i] * 2)), (mod4_pods_mod4_fit$NPREBOT[i]*exp(mod4_pods_mod4_fit$RANC[i] * 5)), (mod4_pods_mod4_fit$NPREBOT[i]*exp(mod4_pods_mod4_fit$RANC[i] * 10)), (mod4_pods_mod4_fit$NPREBOT[i]*exp(mod4_pods_mod4_fit$RANC[i] * 1000))) #haploid
  mod4pod.mod4fit.coord[,4,i] <- mod4pod.mod4fit.coord[,3,i]/2 #diploid
}

# Model 6
# First calculate recent r (R2008) for Model 6 fits
r <- vector(length = nrow(mod4_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod4_pods_mod6_fit$NBOT[i]/mod4_pods_mod6_fit$NPOP08[i])/(mod4_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod4pod.mod6fit.coord <- array(numeric(), c(10,4,nrow(mod4_pods_mod6_fit)))
for (i in 1:nrow(mod4_pods_mod6_fit)) {
  mod4pod.mod6fit.coord[,1,i] <- c(mod4_pods_mod6_fit$TBOT[i]-mod4_pods_mod6_fit$TBOT[i], mod4_pods_mod6_fit$TBOT[i]-3, mod4_pods_mod6_fit$TBOT[i]-2, mod4_pods_mod6_fit$TBOT[i]-1, mod4_pods_mod6_fit$TBOT[i], (mod4_pods_mod6_fit$TBOT[i]+mod4_pods_mod6_fit$TLEN[i]), (mod4_pods_mod6_fit$TBOT[i]+mod4_pods_mod6_fit$TLEN[i]), (mod4_pods_mod6_fit$TBOT[i]+mod4_pods_mod6_fit$TLEN[i]+2), (mod4_pods_mod6_fit$TBOT[i]+mod4_pods_mod6_fit$TLEN[i]+5), (mod4_pods_mod6_fit$TBOT[i]+mod4_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod4pod.mod6fit.coord[,2,i] <- 2008 -2*mod4pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod4pod.mod6fit.coord[,3,i] <- c(mod4_pods_mod6_fit$NPOP08[i], (mod4_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod4_pods_mod6_fit$TBOT[i]-3))), (mod4_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod4_pods_mod6_fit$TBOT[i]-2))), (mod4_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod4_pods_mod6_fit$TBOT[i]-1))), mod4_pods_mod6_fit$NBOT[i], mod4_pods_mod6_fit$NBOT[i], mod4_pods_mod6_fit$NPREBOT[i], (mod4_pods_mod6_fit$NPREBOT[i]*exp(mod4_pods_mod6_fit$RANC[i] * 2)), (mod4_pods_mod6_fit$NPREBOT[i]*exp(mod4_pods_mod6_fit$RANC[i] * 5)), (mod4_pods_mod6_fit$NPREBOT[i]*exp(mod4_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod4pod.mod6fit.coord[,4,i] <- mod4pod.mod6fit.coord[,3,i]/2 #diploid
}

# Model 7
mod4pod.mod7fit.coord <- array(numeric(), c(4,4,nrow(mod4_pods_mod7_fit)))
for (i in 1:nrow(mod4_pods_mod7_fit)) {
  mod4pod.mod7fit.coord[,1,i] <- c(mod4_pods_mod7_fit$TCAR[i]-mod4_pods_mod7_fit$TCAR[i], mod4_pods_mod7_fit$TCAR[i], (mod4_pods_mod7_fit$TCAR[i]+5), (mod4_pods_mod7_fit$TCAR[i]+1000)) #generations going back in time
  mod4pod.mod7fit.coord[,2,i] <- 2008 -2*mod4pod.mod7fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod4pod.mod7fit.coord[,3,i] <- c(mod4_pods_mod7_fit$NPOP08[i], mod4_pods_mod7_fit$NPOP08[i], (mod4_pods_mod7_fit$NANC[i]*exp(mod4_pods_mod7_fit$RANC[i] * 5)), (mod4_pods_mod7_fit$NANC[i]*exp(mod4_pods_mod7_fit$RANC[i] * 1000))) #haploid
  mod4pod.mod7fit.coord[,4,i] <- mod4pod.mod7fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod4$X2, max.mod4$X4, xlab = '', ylab = '', type = 'n', xlim = c(1965,2008), ylim = c(0,70000), las = 1)
lines(max.mod4$X2, max.mod4$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:2) {
  lines(jitter(mod4pod.mod4fit.coord[,2,l], factor = 0.2), mod4pod.mod4fit.coord[,4,l], col = cols) # plots 2 lines: Model 4 fits to Model 4 PODs
}
for (l in 1:4) {
  lines(jitter(mod4pod.mod6fit.coord[,2,l], factor = 0.2), mod4pod.mod6fit.coord[,4,l], col = cols) # plots 4 lines: Model 6 fits to Model 4 PODs
}
for (l in 1:4) {
  lines(jitter(mod4pod.mod7fit.coord[,2,l], factor = 0.2), mod4pod.mod7fit.coord[,4,l], col = cols) # plots 4 lines: Model 7 fits to Model 4 PODs
}

#### Model 5 ####
max.mod5 <- data.frame(matrix(NA, nrow = 10, ncol = 4))
max.mod5[,1] <- c(mod5_ml_parameters$TBOTTWO-8, mod5_ml_parameters$TBOTTWO, mod5_ml_parameters$TBOTTWO, (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO + mod5_ml_parameters$TBOTONE), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO + mod5_ml_parameters$TBOTONE), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO + mod5_ml_parameters$TBOTONE + mod5_ml_parameters$TLENONE), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO + mod5_ml_parameters$TBOTONE + mod5_ml_parameters$TLENONE), (mod5_ml_parameters$TBOTTWO + mod5_ml_parameters$TLENTWO + mod5_ml_parameters$TBOTONE + mod5_ml_parameters$TLENONE + 100)) #generations going back in time
max.mod5[,2] <- 2008 - 2*max.mod5[,1]
max.mod5[,3] <- c(mod5_ml_parameters$NPOP08, mod5_ml_parameters$NPOP08, mod5_ml_parameters$NBOTTWO, mod5_ml_parameters$NBOTTWO, mod5_ml_parameters$NPREBOT, mod5_ml_parameters$NPREBOT, mod5_ml_parameters$NBOTONE, mod5_ml_parameters$NBOTONE, mod5_ml_parameters$NANC, mod5_ml_parameters$NANC) #haploid
max.mod5[,4] <- max.mod5[,3]/2 #diploid

# Model 5 POD fits
mod5pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod5pods", pattern = "*summary.txt", full.names = TRUE)
mod5_files <- lapply(mod5pods_names, function(x)read.table(x, header=F, fill = T))
mod5_files_array <- array(as.numeric(unlist(mod5_files)), dim=c(7, 11, 10))

mod5_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod5_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod5_files_array_sub[1,,i] <- mod5_files_array[1,c(1,4:5,2:3),i]
  mod5_files_array_sub[2,,i] <- mod5_files_array[2,c(1:3,6:7),i]
  mod5_files_array_sub[3,,i] <- mod5_files_array[3,c(1:3,6:7),i]
  mod5_files_array_sub[4,,i] <- mod5_files_array[4,c(1:3,8:9),i]
  mod5_files_array_sub[5,,i] <- mod5_files_array[5,c(1:3,10:11),i]
  mod5_files_array_sub[6,,i] <- mod5_files_array[6,c(1:3,8:9),i]
  mod5_files_array_sub[7,,i] <- mod5_files_array[7,c(1,7:8,5:6),i]
}

mod5.aic <- 2*b-2*(mod5_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod5.ass <- apply(mod5.aic,2,which.min)

# Get ML parameter for each fit POD
mod5_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod5_pods_mlfits[i,] <- mod5_files_array[mod5.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod5_pods_mod2_fit <- mod5_pods_mlfits[which(mod5.ass == 2),-c(8:11)] # models best-fit to model 2
colnames(mod5_pods_mod2_fit) <- colnames(mod2_ml_parameters)
mod5_pods_mod3_fit <- mod5_pods_mlfits[which(mod5.ass == 3),-c(8:11)] # models best-fit to model 3
colnames(mod5_pods_mod3_fit) <- colnames(mod3_ml_parameters)
mod5_pods_mod4_fit <- mod5_pods_mlfits[which(mod5.ass == 4),-c(10:11)] # models best-fit to model 4
colnames(mod5_pods_mod4_fit) <- colnames(mod4_ml_parameters)
mod5_pods_mod5_fit <- mod5_pods_mlfits[which(mod5.ass == 5),] # models best-fit to model 5
colnames(mod5_pods_mod5_fit) <- colnames(mod5_ml_parameters)
mod5_pods_mod6_fit <- mod5_pods_mlfits[which(mod5.ass == 6),-c(10:11)] # models best-fit to model 6
colnames(mod5_pods_mod6_fit) <- colnames(mod6_ml_parameters)

# Make curves for fitted PODs
# Model 2
mod5pod.mod2fit.coord <- array(numeric(), c(10,4,nrow(mod5_pods_mod2_fit)))
for (i in 1:nrow(mod5_pods_mod2_fit)) {
  mod5pod.mod2fit.coord[,1,i] <- c(mod5_pods_mod2_fit$TBOT[i]-mod5_pods_mod2_fit$TBOT[i], mod5_pods_mod2_fit$TBOT[i]-5, mod5_pods_mod2_fit$TBOT[i]-3, mod5_pods_mod2_fit$TBOT[i]-1, mod5_pods_mod2_fit$TBOT[i], mod5_pods_mod2_fit$TBOT[i], (mod5_pods_mod2_fit$TBOT[i]+mod5_pods_mod2_fit$TLEN[i]), (mod5_pods_mod2_fit$TBOT[i]+mod5_pods_mod2_fit$TLEN[i]), (mod5_pods_mod2_fit$TBOT[i]+mod5_pods_mod2_fit$TLEN[i]+2), (mod5_pods_mod2_fit$TBOT[i]+mod5_pods_mod2_fit$TLEN[i]+50)) #generations going back in time
  mod5pod.mod2fit.coord[,2,i] <- 2008 - 2*mod5pod.mod2fit.coord[,1,i]
  mod5pod.mod2fit.coord[,3,i] <- c(mod5_pods_mod2_fit$NPOP08[i], mod5_pods_mod2_fit$NPOP08[i], mod5_pods_mod2_fit$NPOP08[i], mod5_pods_mod2_fit$NPOP08[i], mod5_pods_mod2_fit$NPOP08[i], mod5_pods_mod2_fit$NBOT[i], mod5_pods_mod2_fit$NBOT[i], mod5_pods_mod2_fit$PREBOT[i], mod5_pods_mod2_fit$NPREBOT[i], mod5_pods_mod2_fit$NPREBOT[i], mod5_pods_mod2_fit$NPREBOT[i]) #haploid
  mod5pod.mod2fit.coord[,4,i] <- mod5pod.mod2fit.coord[,3,i]/2 #diploid
}
  
# Model 3
r2008 <- (log(mod5_pods_mod3_fit$NBOT/mod5_pods_mod3_fit$NPOP08)/(mod5_pods_mod3_fit$TBOT)) # check that this is correct; growth rate going back in time

mod5pod.mod3fit.coord <- array(numeric(), c(13,4,nrow(mod5_pods_mod3_fit)))
for (i in 1:nrow(mod5_pods_mod3_fit)) {
  mod5pod.mod3fit.coord[,1,i] <- c(mod5_pods_mod3_fit$TBOT[i]-mod5_pods_mod3_fit$TBOT[i], mod5_pods_mod3_fit$TBOT[i]-6, mod5_pods_mod3_fit$TBOT[i]-5, mod5_pods_mod3_fit$TBOT[i]-4, mod5_pods_mod3_fit$TBOT[i]-3, mod5_pods_mod3_fit$TBOT[i]-2, mod5_pods_mod3_fit$TBOT[i]-1, mod5_pods_mod3_fit$TBOT[i], (mod5_pods_mod3_fit$TBOT[i]+mod5_pods_mod3_fit$TLEN[i]), (mod5_pods_mod3_fit$TBOT[i]+mod5_pods_mod3_fit$TLEN[i]), (mod5_pods_mod3_fit$TBOT[i]+mod5_pods_mod3_fit$TLEN[i]+2), (mod5_pods_mod3_fit$TBOT[i]+mod5_pods_mod3_fit$TLEN[i]+5), (mod5_pods_mod3_fit$TBOT[i]+mod5_pods_mod3_fit$TLEN[i]+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
  mod5pod.mod3fit.coord[,2,i] <- 2008 - 2*mod5pod.mod3fit.coord[,1,i]
  mod5pod.mod3fit.coord[,3,i] <- c(mod5_pods_mod3_fit$NPOP08[i], (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-6))), (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-5))), (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-4))), (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-3))), (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-2))), (mod5_pods_mod3_fit$NPOP08[i]*exp(r2008[i] * (mod5_pods_mod3_fit$TBOT[i]-1))), mod5_pods_mod3_fit$NBOT[i], mod5_pods_mod3_fit$NBOT[i], mod5_pods_mod3_fit$NPREBOT[i], mod5_pods_mod3_fit$NPREBOT[i], mod5_pods_mod3_fit$NPREBOT[i], mod5_pods_mod3_fit$NPREBOT[i]) #haploid
  mod5pod.mod3fit.coord[,4,i] <- mod5pod.mod3fit.coord[,3,i]/2 #diploid
}

# Model 4
mod5pod.mod4fit.coord <- array(numeric(), c(12,4,nrow(mod5_pods_mod4_fit)))
for (i in 1:nrow(mod5_pods_mod4_fit)) {
  mod5pod.mod4fit.coord[,1,i] <- c(mod5_pods_mod4_fit$TBOT[i]-mod5_pods_mod4_fit$TBOT[i], mod5_pods_mod4_fit$TBOT[i]-3, mod5_pods_mod4_fit$TBOT[i]-2, mod5_pods_mod4_fit$TBOT[i]-1, mod5_pods_mod4_fit$TBOT[i], mod5_pods_mod4_fit$TBOT[i], (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]), (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]), (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]+2), (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]+5), (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]+10), (mod5_pods_mod4_fit$TBOT[i]+mod5_pods_mod4_fit$TLEN[i]+1000)) #generations going back in time
  mod5pod.mod4fit.coord[,2,i] <- 2008 -2*mod5pod.mod4fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod5pod.mod4fit.coord[,3,i] <- c(mod5_pods_mod4_fit$NPOP08[i], mod5_pods_mod4_fit$NPOP08[i], mod5_pods_mod4_fit$NPOP08[i], mod5_pods_mod4_fit$NPOP08[i], mod5_pods_mod4_fit$NPOP08[i], mod5_pods_mod4_fit$NBOT[i], mod5_pods_mod4_fit$NBOT[i], mod5_pods_mod4_fit$NPREBOT[i], (mod5_pods_mod4_fit$NPREBOT[i]*exp(mod5_pods_mod4_fit$RANC[i] * 2)), (mod5_pods_mod4_fit$NPREBOT[i]*exp(mod5_pods_mod4_fit$RANC[i] * 5)), (mod5_pods_mod4_fit$NPREBOT[i]*exp(mod5_pods_mod4_fit$RANC[i] * 10)), (mod5_pods_mod4_fit$NPREBOT[i]*exp(mod5_pods_mod4_fit$RANC[i] * 1000))) #haploid
  mod5pod.mod4fit.coord[,4,i] <- mod5pod.mod4fit.coord[,3,i]/2 #diploid
}

# Model 5
mod5pod.mod5fit.coord <- array(numeric(), c(10,4,nrow(mod5_pods_mod5_fit)))
for (i in 1:nrow(mod5_pods_mod5_fit)) {
  mod5pod.mod5fit.coord[,1,i] <- c((mod5_pods_mod5_fit$TBOTTWO[i]-mod5_pods_mod5_fit$TBOTTWO[i]), mod5_pods_mod5_fit$TBOTTWO[i], mod5_pods_mod5_fit$TBOTTWO[i], (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i] + mod5_pods_mod5_fit$TBOTONE[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i] + mod5_pods_mod5_fit$TBOTONE[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i] + mod5_pods_mod5_fit$TBOTONE[i] + mod5_pods_mod5_fit$TLENONE[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i] + mod5_pods_mod5_fit$TBOTONE[i] + mod5_pods_mod5_fit$TLENONE[i]), (mod5_pods_mod5_fit$TBOTTWO[i] + mod5_pods_mod5_fit$TLENTWO[i] + mod5_pods_mod5_fit$TBOTONE[i] + mod5_pods_mod5_fit$TLENONE[i] + 100)) #generations going back in time
  mod5pod.mod5fit.coord[,2,i] <- 2008 - 2*mod5pod.mod5fit.coord[,1,i]
  mod5pod.mod5fit.coord[,3,i] <- c(mod5_pods_mod5_fit$NPOP08[i], mod5_pods_mod5_fit$NPOP08[i], mod5_pods_mod5_fit$NBOTTWO[i], mod5_pods_mod5_fit$NBOTTWO[i], mod5_pods_mod5_fit$NPREBOT[i], mod5_pods_mod5_fit$NPREBOT[i], mod5_pods_mod5_fit$NBOTONE[i], mod5_pods_mod5_fit$NBOTONE[i], mod5_pods_mod5_fit$NANC[i], mod5_pods_mod5_fit$NANC[i]) #haploid
  mod5pod.mod5fit.coord[,4,i] <- mod5pod.mod5fit.coord[,3,i]/2 #diploid
}

# Model 6
# First calculate recent r (R2008) for each SFS
r <- vector(length = nrow(mod5_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod5_pods_mod6_fit$NBOT[i]/mod5_pods_mod6_fit$NPOP08[i])/(mod5_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod5pod.mod6fit.coord <- array(numeric(), c(10,4,nrow(mod5_pods_mod6_fit)))
for (i in 1:nrow(mod5_pods_mod6_fit)) {
  mod5pod.mod6fit.coord[,1,i] <- c(0, mod5_pods_mod6_fit$TBOT[i]-3, mod5_pods_mod6_fit$TBOT[i]-2, mod5_pods_mod6_fit$TBOT[i]-1, mod5_pods_mod6_fit$TBOT[i], (mod5_pods_mod6_fit$TBOT[i]+mod5_pods_mod6_fit$TLEN[i]), (mod5_pods_mod6_fit$TBOT[i]+mod5_pods_mod6_fit$TLEN[i]), (mod5_pods_mod6_fit$TBOT[i]+mod5_pods_mod6_fit$TLEN[i]+2), (mod5_pods_mod6_fit$TBOT[i]+mod5_pods_mod6_fit$TLEN[i]+5), (mod5_pods_mod6_fit$TBOT[i]+mod5_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod5pod.mod6fit.coord[,2,i] <- 2008 -2*mod5pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod5pod.mod6fit.coord[,3,i] <- c(mod5_pods_mod6_fit$NPOP08[i], (mod5_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod5_pods_mod6_fit$TBOT[i]-3))), (mod5_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod5_pods_mod6_fit$TBOT[i]-2))), (mod5_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod5_pods_mod6_fit$TBOT[i]-1))), mod5_pods_mod6_fit$NBOT[i], mod5_pods_mod6_fit$NBOT[i], mod5_pods_mod6_fit$NPREBOT[i], (mod5_pods_mod6_fit$NPREBOT[i]*exp(mod5_pods_mod6_fit$RANC[i] * 2)), (mod5_pods_mod6_fit$NPREBOT[i]*exp(mod5_pods_mod6_fit$RANC[i] * 5)), (mod5_pods_mod6_fit$NPREBOT[i]*exp(mod5_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod5pod.mod6fit.coord[,4,i] <- mod5pod.mod6fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod5$X2, max.mod5$X4, xlab = '', ylab = '', type = 'n', xlim = c(1976,2008), ylim = c(0,32000), las = 1)
lines(max.mod5$X2, max.mod5$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:2) {
  lines(jitter(mod5pod.mod2fit.coord[,2,l], factor = 0.2), mod5pod.mod2fit.coord[,4,l], col = cols) # plots 4 lines: Model 2 fits to Model 5 PODs
}
for (l in 1:4) {
  lines(jitter(mod5pod.mod3fit.coord[,2,l], factor = 0.2), mod5pod.mod3fit.coord[,4,l], col = cols) # plots 4 lines: Model 3 fits to Model 5 PODs
}
for (l in 1:1) {
  lines(jitter(mod5pod.mod4fit.coord[,2,l], factor = 0.2), mod5pod.mod4fit.coord[,4,l], col = cols) # plots 1 line: Model 4 fits to Model 5 PODs
}
for (l in 1:2) {
  lines(mod5pod.mod5fit.coord[,2,l], mod5pod.mod5fit.coord[,4,l], col = cols) # plots 2 line: Model 5 fits to Model 5 PODs
}
for (l in 1:1) {
  lines(jitter(mod5pod.mod6fit.coord[,2,l], factor = 0.2), mod5pod.mod6fit.coord[,4,l], col = cols) # plots 1 line: Model 6 fits to Model 5 PODs
}

#### Model 6 ####
r2008 <- (log(mod6_ml_parameters$NBOT/mod6_ml_parameters$NPOP08)/(mod6_ml_parameters$TBOT)) # check that this is correct; growth rate going back in time

max.mod6 <- data.frame(matrix(NA, nrow = 12, ncol = 4))
max.mod6[,1] <- c(0, mod6_ml_parameters$TBOT-8, mod6_ml_parameters$TBOT-6, mod6_ml_parameters$TBOT-4, mod6_ml_parameters$TBOT-3, mod6_ml_parameters$TBOT-1, mod6_ml_parameters$TBOT, (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+2), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+5), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
max.mod6[,2] <- 2008 - 2*max.mod6[,1]
max.mod6[,3] <- c(mod6_ml_parameters$NPOP08, (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-8))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-6))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-4))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-3))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-1))), mod6_ml_parameters$NBOT, mod6_ml_parameters$NBOT, mod6_ml_parameters$NPREBOT, (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 2)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 5)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 1000))) #haploid; when t is replaced by 9000 generations, NANC = 2877
max.mod6[,4] <- max.mod6[,3]/2 #diploid

# Model 6 POD fits
mod6pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod6pods", pattern = "*summary.txt", full.names = TRUE)
mod6_files <- lapply(mod6pods_names, function(x)read.table(x, header=F, fill = T))
mod6_files_array <- array(as.numeric(unlist(mod6_files)), dim=c(7, 11, 10))

mod6_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod6_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod6_files_array_sub[1,,i] <- mod6_files_array[1,c(1,4:5,2:3),i]
  mod6_files_array_sub[2,,i] <- mod6_files_array[2,c(1:3,6:7),i]
  mod6_files_array_sub[3,,i] <- mod6_files_array[3,c(1:3,6:7),i]
  mod6_files_array_sub[4,,i] <- mod6_files_array[4,c(1:3,8:9),i]
  mod6_files_array_sub[5,,i] <- mod6_files_array[5,c(1:3,10:11),i]
  mod6_files_array_sub[6,,i] <- mod6_files_array[6,c(1:3,8:9),i]
  mod6_files_array_sub[7,,i] <- mod6_files_array[7,c(1,7:8,5:6),i]
}

mod6.aic <- 2*b-2*(mod6_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod6.ass <- apply(mod6.aic,2,which.min)

# Get ML parameter for each fit POD
mod6_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod6_pods_mlfits[i,] <- mod6_files_array[mod6.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod6_pods_mod6_fit <- mod6_pods_mlfits[which(mod6.ass == 6),-c(10:11)]
colnames(mod6_pods_mod6_fit) <- colnames(mod6_ml_parameters)
mod6_pods_mod7_fit <-mod6_pods_mlfits[which(mod6.ass == 7),-c(7:11)]
colnames(mod6_pods_mod7_fit) <- colnames(mod7_ml_parameters)

# Make curves for fitted PODs
# First calculate recent r (R2008) for each SFS
# Model 6
r <- vector(length = nrow(mod6_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod6_pods_mod6_fit$NBOT[i]/mod6_pods_mod6_fit$NPOP08[i])/(mod6_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod6pod.mod6fit.coord <- array(numeric(), c(10,4,nrow(mod6_pods_mod6_fit)))
for (i in 1:nrow(mod6_pods_mod6_fit)) {
  mod6pod.mod6fit.coord[,1,i] <- c(mod6_pods_mod6_fit$TBOT[i]-mod6_pods_mod6_fit$TBOT[i], mod6_pods_mod6_fit$TBOT[i]-7, mod6_pods_mod6_fit$TBOT[i]-6, mod6_pods_mod6_fit$TBOT[i]-5, mod6_pods_mod6_fit$TBOT[i], (mod6_pods_mod6_fit$TBOT[i]+mod6_pods_mod6_fit$TLEN[i]), (mod6_pods_mod6_fit$TBOT[i]+mod6_pods_mod6_fit$TLEN[i]), (mod6_pods_mod6_fit$TBOT[i]+mod6_pods_mod6_fit$TLEN[i]+2), (mod6_pods_mod6_fit$TBOT[i]+mod6_pods_mod6_fit$TLEN[i]+5), (mod6_pods_mod6_fit$TBOT[i]+mod6_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod6pod.mod6fit.coord[,2,i] <- 2008 -2*mod6pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod6pod.mod6fit.coord[,3,i] <- c(mod6_pods_mod6_fit$NPOP08[i], (mod6_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod6_pods_mod6_fit$TBOT[i]-7))), (mod6_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod6_pods_mod6_fit$TBOT[i]-6))), (mod6_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod6_pods_mod6_fit$TBOT[i]-5))), mod6_pods_mod6_fit$NBOT[i], mod6_pods_mod6_fit$NBOT[i], mod6_pods_mod6_fit$NPREBOT[i], (mod6_pods_mod6_fit$NPREBOT[i]*exp(mod6_pods_mod6_fit$RANC[i] * 2)), (mod6_pods_mod6_fit$NPREBOT[i]*exp(mod6_pods_mod6_fit$RANC[i] * 5)), (mod6_pods_mod6_fit$NPREBOT[i]*exp(mod6_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod6pod.mod6fit.coord[,4,i] <- mod6pod.mod6fit.coord[,3,i]/2 #diploid
}

# Model 7
mod6pod.mod7fit.coord <- array(numeric(), c(4,4,nrow(mod6_pods_mod7_fit)))
for (i in 1:nrow(mod6_pods_mod7_fit)) {
  mod6pod.mod7fit.coord[,1,i] <- c(mod6_pods_mod7_fit$TCAR[i]-mod6_pods_mod7_fit$TCAR[i], mod6_pods_mod7_fit$TCAR[i], (mod6_pods_mod7_fit$TCAR[i]+5), (mod6_pods_mod7_fit$TCAR[i]+1000)) #generations going back in time
  mod6pod.mod7fit.coord[,2,i] <- 2008 -2*mod6pod.mod7fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod6pod.mod7fit.coord[,3,i] <- c(mod6_pods_mod7_fit$NPOP08[i], mod6_pods_mod7_fit$NPOP08[i], (mod6_pods_mod7_fit$NANC[i]*exp(mod6_pods_mod7_fit$RANC[i] * 5)), (mod6_pods_mod7_fit$NANC[i]*exp(mod6_pods_mod7_fit$RANC[i] * 1000))) #haploid
  mod6pod.mod7fit.coord[,4,i] <- mod6pod.mod7fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod6$X2, max.mod6$X4, xlab = '', ylab = '', type = 'n', xlim = c(1976,2008), ylim = c(0,70000), las = 1)
lines(max.mod6$X2, max.mod6$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:9) {
  lines(jitter(mod6pod.mod6fit.coord[,2,l], factor = 0.2), mod6pod.mod6fit.coord[,4,l], col = cols) # plots 9 lines: Model 6 fits to Model 6 PODs
}
for (l in 1:1){
  lines(jitter(mod6pod.mod7fit.coord[,2,l], factor = 0.2), mod6pod.mod7fit.coord[,4,l], col = cols) # plots 1 line: Model 7 fit to Model 6 POD
}

#### Model 7 ####
max.mod7 <- data.frame(matrix(NA, nrow = 4, ncol = 4))
max.mod7[,1] <- c(mod7_ml_parameters$TCAR-10, mod7_ml_parameters$TCAR, (mod7_ml_parameters$TCAR+5), (mod7_ml_parameters$TCAR+1000)) #generations going back in time
max.mod7[,2] <- 2008 - 2*max.mod7[,1]
max.mod7[,3] <- c(mod7_ml_parameters$NPOP08, mod7_ml_parameters$NPOP08, (mod7_ml_parameters$NANC*exp(mod7_ml_parameters$RANC * 5)), (mod7_ml_parameters$NANC*exp(mod7_ml_parameters$RANC * 1000))) #haploid
max.mod7[,4] <- max.mod7[,3]/2 #diploid

# Model 7 POD fits
# Model 7
mod7pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod7pods", pattern = "*summary.txt", full.names = TRUE)
mod7_files <- lapply(mod7pods_names, function(x)read.table(x, header=F, fill = T))
mod7_files_array <- array(as.numeric(unlist(mod7_files)), dim=c(7, 11, 10))

mod7_files_array_sub <- array(0, dim = c(7,5,10))
colnames(mod7_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod7_files_array_sub[1,,i] <- mod7_files_array[1,c(1,4:5,2:3),i]
  mod7_files_array_sub[2,,i] <- mod7_files_array[2,c(1:3,6:7),i]
  mod7_files_array_sub[3,,i] <- mod7_files_array[3,c(1:3,6:7),i]
  mod7_files_array_sub[4,,i] <- mod7_files_array[4,c(1:3,8:9),i]
  mod7_files_array_sub[5,,i] <- mod7_files_array[5,c(1:3,10:11),i]
  mod7_files_array_sub[6,,i] <- mod7_files_array[6,c(1:3,8:9),i]
  mod7_files_array_sub[7,,i] <- mod7_files_array[7,c(1,7:8,5:6),i]
}

mod7.aic <- 2*b-2*(mod7_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod7.ass <- apply(mod7.aic,2,which.min)

# Get ML parameter for each fit POD
mod7_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 10))
for (i in 1:10) {
  mod7_pods_mlfits[i,] <- mod7_files_array[mod7.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod7_pods_mod6_fit <- mod7_pods_mlfits[which(mod7.ass == 6),-c(10:11)]
colnames(mod7_pods_mod6_fit) <- colnames(mod6_ml_parameters)
mod7_pods_mod7_fit <-mod7_pods_mlfits[which(mod7.ass == 7),-c(7:11)]
colnames(mod7_pods_mod7_fit) <- colnames(mod7_ml_parameters)

# Make curves for fitted PODs
# Model 6
# First calculate recent r (R2008) for each SFS
r <- vector(length = nrow(mod7_pods_mod6_fit))
for (i in 1:length(r)) {
  r[i] <- (log(mod7_pods_mod6_fit$NBOT[i]/mod7_pods_mod6_fit$NPOP08[i])/(mod7_pods_mod6_fit$TBOT[i])) # check a few to make sure they're correct
}

mod7pod.mod6fit.coord <- array(numeric(), c(12,4,nrow(mod7_pods_mod6_fit)))
for (i in 1:nrow(mod7_pods_mod6_fit)) {
  mod7pod.mod6fit.coord[,1,i] <- c(mod7_pods_mod6_fit$TBOT[i]-mod7_pods_mod6_fit$TBOT[i], mod7_pods_mod6_fit$TBOT[i]-10, mod7_pods_mod6_fit$TBOT[i]-8, mod7_pods_mod6_fit$TBOT[i]-6, mod7_pods_mod6_fit$TBOT[i]-4, mod7_pods_mod6_fit$TBOT[i]-1, mod7_pods_mod6_fit$TBOT[i], (mod7_pods_mod6_fit$TBOT[i]+mod7_pods_mod6_fit$TLEN[i]), (mod7_pods_mod6_fit$TBOT[i]+mod7_pods_mod6_fit$TLEN[i]), (mod7_pods_mod6_fit$TBOT[i]+mod7_pods_mod6_fit$TLEN[i]+2), (mod7_pods_mod6_fit$TBOT[i]+mod7_pods_mod6_fit$TLEN[i]+5), (mod7_pods_mod6_fit$TBOT[i]+mod7_pods_mod6_fit$TLEN[i]+1000)) #generations going back in time
  mod7pod.mod6fit.coord[,2,i] <- 2008 -2*mod7pod.mod6fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod7pod.mod6fit.coord[,3,i] <- c(mod7_pods_mod6_fit$NPOP08[i], (mod7_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod7_pods_mod6_fit$TBOT[i]-10))), (mod7_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod7_pods_mod6_fit$TBOT[i]-8))), (mod7_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod7_pods_mod6_fit$TBOT[i]-6))), (mod7_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod7_pods_mod6_fit$TBOT[i]-4))), (mod7_pods_mod6_fit$NPOP08[i]*exp(r[i] * (mod7_pods_mod6_fit$TBOT[i]-1))), mod7_pods_mod6_fit$NBOT[i], mod7_pods_mod6_fit$NBOT[i], mod7_pods_mod6_fit$NPREBOT[i], (mod7_pods_mod6_fit$NPREBOT[i]*exp(mod7_pods_mod6_fit$RANC[i] * 2)), (mod7_pods_mod6_fit$NPREBOT[i]*exp(mod7_pods_mod6_fit$RANC[i] * 5)), (mod7_pods_mod6_fit$NPREBOT[i]*exp(mod7_pods_mod6_fit$RANC[i] * 1000))) #haploid
  mod7pod.mod6fit.coord[,4,i] <- mod7pod.mod6fit.coord[,3,i]/2 #diploid
}

# Model 7
mod7pod.mod7fit.coord <- array(numeric(), c(4,4,nrow(mod7_pods_mod7_fit)))
for (i in 1:nrow(mod7_pods_mod7_fit)) {
  mod7pod.mod7fit.coord[,1,i] <- c(mod7_pods_mod7_fit$TCAR[i]-mod7_pods_mod7_fit$TCAR[i], mod7_pods_mod7_fit$TCAR[i], (mod7_pods_mod7_fit$TCAR[i]+5), (mod7_pods_mod7_fit$TCAR[i]+1000)) #generations going back in time
  mod7pod.mod7fit.coord[,2,i] <- 2008 -2*mod7pod.mod7fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod7pod.mod7fit.coord[,3,i] <- c(mod7_pods_mod7_fit$NPOP08[i], mod7_pods_mod7_fit$NPOP08[i], (mod7_pods_mod7_fit$NANC[i]*exp(mod7_pods_mod7_fit$RANC[i] * 5)), (mod7_pods_mod7_fit$NANC[i]*exp(mod7_pods_mod7_fit$RANC[i] * 1000))) #haploid
  mod7pod.mod7fit.coord[,4,i] <- mod7pod.mod7fit.coord[,3,i]/2 #diploid
}

# Plot true (ML model) & inferred model for each of 10 PODs
plot(max.mod7$X2, max.mod7$X4, xlab = '', ylab = '', type = 'n', xlim = c(1950,2008), ylim = c(0,70000), las = 1)
lines(max.mod7$X2, max.mod7$X4, lwd = 1.8)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:2) {
  lines(jitter(mod7pod.mod6fit.coord[,2,l], factor = 0.2), mod7pod.mod6fit.coord[,4,l], col = cols) # plots 2 lines: Model 6 fits to Model 7 PODs
}
for (l in 1:8) {
lines(jitter(mod7pod.mod7fit.coord[,2,l], factor = 0.2), mod7pod.mod7fit.coord[,4,l], col = cols) # plots 8 line: Model 7 fits to Model 7 POD
}  
mtext('Year', 1, 2.5, cex = 1.2)


dev.off()

