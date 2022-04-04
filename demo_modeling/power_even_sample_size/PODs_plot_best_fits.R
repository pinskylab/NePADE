#### This script plots the best-fit model for Model 6 with population size as a line. This best fit model was then used to generate 50 Model 6 PODs, and Models 1-7 were fit to each POD. Then line plots of the best-fit model for each Model 6 POD were plotted ####
# 1. Reads in the parameters for Model 6 resulting from 50 fastsimcoal runs
# 2. Determines best fit model and creates a data frame of points to plot as a line
# 3. Reads in the parameters resulting from fitting Models 1-7 to the Model 6 PODS
# 4. Determines the best-fit model for each Model 6 POD
# 5. Creates an array of points to plot each best-fit model to a Model 6 POD as a line
# 6. Plots points from steps 2 and 5, which shows the 'true' Model 6 model and the histories inferred from 50 Model 6 PODs

# Read in fsc data to select ML parameters for Model 6 and plot
# Read in the estimated parameters for each model
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.bestlhoods.summary.txt", header = TRUE)

# Now find the ML run for each model
# Model 6, exponential change in pop size before and after bottleneck
mod6_ml <- max(mod6$MaxEstLhood)
mod6_ml_parameters <- mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]

#### Determine best-fit model for each POD ####
b <- c(1, 5, 5, 6, 9, 6, 3) # number of parameters in each of the models going from model 1 to model 7

#### Plot ML parameters for Model 6 as a line. Then plot the ML model for each of the 50 Model 6 PODS ####
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power_even_sample_size/PODs_fit_even_sampling.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 6.5, 1.5, 1.5), # panel margin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

#### Model 6 ####
r2008 <- (log(mod6_ml_parameters$NBOT/mod6_ml_parameters$NPOP08)/(mod6_ml_parameters$TBOT)) # check that this is correct; growth rate going back in time

max.mod6 <- data.frame(matrix(NA, nrow = 12, ncol = 4))
max.mod6[,1] <- c(0, mod6_ml_parameters$TBOT-8, mod6_ml_parameters$TBOT-6, mod6_ml_parameters$TBOT-4, mod6_ml_parameters$TBOT-3, mod6_ml_parameters$TBOT-1, mod6_ml_parameters$TBOT, (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+2), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+5), (mod6_ml_parameters$TBOT+mod6_ml_parameters$TLEN+1000)) #generations going back in time, 5 and 1000 generations pre-bottleneck were chosen arbitrarily for plotting
max.mod6[,2] <- 2008 - 2*max.mod6[,1]
max.mod6[,3] <- c(mod6_ml_parameters$NPOP08, (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-8))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-6))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-4))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-3))), (mod6_ml_parameters$NPOP08*exp(r2008 * (mod6_ml_parameters$TBOT-1))), mod6_ml_parameters$NBOT, mod6_ml_parameters$NBOT, mod6_ml_parameters$NPREBOT, (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 2)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 5)), (mod6_ml_parameters$NPREBOT*exp(mod6_ml_parameters$RANC * 1000))) #haploid; when t is replaced by 9000 generations, NANC = 2877
max.mod6[,4] <- max.mod6[,3]/2 #diploid

# Model 6 POD fits
mod6pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power_even_sample_size/", pattern = "*summary.txt", full.names = TRUE)
mod6_files <- lapply(mod6pods_names, function(x)read.table(x, header=F, fill = T))
mod6_files_array <- array(as.numeric(unlist(mod6_files)), dim=c(7, 11, 50))

mod6_files_array_sub <- array(0, dim = c(7,5,50))
colnames(mod6_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:50) {
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
mod6_pods_mlfits <- data.frame(matrix(NA, ncol = 11, nrow = 50))
for (i in 1:50) {
  mod6_pods_mlfits[i,] <- mod6_files_array[mod6.ass[i],,i]
}

# Separate data frame best-fit model to POD
mod6_pods_mod1_fit <- mod6_pods_mlfits[which(mod6.ass == 1),-c(4:11)]
colnames(mod6_pods_mod1_fit) <- c("NPOP08", "MaxEstLhood", "MaxObsLhood")
mod6_pods_mod4_fit <-mod6_pods_mlfits[which(mod6.ass == 4),-c(10:11)]
colnames(mod6_pods_mod4_fit) <- c("NPOP08", "NPREBOT", "NBOT", "TBOT", "TLEN", "NANC", "RANC", "MaxEstLhood", "MaxObsLhood")
mod6_pods_mod6_fit <- mod6_pods_mlfits[which(mod6.ass == 6),-c(10:11)]
colnames(mod6_pods_mod6_fit) <- colnames(mod6_ml_parameters)
mod6_pods_mod7_fit <-mod6_pods_mlfits[which(mod6.ass == 7),-c(7:11)]
colnames(mod6_pods_mod7_fit) <- c("NPOP08", "TCAR", "NANC", "RANC", "MaxEstLhood", "MaxObsLhood")

# Make curves for fitted PODs
# Model 4
mod6pod.mod4fit.coord <- array(numeric(), c(12,4,nrow(mod6_pods_mod4_fit)))
for (i in 1:nrow(mod6_pods_mod4_fit)) {
  mod6pod.mod4fit.coord[,1,i] <- c(mod6_pods_mod4_fit$TBOT[i]-mod6_pods_mod4_fit$TBOT[i], mod6_pods_mod4_fit$TBOT[i]-3, mod6_pods_mod4_fit$TBOT[i]-2, mod6_pods_mod4_fit$TBOT[i]-1, mod6_pods_mod4_fit$TBOT[i], mod6_pods_mod4_fit$TBOT[i], (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]), (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]), (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]+2), (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]+5), (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]+10), (mod6_pods_mod4_fit$TBOT[i]+mod6_pods_mod4_fit$TLEN[i]+1000)) #generations going back in time
  mod6pod.mod4fit.coord[,2,i] <- 2008 -2*mod6pod.mod4fit.coord[,1,i] #convert to years assuming summer flounder generation time is 2 years
  mod6pod.mod4fit.coord[,3,i] <- c(mod6_pods_mod4_fit$NPOP08[i], mod6_pods_mod4_fit$NPOP08[i], mod6_pods_mod4_fit$NPOP08[i], mod6_pods_mod4_fit$NPOP08[i], mod6_pods_mod4_fit$NPOP08[i], mod6_pods_mod4_fit$NBOT[i], mod6_pods_mod4_fit$NBOT[i], mod6_pods_mod4_fit$NPREBOT[i], (mod6_pods_mod4_fit$NPREBOT[i]*exp(mod6_pods_mod4_fit$RANC[i] * 2)), (mod6_pods_mod4_fit$NPREBOT[i]*exp(mod6_pods_mod4_fit$RANC[i] * 5)), (mod6_pods_mod4_fit$NPREBOT[i]*exp(mod6_pods_mod4_fit$RANC[i] * 10)), (mod6_pods_mod4_fit$NPREBOT[i]*exp(mod6_pods_mod4_fit$RANC[i] * 1000))) #haploid
  mod6pod.mod4fit.coord[,4,i] <- mod6pod.mod4fit.coord[,3,i]/2 #diploid
}

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

# Plot true (ML model) & inferred model for each of 50 PODs
# 50 PODs first
plot(max.mod6$X2, max.mod6$X4, xlab = '', ylab = '', type = 'n', xlim = c(1940,2008), ylim = c(0,70000), las = 1)
mtext('Year', 1, 2.5, cex = 1.2)
mtext(expression(italic('N'[e])), 2, 3.7, cex = 1.2)

cols <- adjustcolor('gray70', alpha.f = 0.5)
for (l in 1:1) {
  lines(jitter(mod6pod.mod1fit.coord[,2,l], factor = 0.2), mod6pod.mod1fit.coord[,4,l], col = cols) # plots 1 line: Model 1 fits to Model 6 PODs
}
for (l in 2:2) {
  lines(jitter(mod6pod.mod4fit.coord[,2,l], factor = 0.2), mod6pod.mod4fit.coord[,4,l], col = cols) # plots 2 lines: Model 4 fits to Model 6 PODs
}
for (l in 1:43) {
  lines(jitter(mod6pod.mod6fit.coord[,2,l], factor = 0.2), jitter(mod6pod.mod6fit.coord[,4,l]), col = cols) # plots 43 lines: Model 6 fits to Model 6 PODs
}
for (l in 1:4){
  lines(jitter(mod6pod.mod7fit.coord[,2,l], factor = 0.2), mod6pod.mod7fit.coord[,4,l], col = cols) # plots 1 line: Model 7 fit to Model 6 POD
}

# True model on top
lines(max.mod6$X2, max.mod6$X4, lwd = 1.8)

dev.off()

