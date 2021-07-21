library(RColorBrewer)

#### Read in some of the estimated parameters and ML for each SFS simulated under each model ####

# Model 1
mod1pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod1pods", pattern = "*summary.txt", full.names = TRUE)
mod1_files <- lapply(mod1pods_names, function(x)read.table(x, header=F, fill = T))
mod1_files_array <- array(as.numeric(unlist(mod1_files)), dim=c(6, 11, 10))

mod1_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod1_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod1_files_array_sub[1,,i] <- mod1_files_array[1,c(1,4:5,2:3),i]
  mod1_files_array_sub[2,,i] <- mod1_files_array[2,c(1:3,6:7),i]
  mod1_files_array_sub[3,,i] <- mod1_files_array[3,c(1:3,6:7),i]
  mod1_files_array_sub[4,,i] <- mod1_files_array[4,c(1:3,8:9),i]
  mod1_files_array_sub[5,,i] <- mod1_files_array[5,c(1:3,10:11),i]
  mod1_files_array_sub[6,,i] <- mod1_files_array[6,c(1:3,8:9),i]
}

# Model 2
mod2pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod2pods", pattern = "*summary.txt", full.names = TRUE)
mod2_files <- lapply(mod2pods_names, function(x)read.table(x, header=F, fill = T))
mod2_files_array <- array(as.numeric(unlist(mod2_files)), dim=c(6, 11, 10))

mod2_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod2_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod2_files_array_sub[1,,i] <- mod2_files_array[1,c(1,4:5,2:3),i]
  mod2_files_array_sub[2,,i] <- mod2_files_array[2,c(1:3,6:7),i]
  mod2_files_array_sub[3,,i] <- mod2_files_array[3,c(1:3,6:7),i]
  mod2_files_array_sub[4,,i] <- mod2_files_array[4,c(1:3,8:9),i]
  mod2_files_array_sub[5,,i] <- mod2_files_array[5,c(1:3,10:11),i]
  mod2_files_array_sub[6,,i] <- mod2_files_array[6,c(1:3,8:9),i]
}

# Model 3
mod3pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod3pods", pattern = "*summary.txt", full.names = TRUE)
mod3_files <- lapply(mod3pods_names, function(x)read.table(x, header=F, fill = T))
mod3_files_array <- array(as.numeric(unlist(mod3_files)), dim=c(6, 11, 10))

mod3_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod3_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod3_files_array_sub[1,,i] <- mod3_files_array[1,c(1,4:5,2:3),i]
  mod3_files_array_sub[2,,i] <- mod3_files_array[2,c(1:3,6:7),i]
  mod3_files_array_sub[3,,i] <- mod3_files_array[3,c(1:3,6:7),i]
  mod3_files_array_sub[4,,i] <- mod3_files_array[4,c(1:3,8:9),i]
  mod3_files_array_sub[5,,i] <- mod3_files_array[5,c(1:3,10:11),i]
  mod3_files_array_sub[6,,i] <- mod3_files_array[6,c(1:3,8:9),i]
}

# Model 4
mod4pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod4pods", pattern = "*summary.txt", full.names = TRUE)
mod4_files <- lapply(mod4pods_names, function(x)read.table(x, header=F, fill = T))
mod4_files_array <- array(as.numeric(unlist(mod4_files)), dim=c(6, 11, 10))

mod4_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod4_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod4_files_array_sub[1,,i] <- mod4_files_array[1,c(1,4:5,2:3),i]
  mod4_files_array_sub[2,,i] <- mod4_files_array[2,c(1:3,6:7),i]
  mod4_files_array_sub[3,,i] <- mod4_files_array[3,c(1:3,6:7),i]
  mod4_files_array_sub[4,,i] <- mod4_files_array[4,c(1:3,8:9),i]
  mod4_files_array_sub[5,,i] <- mod4_files_array[5,c(1:3,10:11),i]
  mod4_files_array_sub[6,,i] <- mod4_files_array[6,c(1:3,8:9),i]
}

# Model 5
mod5pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod5pods", pattern = "*summary.txt", full.names = TRUE)
mod5_files <- lapply(mod5pods_names, function(x)read.table(x, header=F, fill = T))
mod5_files_array <- array(as.numeric(unlist(mod5_files)), dim=c(6, 11, 10))

mod5_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod5_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod5_files_array_sub[1,,i] <- mod5_files_array[1,c(1,4:5,2:3),i]
  mod5_files_array_sub[2,,i] <- mod5_files_array[2,c(1:3,6:7),i]
  mod5_files_array_sub[3,,i] <- mod5_files_array[3,c(1:3,6:7),i]
  mod5_files_array_sub[4,,i] <- mod5_files_array[4,c(1:3,8:9),i]
  mod5_files_array_sub[5,,i] <- mod5_files_array[5,c(1:3,10:11),i]
  mod5_files_array_sub[6,,i] <- mod5_files_array[6,c(1:3,8:9),i]
}

# Model 6
mod6pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/mod6pods", pattern = "*summary.txt", full.names = TRUE)
mod6_files <- lapply(mod6pods_names, function(x)read.table(x, header=F, fill = T))
mod6_files_array <- array(as.numeric(unlist(mod6_files)), dim=c(6, 11, 10))

mod6_files_array_sub <- array(0, dim = c(6,5,10))
colnames(mod6_files_array_sub) <- c('NPOP08', 'NPREBOT', 'NBOT','MaxEstLhood', 'MaxObsLhood')
for (i in 1:10) {
  mod6_files_array_sub[1,,i] <- mod6_files_array[1,c(1,4:5,2:3),i]
  mod6_files_array_sub[2,,i] <- mod6_files_array[2,c(1:3,6:7),i]
  mod6_files_array_sub[3,,i] <- mod6_files_array[3,c(1:3,6:7),i]
  mod6_files_array_sub[4,,i] <- mod6_files_array[4,c(1:3,8:9),i]
  mod6_files_array_sub[5,,i] <- mod6_files_array[5,c(1:3,10:11),i]
  mod6_files_array_sub[6,,i] <- mod6_files_array[6,c(1:3,8:9),i]
}

##### Calculate AIC for each set of PODs ####
b <- c(1, 5, 5, 6, 9, 6) # number of parameters in each of the models going from model 1 to model 6

# Model 1
mod1.aic <- 2*b-2*(mod1_files_array_sub[,4,] *2.303) # Convert from log10 to ln, first
mod1.ass <- apply(mod1.aic,2,which.min)
mod1.ass.count <- table(mod1.ass)
mod1.ass.count6 <- append(mod1.ass.count, c('4' = 0, '5' = 0, '6'= 0), 3)

# Model 2
mod2.aic <- 2*b-2*(mod2_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod2.ass <- apply(mod2.aic,2,which.min)
mod2.ass.count <- table(mod2.ass)
mod2.ass.count6 <- append(mod2.ass.count, c('1' = 0, '3' = 0, '6'= 0), 3)
mod2.ass.count6 <- mod2.ass.count6[order(names(mod2.ass.count6))]

# Model 3
mod3.aic <- 2*b-2*(mod3_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod3.ass <- apply(mod3.aic,2,which.min)
mod3.ass.count <- table(mod3.ass)
mod3.ass.count6 <- append(mod3.ass.count, c('1' = 0), 0)

# Model 4
mod4.aic <- 2*b-2*(mod4_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod4.ass <- apply(mod4.aic,2,which.min)
mod4.ass.count <- table(mod4.ass)
mod4.ass.count6 <- append(mod4.ass.count, c('1' = 0, '2' = 0, '3'= 0, '5' = 0), 2)
mod4.ass.count6 <- mod4.ass.count6[order(names(mod4.ass.count6))]

# Model 5
mod5.aic <- 2*b-2*(mod5_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod5.ass <- apply(mod5.aic,2,which.min)
mod5.ass.count <- table(mod5.ass)
mod5.ass.count6 <- append(mod5.ass.count, c('1' = 0, '3' = 0, '6'= 0), 3)
mod5.ass.count6 <- mod5.ass.count6[order(names(mod5.ass.count6))]

# Model 6
mod6.aic <- 2*b-2*(mod6_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod6.ass <- apply(mod6.aic,2,which.min)
mod6.ass.count <- table(mod6.ass)
mod6.ass.count6 <- append(mod6.ass.count, c('1' = 0, '2' = 0, '3'= 0, '5' = 0), 2)
mod6.ass.count6 <- mod6.ass.count6[order(names(mod6.ass.count6))]

#### Assemble confusion matrix ####
mat <- rbind(mod1.ass.count6, mod2.ass.count6, mod3.ass.count6, mod4.ass.count6, mod5.ass.count6, mod6.ass.count6)
rownames(mat) <- c('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6')

nocol <- 9
display.brewer.pal(n = nocol, name = 'GnBu')
color <- brewer.pal(n = nocol, name = "GnBu")
ColorLevels <- seq(0, 8, length=length(color)) # this makes the color bar nicer later on

#### Plot ####
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power/model_confusion_matrix.png", width=7, height=5, res=300, units="in")

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,5,1.5,1), font = 2)

# Plot it up!
image(1:ncol(mat), 1:nrow(mat), mat,
      col=color, xlab="True model", ylab="Estimated model",
      axes=FALSE, zlim = c(0, 8),
      main= NA, xlim = c(0.5, 6.5), ylim = c(0.5, 6.5))
box()
abline(0,1)

# Now annotate the plot
axis(side = 1, at=seq(1,6,1), labels=c(1:6),
     cex.axis=1.0)
axis(side = 2, at=seq(1,6,1), labels=c(1:6), las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(6,2.5,3.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=color,xlab="",ylab="",xaxt="n", las = 1)
mtext('Count', side = 3, cex = 1, line = 0.5, adj = 1)

dev.off()

#### Probability of recovery for the Model 6 PODs correctly ID'ed as Model 6 ####
# Model 6
mod6.recov <- vector()
for (k in 1:10) {
  if (mod6_files_array_sub[as.numeric(paste0(mod6.ass[k])),1,k]/mod6_files_array_sub[as.numeric(paste0(mod6.ass[k])),2,k] >= 1) { # if NPOP08/NPREPOT > 1, then TRUE (fully recovered); point estimate for Model 6 NPOP08/NPREBOT = 0.927
    mod6.recov[k] <- 'T'
  } else {
    mod6.recov[k] <- 'F'
  }
}


mod6.recov[which(mod6.ass == 6)]

table(mod6.recov)
table(mod6.recov[which(mod6.ass == 6)])

# Estimate of the error around that true recovery from the ML estimates of Model 6 fit to the PODS
prop.recovery <- mod6_files_array_sub[6,1,]/mod6_files_array_sub[6,2,] # NPOP08/NPREBOT

quantile(prop.recovery, c(0.025, 0.975)) # 95% CIs
quantile(prop.recovery, c(0.05, 0.95)) # 90% CIs
quantile(prop.recovery, c(0.10, 0.90)) # 80% CIs