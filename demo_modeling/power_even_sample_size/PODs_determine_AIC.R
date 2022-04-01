library(RColorBrewer)

#### Read in the ML estimated parameters for each Model 6 POD (50 SFSs simulated using ML Model 6 parameters ####
# Model 6 PODs
mod6pods_names <- list.files("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/power_even_sample_size/", pattern = "*summary.txt", full.names = TRUE)
mod6_files <- lapply(mod6pods_names, function(x)read.table(x, header=F, fill = T))
mod6_files_array <- array(as.numeric(unlist(mod6_files)), dim=c(7, 11, 50))

mod6_files_array_sub <- array(0, dim = c(7,5,50)) # array depth is equal to number of PODs
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

##### Calculate AIC for each Model 6 POD ####
b <- c(1, 5, 5, 6, 9, 6, 3) # number of parameters in each of the models going from model 1 to model 7

# Model 6 PODS
mod6.aic <- 2*b-2*(mod6_files_array_sub[,4,]*2.303) # Convert from log10 to ln, first
mod6.ass <- apply(mod6.aic,2,which.min) # Determine best fit model for each POD
mod6.ass.count <- table(mod6.ass)
mod6.ass.count7 <- append(mod6.ass.count, c('2' = 0, '3'= 0, '5' = 0), 2)
mod6.ass.count7 <- mod6.ass.count7[order(names(mod6.ass.count7))]

# Probability of inferring Model 6 when the true model is Model 6
length(which(mod6.ass == '6'))/length(mod6.ass) # 86%
