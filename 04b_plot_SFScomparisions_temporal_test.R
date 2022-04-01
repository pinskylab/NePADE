#### This script is similar to 03b_plot_observedSFS.R, but compares the observed SFSs to those simulated under the best-fit models ####

library(data.table)
library(wesanderson)

#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 10K loci ####
# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6.60fish_10Kloci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_10Kloci/mod6_pop0_sfs_summary.txt') #100x41
pop97.mod6.60fish_10Kloci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_10Kloci/mod6_pop1_sfs_summary.txt') #100x41
pop94.mod6.60fish_10Kloci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_10Kloci/mod6_pop2_sfs_summary.txt') #100x41

# Take means across each bin of # of minor alleles
pop08.mod6.60fish_10Kloci.avg <- colMeans(pop08.mod6.60fish_10Kloci)
pop97.mod6.60fish_10Kloci.avg <- colMeans(pop97.mod6.60fish_10Kloci)
pop94.mod6.60fish_10Kloci.avg <- colMeans(pop94.mod6.60fish_10Kloci)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod6.60fish_10Kloci.avg.poly <- mean(rowSums(pop08.mod6.60fish_10Kloci[-1])) #avg number of polymorphic snps for pop08 across all simulated SFS
pop97.mod6.60fish_10Kloci.avg.poly <- mean(rowSums(pop97.mod6.60fish_10Kloci[-1])) #avg number of polymorphic snps for pop97 across all simulated SFS
pop94.mod6.60fish_10Kloci.avg.poly <- mean(rowSums(pop94.mod6.60fish_10Kloci[-1])) #avg number of polymorphic snps for pop94 across all simulated SFS

pop08.mod6.60fish_10Kloci.avg.poly.sd <- sd(rowSums(pop08.mod6.60fish_10Kloci[-1]))
pop97.mod6.60fish_10Kloci.avg.poly.sd <- sd(rowSums(pop97.mod6.60fish_10Kloci[-1]))
pop94.mod6.60fish_10Kloci.avg.poly.sd <- sd(rowSums(pop94.mod6.60fish_10Kloci[-1]))

pop08.mod6.60fish_10Kloci.prop <- t(pop08.mod6.60fish_10Kloci.avg[-1]/pop08.mod6.60fish_10Kloci.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.60fish_10Kloci.prop <- t(pop97.mod6.60fish_10Kloci.avg[-1]/pop97.mod6.60fish_10Kloci.avg.poly)
pop94.mod6.60fish_10Kloci.prop <- t(pop94.mod6.60fish_10Kloci.avg[-1]/pop94.mod6.60fish_10Kloci.avg.poly)

n <- max(length(pop08.mod6.60fish_10Kloci.avg[-1]), length(pop97.mod6.60fish_10Kloci.avg[-1]), length(pop94.mod6.60fish_10Kloci.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.60fish_10Kloci.prop <- as.numeric(pop08.mod6.60fish_10Kloci.prop)
length(pop97.mod6.60fish_10Kloci.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.60fish_10Kloci.prop) <- n # adds NAs to the end of the vector

pop97.mod6.60fish_10Kloci.prop[is.na(pop97.mod6.60fish_10Kloci.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.60fish_10Kloci.prop[is.na(pop94.mod6.60fish_10Kloci.prop)] <- 0

m <- rbind(pop94.mod6.60fish_10Kloci.prop, pop97.mod6.60fish_10Kloci.prop, pop08.mod6.60fish_10Kloci.prop)
colnames(m) <- 1:40

#### Plots averaged SFSs simulated using ML parameters from Model 6 (best model) across sampling year/cohort
# Need to first read in the simulated SFSs for Model 6 and calculate the proportion of polymorphic SNPs for each bin (starts at line 63)
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/exp_sfs_model6_temporal_test_60fish_10Kloci.png", width=11, height=3, res=300, units="in")

par(
  mar=c(4.5, 5, 1.5, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12
)

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = 'Averaged SFS based on Model 6 ML parameter values', xlim = c(0, 87), col = col.palette[1:3])

dev.off()

#### 1068 loci and 20 diploids ####
#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 1068 loci ####
# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6.60fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_1068loci/mod6_pop0_sfs_summary.txt') #100x41
pop97.mod6.60fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_1068loci/mod6_pop1_sfs_summary.txt') #100x41
pop94.mod6.60fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/60fish_1068loci/mod6_pop2_sfs_summary.txt') #100x41

# Take means across each bin of # of minor alleles
pop08.mod6.60fish_1068loci.avg <- colMeans(pop08.mod6.60fish_1068loci)
pop97.mod6.60fish_1068loci.avg <- colMeans(pop97.mod6.60fish_1068loci)
pop94.mod6.60fish_1068loci.avg <- colMeans(pop94.mod6.60fish_1068loci)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod6.60fish_1068loci.avg.poly <- mean(rowSums(pop08.mod6.60fish_1068loci[-1])) #avg number of polymorphic snps for pop08 across all simulated SFS
pop97.mod6.60fish_1068loci.avg.poly <- mean(rowSums(pop97.mod6.60fish_1068loci[-1])) #avg number of polymorphic snps for pop97 across all simulated SFS
pop94.mod6.60fish_1068loci.avg.poly <- mean(rowSums(pop94.mod6.60fish_1068loci[-1])) #avg number of polymorphic snps for pop94 across all simulated SFS

pop08.mod6.60fish_1068loci.avg.poly.sd <- sd(rowSums(pop08.mod6.60fish_1068loci[-1]))
pop97.mod6.60fish_1068loci.avg.poly.sd <- sd(rowSums(pop97.mod6.60fish_1068loci[-1]))
pop94.mod6.60fish_1068loci.avg.poly.sd <- sd(rowSums(pop94.mod6.60fish_1068loci[-1]))

pop08.mod6.60fish_1068loci.prop <- t(pop08.mod6.60fish_1068loci.avg[-1]/pop08.mod6.60fish_1068loci.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.60fish_1068loci.prop <- t(pop97.mod6.60fish_1068loci.avg[-1]/pop97.mod6.60fish_1068loci.avg.poly)
pop94.mod6.60fish_1068loci.prop <- t(pop94.mod6.60fish_1068loci.avg[-1]/pop94.mod6.60fish_1068loci.avg.poly)

n <- max(length(pop08.mod6.60fish_1068loci.avg[-1]), length(pop97.mod6.60fish_1068loci.avg[-1]), length(pop94.mod6.60fish_1068loci.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.60fish_1068loci.prop <- as.numeric(pop08.mod6.60fish_1068loci.prop)
length(pop97.mod6.60fish_1068loci.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.60fish_1068loci.prop) <- n # adds NAs to the end of the vector

pop97.mod6.60fish_1068loci.prop[is.na(pop97.mod6.60fish_1068loci.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.60fish_1068loci.prop[is.na(pop94.mod6.60fish_1068loci.prop)] <- 0

m <- rbind(pop94.mod6.60fish_1068loci.prop, pop97.mod6.60fish_1068loci.prop, pop08.mod6.60fish_1068loci.prop)
colnames(m) <- 1:40

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = 'Averaged SFS based on Model 6 ML parameter values', xlim = c(0, 88), col = col.palette[1:3])



plot(c(53,207,301), c(pop94.mod6.avg.poly, pop97.mod6.avg.poly, pop08.mod6.avg.poly), xlab= 'Number of alleles sampled', ylab = 'Avg number of polymorphic sites', xlim = c(20,310))
points(c(40,40,40), c(pop94.mod6.60fish_1068loci.avg.poly, pop97.mod6.60fish_1068loci.avg.poly, pop08.mod6.60fish_1068loci.avg.poly), col = 'tomato')

plot(c(53,207,301), c(pop94.mod6.avg[2], pop97.mod6.avg[2], pop08.mod6.avg[2]), xlab= 'Number of alleles sampled', ylab = 'Avg minor allele count 1', xlim = c(20,310), ylim = c(180,400))
points(c(40,40,40), c(pop94.mod6.60fish_1068loci.avg[2], pop97.mod6.60fish_1068loci.avg[2], pop08.mod6.60fish_1068loci.avg[2]), col = 'tomato')

#### 1068 loci and 200 diploids ####
#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 1068 loci ####
# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6.600fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/600fish_1068loci/mod6_pop0_sfs_summary.txt') #100x401
pop97.mod6.600fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/600fish_1068loci/mod6_pop1_sfs_summary.txt') #100x401
pop94.mod6.600fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/600fish_1068loci/mod6_pop2_sfs_summary.txt') #100x401

# Take means across each bin of # of minor alleles
pop08.mod6.600fish_1068loci.avg <- colMeans(pop08.mod6.600fish_1068loci)
pop97.mod6.600fish_1068loci.avg <- colMeans(pop97.mod6.600fish_1068loci)
pop94.mod6.600fish_1068loci.avg <- colMeans(pop94.mod6.600fish_1068loci)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod6.600fish_1068loci.avg.poly <- mean(rowSums(pop08.mod6.600fish_1068loci[-1])) #avg number of polymorphic snps for pop08 across all simulated SFS
pop97.mod6.600fish_1068loci.avg.poly <- mean(rowSums(pop97.mod6.600fish_1068loci[-1])) #avg number of polymorphic snps for pop97 across all simulated SFS
pop94.mod6.600fish_1068loci.avg.poly <- mean(rowSums(pop94.mod6.600fish_1068loci[-1])) #avg number of polymorphic snps for pop94 across all simulated SFS

pop08.mod6.600fish_1068loci.avg.poly.sd <- sd(rowSums(pop08.mod6.600fish_1068loci[-1]))
pop97.mod6.600fish_1068loci.avg.poly.sd <- sd(rowSums(pop97.mod6.600fish_1068loci[-1]))
pop94.mod6.600fish_1068loci.avg.poly.sd <- sd(rowSums(pop94.mod6.600fish_1068loci[-1]))

pop08.mod6.600fish_1068loci.prop <- t(pop08.mod6.600fish_1068loci.avg[-1]/pop08.mod6.600fish_1068loci.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.600fish_1068loci.prop <- t(pop97.mod6.600fish_1068loci.avg[-1]/pop97.mod6.600fish_1068loci.avg.poly)
pop94.mod6.600fish_1068loci.prop <- t(pop94.mod6.600fish_1068loci.avg[-1]/pop94.mod6.600fish_1068loci.avg.poly)

n <- max(length(pop08.mod6.600fish_1068loci.avg[-1]), length(pop97.mod6.600fish_1068loci.avg[-1]), length(pop94.mod6.600fish_1068loci.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.600fish_1068loci.prop <- as.numeric(pop08.mod6.600fish_1068loci.prop)
length(pop97.mod6.600fish_1068loci.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.600fish_1068loci.prop) <- n # adds NAs to the end of the vector

pop97.mod6.600fish_1068loci.prop[is.na(pop97.mod6.600fish_1068loci.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.600fish_1068loci.prop[is.na(pop94.mod6.600fish_1068loci.prop)] <- 0

m <- rbind(pop94.mod6.600fish_1068loci.prop, pop97.mod6.600fish_1068loci.prop, pop08.mod6.600fish_1068loci.prop)
colnames(m) <- 1:400

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = 'Averaged SFS based on Model 6 ML parameter values', xlim = c(0, 88), col = col.palette[1:3])

#### 1068 loci and 100 diploids ####
#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 1068 loci ####
# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6.300fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/300fish_1068loci/mod6_pop0_sfs_summary.txt') #100x401
pop97.mod6.300fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/300fish_1068loci/mod6_pop1_sfs_summary.txt') #100x401
pop94.mod6.300fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/300fish_1068loci/mod6_pop2_sfs_summary.txt') #100x401

# Take means across each bin of # of minor alleles
pop08.mod6.300fish_1068loci.avg <- colMeans(pop08.mod6.300fish_1068loci)
pop97.mod6.300fish_1068loci.avg <- colMeans(pop97.mod6.300fish_1068loci)
pop94.mod6.300fish_1068loci.avg <- colMeans(pop94.mod6.300fish_1068loci)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod6.300fish_1068loci.avg.poly <- mean(rowSums(pop08.mod6.300fish_1068loci[-1])) #avg number of polymorphic snps for pop08 across all simulated SFS
pop97.mod6.300fish_1068loci.avg.poly <- mean(rowSums(pop97.mod6.300fish_1068loci[-1])) #avg number of polymorphic snps for pop97 across all simulated SFS
pop94.mod6.300fish_1068loci.avg.poly <- mean(rowSums(pop94.mod6.300fish_1068loci[-1])) #avg number of polymorphic snps for pop94 across all simulated SFS

pop08.mod6.300fish_1068loci.avg.poly.sd <- sd(rowSums(pop08.mod6.300fish_1068loci[-1]))
pop97.mod6.300fish_1068loci.avg.poly.sd <- sd(rowSums(pop97.mod6.300fish_1068loci[-1]))
pop94.mod6.300fish_1068loci.avg.poly.sd <- sd(rowSums(pop94.mod6.300fish_1068loci[-1]))

pop08.mod6.300fish_1068loci.prop <- t(pop08.mod6.300fish_1068loci.avg[-1]/pop08.mod6.300fish_1068loci.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.300fish_1068loci.prop <- t(pop97.mod6.300fish_1068loci.avg[-1]/pop97.mod6.300fish_1068loci.avg.poly)
pop94.mod6.300fish_1068loci.prop <- t(pop94.mod6.300fish_1068loci.avg[-1]/pop94.mod6.300fish_1068loci.avg.poly)

n <- max(length(pop08.mod6.300fish_1068loci.avg[-1]), length(pop97.mod6.300fish_1068loci.avg[-1]), length(pop94.mod6.300fish_1068loci.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.300fish_1068loci.prop <- as.numeric(pop08.mod6.300fish_1068loci.prop)
length(pop97.mod6.300fish_1068loci.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.300fish_1068loci.prop) <- n # adds NAs to the end of the vector

pop97.mod6.300fish_1068loci.prop[is.na(pop97.mod6.300fish_1068loci.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.300fish_1068loci.prop[is.na(pop94.mod6.300fish_1068loci.prop)] <- 0

m <- rbind(pop94.mod6.300fish_1068loci.prop, pop97.mod6.300fish_1068loci.prop, pop08.mod6.300fish_1068loci.prop)
colnames(m) <- 1:200

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = 'Averaged SFS based on Model 6 ML parameter values', xlim = c(0, 88), col = col.palette[1:3])




plot(c(53,207,301), c(pop94.mod6.avg.poly, pop97.mod6.avg.poly, pop08.mod6.avg.poly), xlab= 'Number of alleles sampled', ylab = 'Avg number of polymorphic sites', xlim = c(20,410))
points(c(40,40,40), c(pop94.mod6.60fish_1068loci.avg.poly, pop97.mod6.60fish_1068loci.avg.poly, pop08.mod6.60fish_1068loci.avg.poly), col = 'tomato')
points(c(400,400,400), c(pop94.mod6.600fish_1068loci.avg.poly, pop97.mod6.600fish_1068loci.avg.poly, pop08.mod6.600fish_1068loci.avg.poly), col = 'gold')
points(c(200,200,200), c(pop94.mod6.300fish_1068loci.avg.poly, pop97.mod6.300fish_1068loci.avg.poly, pop08.mod6.300fish_1068loci.avg.poly), col = 'purple')


plot(c(53,207,301), c(pop94.mod6.avg[2], pop97.mod6.avg[2], pop08.mod6.avg[2]), xlab= 'Number of alleles sampled', ylab = 'Avg minor allele count 1', xlim = c(20,410), ylim = c(130,350))
points(c(40,40,40), c(pop94.mod6.60fish_1068loci.avg[2], pop97.mod6.60fish_1068loci.avg[2], pop08.mod6.60fish_1068loci.avg[2]), col = 'tomato')
points(c(400,400,400), c(pop94.mod6.600fish_1068loci.avg[2], pop97.mod6.600fish_1068loci.avg[2], pop08.mod6.600fish_1068loci.avg[2]), col = 'gold')
points(c(200,200,200), c(pop94.mod6.300fish_1068loci.avg[2], pop97.mod6.300fish_1068loci.avg[2], pop08.mod6.300fish_1068loci.avg[2]), col = 'purple')



###########











#### Read in observed SFS for each larval cohort ####
# Okay, now read in single population observed SFSs for plotting
# pop08.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop0.obs', skip = 1, header = TRUE)
# pop97.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop1.obs', skip = 1, header = TRUE)
# pop94.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop2.obs', skip = 1, header = TRUE)
# 
# hist(as.numeric(pop08.obs)) # 981 polymorphic snps
# hist(as.numeric(pop97.obs)) # 941 polymorphic snps
# hist(as.numeric(pop94.obs)) # 551 polymorphic snps
# 
# # All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns. This makes plotting easier.
# pop08.obs.prop <- t(pop08.obs[-1]/sum(pop08.obs[-1])) # denominator is number of SNPs (only polymorphic sites)
# pop97.obs.prop <- t(pop97.obs[-1]/sum(pop97.obs[-1]))
# pop94.obs.prop <- t(pop94.obs[-1]/sum(pop94.obs[-1]))
# 
# n <- max(length(pop08.obs[-1]), length(pop97.obs[-1]), length(pop94.obs[-1])) # determines max vector length (300) and then makes all shorter vectors 300 when including only polymorphic sites
# pop08.obs.prop <- as.numeric(pop08.obs.prop)
# length(pop97.obs.prop) <- n # adds NAs to the end of the vector
# length(pop94.obs.prop) <- n # adds NAs to the end of the vector
# 
# pop97.obs.prop[is.na(pop97.obs.prop)] <- 0 # replaces NAs with zeros
# pop94.obs.prop[is.na(pop94.obs.prop)] <- 0
# 
# m <- rbind(pop94.obs.prop, pop97.obs.prop, pop08.obs.prop)
# colnames(m) <- 1:n # only polymorphic sites

#### Plots to compare demographic scenarios with observed ####
# col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
# palette(col.palette)
# 
# early <- rbind(pop94.mod4.prop, pop94.mod6.prop, pop94.mod7.prop, t(pop94.obs.prop))
# colnames(early) <- 1:300 #1:52
# mid <- rbind(pop97.mod4.prop, pop97.mod6.prop, pop97.mod7.prop, t(pop97.obs.prop))
# colnames(mid) <- 1:300 #1:206
# late <- rbind(pop08.mod4.prop, pop08.mod6.prop, pop08.mod7.prop, t(pop08.obs.prop))
# colnames(late) <- 1:300
# 
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/sfs.comparison.top3.png",width=12, height=9, res=300, units="in")
# par(mfrow = c(3,1),
#     mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=12
# )
# 
# # Model 4, 6 & 7 plus observed
# barplot(early, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1994-1995 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.5), las = 1)
# barplot(mid, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1997-1998 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.30), las = 1)
# barplot(late, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '2008-2009 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.20), las = 1)
# 
# dev.off()