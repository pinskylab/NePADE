#### This script is similar to 03b_plot_observedSFS.R, but compares the observed SFSs to those simulated under the best-fit models ####

library(data.table)
library(wesanderson)

#### Read in observed SFS for each larval cohort ####
# Okay, now read in single population observed SFSs for plotting
pop08.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop0.obs', skip = 1, header = TRUE)
pop97.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop1.obs', skip = 1, header = TRUE)
pop94.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne279_1068loci.res/Ne279_1068loci_MAFpop2.obs', skip = 1, header = TRUE)

hist(as.numeric(pop08.obs)) # 981 polymorphic snps
hist(as.numeric(pop97.obs)) # 941 polymorphic snps
hist(as.numeric(pop94.obs)) # 551 polymorphic snps

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns. This makes plotting easier.
pop08.obs.prop <- t(pop08.obs[-1]/sum(pop08.obs[-1])) # denominator is number of SNPs (only polymorphic sites)
pop97.obs.prop <- t(pop97.obs[-1]/sum(pop97.obs[-1]))
pop94.obs.prop <- t(pop94.obs[-1]/sum(pop94.obs[-1]))

n <- max(length(pop08.obs[-1]), length(pop97.obs[-1]), length(pop94.obs[-1])) # determines max vector length (300) and then makes all shorter vectors 300 when including only polymorphic sites
pop08.obs.prop <- as.numeric(pop08.obs.prop)
length(pop97.obs.prop) <- n # adds NAs to the end of the vector
length(pop94.obs.prop) <- n # adds NAs to the end of the vector

pop97.obs.prop[is.na(pop97.obs.prop)] <- 0 # replaces NAs with zeros
pop94.obs.prop[is.na(pop94.obs.prop)] <- 0

m <- rbind(pop94.obs.prop, pop97.obs.prop, pop08.obs.prop)
colnames(m) <- 1:n # only polymorphic sites

#### Comparing averaged SFSs under different demographic scenarios to the observed SFS ####
# Reading in individual population SFSs simulated for Model 4 (exponential growth with instantaneous bottleneck and recovery)
pop08.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop0_sfs_summary.txt') #100x301
pop97.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop1_sfs_summary.txt') #100x207
pop94.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop2_sfs_summary.txt') #100x53

# Take means across each bin of # of minor alleles (column-wise)
pop08.mod4.avg <- colMeans(pop08.mod4)
pop97.mod4.avg <- colMeans(pop97.mod4)
pop94.mod4.avg <- colMeans(pop94.mod4)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod4.avg.poly <- mean(rowSums(pop08.mod4[-1])) #avg number of polymorphic snps for pop08
pop97.mod4.avg.poly <- mean(rowSums(pop97.mod4[-1])) #avg number of polymorphic snps for pop97
pop94.mod4.avg.poly <- mean(rowSums(pop94.mod4[-1])) #avg number of polymorphic snps for pop94

pop08.mod4.prop <- t(pop08.mod4.avg[-1]/pop08.mod4.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod4.prop <- t(pop97.mod4.avg[-1]/pop97.mod4.avg.poly)
pop94.mod4.prop <- t(pop94.mod4.avg[-1]/pop94.mod4.avg.poly)

n <- max(length(pop08.mod4.avg[-1]), length(pop97.mod4.avg[-1]), length(pop94.mod4.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod4.prop <- as.numeric(pop08.mod4.prop)
length(pop97.mod4.prop) <- n # adds NAs to the end of the vector
length(pop94.mod4.prop) <- n # adds NAs to the end of the vector

pop97.mod4.prop[is.na(pop97.mod4.prop)] <- 0 # Replaces NAs with 0
pop94.mod4.prop[is.na(pop94.mod4.prop)] <- 0

m <- rbind(pop94.mod4.prop, pop97.mod4.prop, pop08.mod4.prop)
colnames(m) <- 1:300

# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop0_sfs_summary.txt') #100x301
pop97.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop1_sfs_summary.txt') #100x207
pop94.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop2_sfs_summary.txt') #100x53

# Take means across each bin of # of minor alleles
pop08.mod6.avg <- colMeans(pop08.mod6)
pop97.mod6.avg <- colMeans(pop97.mod6)
pop94.mod6.avg <- colMeans(pop94.mod6)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod6.avg.poly <- mean(rowSums(pop08.mod6[-1])) #avg number of polymorphic snps for pop08
pop97.mod6.avg.poly <- mean(rowSums(pop97.mod6[-1])) #avg number of polymorphic snps for pop97
pop94.mod6.avg.poly <- mean(rowSums(pop94.mod6[-1])) #avg number of polymorphic snps for pop94

pop08.mod6.prop <- t(pop08.mod6.avg[-1]/pop08.mod6.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.prop <- t(pop97.mod6.avg[-1]/pop97.mod6.avg.poly)
pop94.mod6.prop <- t(pop94.mod6.avg[-1]/pop94.mod6.avg.poly)

n <- max(length(pop08.mod6.avg[-1]), length(pop97.mod6.avg[-1]), length(pop94.mod6.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.prop <- as.numeric(pop08.mod6.prop)
length(pop97.mod6.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.prop) <- n # adds NAs to the end of the vector

pop97.mod6.prop[is.na(pop97.mod6.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.prop[is.na(pop94.mod6.prop)] <- 0

m <- rbind(pop94.mod6.prop, pop97.mod6.prop, pop08.mod6.prop)
colnames(m) <- 1:300

# Reading in individual population SFSs simulated under Model 7 (ancestral growth before leveling off to some carrying capacity)
pop08.mod7 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 7/mod7_pop0_sfs_summary.txt') #100x301
pop97.mod7 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 7/mod7_pop1_sfs_summary.txt') #100x207
pop94.mod7 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 7/mod7_pop2_sfs_summary.txt') #100x53

# Take means across each bin of # of minor alleles
pop08.mod7.avg <- colMeans(pop08.mod7)
pop97.mod7.avg <- colMeans(pop97.mod7)
pop94.mod7.avg <- colMeans(pop94.mod7)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod7.avg.poly <- mean(rowSums(pop08.mod7[-1])) #avg number of polymorphic snps for pop08
pop97.mod7.avg.poly <- mean(rowSums(pop97.mod7[-1])) #avg number of polymorphic snps for pop97
pop94.mod7.avg.poly <- mean(rowSums(pop94.mod7[-1])) #avg number of polymorphic snps for pop94

pop08.mod7.prop <- t(pop08.mod7.avg[-1]/pop08.mod7.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod7.prop <- t(pop97.mod7.avg[-1]/pop97.mod7.avg.poly)
pop94.mod7.prop <- t(pop94.mod7.avg[-1]/pop94.mod7.avg.poly)

n <- max(length(pop08.mod7.avg[-1]), length(pop97.mod7.avg[-1]), length(pop94.mod7.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod7.prop <- as.numeric(pop08.mod7.prop)
length(pop97.mod7.prop) <- n # adds NAs to the end of the vector
length(pop94.mod7.prop) <- n # adds NAs to the end of the vector

pop97.mod7.prop[is.na(pop97.mod7.prop)] <- 0 # Replaces NAs with 0
pop94.mod7.prop[is.na(pop94.mod7.prop)] <- 0

m <- rbind(pop94.mod7.prop, pop97.mod7.prop, pop08.mod7.prop)
colnames(m) <- 1:300


#### Plots to compare demographic scenarios with observed ####
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

early <- rbind(pop94.mod4.prop, pop94.mod6.prop, pop94.mod7.prop, t(pop94.obs.prop))
colnames(early) <- 1:300 #1:52
mid <- rbind(pop97.mod4.prop, pop97.mod6.prop, pop97.mod7.prop, t(pop97.obs.prop))
colnames(mid) <- 1:300 #1:206
late <- rbind(pop08.mod4.prop, pop08.mod6.prop, pop08.mod7.prop, t(pop08.obs.prop))
colnames(late) <- 1:300

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/sfs.comparison.top3.png",width=12, height=9, res=300, units="in")
par(mfrow = c(3,1),
    mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Model 4, 6 & 7 plus observed
barplot(early, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1994-1995 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.5), las = 1)
barplot(mid, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1997-1998 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.30), las = 1)
barplot(late, legend = c('Model 4', 'Model 6', 'Model 7', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '2008-2009 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.20), las = 1)

dev.off()
