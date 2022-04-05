#### This script is similar to 04a_plot_SFScomparisions.R, but compares the expected SFSs across cohorts when equal numbers of indivdiuals are sampled ####

library(data.table)
library(wesanderson)

#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 1068 loci ####
#### 1068 loci and 240 fish (80 fish or 160 alleles per cohort) ####
#### Comparing averaged SFSs under Model 6 with even sampling across cohorts and 1068 loci ####
# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6.240fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/240fish_1068loci/mod6_pop0_sfs_summary.txt') #100x161
pop97.mod6.240fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/240fish_1068loci/mod6_pop1_sfs_summary.txt') #100x161
pop94.mod6.240fish_1068loci <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model6_temporal_test/240fish_1068loci/mod6_pop2_sfs_summary.txt') #100x161

dim(pop08.mod6.240fish_1068loci)
dim(pop97.mod6.240fish_1068loci)
dim(pop94.mod6.240fish_1068loci)

# Take means across each bin of # of minor alleles
pop08.mod6.240fish_1068loci.avg <- colMeans(pop08.mod6.240fish_1068loci)
pop97.mod6.240fish_1068loci.avg <- colMeans(pop97.mod6.240fish_1068loci)
pop94.mod6.240fish_1068loci.avg <- colMeans(pop94.mod6.240fish_1068loci)

# Code holdover from when the populations were different sizes, but still converts to proportion of SNPs (& then add zeros so that all cohorts have same number of columns, even though this is already true in this case)
pop08.mod6.240fish_1068loci.avg.poly <- mean(rowSums(pop08.mod6.240fish_1068loci[-1])) #avg number of polymorphic snps for pop08 across all simulated SFS
pop97.mod6.240fish_1068loci.avg.poly <- mean(rowSums(pop97.mod6.240fish_1068loci[-1])) #avg number of polymorphic snps for pop97 across all simulated SFS
pop94.mod6.240fish_1068loci.avg.poly <- mean(rowSums(pop94.mod6.240fish_1068loci[-1])) #avg number of polymorphic snps for pop94 across all simulated SFS

pop08.mod6.240fish_1068loci.avg.poly.sd <- sd(rowSums(pop08.mod6.240fish_1068loci[-1]))
pop97.mod6.240fish_1068loci.avg.poly.sd <- sd(rowSums(pop97.mod6.240fish_1068loci[-1]))
pop94.mod6.240fish_1068loci.avg.poly.sd <- sd(rowSums(pop94.mod6.240fish_1068loci[-1]))

pop08.mod6.240fish_1068loci.prop <- t(pop08.mod6.240fish_1068loci.avg[-1]/pop08.mod6.240fish_1068loci.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod6.240fish_1068loci.prop <- t(pop97.mod6.240fish_1068loci.avg[-1]/pop97.mod6.240fish_1068loci.avg.poly)
pop94.mod6.240fish_1068loci.prop <- t(pop94.mod6.240fish_1068loci.avg[-1]/pop94.mod6.240fish_1068loci.avg.poly)

n <- max(length(pop08.mod6.240fish_1068loci.avg[-1]), length(pop97.mod6.240fish_1068loci.avg[-1]), length(pop94.mod6.240fish_1068loci.avg[-1])) # determines max vector length of only polymorphic sites (300) and then makes all shorter vectors 300 for easier plotting
pop08.mod6.240fish_1068loci.prop <- as.numeric(pop08.mod6.240fish_1068loci.prop)
length(pop97.mod6.240fish_1068loci.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.240fish_1068loci.prop) <- n # adds NAs to the end of the vector

pop97.mod6.240fish_1068loci.prop[is.na(pop97.mod6.240fish_1068loci.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.240fish_1068loci.prop[is.na(pop94.mod6.240fish_1068loci.prop)] <- 0

m <- rbind(pop94.mod6.240fish_1068loci.prop, pop97.mod6.240fish_1068loci.prop, pop08.mod6.240fish_1068loci.prop)
colnames(m) <- 1:160

col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/exp_sfs_model6_temporal_test_240fish_1068loci.png", width=11, height=3, res=300, units="in")

par(
  mar=c(4.5, 5, 1.5, 1), # panel margin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12
)

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = 'Averaged SFS based on Model 6 ML parameter values', xlim = c(0, 106), col = col.palette[1:3])

dev.off()

