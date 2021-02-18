library(ade4)
library(adegenet)
library(hierfstat)
library(pegas)
library(zvau) # for function that writes genind object to GENEPOP format, or if having trouble installing, load necessary function
library(data.table)
library(wesanderson)

# Needed function
writeGenPop <- function(gi, file.name, comment) {
  
  if (is.list(gi)) {
    # do all genind objects have the same number of loci?
    if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
    gi.char <- gi[[1]]
    loc.names <- locNames(gi[[1]])
  } else {
    gi.char <- gi
    loc.names <- locNames(gi)
  }
  
  # Calculate the length of two alleles.
  lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
  lng <- unique(nchar(lng))
  
  stopifnot(length(lng) == 1)
  
  cat(paste(comment, "\n"), file = file.name)
  cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
  
  if (is.list(gi)) {
    pop.names <- seq_len(length(gi))
  } else {
    pop.names <- popNames(gi)
  }
  
  for (i in pop.names) {
    cat("pop\n", file = file.name, append = TRUE)
    if (is.list(gi)) {
      intm <- gi[[i]]
      loc.names <- locNames(gi[[i]])
    } else {
      intm <- gi[pop(gi) == i, drop = FALSE]
    }
    ind.names <- indNames(intm)
    intm <- genind2df(intm, sep = "")
    intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
    out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
    write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
  
  return(NULL)
}

#### No MAF filter has been applied to these data: 280 larvae & 3821 loci ####
ne_280fish_3821_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.140trimmed.280fish.firstsnp.genepop.gen", ncode = 3L) # 7 PADE from Fall 2009 removed, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140

# Remove loci with 3 or 4 alleles because I think it's causing a problem
ToKeep <- which(ne_280fish_3821_nomaf@loc.n.all == 2) #3681 loci
ne_280fish_3681_nomaf_biallelic <- ne_280fish_3821_nomaf[loc = ToKeep]

pops <- as.data.frame(ne_280fish_3681_nomaf_biallelic@pop)
data <- as.data.frame(ne_280fish_3681_nomaf_biallelic@tab)

na.count <- sapply(data, function(y) sum(length(which(is.na(y)))))
na.count <- data.frame(na.count)
na.count

# Histogram of SNPs with missing data
counts <- data.matrix(na.count)
hist(counts)

na.count$names <- rownames(na.count)

# Creating a list of SNP names that have no missing data
keeps <- which(na.count$na.count == '0')
keepers <- na.count[keeps,]

# Now subsetting the data to only allels with SNP names in the keep list
data_sub <- data[,keepers$names]
dim(data_sub) # 280 x 2392

# Check to see if all SNPs having no missing data & write str file 
test <- sapply(data_sub, function(y) sum(length(which(is.na(y)))))
test <- data.frame(test) # all zeros
test <- data.matrix(test)
summary(test)

data_sub <- as.genind(data_sub)
data_sub@pop <- as.factor(pops[,1])

writeGenPop(data_sub, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1196loci_complete.gen", comment = '1196 loci with no missing data across 280 PADE, no MAF')

#### Plot observed SFS ####
nomaf.msfs <- fread("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne280_1196loci_nomaf.res/Ne280_1196loci_nomaf_MSFS.obs", skip = 2) # Read in MSFS and check that there are 1196 loci
dim(nomaf.msfs)
sum(nomaf.msfs[1,]) #1196 snps, yes
table(t(nomaf.msfs)) #the number of snps in each of the categories

# Okay, now read in single population observed SFSs for plotting
pop08.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne280_1196loci_nomaf.res/Ne280_1196loci_nomaf_MAFpop0.obs', skip = 1, header = TRUE)
pop97.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne280_1196loci_nomaf.res/Ne280_1196loci_nomaf_MAFpop1.obs', skip = 1, header = TRUE)
pop94.obs <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne280_1196loci_nomaf.res/Ne280_1196loci_nomaf_MAFpop2.obs', skip = 1, header = TRUE)

hist(as.numeric(pop08.obs)) # 1161 polymorphic snps
hist(as.numeric(pop97.obs)) # 1088 polymorphic snps
hist(as.numeric(pop94.obs)) # 685 polymorphic snps

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.obs.prop <- t(pop08.obs[-1]/sum(pop08.obs[-1])) # denominator is number of SNPs (only polymorphic sites)
pop97.obs.prop <- t(pop97.obs[-1]/sum(pop97.obs[-1]))
pop94.obs.prop <- t(pop94.obs[-1]/sum(pop94.obs[-1]))

# n <- max(length(pop08), length(pop97), length(pop94)) # determines max vector length (307) and then makes all shorter vectors 307 when including nonpolymorphic sites
n <- max(length(pop08.obs[-1]), length(pop97.obs[-1]), length(pop94.obs[-1])) # determines max vector length (306) and then makes all shorter vectors 306 when including only polymorphic sites
pop08.obs.prop <- as.numeric(pop08.obs.prop)
length(pop97.obs.prop) <- n # adds NAs to the end of the vector
length(pop94.obs.prop) <- n # adds NAs to the end of the vector

pop97.obs.prop[is.na(pop97.obs.prop)] <- 0
pop94.obs.prop[is.na(pop94.obs.prop)] <- 0

m <- rbind(pop94.obs.prop, pop97.obs.prop, pop08.obs.prop)
colnames(m) <- 1:306 # only polymorphic sites

# Plot
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/obs_sfs_polyonly.png", width=11, height=3, res=300, units="in")

par(
  mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12
)

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of SNPs', main = 'Observed SFS', xlim = c(0, 105), col = col.palette[1:3])

dev.off()

#### Comparing averaged SFSs under different demographic scenarios to the observed SFS ####
# Reading in individual population SFSs simulated under Model 2 (instantaneous bottleneck)
pop08.mod2 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 2/mod2_pop0_sfs_summary.txt') #100x307
pop97.mod2 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 2/mod2_pop1_sfs_summary.txt') #100x207
pop94.mod2 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 2/mod2_pop2_sfs_summary.txt') #100x49

# Take means across each bin of # of minor alleles
pop08.mod2.avg <- colMeans(pop08.mod2)
pop97.mod2.avg <- colMeans(pop97.mod2)
pop94.mod2.avg <- colMeans(pop94.mod2)

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.mod2.avg.poly <- mean(rowSums(pop08.mod2[-1])) #avg number of polymorphic snps for pop08
pop97.mod2.avg.poly <- mean(rowSums(pop97.mod2[-1])) #avg number of polymorphic snps for pop97
pop94.mod2.avg.poly <- mean(rowSums(pop94.mod2[-1])) #avg number of polymorphic snps for pop94

pop08.mod2.prop <- t(pop08.mod2.avg[-1]/pop08.mod2.avg.poly) #prop of polymorphic snps vs # of minor alleles
pop97.mod2.prop <- t(pop97.mod2.avg[-1]/pop97.mod2.avg.poly)
pop94.mod2.prop <- t(pop94.mod2.avg[-1]/pop94.mod2.avg.poly)

n <- max(length(pop08.mod2.avg[-1]), length(pop97.mod2.avg[-1]), length(pop94.mod2.avg[-1])) # determines max vector length of only polymorphic sites (306) and then makes all shorter vectors 306
pop08.mod2.prop <- as.numeric(pop08.mod2.prop)
length(pop97.mod2.prop) <- n # adds NAs to the end of the vector
length(pop94.mod2.prop) <- n # adds NAs to the end of the vector

pop97.mod2.prop[is.na(pop97.mod2.prop)] <- 0 # Replaces NAs with 0
pop94.mod2.prop[is.na(pop94.mod2.prop)] <- 0

m <- rbind(pop94.mod2.prop, pop97.mod2.prop, pop08.mod2.prop)
colnames(m) <- 1:306

# Reading in individual population SFSs simulated for Model 4 (exponential growth with instantaneous bottleneck and recovery)
pop08.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop0_sfs_summary.txt') #100x307
pop97.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop1_sfs_summary.txt') #100x207
pop94.mod4 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 4/mod4_pop2_sfs_summary.txt') #100x49

# Take means across each bin of # of minor alleles
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

n <- max(length(pop08.mod4.avg[-1]), length(pop97.mod4.avg[-1]), length(pop94.mod4.avg[-1])) # determines max vector length of only polymorphic sites (306) and then makes all shorter vectors 306
pop08.mod4.prop <- as.numeric(pop08.mod4.prop)
length(pop97.mod4.prop) <- n # adds NAs to the end of the vector
length(pop94.mod4.prop) <- n # adds NAs to the end of the vector

pop97.mod4.prop[is.na(pop97.mod4.prop)] <- 0 # Replaces NAs with 0
pop94.mod4.prop[is.na(pop94.mod4.prop)] <- 0

m <- rbind(pop94.mod4.prop, pop97.mod4.prop, pop08.mod4.prop)
colnames(m) <- 1:306

# Reading in individual population SFSs simulated for Model 6 (exponential change in population size before and after bottleneck)
pop08.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop0_sfs_summary.txt') #100x307
pop97.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop1_sfs_summary.txt') #100x207
pop94.mod6 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/model 6/mod6_pop2_sfs_summary.txt') #100x49

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

n <- max(length(pop08.mod6.avg[-1]), length(pop97.mod6.avg[-1]), length(pop94.mod6.avg[-1])) # determines max vector length of only polymorphic sites (306) and then makes all shorter vectors 306
pop08.mod6.prop <- as.numeric(pop08.mod6.prop)
length(pop97.mod6.prop) <- n # adds NAs to the end of the vector
length(pop94.mod6.prop) <- n # adds NAs to the end of the vector

pop97.mod6.prop[is.na(pop97.mod6.prop)] <- 0 # Replaces NAs with 0
pop94.mod6.prop[is.na(pop94.mod6.prop)] <- 0

m <- rbind(pop94.mod6.prop, pop97.mod6.prop, pop08.mod6.prop)
colnames(m) <- 1:306

#### Plots to compare demographic scenarios with observed ####
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

# early <- rbind(pop94.mod2.prop, pop94.mod4.prop, t(pop94.obs.prop))
# colnames(early) <- 1:306 #1:48
# mid <- rbind(pop97.mod2.prop, pop97.mod4.prop, t(pop97.obs.prop))
# colnames(mid) <- 1:306 #1:206
# late <- rbind(pop08.mod2.prop, pop08.mod4.prop, t(pop08.obs.prop))
# colnames(late) <- 1:306

early <- rbind(pop94.mod2.prop, pop94.mod4.prop, pop94.mod6.prop, t(pop94.obs.prop))
colnames(early) <- 1:306 #1:48
mid <- rbind(pop97.mod2.prop, pop97.mod4.prop, pop97.mod6.prop, t(pop97.obs.prop))
colnames(mid) <- 1:306 #1:206
late <- rbind(pop08.mod2.prop, pop08.mod4.prop, pop08.mod6.prop, t(pop08.obs.prop))
colnames(late) <- 1:306


png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/sim_sfs/sfs.comparison.top3.png",width=12, height=9, res=300, units="in")
par(mfrow = c(3,1),
    mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# Model 2 & 4, plus observed
# barplot(early, legend = c('Model 2', 'Model 4', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1994-1995 cohort', xlim = c(0, 106), col = col.palette[1:3], ylim = c(0,0.5))
# barplot(mid, legend = c('Model 2','Model 4', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1997-1998 cohort', xlim = c(0, 106), col = col.palette[1:3], ylim = c(0,0.25))
# barplot(late, legend = c('Model 2','Model 4', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '2008-2009 cohort', xlim = c(0, 106), col = col.palette[1:3], ylim = c(0,0.20))
# dev.off()

# Model 2, 4 & 6, plus observed
barplot(early, legend = c('Model 2', 'Model 4', 'Model 6', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1994-1995 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.5), las = 2)
barplot(mid, legend = c('Model 2','Model 4', 'Model 6', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '1997-1998 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.25), las = 2)
barplot(late, legend = c('Model 2','Model 4', 'Model 6', 'Observed'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of polymorphic SNPs', main = '2008-2009 cohort', xlim = c(0, 133), col = col.palette[c(1:3,5)], ylim = c(0,0.20), las = 2)
dev.off()
