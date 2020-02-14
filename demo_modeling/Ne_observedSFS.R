library(ade4)
library(adegenet)
library(hierfstat)
library(pegas)
library(zvau) # for function that writes genind object to GENEPOP format

#### A filter of MAF > 0.05 has been applied to these data ####
ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140

pops <- as.data.frame(ne_data@pop)
data <- as.data.frame(ne_data@tab)

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
dim(data_sub) # 285 x 2056

# Check to see if all SNPs having no missing data & write str file 
test <- sapply(data_sub, function(y) sum(length(which(is.na(y)))))
test <- data.frame(test) # all zeros
test <- data.matrix(test)
summary(test)

# Remove single locus with 3 alleles because I think it's causing a problem
data_sub <- data_sub[, -c(141:143, 1230:1232, 1299:1301, 1556:1558, 1809:1811, 1932:1934)] # Removes 3 alleles for SNP86, 833, 1202, 1252, 1358, 1433,1467,1531,1534,1639,1871,1911,1961,1989,2006,2040,2078,2163,2252,2256,2311,2369,2386,2464,2470,2524,2737,2794,2853,2939,3012,3131 
dim(data_sub) # 285 x 2038, 1019 loci

data_sub <- as.genind(data_sub)
data_sub@pop <- as.factor(pops[,1])

writeGenPop(data_sub, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1019loci_complete.gen", comment = '1019 loci with no missing data across 285 PADE')

#### No MAF filter has been applied to these data ####
ne_data_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140

pops <- as.data.frame(ne_data_nomaf@pop)
data <- as.data.frame(ne_data_nomaf@tab)

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
dim(data_sub) # 285 x 2632

# Check to see if all SNPs having no missing data & write str file 
test <- sapply(data_sub, function(y) sum(length(which(is.na(y)))))
test <- data.frame(test) # all zeros
test <- data.matrix(test)
summary(test)

# Remove single loci with 3 alleles because I think it's causing a problem
data_sub <- data_sub[, -c(57:59, 338:340, 349:351, 448:450, 589:591, 690:692, 751:753, 800:802, 805:807, 890:892, 893:895, 932:934, 1119:1121, 1192:1194, 1253:1255, 1322:1324, 1531:1533, 1538:1540, 1551:1553, 1556:1558, 1583:1585, 1620:1622, 1627:1629, 1726:1728, 1795:1797, 1830:1832, 1845:1847, 1914:1916, 1919:1921, 1994:1996, 2075:2077, 2122:2124, 2129:2131, 2136:2138, 2143:2145, 2268:2270, 2271:2273, 2432:2434, 2511:2513, 2570:2572)] # Removes loci with 3 alleles 
dim(data_sub) # 285 x 2512, 1256 loci

data_sub <- as.genind(data_sub)
data_sub@pop <- as.factor(pops[,1])

writeGenPop(data_sub, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1256loci_complete.gen", comment = '1256 loci with no missing data across 285 PADE, no MAF')

#### Opportunities to read in and visualize single, pairwise or multiSFSs ####
# Read in multiSFS file
library(data.table)

maf05.msfs <- fread("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing.res/stable_fixing_MSFS.obs", skip = 2)
dim(maf05.msfs)
sum(maf05.msfs[1,])
table(t(maf05.msfs))

nomaf.msfs <- fread("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing_nomaf2.res/stable_fixing_nomaf2_MSFS.obs", skip = 2)
dim(nomaf.msfs)
sum(nomaf.msfs[1,])
table(t(nomaf.msfs))
which(nomaf.msfs == 40) # maf = 0.0070
which(nomaf.msfs == 53) # maf = 0.0508

# Read in single population SFSs
pop08 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing_nomaf2.res/stable_fixing_nomaf2_MAFpop0.obs', skip = 1, header = TRUE)
pop97 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing_nomaf2.res/stable_fixing_nomaf2_MAFpop1.obs', skip = 1, header = TRUE)
pop94 <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/stable_fixing_nomaf2.res/stable_fixing_nomaf2_MAFpop2.obs', skip = 1, header = TRUE)

hist(as.numeric(pop08))
hist(as.numeric(pop97))
hist(as.numeric(pop94))

# All the populations are different sizes, so need to convert to proportion of SNPs & then add zeros so that all cohorts have same number of columns
pop08.prop <- t(pop08/sum(pop08)) # sum of all these is 1256, the number of SNPs
pop97.prop <- t(pop97/sum(pop97))
pop94.prop <- t(pop94/sum(pop94))

n <- max(length(pop08), length(pop97), length(pop94)) # determines max vector length (317) and then makes all shorter vectors 317
pop08.prop <- as.numeric(pop08.prop)
length(pop97.prop) <- n # adds NAs to the end of the vector
length(pop94.prop) <- n # adds NAs to the end of the vector

pop97.prop[is.na(pop97.prop)] <- 0
pop94.prop[is.na(pop94.prop)] <- 0

m <- rbind(pop94.prop, pop97.prop, pop08.prop)
colnames(m) <- 0:316

# Plot
library(wesanderson)
col.palette <- wes_palette("FantasticFox1", 5, type = "discrete")
palette(col.palette)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/obs_sfs.png", width=11, height=3, res=300, units="in")

par(
  mar=c(4.5, 4, 1.5, 1), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12
)

barplot(m, legend = c('1994-1995 cohort', '1997-1998 cohort','2008-2009 cohort'), beside = TRUE, xlab = 'Number of minor alleles', ylab = 'Proportion of SNPs', xlim = c(0, 590), col = col.palette[1:3])

dev.off()


# This doesn't actually work because the STRUCTURE file is counting the number of alleles (0, 1 or 2)
# Write as a STRUCTURE file using the function below
# obj: genind object
# file: file name to write
# pops: whether to include population info in the file
# Function is flexible with regards to ploidy, although genotypes are
# considered to be unambiguous.
# Missing data must be recorded as NA in obj@tab.

  # example use: 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)

genind2structure <- function(obj, file="", pops=TRUE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

# Now actually make/write the STRUCTURE file
genind2structure(data_sub, file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1256loci_complete.txt", pops=TRUE)
