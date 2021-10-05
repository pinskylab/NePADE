#### This script plots the genepop file resulting from 02_additional_filters.R in various ways & calculates population genetic statistics for data set composed of 3752 loci across 278 larvae (missing data allowed) ####

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)
library(diveRsity)
library(wesanderson)
library(tidyr)
library(readxl)
library(RColorBrewer)
library(plotly)
library(stringr)
library(boot)

# Read in genepop file  with population identifiers for each of the three time periods
ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_278PADE_3752loci_missingallowed.gen", ncode = 3L) # 3752 loci across 278 larvae

# Need to modify names in allele count data (genind object) so that I can match them up with names in the larval database
freq_names <- as.vector(rownames(ne_data@tab))
freq_names_split <- do.call(rbind, strsplit(as.character(freq_names), 'L'))
freq_names_split2 <- separate(as.data.frame(freq_names_split), V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
freq.newnames <- paste(freq_names_split2$name1, freq_names_split2$name3, freq_names_split2$name2, freq_names_split2$name4, sep = '')
rownames(ne_data@tab) <- freq.newnames # replaces rownames in genind object with PADEXX_XXX formatting

# Read in and create complete meta data by adding Ne fish to connectivity metadata
meta <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt", header = TRUE) # read in connectivity metadata
meta2 <- read_excel("~/Documents/Graduate School/Rutgers/Summer Flounder/Larvae/Ne/Ne_larvs.xlsx", sheet = 1) # read in Ne metadata

ne_meta <- merge(meta, as.data.frame(meta2), by = c('PinskyID', 'Sampling.Date', 'Month', 'Date', 'Year', 'Place'), all.x = TRUE, all.y = TRUE)

# Subset metadata to fish being used in Ne project
meta_sub <- ne_meta[ne_meta$PinskyID %in% rownames(ne_data$tab),]
meta_sub2 <- meta_sub[,-c(16:46)] # Get rid of a bunch of irrelevent columns

# Now order them so the names of the genind data are the same as the associated metadata
ordered_meta_sub2 <- meta_sub2[match(rownames(ne_data@tab), meta_sub2$PinskyID),] # order them so they are the same as the genind object

ordered_meta_sub2$PinskyID == rownames(ne_data@tab) # Double check names are the same

# Create a few new strata so I can color the PCA by year
pop_strata <- data.frame(cbind(ordered_meta_sub2$Year, ordered_meta_sub2$Place)) #1 = NC; 3 = NJ
strata(ne_data) <- pop_strata

#### Do PCA & plot ####
sum(is.na(ne_data$tab)) #25520
X <- scaleGen(ne_data, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plot PCA based on three time periods
col <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col)
s.class(pca1$li, pop(ne_data), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-10,50), ylim = c(-120,100))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = -2)
axis(2, at=seq(-70,80, by = 10), labels=seq(-70,80, by= 10), line = -6, las = 2)
mtext("PC1 (0.71%)", side = 1, line = 1.25)
mtext("PC2 (0.70%)", side = 2, line = -3)

legend(40, 40,
       legend=c("2008-2009 (n = 149)", "1994-1995 (n = 26)", "1997-1998 (n = 103)"),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1,
       cex = 0.8)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Plot PCA by Year
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X1, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-10,50), ylim = c(-120,100))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = -2)
axis(2, at=seq(-70,80, by = 10), labels=seq(-70,80, by= 10), line = -6, las = 2)
mtext("PC1 (0.71%)", side = 1, line = 1.25)
mtext("PC2 (0.70%)", side = 2, line = -3)

legend(40, 40,
       legend=levels(ne_data@strata$X1),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1,
       cex = 0.65)

# Plot PCA by capture Location
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X2, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-10,50), ylim = c(-120,100))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = -2)
axis(2, at=seq(-70,80, by = 10), labels=seq(-70,80, by= 10), line = -6, las = 2)
mtext("PC1 (0.71%)", side = 1, line = 1.25)
mtext("PC2 (0.70%)", side = 2, line = -3)

legend(20, 40,
       legend=c('Little Egg Inlet, NJ', 'Beaufort, NC'),
       pch=19,
       col = col,
       bty = "n",
       y.intersp = 1,
       cex = 0.7)

# Examine PC loadings. Which loci contribute most to PCs?
s.arrow(pca1$c1)
loadingplot(pca1$c1^2)

pc_names <- names(which(rowSums(pca1$c1^2) > .005)) # these are the names of the alleles that contribute the most to PCs


#######################################################################
#### Diversity metrics ####
ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_278PADE_3752loci_missingallowed.gen", ncode = 3L) # 3752 loci across 278 larvae

# Test SNPs for HWE
hwe <- hw.test(ne_data,res="matrix")
pval <- hwe[hwe[,"Pr.exact"] < 0.001,] # p<0.001, exact test
length(pval[,"Pr.exact"]) 

#### Now get genind object ready to remove SNPs not in HWE ####
cols <- colnames(ne_data@tab)
cols.split <- data.frame(do.call('rbind', strsplit(as.character(cols),'.',fixed=TRUE)))

# Make column names for loci in HWE only
cols.split.hwe <- cols.split[-which(cols.split$X1 %in% rownames(pval)),] 
cols.hwe.joined <- paste(cols.split.hwe$X1, cols.split.hwe$X2, sep = '.')

# Replace genind object column names with new names
colnames(ne_data@tab) <- cols.split$X1

# Remove loci not in HWE, add new column names & make into genind object
ne_data_hwe <- ne_data@tab[, -which(colnames(ne_data@tab) %in% rownames(pval))] 
dim(ne_data_hwe) # 284 x 7737
colnames(ne_data_hwe) <- cols.hwe.joined

ne_data_hwe_genind <- as.genind(ne_data_hwe) #284 x 3752 loci
pop(ne_data_hwe_genind) <- pop(ne_data)

# diversity table using poppr
poppr(ne_data) # He uses Nei's gene diversity

# expected heterozygosity
Hs(ne_data) # pretty similar between 3 time periods, but slightly different than using poppr function above
Hs(ne_data, ne_data@strata$X1)

# Heterozygosity across all individuals
div <- summary(ne_data)

# subset the genind object by population, so heterozygosity can be calculated by population
early.genind <- popsub(ne_data, sublist = 'PADE_95011L2524')
mid.genind <- popsub(ne_data, sublist = 'PADE_98027L2146')
late.genind <- popsub(ne_data, sublist = 'PADE_09151L2330')

early.genind <- popsub(ne_data_hwe_genind, sublist = 'PADE_95011L2524')
mid.genind <- popsub(ne_data_hwe_genind, sublist = 'PADE_98027L2146')
late.genind <- popsub(ne_data_hwe_genind, sublist = 'PADE_09151L2330')

early.he <- summary(early.genind)
mid.he <- summary(mid.genind)
late.he <- summary(late.genind)

mean(early.he$Hobs)
mean(early.he$Hexp)
mean(mid.he$Hobs)
mean(mid.he$Hexp)
mean(late.he$Hobs)
mean(late.he$Hexp)

bartlett.test(list(early.he$Hexp, early.he$Hobs)) #test for differences between Hexp and Hobs; different for earliest cohort
bartlett.test(list(mid.he$Hexp, mid.he$Hobs)) # different for middle cohort
bartlett.test(list(late.he$Hexp, late.he$Hobs)) # different for late time period

plot(div$Hobs, xlab = 'Loci number', ylab = 'Observed Heterozygosity', main = 'Observed heterozygosity per locus')
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")
bartlett.test(list(div$Hexp, div$Hobs))

# FIS by cohort
early.hier <- genind2hierfstat(early.genind) # convert to hierfstat format
early.hier$pop <- as.numeric(early.hier$pop) # convert pop factor to numeric so that bastic.stats function works
mid.hier <- genind2hierfstat(mid.genind)
mid.hier$pop <- as.numeric(mid.hier$pop)
late.hier <- genind2hierfstat(late.genind)
late.hier$pop <- as.numeric(late.hier$pop)

stats.early <- basic.stats(early.hier)
stats.mid <- basic.stats(mid.hier)
stats.late <- basic.stats(late.hier)

#### Pi calculations for 3821 loci ####
# Read in window pi calculations from vcftools on Amarel
early.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/namess94_95.windowed.pi", header = TRUE) #2090 x 5
mid.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/namess97_98.windowed.pi", header = TRUE) #3348 x 5
late.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/namess08_09.windowed.pi", header = TRUE) #3322 x 5

# This was only for trouble shooting which loci were missing from the pi calculations. Code to acutally calculate pi is below.
# Read in CHOM column from SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.vcf so I can figure out which ones are missing in the site pi files
contig_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/all_contigs3821.txt")
contig_names$index <- rownames(contig_names)

early.win.missing <- contig_names[!contig_names[,1] %in% early.win.pi$CHROM,] #names and indices of snps missing for window pi calculations
early.win.missing$index <- paste('SNP', early.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object
mid.win.missing <- contig_names[!contig_names[,1] %in% mid.win.pi$CHROM,]
mid.win.missing$index <- paste('SNP', mid.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object
late.win.missing <- contig_names[!contig_names[,1] %in% late.win.pi$CHROM,]
late.win.missing$index <- paste('SNP', late.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object

# Subset genind objects to loci for which pi wasn't calculated
# Early
early.names <- colnames(early.genind@tab) # prep names
early.names.split <- data.frame(do.call('rbind', strsplit(as.character(early.names),'.',fixed=TRUE)))
colnames(early.genind@tab) <- early.names.split$X1 # replace SNPX_X.YYY with only SNP_X

early.missing.genind <- early.genind@tab[,colnames(early.genind@tab) %in% early.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
early.missing.genind.unique <- early.missing.genind[, !duplicated(colnames(early.missing.genind))]

early.pi <- early.genind@tab[,!colnames(early.genind@tab) %in% early.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
early.pi.genind.unique <- early.pi[, !duplicated(colnames(early.pi))]
early.source.table <- matrix(NA, nrow = 2153, ncol = 4, byrow = TRUE, dimnames = list(colnames(early.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(early.pi.genind.unique)){
  early.source.table[i,] <- as.matrix(t(table(factor(early.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
early.source.table.pi <- cbind(early.source.table, early.win.pi$PI) # cbind pi values with allele counts
early.source.table.pi.ordered <- early.source.table.pi[order(-early.source.table.pi[,'2']),]
write.table(early.source.table.pi.ordered, '~/Desktop/early.source.table.pi.ordered.txt')

early.lookup.table <- matrix(NA, nrow = 1826, ncol = 4, byrow = TRUE, dimnames = list(colnames(early.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(early.missing.genind.unique)){
  early.lookup.table[i,] <- as.matrix(t(table(factor(early.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
early.lookup.table.ordered <- early.lookup.table[order(-early.lookup.table[,'2']),]
write.table(early.lookup.table.ordered, '~/Desktop/early.lookup.table.ordered.txt')

# Mid
mid.names <- colnames(mid.genind@tab) # prep names
mid.names.split <- data.frame(do.call('rbind', strsplit(as.character(mid.names),'.',fixed=TRUE)))
colnames(mid.genind@tab) <- mid.names.split$X1 # replace SNPX_X.YYY with only SNP_X

mid.missing.genind <- mid.genind@tab[,colnames(mid.genind@tab) %in% mid.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
mid.missing.genind.unique <- mid.missing.genind[, !duplicated(colnames(mid.missing.genind))]

mid.pi <- mid.genind@tab[,!colnames(mid.genind@tab) %in% mid.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
mid.pi.genind.unique <- mid.pi[, !duplicated(colnames(mid.pi))]
mid.source.table <- matrix(NA, nrow = 3475, ncol = 4, byrow = TRUE, dimnames = list(colnames(mid.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(mid.pi.genind.unique)){
  mid.source.table[i,] <- as.matrix(t(table(factor(mid.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
mid.source.table.pi <- cbind(mid.source.table, mid.win.pi$PI) # cbind pi values with allele counts
mid.source.table.pi.ordered <- mid.source.table.pi[order(-mid.source.table.pi[,'2']),]
write.table(mid.source.table.pi.ordered, '~/Desktop/mid.source.table.pi.ordered.txt')

mid.lookup.table <- matrix(NA, nrow = 504, ncol = 4, byrow = TRUE, dimnames = list(colnames(mid.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(mid.missing.genind.unique)){
  mid.lookup.table[i,] <- as.matrix(t(table(factor(mid.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
mid.lookup.table.ordered <- mid.lookup.table[order(-mid.lookup.table[,'2']),]
write.table(mid.lookup.table.ordered, '~/Desktop/mid.lookup.table.ordered.txt')

for (i in 1:ncol(mid.missing.genind.unique)){
  print(table(mid.missing.genind.unique[,i], useNA = 'ifany'))
}

# Late
late.names <- colnames(late.genind@tab) # prep names
late.names.split <- data.frame(do.call('rbind', strsplit(as.character(late.names),'.',fixed=TRUE)))
colnames(late.genind@tab) <- late.names.split$X1 # replace SNPX_X.YYY with only SNP_X

late.missing.genind <- late.genind@tab[,colnames(late.genind@tab) %in% late.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
late.missing.genind.unique <- late.missing.genind[, !duplicated(colnames(late.missing.genind))]

late.pi <- late.genind@tab[,!colnames(late.genind@tab) %in% late.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
late.pi.genind.unique <- late.pi[, !duplicated(colnames(late.pi))]
late.source.table <- matrix(NA, nrow = 3478, ncol = 4, byrow = TRUE, dimnames = list(colnames(late.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(late.pi.genind.unique)){
  late.source.table[i,] <- as.matrix(t(table(factor(late.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
late.source.table.pi <- cbind(late.source.table, late.win.pi$PI) # cbind pi values with allele counts
late.source.table.pi.ordered <- late.source.table.pi[order(-late.source.table.pi[,'2']),]

late.lookup.table <- matrix(NA, nrow = 501, ncol = 4, byrow = TRUE, dimnames = list(colnames(late.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(late.missing.genind.unique)){
  late.lookup.table[i,] <- as.matrix(t(table(factor(late.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
late.lookup.table.ordered <- late.lookup.table[order(-late.lookup.table[,'2']),]

for (i in 1:ncol(late.missing.genind.unique)){
  print(table(late.missing.genind.unique[,i], useNA = 'ifany'))
}

which(apply(early.genind@tab, 2, function(x) length(unique(x))) == 1) # these are all loci with one unique value
which(apply(mid.genind@tab, 2, function(x) length(unique(x))) == 1) 
which(apply(late.genind@tab, 2, function(x) length(unique(x))) == 1)

# Mean pi across 140bp windows
mean(early.win.pi$PI)
mean(mid.win.pi$PI)
mean(late.win.pi$PI)

# add zeros for all the loci that got dropped, then take mean
early.pi3821 <- c(early.win.pi$PI, rep(0,1731))
mid.pi3821 <- c(mid.win.pi$PI, rep(0,473))
late.pi3821 <- c(late.win.pi$PI, rep(0,499))

mean(early.pi3821)
mean(mid.pi3821)
mean(late.pi3821)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

early.boot <- boot(early.pi3821, samplemean, 1000)
plot(early.boot)
early.ci <- boot.ci(early.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
early.ci$normal # view 95% CIs

mid.boot <- boot(mid.pi3821, samplemean, 1000)
plot(mid.boot)
mid.ci <- boot.ci(mid.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.ci$normal # view 95% CIs

late.boot <- boot(late.pi3821, samplemean, 1000)
plot(late.boot)
late.ci <- boot.ci(late.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.ci$normal # view 95% CIs

# Plot average pi
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} # error bar function

options(scipen = 5)
barplot(c(early.boot$t0, mid.boot$t0, late.boot$t0), ylim = c(0,0.0006), xlab = 'Larval cohort', ylab = 'Average nucleotide diversity (π)')
error.bar(c(0.7,1.9,3.1), c(early.boot$t0, mid.boot$t0, late.boot$t0), c(early.ci$normal[3]-early.boot$t0, mid.ci$normal[3]-mid.boot$t0, late.ci$normal[3]-late.boot$t0), c(early.boot$t0-early.ci$normal[2], mid.boot$t0-mid.ci$normal[2], late.boot$t0-late.ci$normal[2]))
axis(1, at=c(0.7,1.9,3.1), labels = c('1994', '1997', '2008'))

#### Pi calculations for all SNPs on contigs ####
# Read in data containing all loci across 280 fish
data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.140trimmed.280fish.gen", ncode = 3L)

contig_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/all_contigs3821.txt")
contig_names$index <- rownames(contig_names)

# subset the genind object by population, so I can see if loci for which pi wasn't calculated are monomorphic
early.genind.all <- popsub(data, sublist = 'PADE_95011L2524') # early
mid.genind.all <- popsub(data, sublist = 'PADE_98027L2146') # mid
late.genind.all <- popsub(data, sublist = 'PADE_09151L2330') # late

# Read in window pi calculations from vcftools on Amarel
early.win.pi.all <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names94_95_allsites.windowed.pi", header = TRUE) #3685 x 5
mid.win.pi.all <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_allsites.windowed.pi", header = TRUE) #3997 x 5
late.win.pi.all <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_allsites.windowed.pi", header = TRUE) #3986 x 5

# VCFtools only outputs pi when there is at least one SNP, so sites that are invariant will get dropped from the file. Which contigs are not in vcftools pi calc file? Spot check in the cohort genind object to make sure that site is invariant
early.missing <- contig_names$V1[!contig_names$V1 %in% early.win.pi.all$CHROM] #287 contigs missing
mid.missing <- contig_names$V1[!contig_names$V1 %in% mid.win.pi.all$CHROM] #59
late.missing <- contig_names$V1[!contig_names$V1 %in% late.win.pi.all$CHROM] #75

# Keeping only one entry per contig
early.win.pi.all.140 <- early.win.pi.all[!duplicated(early.win.pi.all[,1]),] #3534 x 5
mid.win.pi.all.140 <- mid.win.pi.all[!duplicated(mid.win.pi.all[,1]),] #3762 x 5
late.win.pi.all.140 <- late.win.pi.all[!duplicated(late.win.pi.all[,1]),] #3746 x 5

# add zeros for all the loci that got dropped, then take mean
early.pi.all <- c(early.win.pi.all.140$PI, rep(0,length(contig_names$V1)-length(early.win.pi.all.140$CHROM)))
mid.pi.all <- c(mid.win.pi.all.140$PI, rep(0,length(contig_names$V1)-length(mid.win.pi.all.140$CHROM)))
late.pi.all <- c(late.win.pi.all.140$PI, rep(0,length(contig_names$V1)-length(late.win.pi.all.140$CHROM)))

mean(early.pi.all)
mean(mid.pi.all)
mean(late.pi.all)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

early.boot <- boot(early.pi.all, samplemean, 1000)
plot(early.boot)
early.ci <- boot.ci(early.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
early.ci$normal # view 95% CIs

mid.boot <- boot(mid.pi.all, samplemean, 1000)
plot(mid.boot)
mid.ci <- boot.ci(mid.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.ci$normal # view 95% CIs

late.boot <- boot(late.pi.all, samplemean, 1000)
plot(late.boot)
late.ci <- boot.ci(late.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.ci$normal # view 95% CIs

# How different is pi if all 140 bp windows are kept?
# add zeros for all the loci that got dropped, then take mean
early.pi.all140 <- c(early.win.pi.all$PI, rep(0,length(early.missing))) #3972
mid.pi.all140 <- c(mid.win.pi.all$PI, rep(0,length(mid.missing))) #4056
late.pi.all140 <- c(late.win.pi.all$PI, rep(0,length(late.missing))) #4061

mean(early.pi.all140)
mean(mid.pi.all140)
mean(late.pi.all140)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

early.boot <- boot(early.pi.all140, samplemean, 1000)
plot(early.boot)
early.ci <- boot.ci(early.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
early.ci$normal # view 95% CIs

mid.boot <- boot(mid.pi.all140, samplemean, 1000)
plot(mid.boot)
mid.ci <- boot.ci(mid.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.ci$normal # view 95% CIs

late.boot <- boot(late.pi.all140, samplemean, 1000)
plot(late.boot)
late.ci <- boot.ci(late.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.ci$normal # view 95% CIs


png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/pi1296_barplots.png",width=6, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

options(scipen = 5)
barplot(c(early.boot$t0, mid.boot$t0, late.boot$t0), ylim = c(0,0.005), xlab = 'Larval cohort', ylab = 'Average nucleotide diversity (π)')
error.bar(c(0.7,1.9,3.1), c(early.boot$t0, mid.boot$t0, late.boot$t0), c(early.ci$normal[3]-early.boot$t0, mid.ci$normal[3]-mid.boot$t0, late.ci$normal[3]-late.boot$t0), c(early.boot$t0-early.ci$normal[2], mid.boot$t0-mid.ci$normal[2], late.boot$t0-late.ci$normal[2]))
axis(1, at=c(0.7,1.9,3.1), labels = c('1994', '1997', '2008'))

dev.off()

####################
# FST
ne_data
larvs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref_alltrimmed140/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.txt", sep="\t", header = TRUE) # STRUCTURE formatted input

even_indexes<-seq(2,570,2)
odd_indexes<-seq(1,569,2)

odds <- data.frame(larvs[odd_indexes,]) # 285 x 3135
odds2 <- odds[,-c(1:2)] # 285 x 3133
evens <- data.frame(larvs[even_indexes,]) # 285 x 3135
evens2 <- evens[,-c(1:2)] # 285 x 3133

s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 285, ncol = 3133))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 285 x 3133

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop)) #check to make sure correct number of individuals in each cohort

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 285 x 3134

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')

options("scipen"=100, "digits"=4) # forces output to not be in scientific notation

# Tajima's D

# Nucleotide diversity
early <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref_alltrimmed140/names94_95_pi.sites.pi", header = TRUE)
mid <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref_alltrimmed140/names97_98_pi.sites.pi", header = TRUE)
late <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref_alltrimmed140/names08_09_pi.sites.pi", header = TRUE)

contigs <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref_alltrimmed140/contig_names.txt", header = FALSE)

early.excluded <- contigs[!contigs[,1] %in% as.character(early$CHROM),] # contigs/SNPs for which site-pi was not calculated
mid.excluded <- contigs[!contigs[,1] %in% as.character(mid$CHROM),]
late.excluded <- contigs[!contigs[,1] %in% as.character(late$CHROM),]


hist(early$PI, col=rgb(1,0,0,0.5), main = "", xlab = "Log per-site pi (π)", ylim = c(0,900))
hist(mid$PI, add = TRUE, col=rgb(0,0,1,0.5))
hist(late$PI, add = TRUE, col=rgb(0,1,0,0.5))
legend('topright', legend = c('1994-1995', '1997-1998', '2008-2009'), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(0,1,0,0.5)), pch = 15, cex = 0.8, title = 'Cohort', text.font = 1)

#### Sort fish based on site of capture for pi calculation sensitivity analysis ####
ordered_meta_sub2 # from top of script
data <- ordered_meta_sub2[, c("PinskyID","Year","Place")]

# Subset data to figure out IDs of NJ and NC fish
names94_95_NJ <- data[which(data$Place == 'Little Egg Inlet, NJ' & data$Year <= 1995),] #24 fish
names94_95_NC <- data[which(data$Place == 'Beaufort, NC' & data$Year <= 1995),] #0

names97_98_NJ <- data[which(data$Year >= 1997 & data$Year < 2008 & data$Place == 'Little Egg Inlet, NJ'),] # order of the arguments really matters here; 85 fish
names97_98_NC <- data[which(data$Year >= 1997 & data$Year < 2008 & data$Place == 'Beaufort, NC'),] # 18

names08_09_NJ <- data[which(data$Place == 'Little Egg Inlet, NJ' & data$Year >= 2008),] #146
names08_09_NC <- data[which(data$Place == 'Beaufort, NC' & data$Year >= 2008),]

# Read in names that correspond to genetic data
genetic.names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names.txt")

genetic.names_split <- do.call(rbind, strsplit(as.character(genetic.names$V1), 'L'))
genetic.names_split2 <- separate(as.data.frame(genetic.names_split), V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
genetic.newnames <- paste(freq_names_split2$name1, freq_names_split2$name3, freq_names_split2$name2, freq_names_split2$name4, sep = '') #PADEXX_XXX formatting

# Subset genetic names based on which fish came from NJ/NC and when
names97_98_NJ_genetic <- genetic.names[genetic.newnames %in% names97_98_NJ$PinskyID,]
names97_98_NC_genetic <- genetic.names[genetic.newnames %in% names97_98_NC$PinskyID,]

names08_09_NJ_genetic <- genetic.names[genetic.newnames %in% names08_09_NJ$PinskyID,]
names08_09_NC_genetic <- genetic.names[genetic.newnames %in% names08_09_NC$PinskyID,]

# Split the genetic names
write.table(names97_98_NJ_genetic, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_NJ.txt", col.names = FALSE, row.names = FALSE)
write.table(names97_98_NC_genetic, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_NC.txt", col.names = FALSE, row.names = FALSE)
write.table(names08_09_NJ_genetic, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_NJ.txt", col.names = FALSE, row.names = FALSE)
write.table(names08_09_NC_genetic, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_NC.txt", col.names = FALSE, row.names = FALSE)

# Read in window pi calculations for NJ/NC and cohort combos calculated by vcftools on Amarel
mid.nj.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_NJ.windowed.pi", header = TRUE) #3353 x 5
mid.nc.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_NC.windowed.pi", header = TRUE) #2047 x 5
late.nj.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_NJ.windowed.pi", header = TRUE) #3449x 5
late.nc.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_NC.windowed.pi", header = TRUE) #1690 x 5

mean(mid.nj.win.pi$PI)
mean(mid.nc.win.pi$PI)
mean(late.nj.win.pi$PI)
mean(late.nc.win.pi$PI)

# add zeros for all the loci that got dropped, then take mean
mid.nj.win.pi3979 <- c(mid.nj.win.pi$PI, rep(0,626))
mid.nc.win.pi3979 <- c(mid.nc.win.pi$PI, rep(0,1932))
late.nj.win.pi3979 <- c(late.nj.win.pi$PI, rep(0,530))
late.nc.win.pi3979 <- c(late.nc.win.pi$PI, rep(0,2289))

mean(mid.nj.win.pi3979)
mean(mid.nc.win.pi3979)
mean(late.nj.win.pi3979)
mean(late.nc.win.pi3979)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

mid.nj.boot <- boot(mid.nj.win.pi3979, samplemean, 1000)
plot(mid.nj.boot)
mid.nj.ci <- boot.ci(mid.nj.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.nj.ci$normal # view 95% CIs

mid.nc.boot <- boot(mid.nc.win.pi3979, samplemean, 1000)
plot(mid.nc.boot)
mid.nc.ci <- boot.ci(mid.nc.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.nc.ci$normal # view 95% CIs

late.nj.boot <- boot(late.nj.win.pi3979, samplemean, 1000)
plot(late.nj.boot)
late.nj.ci <- boot.ci(late.nj.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.nj.ci$normal # view 95% CIs

late.nc.boot <- boot(late.nc.win.pi3979, samplemean, 1000)
plot(late.nc.boot)
late.nc.ci <- boot.ci(late.nc.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.nc.ci$normal # view 95% CIs


