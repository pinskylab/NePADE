setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE")

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
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # wants .gen extension, file is the same as that with .txt extension; 3 digits used to denote nucleotide
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/SNP.DP3g95maf05.FIL.FIL.recode.trimmed144.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset because the read lengths of the different sequencing runs are different
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140

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
sum(is.na(ne_data$tab)) #20362 
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
s.class(pca1$li, pop(ne_data), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-22,50), ylim = c(-50,65))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = 2)
axis(2, at=seq(-30,50, by = 10), labels=seq(-30,50, by= 10), line = 0, las = 2)
mtext("PC1 (0.80%)", side = 1, line = 4)
mtext("PC2 (0.66%)", side = 2, line = 2.5)

legend(20, 40,
       legend=c("2008-2009 (n = 158)", "1994-1995 (n = 24)", "1997-1998 (n = 103)"),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1,
       cex = 0.65)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Plot PCA by Year
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X1, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-22,50), ylim = c(-50,60))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = 2)
axis(2, at=seq(-30,50, by = 10), labels=seq(-30,50, by= 10), line = 0, las = 2)
mtext("PC1 (0.80%)", side = 1, line = 4)
mtext("PC2 (0.66%)", side = 2, line = 2.5)

legend(40, 40,
       legend=levels(ne_data@strata$X1),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1,
       cex = 0.65)

# Plot PCA by capture Location
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X2, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-22,50), ylim = c(-50,60))
axis(1, at=seq(-20,60, by=10), labels=seq(-20,60, by= 10), line = 2)
axis(2, at=seq(-30,50, by = 10), labels=seq(-30,50, by= 10), line = 0, las = 2)
mtext("PC1 (0.80%)", side = 1, line = 4)
mtext("PC2 (0.66%)", side = 2, line = 2.5)

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

plot(ne_data@tab[,"SNP_32.110"]) # lots of the 98 larvae are clearly different from rest
plot(ne_data@tab[,"SNP_1575.100"])

for (i in 1:length(pc_names)){
  plot(ne_data@tab[, paste0(pc_names[i])])
}

#### If I take out the loci that contribute the most to the PCs, what does PCA look like? ####
ne_data.subloci <- as.genind(ne_data@tab[,!colnames(ne_data@tab) %in% pc_names])
ordered_meta_sub2$PinskyID == rownames(ne_data.subloci@tab) # Double check names are the same
pop_strata <- data.frame(cbind(ordered_meta_sub2$Year, ordered_meta_sub2$Place)) # Same as for PCA using all loci
strata(ne_data.subloci) <- pop_strata # Add strata into genind object

# Redo PCA & plot #
sum(is.na(ne_data.subloci@tab)) #10131 
Y <- scaleGen(ne_data.subloci, NA.method = "mean")
dim(Y)
class (Y)

# make PCA
pca2 <- dudi.pca(Y,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca2$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca2

# Plot PCA based on years
col <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col)
s.class(pca2$li, ne_data.subloci@strata$X1, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-60,20), ylim = c(-35,11))
axis(1, at=seq(-50,10, by=10), labels=seq(-50,10, by= 10), line = 0.5)
axis(2, at=seq(-50,20, by = 10), labels=seq(-50,20, by= 10), line = 0, las = 2)
mtext("PC1 (0.80%)", side = 1, line = 3.5)
mtext("PC2 (0.72%)", side = 2, line = 3)

legend(-30, -20,
       legend=levels(ne_data.subloci@strata$X1),
       pch=19,
       col = col,
       bty = "n",
       y.intersp = 1)

eig_percent <- round((pca2$eig/(sum(pca2$eig)))*100,2)
eig_percent [1:3]

# If I take out the individuals from other SEQ runs, what does the PCA look like?


#######################################################################
#### Diversity metrics ####
# Need to decide if I allow missing data or not in these calculations? Also the no MAF dataset or the 0.05 MAF dataset?
# ne_data_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1256loci_complete.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140; dataset where no MAF filter applied, and only sites with no missing data: 1256 loci
# ne_data_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140; no MAF filter applied and some missing data allowed: 3979 loci
# ne281_3721_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.281fish.firstsnp.genepop.gen", ncode = 3L) # 281 larvae and 3721 loci
ne280_3821_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.140trimmed.280fish.firstsnp.genepop.gen", ncode = 3L) # 280 larvae and 3821 loci

# what fish are different between SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.genepop.gen (285) and SNP.DP3g95nomaf.FIL.FIL.recode.140trimmed.280fish.firstsnp.genepop.gen (280)? I removed 7 fish but there is only a numerical difference of 5 fish
setdiff(rownames(ne_data_nomaf@tab), rownames(ne280_3821_nomaf@tab)) # these are the 7 fish that were removed
setdiff(rownames(ne280_3821_nomaf@tab), rownames(ne_data_nomaf@tab)) # these two fish were added

# Test SNPs for HWE
hwe <- hw.test(ne_data_nomaf,res="matrix")
pval <- hwe[hwe[,"Pr.exact"] < 0.0100,] # p<0.01, exact test
length(pval[,"Pr.exact"]) 

#### Do PCA & plot ####
sum(is.na(ne280_3821_nomaf$tab)) #26050
X <- scaleGen(ne280_3821_nomaf, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1
loadingplot(pca1$c1^2)

# Plot PCA based on three time periods
col <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col)
s.class(pca1$li, pop(ne280_3821_nomaf), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-40,100), ylim = c(-80,100))
axis(1, at=seq(-20,100, by=10), labels=seq(-20,100, by= 10), line = 2)
axis(2, at=seq(-50,60, by = 10), labels=seq(-50,60, by= 10), line = 0, las = 2)
mtext("PC1 (0.84%)", side = 1, line = 4)
mtext("PC2 (0.82%)", side = 2, line = 2)

legend(20, 40,
       legend=c("2008-2009 (n = 153)", "1994-1995 (n = 24)", "1997-1998 (n = 103)"),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 0.8,
       cex = 0.85)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# genotype accumulation curve
# genotype_curve(ne_data, sample = 100, quiet = TRUE)

# diversity table using poppr
poppr(ne280_3821_nomaf) # He uses Nei's gene diversity

# expected heterozygosity
Hs(ne280_3821_nomaf) # pretty similar between 3 time periods, but slightly different than using poppr function above
Hs(ne280_3821_nomaf, ne280_3821_nomaf@strata$X1)

# Heterozygosity across all individuals
div <- summary(ne280_3821_nomaf)

# subset the genind object by population, so heterozygosity can be calculated by population
early.genind <- popsub(ne280_3821_nomaf, sublist = 'PADE_95011L2524')
mid.genind <- popsub(ne280_3821_nomaf, sublist = 'PADE_98027L2146')
late.genind <- popsub(ne280_3821_nomaf, sublist = 'PADE_09151L2330')

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

#### Pi calculations for 1296 loci ####
# How about trying a vcf that contains no missing data: 1296 loci across 285 fish
no.missing <- read.structure("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/no_missing_1296loci_285larvae.str",
                       n.ind = 285, n.loc = 1296, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# subset the genind object by population, so I can see if loci for which pi wasn't calculated are monomorphic
early.genind.no.missing <- popsub(no.missing, sublist = '1') # early
mid.genind.no.missing <- popsub(no.missing, sublist = '2') # mid
late.genind.no.missing <- popsub(no.missing, sublist = '3') # late

# Read in window pi calculations from vcftools on Amarel
early.win.pi.no.missing <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names94_95_no_missing.windowed.pi", header = TRUE) #2153 x 5
mid.win.pi.no.missing <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98_no_missing.windowed.pi", header = TRUE) #3475 x 5
late.win.pi.no.missing <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09_no_missing.windowed.pi", header = TRUE) #3478 x 5

# Read in CHOM column from no_missing.recode.vcf so I can figure out which ones are missing in the site pi files
contig_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/contig_names_1296.txt")
contig_names$index <- rownames(contig_names)

early.win.missing <- contig_names[!contig_names[,1] %in% early.win.pi.no.missing$CHROM,]
early.win.missing$index <- paste('SNP', early.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object
mid.win.missing <- contig_names[!contig_names[,1] %in% mid.win.pi.no.missing$CHROM,]
mid.win.missing$index <- paste('SNP', mid.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object
late.win.missing <- contig_names[!contig_names[,1] %in% late.win.pi.no.missing$CHROM,]
late.win.missing$index <- paste('SNP', late.win.missing$index, sep = '_') # add 'SNP_' to the index number for easier matching with genind object

# Early
early.names <- colnames(early.genind.no.missing@tab) # prep names
early.names.split <- data.frame(do.call('rbind', strsplit(as.character(early.names),'.',fixed=TRUE)))
colnames(early.genind.no.missing@tab) <- early.names.split$X1 # replace SNPX_X.YYY with only SNP_X

early.missing.genind <- early.genind.no.missing@tab[,colnames(early.genind.no.missing@tab) %in% early.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
early.missing.genind.unique <- early.missing.genind[, !duplicated(colnames(early.missing.genind))] # this is actually unnecessary because all sites are monomorphic with no missing data

early.pi <- early.genind.no.missing@tab[,!colnames(early.genind.no.missing@tab) %in% early.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
early.pi.genind.unique <- early.pi[, !duplicated(colnames(early.pi))]
early.source.table <- matrix(NA, nrow = 738, ncol = 4, byrow = TRUE, dimnames = list(colnames(early.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(early.pi.genind.unique)){
  early.source.table[i,] <- as.matrix(t(table(factor(early.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
early.source.table.pi <- cbind(early.source.table, early.win.pi.no.missing$PI) # cbind pi values with allele counts
early.source.table.pi.ordered <- early.source.table.pi[order(-early.source.table.pi[,'2']),]

early.lookup.table <- matrix(NA, nrow = 558, ncol = 4, byrow = TRUE, dimnames = list(colnames(early.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(early.missing.genind.unique)){
  early.lookup.table[i,] <- as.matrix(t(table(factor(early.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
early.lookup.table.ordered <- early.lookup.table[order(-early.lookup.table[,'2']),] #thse are indeed all monomorphic sites!

# Mid
mid.names <- colnames(mid.genind.no.missing@tab) # prep names
mid.names.split <- data.frame(do.call('rbind', strsplit(as.character(mid.names),'.',fixed=TRUE)))
colnames(mid.genind.no.missing@tab) <- mid.names.split$X1 # replace SNPX_X.YYY with only SNP_X

mid.missing.genind <- mid.genind.no.missing@tab[,colnames(mid.genind.no.missing@tab) %in% mid.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
mid.missing.genind.unique <- mid.missing.genind[, !duplicated(colnames(mid.missing.genind))]

mid.pi <- mid.genind.no.missing@tab[,!colnames(mid.genind.no.missing@tab) %in% mid.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
mid.pi.genind.unique <- mid.pi[, !duplicated(colnames(mid.pi))]
mid.source.table <- matrix(NA, nrow = 1180, ncol = 4, byrow = TRUE, dimnames = list(colnames(mid.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(mid.pi.genind.unique)){
  mid.source.table[i,] <- as.matrix(t(table(factor(mid.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
mid.source.table.pi <- cbind(mid.source.table, mid.win.pi.no.missing$PI) # cbind pi values with allele counts
mid.source.table.pi.ordered <- mid.source.table.pi[order(-mid.source.table.pi[,'2']),]

mid.lookup.table <- matrix(NA, nrow = 116, ncol = 4, byrow = TRUE, dimnames = list(colnames(mid.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(mid.missing.genind.unique)){
  mid.lookup.table[i,] <- as.matrix(t(table(factor(mid.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
mid.lookup.table.ordered <- mid.lookup.table[order(-mid.lookup.table[,'2']),] #thse are indeed all monomorphic sites!

# Late
late.names <- colnames(late.genind.no.missing@tab) # prep names
late.names.split <- data.frame(do.call('rbind', strsplit(as.character(late.names),'.',fixed=TRUE)))
colnames(late.genind.no.missing@tab) <- late.names.split$X1 # replace SNPX_X.YYY with only SNP_X

late.missing.genind <- late.genind.no.missing@tab[,colnames(late.genind.no.missing@tab) %in% late.win.missing$index] # these are all the alleles at loci for which pi wasn't calculated
late.missing.genind.unique <- late.missing.genind[, !duplicated(colnames(late.missing.genind))]

late.pi <- late.genind.no.missing@tab[,!colnames(late.genind.no.missing@tab) %in% late.win.missing$index] # these are all the alelles at loci for which pi WAS calculated
late.pi.genind.unique <- late.pi[, !duplicated(colnames(late.pi))]
late.source.table <- matrix(NA, nrow = 1262, ncol = 4, byrow = TRUE, dimnames = list(colnames(late.pi.genind.unique), c('0','1','2', 'NA')))
for (i in 1:ncol(late.pi.genind.unique)){
  late.source.table[i,] <- as.matrix(t(table(factor(late.pi.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,]
}
late.source.table.pi <- cbind(late.source.table, late.win.pi.no.missing$PI) # cbind pi values with allele counts
late.source.table.pi.ordered <- late.source.table.pi[order(-late.source.table.pi[,'2']),]

late.lookup.table <- matrix(NA, nrow = 34, ncol = 4, byrow = TRUE, dimnames = list(colnames(late.missing.genind.unique), c('0','1','2', 'NA'))) # these are the loci I have to assign a pi value to 
for (i in 1:ncol(late.missing.genind.unique)){
  late.lookup.table[i,] <- as.matrix(t(table(factor(late.missing.genind.unique[,i], levels = 0:2), useNA = 'always')))[1,] #0-2 = counts; NA = missing data
}
late.lookup.table.ordered <- late.lookup.table[order(-late.lookup.table[,'2']),] #thse are indeed all monomorphic sites!

# add zeros for all the loci that got dropped, then take mean
early.pi1296 <- c(early.win.pi.no.missing$PI, rep(0,558))
mid.pi1296 <- c(mid.win.pi.no.missing$PI, rep(0,116))
late.pi1296 <- c(late.win.pi.no.missing$PI, rep(0,34))

mean(early.pi1296)
mean(mid.pi1296)
mean(late.pi1296)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

early.boot <- boot(early.pi1296, samplemean, 1000)
plot(early.boot)
early.ci <- boot.ci(early.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
early.ci$normal # view 95% CIs

mid.boot <- boot(mid.pi1296, samplemean, 1000)
plot(mid.boot)
mid.ci <- boot.ci(mid.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.ci$normal # view 95% CIs

late.boot <- boot(late.pi1296, samplemean, 1000)
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
barplot(c(early.boot$t0, mid.boot$t0, late.boot$t0), ylim = c(0,0.0007), xlab = 'Larval cohort', ylab = 'Average nucleotide diversity (π)')
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


