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
library(boot)

# Read in genepop file  with population identifiers for each of the three time periods
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # wants .gen extension, file is the same as that with .txt extension; 3 digits used to denote nucleotide
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/SNP.DP3g95maf05.FIL.FIL.recode.trimmed144.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset because the read lengths of the different sequencing runs are different
# ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/troubleshooting/newref/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference
ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140

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
pop_strata <- data.frame(cbind(ordered_meta_sub2$Year, ordered_meta_sub2$Place))
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
s.class(pca1$li, pop(ne_data), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-22,50), ylim = c(-50,60))
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
ne_data_nomaf <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # troubleshooting the Ne dataset, these are SNPs called from the new larval reference and all mapped reads have been trimmed to 140; no MAF filter applied and some missing data allowed: 3979 loci

# genotype accumulation curve
# genotype_curve(ne_data, sample = 100, quiet = TRUE)

# diversity table using poppr
poppr(ne_data_nomaf) # He uses Nei's gene diversity

# expected heterozygosity
Hs(ne_data_nomaf) # pretty similar between 3 time periods, but slightly different than using poppr function above
Hs(ne_data_nomaf, ne_data_nomaf@strata$X1)

# Heterozygosity across all individuals
div <- summary(ne_data_nomaf)

# subset the genind object by population, so heterozygosity can be calculated by population
early.genind <- popsub(ne_data_nomaf, sublist = 'PADE_95011L2524')
mid.genind <- popsub(ne_data_nomaf, sublist = 'PADE_98027L2146')
late.genind <- popsub(ne_data_nomaf, sublist = 'PADE_09151L2330')

early.he <- summary(early.genind)
mid.he <- summary(mid.genind)
late.he <- summary(late.genind)

mean(early.he$Hobs)
mean(mid.he$Hobs)
mean(late.he$Hobs)

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

# Site pi
# Read in pi calculations from vcftools on Amarel
early.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names94_95.sites.pi", header = TRUE)
mid.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98.sites.pi", header = TRUE)
late.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09.sites.pi", header = TRUE)

# Read in CHOM column from SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.vcf so I can figure out which ones are missing in the site pi files
contig_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/contig_names_3979.txt")

early.missing <- contig_names[!contig_names[,1] %in% early.pi$CHROM,]
mid.missing <- contig_names[!contig_names[,1] %in% mid.pi$CHROM,]
late.missing <- contig_names[!contig_names[,1] %in% late.pi$CHROM,]

# Read in window pi calculations from vcftools on Amarel
early.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names94_95.windowed.pi", header = TRUE) #2153 x 5
mid.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names97_98.windowed.pi", header = TRUE) #3475 x 5
late.win.pi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/names08_09.windowed.pi", header = TRUE) #3478 x 5

# Read in CHOM column from SNP.DP3g95nomaf.FIL.FIL.recode.firstsnp.vcf so I can figure out which ones are missing in the site pi files
contig_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/contig_names_3979.txt")
contig_names$index <- rownames(contig_names)

early.win.missing <- contig_names[!contig_names[,1] %in% early.win.pi$CHROM,] #names and indices of snps missing for window pi calculations
mid.win.missing <- contig_names[!contig_names[,1] %in% mid.win.pi$CHROM,]
late.win.missing <- contig_names[!contig_names[,1] %in% late.win.pi$CHROM,]

which(apply(early.genind@tab, 2, function(x) length(unique(x))) == 1) # these are all loci with one unique value
which(apply(mid.genind@tab, 2, function(x) length(unique(x))) == 1) 
which(apply(late.genind@tab, 2, function(x) length(unique(x))) == 1)

# Mean pi across 140bp windows
mean(early.win.pi$PI)
mean(mid.win.pi$PI)
mean(late.win.pi$PI)

# add zeros for all the loci that got dropped, then take mean
early.pi3979 <- c(early.win.pi$PI, rep(0,1826))
mid.pi3979 <- c(mid.win.pi$PI, rep(0,504))
late.pi3979 <- c(late.win.pi$PI, rep(0,501))

mean(early.pi3979)
mean(mid.pi3979)
mean(late.pi3979)

# bootstrap pi to get CI
samplemean <- function(x, d) {
  return(mean(x[d]))
}

early.boot <- boot(early.pi3979, samplemean, 1000)
plot(early.boot)
early.ci <- boot.ci(early.boot, conf = 0.95, type = 'all')
early.ci$normal # view 95% CIs

mid.boot <- boot(mid.pi3979, samplemean, 1000)
plot(mid.boot)
mid.ci <- boot.ci(mid.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
mid.ci$normal # view 95% CIs

late.boot <- boot(late.pi3979, samplemean, 1000)
plot(late.boot)
late.ci <- boot.ci(late.boot, conf = 0.95, type = c('norm', 'basic', 'perc'))
late.ci$normal # view 95% CIs

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



