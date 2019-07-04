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
library(XLConnect)
library(RColorBrewer)
library(plotly)

# Read in genepop file  with population identifiers for each of the three time periods
ne_data <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/SNP.DP3g95maf05.FIL.FIL.recode.firstsnp.genepop.gen", ncode = 3L) # wants .gen extension, file is the same as that with .txt extension; 3 digits used to denote nucleotide

# Need to modify names in allele count data (genind object) so that I can match them up with names in the larval database
freq_names <- as.vector(rownames(ne_data@tab))
freq_names_split <- do.call(rbind, strsplit(as.character(freq_names), 'L'))
freq_names_split2 <- separate(as.data.frame(freq_names_split), V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
freq.newnames <- paste(freq_names_split2$name1, freq_names_split2$name3, freq_names_split2$name2, freq_names_split2$name4, sep = '')
rownames(ne_data@tab) <- freq.newnames # replaces rownames in genind object with PADEXX_XXX formatting

# Read in and create complete meta data by adding Ne fish to connectivity metadata
meta <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt", header = TRUE) # read in connectivity metadata
meta2 <- readWorksheetFromFile("~/Documents/Graduate School/Rutgers/Summer Flounder/Larvae/Ne/Ne_larvs.xlsx", sheet = 1) # read in Ne metadata

ne_meta <- merge(meta, meta2, by = c('PinskyID', 'Sampling.Date', 'Month', 'Date', 'Year', 'Place'), all.x = TRUE, all.y = TRUE)

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
sum(is.na(ne_data$tab)) #10859 
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
s.class(pca1$li, pop(ne_data), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-32,45), ylim = c(-45,11))
axis(1, at=seq(-20,40, by=10), labels=seq(-20,40, by= 10), line = 0.5)
axis(2, at=seq(-50,10, by = 10), labels=seq(-50,10, by= 10), line = 0, las = 2)
mtext("PC1 (2.29%)", side = 1, line = 3.5)
mtext("PC2 (0.78%)", side = 2, line = 3)

legend(9, -25,
       legend=c("2008-2009 (n = 172)", "1994-1995 (n = 27)", "1997-1998 (n = 103)"),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Plot PCA by Year
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X1, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-20,45), ylim = c(-45,11))
axis(1, at=seq(-10,40, by=10), labels=seq(-10,40, by= 10), line = 0.5)
axis(2, at=seq(-50,10, by = 10), labels=seq(-50,10, by= 10), line = 0, las = 2)
mtext("PC1 (2.29%)", side = 1, line = 3.5)
mtext("PC2 (0.78%)", side = 2, line = 3)

legend(25, -25,
       legend=levels(ne_data@strata$X1),
       pch=c(19, 19, 19),
       col = col,
       bty = "n",
       y.intersp = 1)

# Plot PCA by capture Location
col <- brewer.pal(6, "Paired")
s.class(pca1$li, ne_data@strata$X2, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-20,45), ylim = c(-45,11))
axis(1, at=seq(-10,40, by=10), labels=seq(-10,40, by= 10), line = 0.5)
axis(2, at=seq(-50,10, by = 10), labels=seq(-50,10, by= 10), line = 0, las = 2)
mtext("PC1 (2.29%)", side = 1, line = 3.5)
mtext("PC2 (0.78%)", side = 2, line = 3)

legend(15, -25,
       legend=c('Little Egg Inlet, NJ', 'Beaufort, NC'),
       pch=19,
       col = col,
       bty = "n",
       y.intersp = 1)

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
# genotype accumulation curve
genotype_curve(ne_data, sample = 100, quiet = TRUE)

# diversity table using poppr
poppr(ne_data) # He uses Nei's gene diversity

# expected heterozygosity
Hs(ne_data) # pretty similar between 3 time periods, but slightly different than using poppr function above
Hs(ne_data, ne_data@strata$X1)

div <- summary(ne_data)

plot(div$Hobs, xlab = 'Loci number', ylab = 'Observed Heterozygosity', main = 'Observed heterozygosity per locus')
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")
bartlett.test(list(div$Hexp, div$Hobs))

# FST
ne_data
larvs <- read.table("structure_input_Jan9_2019_293_1904_9pops.txt", sep="\t", header = TRUE)

even_indexes<-seq(2,586,2)
odd_indexes<-seq(1,585,2)

odds <- data.frame(larvs[odd_indexes,]) # 293 x 1906
odds2 <- odds[,-c(1:2)] # 293 x 1904
evens <- data.frame(larvs[even_indexes,]) # 293 x 1906
evens2 <- evens[,-c(1:2)] # 293 x 1904

s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 293, ncol = 1904))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 293 x 1904

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop))

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 293 x 1905

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')

options("scipen"=100, "digits"=4) # forces output to not be in scientific notation

