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



#######################################################################
#### Diversity metrics ####
div <- summary(ne_data)

plot(div$Hobs, xlab = 'Loci number', ylab = 'Observed Heterozygosity', main = 'Observed heterozygosity per locus')
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")
bartlett.test(list(div$Hexp, div$Hobs))
