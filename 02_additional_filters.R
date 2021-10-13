#### This script starts with the first SNP on a contig with missing data allowed (3905 loci across 284 summer flounder). Then removes several fish that represent potential contamination/were previously ID as putative siblings by Colony. Then removes loci not in HWP ####

library(ade4)
library(adegenet)
library(wesanderson)
library(pegas)

nomac_missing <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen", ncode = 3L) #3905 loci x 284 fish

#### Calculate heterozygosity within individuals ####
allele.matrix <- data.frame(nomac_missing@tab[,-4881]) #284 x 8058 # remove one locus that's homozygous across all individuals
num.homo.sites <- vector() # count of the number of homozygous loci in a given individual
num.hetero.sites <- vector() # counts the number of alleles with a count of 1 in a given individual (heterozygous). Divide by 2 to get number of heterozygous loci.

# Make sure number of 0's, 1's, 2's and NAs add up to the number of alleles (8058 alleles)
allele.counts.sums <- vector()
for (g in 1:nrow(allele.matrix)){
        allele.counts.sums[g] <- sum(table(as.numeric(allele.matrix[g,]), useNA = 'always'))
}

# First calculate the number of homozygous and heterozygous sites within each individual
for (i in 1:nrow(allele.matrix)) {
        num.homo.sites[i] <- length(which(allele.matrix[i,] == 2)) # best to use a count of 2 as homozygous (rather than 0) because some loci have more than 2 alleles. Because of this, the number of 2's in an individual will not equal the number of 0's, and number of 0's > number of 2's
        num.hetero.sites[i] <- length(which(allele.matrix[i,] == 1))/2
}

# Calculate the number of genotyped loci for each individual
num.geno.loci <- num.homo.sites + num.hetero.sites

# Now put these together: # of heterozygous sites/# of genotyped loci in an individual
prop.hetero.sites <- num.hetero.sites/num.geno.loci

hist(prop.hetero.sites, main = '', xlab = 'Within individual heterozygosity')
abline(v = mean(prop.hetero.sites), col = 'tomato')
mean(prop.hetero.sites) + 3*sd(prop.hetero.sites) # 0.09632892
filter <- allele.matrix[which(prop.hetero.sites >= 0.09632892),]
rownames(filter) # "PADE_08008L2204" "PADE_08115L2219" "PADE_08127L2235" "PADE_09137L2222" "PADE_09143L2374"

# Now remove fish that were found to be highly heterozygous
# sibs <- c('PADE_08115L2219','PADE_09137L2222','PADE_08127L2235','PADE_09079L2317','PADE_09113L2284','PADE_09143L2374','PADE_08008L2204','PADE_08133L2209','PADE_09138L2213') # these are all the putative sibs identified by Colony
sibs2 <- rownames(filter) # these are especially heterozygous based on the above analysis, plus they were ID'ed as putative sibs by Colony

nomac_missing_removehetfish <- as.genind(nomac_missing@tab[! rownames(nomac_missing@tab) %in% sibs2,]) # 279 fish across 3905 loci
pop(nomac_missing_removehetfish) <- pop(nomac_missing)[-match(sibs2, rownames(nomac_missing@tab))] # removes the appropriate fish from the population section of genind object

#### Next, ID loci that are not in HWE ####
# There are multiple ways to do this
# Using permuted p-values
hwe <- hw.test(nomac_missing_removehetfish,res="matrix")
pval <- hwe[hwe[,"Pr.exact"] < 0.001,] # This filter is used by the HapMap data set
dim(pval) #156 x 4; variable depending on the iteration
write.table(hwe, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_HWP_test.txt")

# Now get genind object ready to remove SNPs not in HWE
cols <- colnames(nomac_missing_removehetfish@tab)
cols.split <- data.frame(do.call('rbind', strsplit(as.character(cols),'.',fixed=TRUE)))

# Make column names for loci in HWE only
cols.split.hwe <- cols.split[-which(cols.split$X1 %in% rownames(pval)),] 
cols.hwe.joined <- paste(cols.split.hwe$X1, cols.split.hwe$X2, sep = '.')

# Replace genind object column names with new names
colnames(nomac_missing_removehetfish@tab) <- cols.split$X1

# Remove loci not in HWE, add new column names & make into genind object
nomac_missing_removehetfish_hwe <- nomac_missing_removehetfish@tab[, -which(colnames(nomac_missing_removehetfish@tab) %in% rownames(pval))] 
dim(nomac_missing_removehetfish_hwe) # 279 x 7725
colnames(nomac_missing_removehetfish_hwe) <- cols.hwe.joined

nomac_missing_removehetfish_hwe_genind <- as.genind(nomac_missing_removehetfish_hwe) #279 x 3749 loci
pop(nomac_missing_removehetfish_hwe_genind) <- pop(nomac_missing_removehetfish)

#### Now let's look at a PCA ####
dim(nomac_missing_removehetfish_hwe_genind@tab)
sum(is.na(nomac_missing_removehetfish_hwe_genind@tab)) 
A <- scaleGen(nomac_missing_removehetfish_hwe_genind, NA.method = "mean")
pcaA <- dudi.pca(A,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pcaA$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

col <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col)

s.class(pcaA$li, nomac_missing_removehetfish_hwe_genind@pop, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-5,0), ylim = c(-120,100))
axis(1, at=seq(-20,30, by=10), labels=seq(-20,30, by= 10), line = -1)
axis(2, at=seq(-80,80, by = 10), labels=seq(-80,80, by= 10), line = -9, las = 2)

eig_percent <- round((pcaA$eig/(sum(pcaA$eig)))*100,2)
eig_percent [1:3]
mtext("PC1 (0.71%)", side = 1, line = 1)
mtext("PC2 (0.70%)", side = 2, line = -6.5)
legend(15, -20,
       legend=levels(nomac_missing_removehetfish_hwe_genind@pop),
       pch=19,
       col = col,
       bty = "n",
       y.intersp = 0.8)

#### Write new genepop file: all loci in HWE, putative sibs/contamination removed ####
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

writeGenPop(nomac_missing_removehetfish_hwe_genind, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_279PADE_3749loci_missingallowed.gen", comment = '3749 loci with no missing data across 279 PADE, no MAF or MAC, putative siblings/contamination removed, loci in HWE only')

