#### This script starts with the first SNP on a contig with missing data allowed (3905 loci across 284 summer flounder) & removes loci not in HWE using an exact test. Then removes 9 fish that represent potential contamination/were previously ID as putative siblings by Colony ####

library(ade4)
library(adegenet)
library(wesanderson)
library(pegas)

nomac_missing <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.firstsnp.genepop.gen", ncode = 3L) #3905 loci x 284 fish

# ID loci that are not in HWE
# There are multiple ways to do this
# Using permuted p-values
hwe <- hw.test(nomac_missing,res="matrix")
pval <- hwe[hwe[,"Pr.exact"] < 0.001,]
dim(pval) #153 x 4

#### Now get genind object ready to remove SNPs not in HWE ####
cols <- colnames(nomac_missing@tab)
cols.split <- data.frame(do.call('rbind', strsplit(as.character(cols),'.',fixed=TRUE)))

# Make column names for loci in HWE only
cols.split.hwe <- cols.split[-which(cols.split$X1 %in% rownames(pval)),] 
cols.hwe.joined <- paste(cols.split.hwe$X1, cols.split.hwe$X2, sep = '.')

# Replace genind object column names with new names
colnames(nomac_missing@tab) <- cols.split$X1

# Remove loci not in HWE, add new column names & make into genind object
nomac_missing_hwe <- nomac_missing@tab[, -which(colnames(nomac_missing@tab) %in% rownames(pval))] 
dim(nomac_missing_hwe) # 284 x 7737
colnames(nomac_missing_hwe) <- cols.hwe.joined

nomac_missing_hwe_genind <- as.genind(nomac_missing_hwe) #284 x 3752 loci
pop(nomac_missing_hwe_genind) <- pop(nomac_missing)

#### Remove fish that were previously found to be highly heterozygous
sibs <- c('PADE_08115L2219','PADE_09137L2222','PADE_08127L2235','PADE_09079L2317','PADE_09113L2284','PADE_09143L2374','PADE_08008L2204','PADE_08133L2209','PADE_09138L2213') # these are all the putative sibs
sibs2 <- c('PADE_08115L2219','PADE_09137L2222','PADE_08127L2235','PADE_09143L2374','PADE_08008L2204','PADE_09138L2213') # these are the 6 that are especially heterozygous

nomac_missing_hwe_genind_remove9fish <- as.genind(nomac_missing_hwe_genind@tab[! rownames(nomac_missing_hwe_genind@tab) %in% sibs2,]) # 275 fish across 3752 loci or 278 fish across 3752 loci
pop(nomac_missing_hwe_genind_remove9fish) <- pop(nomac_missing_hwe_genind)[-match(sibs2, rownames(nomac_missing_hwe_genind@tab))] # removes the appropriate fish from the population section of genind object

#### Now let's look at a PCA ####
dim(nomac_missing_hwe_genind_remove9fish@tab)
sum(is.na(nomac_missing_hwe_genind_remove9fish@tab)) 
A <- scaleGen(nomac_missing_hwe_genind_remove9fish, NA.method = "mean")
pcaA <- dudi.pca(A,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pcaA$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

col <- wes_palette("Darjeeling1", 5, type = "discrete")
palette(col)

s.class(pcaA$li, nomac_missing_hwe_genind_remove9fish@pop, xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, clabel = 0, xlim = c(-10,0), ylim = c(-150,100))
axis(1, at=seq(-20,30, by=10), labels=seq(-20,30, by= 10), line = 1)
axis(2, at=seq(-100,80, by = 10), labels=seq(-100,80, by= 10), line = -7, las = 2)

eig_percent <- round((pcaA$eig/(sum(pcaA$eig)))*100,2)
eig_percent [1:3]
mtext("PC1 (0.72%)", side = 1, line = 3)
mtext("PC2 (0.70%)", side = 2, line = -4.5)
legend(15, -20,
       legend=levels(nomac_missing_hwe_genind_remove9fish@pop),
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

writeGenPop(nomac_missing_hwe_genind_remove9fish, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_278PADE_3752loci_missingallowed.gen", comment = '3752 loci with no missing data across 278 PADE, no MAF or MAC, loci in HWE only, putative siblings/contamination removed')

