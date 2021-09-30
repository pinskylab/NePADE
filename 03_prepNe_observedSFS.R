#### This script reads in the genepop file output by 02_additional_filters.R, keeps only loci with two alleles, and removes loci with missing data across 275 larvae ####

library(ade4)
library(adegenet)
library(hierfstat)
library(pegas)
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
ne_278fish_3752loci <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/Ne_278PADE_3752loci_missingallowed.gen", ncode = 3L)

# Remove loci with 3 or 4 alleles because I think it's causing a problem
ToKeep <- which(ne_278fish_3752loci@loc.n.all == 2) #3520 loci
ne_278fish_3520_biallelic <- ne_278fish_3752loci[loc = ToKeep]

pops <- as.data.frame(ne_278fish_3520_biallelic@pop)
data <- as.data.frame(ne_278fish_3520_biallelic@tab)

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
dim(data_sub) # 278 x 2140

# Check to see if all SNPs having no missing data & write str file 
test <- sapply(data_sub, function(y) sum(length(which(is.na(y)))))
test <- data.frame(test) # all zeros
test <- data.matrix(test)
summary(test)

data_sub <- as.genind(data_sub)
data_sub@pop <- as.factor(pops[,1])

writeGenPop(data_sub, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_278PADE_1070loci_complete.gen", comment = '1070 biallelic loci with no missing data across 278 PADE, no MAF or MAC filters, all loci in HWE, 6 highly heterozygous PADE removed')

