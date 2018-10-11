data <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files/masterPADElarvae.txt", header = TRUE)

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)
library(diveRsity)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

#### Dealing with loci with multiple alleles ####
multiallelic <- names(which(full@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")

multiallelic_data <- full@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(full@tab))] # mostly works, but still includes SNP_145*
multiallelic_data <- multiallelic_data[,-c(59:78)] # gets rid of additional SNP_145* that are actually biallelic but got picked up in the grep above
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #67
biallelic_only <- full.dataframe[,!colnames(full.dataframe) %in% multiallelic_names] # 528x3764 (before: 3831-67 = 3764)

# Separate genind for adults and larvae
adults <- rownames(biallelic_only)[205:439] #528 fish total
larvs <- biallelic_only[!(rownames(biallelic_only) %in% adults),]
# write.table(larvs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/larvs_allelecounts.txt", sep="\t")
# adult_data <- biallelic_only[(rownames(biallelic_only) %in% adults),]
# dim(larvs) #293 x 3764
# larvs <- as.genind(larvs)
# larvs@pop <- full@pop[which(full@pop != 4)]

# Modify larvae names so that I can subset them later
library(tidyr)

gen.fishnames.split <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(larvs)), 'L', fixed = TRUE)))
gen.fishnames.split2 <- separate(gen.fishnames.split, V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
gen.newnames <- paste(gen.fishnames.split2$name1, gen.fishnames.split2$name3, gen.fishnames.split2$name2, gen.fishnames.split2$name4, sep = '')
rownames(larvs) <- gen.newnames

#### Subset master data to those fish with genetic data ####
data2 <- data[as.character(data$PinskyID) %in% rownames(larvs),]
dim(data2)

table(data2$Year)
table(data2$Sampling.Date)

# And now to just NJ
dataNJ <- data2[which(data2$Place == 'Little Egg Inlet, NJ'),]
dim(dataNJ)
table(dataNJ$Year)
table(dataNJ$Sampling.Date)

# Just fish in NC
dataNC <- data2[which(data2$Place == 'Beaufort, NC'),]
dim(dataNC)
table(dataNC$Year)
table(dataNC$Sampling.Date)

#### Read in names of fish that I tried sequencing ####
seq_names <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Effective population size/namelist.txt")

rownames(full@tab)

dropouts <- seq_names[!as.character(seq_names[,1]) %in% rownames(full@tab),] # 151

# Change naming format to match date with ID
dropouts.split <- separate(as.data.frame(dropouts, nrow = 151, ncol = 1), dropouts, c("name1", "name2", "name3", "name4", "name5"),sep = c(4,5,7,10))
dropouts.newnames <- paste(dropouts.split$name1, dropouts.split$name3, dropouts.split$name2, dropouts.split$name4, sep = '')
table(dropouts.split$name3) # 45 adults

# When were the dropouts sampled?
dropouts.data <- data[data$PinskyID %in% dropouts.newnames,] # 106 x 39; makes sense because no adults in here
dropouts.data.early <- dropouts.data[dropouts.data$Year %in% c('1991', '1992'),] # subsample to years: 9 x 39
dropouts.data.2008 <- dropouts.data[dropouts.data$Year %in% c('2008'),] # subsample to years: 10 x 39
dropouts.data.2009 <- dropouts.data[dropouts.data$Year %in% c('2009'),] # 17 x 39

dropouts.data.2008.fall <- dropouts.data.2008[dropouts.data.2008$Month %in% c(10,11,12),] #4 fall 2008 fish
dropouts.data.2009.winter <- dropouts.data.2009[dropouts.data.2009$Month %in% c(1,2,3,4),] #11 winter 2009 fish

#### How many fish from each year did I try to sequence? ####
extractions <- read.csv("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Effective population size/extraction_list.csv", header = TRUE)

# Only keep fish that I tried extracting
extractedyes <- extractions[which(extractions$Extracted. == 'TRUE'),] #428
table(extractedyes$Year)

# Which fish didn't get extracted?
extractedno <- extractions[which(extractions$Extracted. == 'FALSE'),] #45
table(extractedno$Year)

extractedno.2008 <- extractedno[which(extractedno$Year == '2008'),]
extractedno.2009 <- extractedno[which(extractedno$Year == '2009'),]
