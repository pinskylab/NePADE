setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/")

counts.ne <- read.table("countsperindiv2.txt")

sum(counts.ne) # 190802200

mean(counts.ne[,1]) # 576441.7
colMeans(counts.ne)

sd(counts.ne[,1]) # 626768.1
summary(counts.ne)

######################################################
# Coverage across all individuals and reference sites
cov <- read.table("cov.stats2.txt") # this is cov.stats that ddocent outputs. File is in .bed format: 1) contig name, 2) contig start position, 3) contig end position, 4) coverage

(sum(cov$V4)/length(cov$V4))/331 # 24.8 when 331 larvae; avg # reads/contig/individual
mean(cov$V4/331)

# Example to calculate coverage across SNPs
# Coverage across 284 individuals and 3905 SNPs
# read in list of contig names containing 3905 SNPs
all <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/all_contigs3905.txt")

all.cov <- cov[cov$V1%in%all$V1,] # keep 3905 contigs

(sum(all.cov$V4)/length(all.cov$V4))/284 # 61.1 when 284 larvae
mean(all.cov$V4/284) # same as above

