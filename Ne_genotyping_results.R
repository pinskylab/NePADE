setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/")

counts.ne <- read.table("countsperindiv.txt")

sum(counts.ne) #197931570

mean(counts.ne[,1]) #585596.4
colMeans(counts.ne)

sd(counts.ne[,1]) #628558.2
summary(counts.ne)

######################################################
# Coverage across all individuals and reference sites
cov <- read.table("cov.stats.txt") # this is cov.stats that ddocent outputs

(sum(cov$V4)/length(cov$V4))/338 # 25.2X; 338 larval individuals
mean(cov$V4/338)

# Example to calculate coverage across SNPs
# Coverage across 285 individuals and 3979 SNPs
# read in list of contig names containing 3979 SNPs
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/")
all <- read.table("all_contigs3979.txt")

all.cov <- cov[cov$V1%in%all$V1,] # keep 3979 contigs

(sum(all.cov$V4)/length(all.cov$V4))/285 # 62X
mean(all.cov$V4/285) # same as above

