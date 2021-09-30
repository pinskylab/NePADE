#######################################################
#### Keeping first SNP on a contig for Ne analysis ####
#######################################################

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/")

data <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.vcf", skip = 66, sep = "\t", header = FALSE) #based on new reference using larval fish and mapping of all reads trimmed to 140, no maf, 7 fish from Fall 2009 removed, no mac

# Keeping only the first SNP at each contig
vcf_firstsnps <- data[!duplicated(data[,1]),] # these are the first snps at each contig
dim(vcf_firstsnps) # 3905 x 293 SNP.DP3g95nomafnomac.FIL.FIL.recode.140trimmed.284fish.vcf

write.table(vcf_firstsnps, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/newref_alltrimmed140/SNP.DP3g95nomafmac.FIL.FIL.recode.140trimmed.284fish.firstsnp.vcf", sep="\t", col.names = FALSE, row.names = FALSE)

# Manually removed all the " and replaced with nothing. Add header columns back in too prior to converting to genepop format.