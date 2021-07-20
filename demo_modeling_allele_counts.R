library(ade4)
library(adegenet)

# Minor allele count of 3 that resulted in 1196 loci across 280 fish
mac3 <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1196loci_complete.gen", ncode = 3L)

# Allele counts?
dim(mac3@tab) #280 x 2392
mac3.allele.counts <- colSums(mac3@tab)
summary(mac3.allele.counts)
mac3.allele.counts[which(mac3.allele.counts < 3)] #31 alleles with counts of 1 or 2
table(mac3.allele.counts[which(mac3.allele.counts < 3)])
hist(mac3.allele.counts, breaks = 559)

alleles.mac3 <- which(mac3.allele.counts < 3)
for (i in alleles.mac3) {
  print(names(which(mac3@tab[,print(i)] != 0)))
  readline(prompt="Press [enter] to continue")
}

# No minor allele count that resulted in 1084 loci across 284 fish
nomac <- read.genepop("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/Ne_PADE_1084loci_complete.gen", ncode = 3L)

# Allele counts?
dim(nomac@tab) #284 x 2168
nomac.allele.counts <- colSums(nomac@tab)
summary(nomac.allele.counts)
nomac.allele.counts[which(nomac.allele.counts < 3)] #184 alleles with counts of 1 or 2
table(nomac.allele.counts[which(nomac.allele.counts < 3)])
hist(nomac.allele.counts, breaks = 567)

alleles.nomac <- which(nomac.allele.counts < 3)
for (i in alleles.nomac) {
  print(names(which(nomac@tab[,print(i)] != 0)))
  readline(prompt="Press [enter] to continue")
}
