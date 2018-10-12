#### This script uses information from the summer flounder stock assessment to calculate the number of potenital breeders in each year class ####
library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# Read in PADE stock assessment data from Mark Terceiro
data <- read_excel_allsheets("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/PADE_stock_assessment16.xlsx")

N <- as.data.frame(lapply(data$N_at_Age[-1,], as.integer)) # Number of fish at age
colnames(N) <- data$N_at_Age[1,]
rownames(N) <- data$N_at_Age[-1,1]
N <- N[,-1]*1000 # N is in 000's, so multiplying by 1000, to add 000 to each value

F <- data$F_at_Age[-1,] # estimated fishing mortality at age
colnames(F) <- data$F_at_Age[1,]
rownames(F) <- data$F_at_Age[-1,1]
F[,-1] <- round(F[,-1],3)

M <- data$M_at_Age[-1,] # input fishing mortality at age
colnames(M) <- data$M_at_Age[1,]
rownames(M) <- data$M_at_Age[-1,1]
M[,-1] <- round(M[,-1],3)

prop <- data$Prop_Mat_at_Age[-1,] # Proportion mature at age in November (peak spawning) of a given year
colnames(prop) <- data$Prop_Mat_at_Age[1,]
rownames(prop) <- data$Prop_Mat_at_Age[-1,1]
prop <- prop[,-1]

# Calculate N for a given age/month
# Nt+1 = Nt*e^-pZ where p = month and Z = M*F
# p = month/12; Z = M*F
Z <- F[,-1]*M[,-1]
p <- 10 # for example, October

N_atmonth <- vector()
for (i in 1:length(unlist(N))){
  N_atmonth[i] <- unlist(N)[i]*(exp(1))^(-(p/12)*unlist(Z)[i])
}

N_atmonth_matrix <- matrix(N_atmonth, nrow = 34, ncol = 8)
N_atmonth_matrix_year <- rowSums(N_atmonth_matrix)

# Now calculate proportion mature for a given age in November based on N
prop_Nov <- as.vector(N_atmonth * unlist(prop))

# Convert this to matrix of appropriate dimensions
N_spawners <- matrix(prop_Nov, nrow = 34, ncol = 8)
colnames(N_spawners) <- colnames(N)
rownames(N_spawners) <- rownames(N)

# Check to make sure calculations are correct
check <- vector()
for (c in 1:length(unlist(N_spawners))) {
  check[c] <- N_atmonth[c] * unlist(prop)[c] == unlist(N_spawners)[c]
  }

check.matrix <- matrix(check, nrow = 34, ncol = 8) # should be all TRUE

# Sum number of spawners across all age classes
N_spawners_allages <- as.vector(rowSums(N_spawners))
year <- as.vector(rownames(N_spawners))
N_spawners_allages <- as.data.frame(cbind(as.integer(year), as.integer(N_spawners_allages)))

plot(year, N_atmonth_matrix_year, xlab = 'Year', ylab = 'Abundance at peak spawning (Nov)', pch = 19)
lw1 <- loess(N_atmonth_matrix_year ~ year)
lines(predict(lw1), x = N_spawners_allages$V1, col = 'black') # plots loess line for number of fish at peak spawning
points(N_spawners_allages$V1, N_spawners_allages$V2, col = "blue", pch = 19)
lw2 <- loess(N_spawners_allages$V2 ~ N_spawners_allages$V1)
lines(predict(lw2), x = N_spawners_allages$V1, col = 'blue') # plots loess line for mature fish

legend("bottomright",
       legend = c('Estimated abundance', 'Estimated mature'),
       pch=19,
       col = c('black', 'blue'))
