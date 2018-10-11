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
data <- read_excel_allsheets("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Effective population size/PADE_stock_assessment16.xlsx")

N <- as.data.frame(lapply(data$N_at_Age[-1,], as.integer)) # Number of fish at age
colnames(N) <- data$N_at_Age[1,]

F <- data$F_at_Age[-1,]
colnames(F) <- data$F_at_Age[1,]
F[,-1] <- round(F[,-1],3)

M <- data$M_at_Age[-1,]
colnames(M) <- data$M_at_Age[1,]
M[,-1] <- round(M[,-1],3)

# p = month/12; Z = M*F
Z <- F[,-1]*M[,-1]

p <- 10
Nt+1 = Nt*e^-pZ
