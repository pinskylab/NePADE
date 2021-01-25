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

M <- data$M_at_Age[-1,] # input natural mortality at age
colnames(M) <- data$M_at_Age[1,]
rownames(M) <- data$M_at_Age[-1,1]
M[,-1] <- round(M[,-1],3)

prop <- data$Prop_Mat_at_Age[-1,] # Proportion mature at age in November (peak spawning) of a given year
colnames(prop) <- data$Prop_Mat_at_Age[1,]
rownames(prop) <- data$Prop_Mat_at_Age[-1,1]
prop <- prop[,-1]

# Calculate N for a given age/month
# Nt+1 = Nt*e^-pZ where p = month and Z = M+F (natural mortality rate + fishing mortality rate)
# p = month/12; Z = M+F
Z <- F[,-1]+M[,-1]
p <- 1 # for example, January

N_atmonth <- vector()
for (i in 1:length(unlist(N))){
  N_atmonth[i] <- unlist(N)[i]*(exp(1))^(-(p/12)*unlist(Z)[i]) #e=2.718282
}

N_atmonth_matrix <- matrix(N_atmonth, nrow = 34, ncol = 8) # this is Nc
N_atmonth_matrix_year <- rowSums(N_atmonth_matrix)

# Now calculate proportion mature for a given age in November based on N
prop_Jan <- as.vector(N_atmonth * unlist(prop))

# Convert this to matrix of appropriate dimensions
N_spawners <- matrix(prop_Jan, nrow = 34, ncol = 8) # this is Nb
colnames(N_spawners) <- colnames(N)
rownames(N_spawners) <- rownames(N)
N_spawners_year <- rowSums(N_spawners)

# Check to make sure calculations are correct
check <- vector()
for (c in 1:length(unlist(N_spawners))) {
  check[c] <- N_atmonth[c] * unlist(prop)[c] == unlist(N_spawners)[c]
}

check.matrix <- matrix(check, nrow = 34, ncol = 8) # should be all TRUE
table(check.matrix)

# Mean length at age in Jan based on 'Guideliness for Estimating Lengths at Age for 18 Northwest Atlantic Finfish and Shellfish Species'
lengths <- c(12.8,	24.3,	34.1,	42.4,	49.6,	55.6,	60.8, 65.2) # ages 1-8 for females 1975-1988, with age defined by Jan 1, so age 1 here is actually YOY or age 0
lengths_mm <- lengths *10

lengths_mm_matrix <- matrix(rep(lengths_mm, 34),
                            ncol = length(lengths_mm),
                            byrow = T)

# Fecundity per individual at each age based on length at age from Morse 1981
# F = 0.0007975*length^3.402
indiv_fec_by_age <- 0.0007975*lengths_mm_matrix^3.402

# Using individual fecundity at age, calculate fecundity for each age class given number of summer flounder in each age class
fec_by_age <- indiv_fec_by_age * N_spawners
fec_by_age_year <- rowSums(fec_by_age) # Total fecundity for each year

# What proportion of each age class is contributing to fecundity?
prop_fec_by_age <- fec_by_age/fec_by_age_year
rowSums(prop_fec_by_age) # Check that rows sum to 1; yes

# Now translate this back into number of individuals
N_contributing_eggs <- prop_fec_by_age * N_spawners_year
N_contributing_eggs_year <- rowSums(N_contributing_eggs)

# Calculate G for females for each year
ages <- c(0:7)
g_female_year <- rowSums(sweep(N_contributing_eggs, MARGIN = 2, ages, '*'))/N_contributing_eggs_year # average age of female spawners for each year
mean(g_female_year[1:7]) #based on 1982-1988

g_male_year <- rowSums(sweep(N_spawners, MARGIN = 2, ages, '*'))/N_spawners_year # assuming equal sex ratio through time and that all mature males have equal chance of reproducing; in reality, more males are shorter length and more females are longer
mean(g_male_year[1:7]) #based on 1982-1988

mean(c(g_female_year[1:27],g_male_year[1:27])) # average G of males and females from 1982-2008

# Plot g over time
plot(c(1982:2015), g_female_year, xlab = 'Year', ylab = 'Generation length', col = 'tomato', pch = 19, ylim = c(0,5))
points(c(1982:2015), g_male_year, col = 'blue', pch = 19)
lines(c(1982:2015), g_female_year, col = 'tomato')
lines(c(1982:2015), g_male_year, col = 'blue')
legend('bottomright',
       legend = c('Females', 'Males'),
       col = c('tomato', 'blue'),
       pch = 19)

# Plot N contributing eggs in each age class for each year
for (i in 1:nrow(N_contributing_eggs)) {
plot(ages, N_contributing_eggs[i,], xlab = 'Age', ylab = 'Number of females contributing eggs', main = paste0(row.names(N_contributing_eggs)[i]))
  readline(prompt="Press [enter] to continue")
}

# Plot N spawners in each age class for each year
for (i in 1:nrow(N_spawners)) {
  plot(ages, N_spawners[i,], xlab = 'Age', ylab = 'Number of mature fish', main = paste0(row.names(N_spawners)[i]))
  readline(prompt="Press [enter] to continue")
}


row.names(N_contributing_eggs)[i]


