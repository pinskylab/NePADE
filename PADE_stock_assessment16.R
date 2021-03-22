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
p <- 10 # for example, October

N_atmonth <- vector()
for (i in 1:length(unlist(N))){
  N_atmonth[i] <- unlist(N)[i]*(exp(1))^(-(p/12)*unlist(Z)[i]) #e=2.718282
}

N_atmonth_matrix <- matrix(N_atmonth, nrow = 34, ncol = 8) # this is Nc
N_atmonth_matrix_year <- rowSums(N_atmonth_matrix)

# Now calculate proportion mature for a given age in November based on N
prop_Nov <- as.vector(N_atmonth * unlist(prop))

# Convert this to matrix of appropriate dimensions
N_spawners <- matrix(prop_Nov, nrow = 34, ncol = 8) # this is Nb
colnames(N_spawners) <- colnames(N)
rownames(N_spawners) <- rownames(N)

# Check to make sure calculations are correct
check <- vector()
for (c in 1:length(unlist(N_spawners))) {
  check[c] <- N_atmonth[c] * unlist(prop)[c] == unlist(N_spawners)[c]
  }

check.matrix <- matrix(check, nrow = 34, ncol = 8) # should be all TRUE
table(check.matrix)

# Sum number of spawners across all age classes
N_spawners_allages <- as.vector(rowSums(N_spawners))
year <- as.vector(rownames(N_spawners))
N_spawners_allages <- as.data.frame(cbind(as.integer(year), as.integer(N_spawners_allages)))

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/Nc_and_Nb_overtime.png",width=7, height=5, res=300, units="in")
par(mar=c(4.5, 5, 1.5, 1), # panel magin size in "line number" units
    mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
    tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
    cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
    ps=12
)

# options(scipen = 5)

plot(year, N_atmonth_matrix_year/1000000, xlab = 'Year', ylab = 'Abundance (in millions) at peak spawning (Nov)', pch = 19, ylim = c(10,110))
lines(year, N_atmonth_matrix_year/1000000, col = 'black', pch = 19)
# lw1 <- loess(N_atmonth_matrix_year ~ year)
# lines(predict(lw1), x = N_spawners_allages$V1, col = 'black') # plots loess line for number of fish at peak spawning
points(N_spawners_allages$V1, N_spawners_allages$V2/1000000, col = "gray70", pch = 19)
lines(N_spawners_allages$V1, N_spawners_allages$V2/1000000, col = "gray70", pch = 19)
# lw2 <- loess(N_spawners_allages$V2 ~ N_spawners_allages$V1)
# lines(predict(lw2), x = N_spawners_allages$V1, col = 'gray70') # plots loess line for mature fish

legend("bottomright",
       legend = c('Estimated total abundance', 'Estimated mature spawners'),
       pch=19,
       col = c('black', 'gray70'))

# legend("bottomright",
#        legend = c(as.expression(bquote('Estimated total abundance (N'['C']*')')), as.expression(bquote('Estimated mature (N'['B']*')'))),
#        pch=19,
#        col = c('black', 'gray70'))

dev.off()

#### Calculating Ne/Nc
N_spawners_allages # This is the number of breeders of all age classes in a year

# Read in Ne data
mod6 <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/NePADE/demo_modeling/model_results/model6.bestlhoods.summary.txt", header = TRUE) # Read in ML parameters from Model 6

# Find maximum from best model
# Exponential growth of NANC, then bottleneck & instantaneous recovery (Model 4)
max(mod6$MaxEstLhood)
mod6_ml_parameters <- mod6[which(mod6$MaxEstLhood == max(mod6$MaxEstLhood)),]

# 1988
(mod6_ml_parameters$NANC/2)/N_spawners_allages$V2[7] #calculate year by 2008 - (2*TBOT) - (2*TLEN)

# 1990
(mod6_ml_parameters$NBOT/2)/N_spawners_allages$V2[9] #calculate year by 2008 - (2*TBOT)

# 2008
(mod6_ml_parameters$NPOP08/2)/N_spawners_allages$V2[27]

#### Calculating r for summer flounder using estimates of N ####
N
count <- rowSums(N)
year <- as.numeric(rownames(N))

plot(year, count)

# Calculating r going back in time, so negative r mean population expansion
# For a single year
r <- vector()
for (i in 1:length(count)) {
  r[i] <- (log(count[i]/count[i+1]))/abs((year[i]-year[i+1]))
}

# Across 2 years
r <- vector()
for (i in 1:length(count)) {
  r[i] <- (log(count[i]/count[i+2]))/abs((year[i]-year[i+2]))
}

# Across 3 years
r <- vector()
for (i in 1:length(count)) {
  r[i] <- (log(count[i]/count[i+3]))/abs((year[i]-year[i+3]))
}

# Across 4 years
r <- vector()
for (i in 1:length(count)) {
  r[i] <- (log(count[i]/count[i+4]))/abs((year[i]-year[i+4]))
}
