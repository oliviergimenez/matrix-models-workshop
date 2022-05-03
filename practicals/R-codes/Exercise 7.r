# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# Exercise 7: Management of overabundant cormorant

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Density-Independent matrix model for cormorants in Northern Europe           #
################################################################################
library(MASS) # an R package needed for some matrix functions
library(popbio)

# Define demographic parameters
cf <- 1         # paramter to change clutch size by a proportionate amount
cs <- 1         # paramter to change survival by a proportionate amount
a3 <- 0.1       # proportion females breeding at age 3
a4 <- 0.6       # proportion females breeding at age 4
c <- 3.5*cf     # clutch size times management-induced change
h <- 0.7        # nestling survival
s0 <- 0.64*cs   # first-year survival times management-induced change
s1 <- 0.86*cs   # surviavl from age 1 to 2
sa <- 0.89*cs   # survival of older birds

# Next, create the 5x5 density-independent matrix model for cormorants
A <- matrix(c(
  0, 0, s0*a3*h*c/2, s0*a4*h*c/2, s0*h*c/2,
  s1, 0, 0, 0, 0,
  0, sa, 0, 0, 0,
  0, 0, sa, 0, 0,
  0, 0, 0, sa, sa),nrow=5, byrow=T)

# call the 'lambda' function from the 'popbio' package to find the long-term
# population growth rate
lambda(A)

################################################################################
# Density-Dependent matrix model for cormorants in Northern Europe             #
################################################################################
### For question 3
Bseq <- seq(12800,13200,by=1)  # sequence of DD coefficients B to loop over
eqnbc <- matrix(0,length(Bseq),1) # storage for equilibrium number of breeding females

for (i in 1:length(Bseq)){
  n <- rep(1,5)                 # initial vector of abundance
  nbc <- (a3*n[3]+a4*n[4]+n[5]) # number of breeding females
  B <- Bseq[i]                  # coefficient of density dependence
  for (t in 1:100){
    h <- 0.7*exp(-nbc/B)
    A <- matrix(c(
      0, 0, s0*a3*h*c/2, s0*a4*h*c/2, s0*h*c/2,
      s1, 0, 0, 0, 0,
      0, sa, 0, 0, 0,
      0, 0, sa, 0, 0,
      0, 0, 0, sa, sa),nrow=5, byrow=T)
    n <- A%*%n
    nbc <- (a3*n[3]+a4*n[4]+n[5])
  }
  eqnbc[i] <- nbc
}

par(mar = c(5, 6, 4, 2),mgp = c(3.5, 1, 0))
plot(Bseq, eqnbc, xlab = list("DD Coefficient B", cex = 1.5),
  ylab = list("Equilibrium Abundanc of Breeders", cex = 1.5), col = "black",
  type = "l", lwd = 2, cex.axis = 1.1)
abline(h = 20000, col = "darkgray", lty = 3, lwd = 2)
      

### For question 4
rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful
# Define demographic parameters
cf <- 1         # paramter to change clutch size by a proportionate amount
cs <- 1         # paramter to change survival by a proportionate amount
a3 <- 0.1       # proportion females breeding at age 3
a4 <- 0.6       # proportion females breeding at age 4
c <- 3.5*cf     # clutch size times management-induced change
s0 <- 0.64*cs   # first-year survival times management-induced change
s1 <- 0.86*cs   # surviavl from age 1 to 2
sa <- 0.89*cs   # survival of older birds

B <- 12921
n <- rep(1,5)                 # initial vector of abundance
nbc <- (a3*n[3]+a4*n[4]+n[5]) # number of breeding females
for (t in 1:100){
  h <- 0.7*exp(-nbc/B)
  A <- matrix(c(
    0, 0, s0*a3*h*c/2, s0*a4*h*c/2, s0*h*c/2,
    s1, 0, 0, 0, 0,
    0, sa, 0, 0, 0,
    0, 0, sa, 0, 0,
    0, 0, 0, sa, sa),nrow=5, byrow=T)
  n <- A%*%n
  nbc <- (a3*n[3]+a4*n[4]+n[5])
}

nbc


