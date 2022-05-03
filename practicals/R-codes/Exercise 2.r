# Matrix Models for Population Management & Conservation
# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# EXERCISE 2: Perturbation analysis of white stork model

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Piece 1 of code                                                              #
################################################################################
library(MASS) # an R package needed for some matrix functions
library(popbio)
# It turns out that the 'popbio' package contains a built-in function that can 
# be called with just a single line of R code to calculate the long-term 
# geometric rate of population growth for a matrix population model. 
# Lets gain some practice using this package for analysis of matrix models.

# First, define the demographic parameters
u <- 0.45
r <- 0.818
b <- 2.9
s0 <- 0.482
s1 <- 0.75

# Next, create the pre birth-pulse white stork matrix population model
A <- matrix(c(
  0, 0, s0*u*r*b/2, s0*r*b/2,
  s1, 0, 0, 0,
  0, s1, 0, 0,
  0, 0, s1, s1), nrow = 4, byrow = TRUE)

# Then use the following popbio function; quite easy!
lambda(A)

# Perhaps too easy. The following code demonstrates the matrix algebra that 
# is used 'behind the scences' of the lambda function in the popbio package.
rows <- dim(A)[1]
cols <- dim(A)[2]
eig <- eigen(A)           # eigenvalues of A
EigVecs <- eig$vectors    # eigenvectors of A
Lambdas <- Re(eig$values) # real number components of eigenvalues
Lambda <- max(Lambdas)    # long-term geometric rate of population growth
Lambda

################################################################################
# Piece 2 of code                                                              #
################################################################################
pos <- which.max(Lambdas)  # finding the position of the dominant eigenvalue
w <- Re(eig$vectors[1:rows,pos]) # its associated right eigenvector
sad <- w/(sum(w))
sad <- round(sad,3) # scaled dominant right eigenvector: Stable Age Distribution
# In the following, the ginv function inverts a matrix by calculating the 
# Moore-Penrose generalized inverse of a matrix.
V <- Conj(ginv(EigVecs))   # left eigenvector; NOTE this notation from H Caswell
v <- Re(t(t(V[pos,])))     # dominant left eigenvector
rv <- v/(sum(v))
rv <- round(rv,3)          # scaled to provide proportional Reproductive Values

# Sensitivity analysis of Lambda to unit changes in matrix elements; note that
# %*% in R implies matrix multiplication, and * implies scalar multiplication
# or element-by-element multiplication of matrix entries.
senmat <- v%*%t(w)         # raw sensitivity matrix
senmat[A==0] <- 0          # puts 0s in locations where vital rate does not exist
senmat <- round(senmat,3)
senmat

# Elasticity analysis of Lambda to proportional changes in matrix elements
emat <- A/Lambda*(v%*%t(w))
emat <- round(emat,3)
emat

# To calculate sensitivities and elasticities for the lower-level vital rates
# we need some elaborate code to implement symbolic algebra. See a nice
# paper by Ezard et al. 2010 in J. Appl. Ecol. with such code.
# It turns out that the 'popbio' package already has this implemented. So lets
# gain some practice using this package to calculate lower-level sensitivities
# and elasticities

# Just put the vital rates in a list, and write the matrix as an expression
stork.vr <- list(u=0.45,r=0.818,b=2.9,s0=0.482,s1=0.75)
stork.A <- expression(0, 0, s0*u*r*b/2, s0*r*b/2,
  s1, 0, 0, 0,
  0, s1, 0, 0,
  0, 0, s1, s1)

# then apply the following popbio function
llsenselas <- vitalsens(stork.A,stork.vr)
llsenselas