# Matrix Models for Population Management & Conservation
# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# EXERCISE 5: Herbivory and environmental variation in the common kidney vetch

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Grazed population                                                            #
################################################################################
library(MASS) # an R package needed for some matrix functions
# Finding the frequency of Poor years that can be sustained
tspan <- 100000  # following ergodic theory, we just need 1 long simulation to
                 # estimate the stochastic growth rate
# Define demographic parameters that do not vary over time
CS <- 0.40
C1 <- 0.19
C2 <- 0.40
CJ <- 0.10
C3 <- 0.84
C3p<- 0.84
C4 <- 0.50
C4p<- 0.50

# generate sequence of probabilities of Poor years to iterate over
P <- seq(0.8, 0.83, by = 0.005)
# storage place for stochastic growth rates across different values of P
Lambda_s <- matrix(0,length(P),1)

for (p in 1:length(P)){
  # Define vector of initial abundance
  n <- c(1/6,1/6,1/6,1/6,1/6,1/6)
  # Temporary storage place for per time step growth rates for eventual
  # calculation of the stochastic growth rates
  gr <- matrix(0,tspan-1,1)
  for (t in 2:tspan){
    X <- rbinom(1, 1, P[p])
    f1 <- X*0 + (1-X)*54
    f2 <- X*0.4 + (1-X)*82
    f2p <- X*0.7 + (1-X)*130
    f3 <- f2
    A <- matrix(c(
    CS, f1*C2, f2*C3, f2p*C3p, f3*C4, f3*C4p,
    C1, 0, 0, 0, 0, 0,
    0, C2, 0, 0, 0, 0,
    0, CJ, 0, 0, 0, 0,
    0, 0, 0, 0, C3p, 0,
    0, 0, 0, C3, 0, 0),nrow=6, byrow=T)
    n <- A%*%n
    gr[t-1] <- sum(n)
    n <- n/gr[t-1]    # this is a trick to avoid numerical problems of very
                      # large or small abundance values and it does not affect
                      # the result: Caswell 2001
  }
ln_lambda_s <- mean(log(gr))
Lambda_s[p] <- exp(ln_lambda_s)
}

plot(P,Lambda_s,xlab="Probability of Poor Year",ylab="Stochastic Growth Rate")


################################################################################
# Ungrazed population                                                          #
################################################################################
# Finding the frequency of Poor years that can be sustained
tspan <- 100000  # following ergodic theory, we just need 1 long simulation to
                 # estimate the stochastic growth rate
# Define demographic parameters that do not vary over time
CS <- 0.43
C1 <- 0.14
C2 <- 0.50
CJ <- 0.17
C3 <- 0.74
C3p <- 0.47
C4 <- 0.48
C4p <- 0.32

# generate sequence of probabilities of Poor years to iterate over
P <- seq(0.85, 0.88, by = 0.005)
# storage place for stochastic growth rates across different values of P
Lambda_s <- matrix(0,length(P),1)

for (p in 1:length(P)){
  # Define vector of initial abundance
  n <- c(1/6,1/6,1/6,1/6,1/6,1/6)
  # Temporary storage place for per time step growth rates for eventual
  # calculation of the stochastic growth rates
  gr <- matrix(0,tspan-1,1)
  for (t in 2:tspan){
    X <- rbinom(1, 1, P[p])
    f1 <- X*0 + (1-X)*72
    f2 <- X*1.0 + (1-X)*246
    f2p <- X*1.6 + (1-X)*370
    f3 <- f2
    A <- matrix(c(
    CS, f1*C2, f2*C3, f2p*C3p, f3*C4, f3*C4p,
    C1, 0, 0, 0, 0, 0,
    0, C2, 0, 0, 0, 0,
    0, CJ, 0, 0, 0, 0,
    0, 0, 0, 0, C3p, 0,
    0, 0, 0, C3, 0, 0),nrow=6, byrow=T)
    n <- A%*%n
    gr[t-1] <- sum(n)
    n <- n/gr[t-1]    # this is a trick to avoid numerical problems of very
                      # large or small abundance values and it does not affect
                      # the result: Caswell 2001
  }
ln_lambda_s <- mean(log(gr))
Lambda_s[p] <- exp(ln_lambda_s)
}

plot(P,Lambda_s,xlab="Probability of Poor Year",ylab="Stochastic Growth Rate")

################################################################################
# Reproductive value, stable stage distribution, and elasticity in a           #
# stochastic environment; example with ungrazed kidney vetch.                  #
################################################################################
rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful
rows <- 6
cols <- rows

tspan <- 100000

# Define demographic parameters that do not vary over time
CS <- 0.43
C1 <- 0.14
C2 <- 0.50
CJ <- 0.17
C3 <- 0.74
C3p <- 0.47
C4 <- 0.48
C4p <- 0.32

p <- 0.87

# Specify an arbitrary initial reproductive value vector and call it vvec.
# Specify an arbitrary initial population structure vector and call it wvec.
# (The long-term stochastic results are insensitive to initial population structure)
v <- matrix(0,rows,(tspan+1))
w <- matrix(0,rows,(tspan+1))
vvec <- (matrix(1,1,cols))/cols
v[,(tspan+1)] <- t(vvec)
wvec <- (matrix(1,rows,1))/rows
w[,1] <- wvec

# generate and store random values of fecundities
f1 <- f2 <- f2p <- f3 <- matrix(0,tspan,1)
for (i in tspan:1) {
    X <- rbinom(1, 1, p)
    f1[i] <- X*0 + (1-X)*72
    f2[i] <- X*1.0 + (1-X)*246
    f2p[i] <- X*1.6 + (1-X)*370
    f3[i] <- f2[i]
}

# generate and store the sequence of reproductive value vectors in reverse time
# ...strange but it works this way because r.v. is a measure of current and
# 'future' contribution to the population. One must start with a 'future'
# quantity and work backwards until stability is reached.
for (i in tspan:1) {
    A <- matrix(c(
    CS, f1[i]*C2, f2[i]*C3, f2p[i]*C3p, f3[i]*C4, f3[i]*C4p,
    C1, 0, 0, 0, 0, 0,
    0, C2, 0, 0, 0, 0,
    0, CJ, 0, 0, 0, 0,
    0, 0, 0, 0, C3p, 0,
    0, 0, 0, C3, 0, 0),nrow=6, byrow=T)
    vvec <- vvec%*%A
    vvec <- vvec/sum(vvec)
    v[,i] <- t(vvec)
}
## the reproductive values in a stochastic environment ##
RV <- c(mean(v[1,]),mean(v[2,]),mean(v[3,]),mean(v[4,]),mean(v[5,]),mean(v[6,]))
RV

# begin looping forward through time to calculate the stable stage distribution
# and population growth rate in a stochastic environment
gr <- matrix(0,tspan,1)

for (i in 1:tspan) {
    A <- matrix(c(
    CS, f1[i]*C2, f2[i]*C3, f2p[i]*C3p, f3[i]*C4, f3[i]*C4p,
    C1, 0, 0, 0, 0, 0,
    0, C2, 0, 0, 0, 0,
    0, CJ, 0, 0, 0, 0,
    0, 0, 0, 0, C3p, 0,
    0, 0, 0, C3, 0, 0),nrow=6, byrow=T)
    wvec <- A%*%wvec
    gr[i] <- sum(wvec)           # per time step population growth rate
    wvec <- wvec/gr[i]           # per time step measures of age distribution
    w[,(i+1)] <- wvec
}
# Stable stage distribution
SSD <- c(mean(w[1,]),mean(w[2,]),mean(w[3,]),mean(w[4,]),mean(w[5,]),mean(w[6,]))
ln_lambda_s = mean(log(gr))       # the stochastic population growth rate
Lambda_s = exp(ln_lambda_s)       # the stochastic population growth rate on the
                                  # non-logged scale
SSD
Lambda_s

# Elasticity analysis of Lambda_s to simultaneous proportional changes in
# both the mean and variance of matrix-entry vital rates
EMAT = matrix(0,rows,cols)        # initialize the elasticity matrix
# calculate summations
for (i in 1:tspan) {
    A <- matrix(c(
    CS, f1[i]*C2, f2[i]*C3, f2p[i]*C3p, f3[i]*C4, f3[i]*C4p,
    C1, 0, 0, 0, 0, 0,
    0, C2, 0, 0, 0, 0,
    0, CJ, 0, 0, 0, 0,
    0, 0, 0, 0, C3p, 0,
    0, 0, 0, C3, 0, 0),nrow=6, byrow=T)
    EMAT <- EMAT +(v[,(i+1)]%*%t(w[,i])*A)/as.numeric(gr[i]*t(v[,(i+1)])%*%w[,(i+1)])
    }
# divide summation by the time span
EMAT <- round(EMAT/(tspan+1),3)
EMAT

# a check to make sure calculations are correct
sum(sum(EMAT))