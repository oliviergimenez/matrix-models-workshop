# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# Exercise 6: Exploitation in a matrix model with two sexes

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Find the harvest rate of males in a polygynous red deer population that      #
# leads to stationary population growth.                                       #
################################################################################
library(MASS) # an R package needed for some matrix functions

tspan <- 1000    # We have to project the population to measure the population
                 # growth rate because of the frequency-dependent
                 # breeding probability

# Define demographic parameters
fec <- 0.8    # typical fecundity
a2 <- 0.6     # proportion of breeders at age 2
sr <- 0.6     # sex ratio: % females at birth
s0 <- 0.6     # offspring survival to first birthday
s <- 0.92     # survival of fully grown females
q <- 0.85     # natural survival of prime-age males in absence of hunting
d <- 0.75     # natural survival of old males in absence of hunting
f <- fec*sr*s0      # net fertility of females produced per female
m <- fec*(1-sr)*s0  # net fertility of males produced per female
b <- 1        # parameter for frequency-dependent breeding probability

# vector of male harvest rates to loop over
hr <- seq(0,0.4,by=0.01)

# storage place for asymptotic population growth rates
lambda <- matrix(0,length(hr),1)
lambda_add <- matrix(0,length(hr),1)

for (i in 1:length(hr)){
  h <- hr[i]
  # Define vector of initial abundance
  n <- rep(1/24,24)
  # Temporary storage place for per time step growth rates
  gr <- matrix(0,tspan-1,1)
  for(t in 2:tspan){
    # Define number of breeding females and males
    nf <- sum(n[2:12])
    nm <- sum(n[16:24])
    # Define frequency-dependent breeding probability
    p <- 1/(1 + exp(-b*nm/nf))
    A <- matrix(c(
    0,a2*p*f,p*f,p*f,p*f,p*f,p*f,p*f,p*f,p*f,p*f,p*f,0,0,0,0,0,0,0,0,0,0,0,0,
    s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,s,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,a2*p*m,p*m,p*m,p*m,p*m,p*m,p*m,p*m,p*m,p*m,p*m,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,q*(1-h),0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q*(1-h),0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q*(1-h),0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q*(1-h),0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q*(1-h),0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,d*(1-h),0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,d*(1-h),0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,d*(1-h),0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,d*(1-h),0),nrow=24, byrow=T)
    n <- A%*%n
    gr[t-1] <- sum(n)
    n <- n/gr[t-1]    # this is a trick to avoid numerical problems of very
                      # large or small abundance values and it does not affect
                      # the result: Caswell 2001
  }
  # the population growth rates stabilize long before 'tspan' years so the
  # last value is definitely the asymptotic population growth rate
  lambda[i] <- gr[tspan-1]
  # calculation of hypothetical population growth rate if harvest of males
  # were completely additive at the population level
  lambda_add[i] <- 1.0111253*(1-h)
}

par(mar = c(5, 6, 4, 2),mgp = c(3.5, 1, 0))
plot(hr, lambda, xlab = list("Male Harvest Rate", cex = 2), ylim=c(0.96,1.02),
  ylab = list("Population Growth Rate", cex = 2), col = "black", type = "l",
  lwd = 2, cex.axis = 1.5)
lines(hr, lambda_add, col="red",lwd = 2)
abline(h = 1, col = "darkgray", lty = 3, lwd = 2)
legend(x = 0.2, y = 0.98, legend = c("actual","additive","stationary"),
  lty = c(1,1,3),lwd = c(2,2,2), col = c("black","red","darkgray"),
  bty = "n", cex = 1.5)
  
lambda
hr[24]
