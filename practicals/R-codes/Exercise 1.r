# Matrix Models for Population Management & Conservation
# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons 

# EXERCISE 1: barn swallows

################################################################################
# step 3: simple matrix model projections                                      #
################################################################################

# Define parameters
f1 <- 1.5
f2 <- 3
s0 <- 0.2
s1 <- 0.5
s2 <- 0.65

# Create the pre birth-pulse swallow matrix population model
A <- matrix(c(
  s0*f1, s0*f2,
  s1,    s2), nrow = 2, byrow = TRUE)
  
tspan <- 20                         # time span for projections
rows <- dim(A)[1]

# Build some matrices for storing eventual output
n <- matrix(0,rows,tspan)     # storage of age-specific abundance
N <- matrix(0,tspan,1)        # storage of total abundance
gr <- matrix(0,tspan-1,1)     # storage of time-specific population growth rates

n[,1] <- c(10,0)              # initial population abundance in each age class

# Project population forward and store output; note that
# %*% in R implies matrix multiplication, and * implies scalar multiplication
# or element-by-element multiplication of matrix entries.
for (t in 1:(tspan-1)) {
    n[,t+1] <- A%*%n[,t]      # %*% = matrix multiplication in R
    N[t+1] <- sum(n[,t+1])
    gr[t] <- sum(n[,t+1])/sum(n[,t])  # per time step population growth rate                                        
}

n
N
gr  

################################################################################
# step 4: projection over different initial conditions                         #
################################################################################

rm(list=ls(all=TRUE))   # clears the R memory; to give us a fresh start

# Define parameters
f1 <- 1.5
f2 <- 3
s0 <- 0.2
s1 <- 0.5
s2 <- 0.65

# Create the pre birth-pulse swallow matrix population model
A <- matrix(c(
  s0*f1, s0*f2,
  s1,    s2), nrow = 2, byrow = TRUE)
  
tspan <- 20                         # time span for projections
rows <- dim(A)[1]
cols <- dim(A)[2]

# Build some matrices for storing eventual output
n <- matrix(0,rows*2,tspan)   # storage of age-specific abundances
N <- matrix(0,tspan,2)     # storage of total abundances
gr <- matrix(0,tspan-1,2)  # storage of time-specific population growth rates

n[,1] <- c(10,0,0,10)         # initial population abundances in each age class
                              # for two different initial conditions

# Project population forward for each initial condition and store output
for (j in 1:2) {
  N[1,j] <- sum(n[(j*2-1):(j*2),1])
  for (t in 1:(tspan-1)) {
    n[(j*2-1):(j*2),t+1] <- A%*%n[(j*2-1):(j*2),t]      
    N[t+1,j] <- sum(n[(j*2-1):(j*2),t+1])
    gr[t,j] <- sum(n[(j*2-1):(j*2),t+1])/sum(n[(j*2-1):(j*2),t])                                          
  }
}

par(mfrow=c(3, 1))                # Set graphics window to 3 rows with 1 column
plot(1:tspan,n[1,],type="l",xlab="Time",ylab="age-specific abundance",
  ylim=c(0,15))
lines(1:tspan,n[2,],lty=3)
lines(1:tspan,n[3,],col="red")
lines(1:tspan,n[4,],lty=3,col="red")
plot(1:tspan,N[,1],type="l",xlab="Time",ylab="Total Abundance",ylim=c(0,30))
lines(1:tspan,N[,2],col="red")
plot(1:(tspan-1),gr[,1],type="l",xlab="Time",ylab="per time step growth rate",
  ylim=c(0.8,1.25))
lines(1:(tspan-1),gr[,2],col="red") 

################################################################################
# steps 5-6: on your own                                                       #
################################################################################
