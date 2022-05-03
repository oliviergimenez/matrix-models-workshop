# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# Exercise 9: Demographic stochasticity & DD for house sparrows

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Part 1 of code: visualizing demographic stochasticity in a density-dependent #
# population of house sparrows.                                                #
################################################################################
# define parameters
s0low <- 0.2     # first-year survival at low population density
s1low <- 0.5     # survival of older sparrows at low population density
b <- 0.01        # strength of density dependence
f <- 6/2         # fecundity: females produced per female

N <- matrix(0,100,1)
N[1] <- 10
for(t in 2:100){
  s0 <- s0low*exp(-b*N[t-1])
  s1 <- s1low*exp(-b*N[t-1])
  N[t] <- rpois(1,N[t-1]*s0*f)+ rbinom(1,N[t-1],s1)  #Dem stoch w/ branch proc
}

plot(1:100, N, xlab = list("Time",cex=1.5),ylab = list("Abundance",cex=1.5),
  col = "black",type = "l",lwd = 2,cex.axis = 1.1)
  
################################################################################
# Part 2 of code: extinction probabilities, mean extinction time, and the      #
# quasi-stationary distribution in the presence of demographic stochasticity   #
# and density-dependence.                                                      #
################################################################################
# Specify the number of simulations, time steps for each simulation, and the
# pseudo-extinction threshold
sims <- 10000
tspans <- c(100,200,300,400)
threshold <- 1

# Storage for indicators on each simulation determining whether or not the
# population ever dropped below the pseudo-extinction threshold.
ext_ind <- matrix(0,sims,length(tspans))

# Storage for times when simulations drop below the extinction threshold
ext_tm <- matrix(0,sims,length(tspans))

# Storage for quasi-stationary distribution of abundance for simulations that
# do not go extinct
QSD_n <- matrix(0,sims,length(tspans))

# Storage for summary output. Column 1 = extinction probabilities, Column 2 =
# extinction times, Column 3 = quasi-stationary distributions, Rows = time
# spans examined
sumtable <- matrix(0,4,3)

# define parameters
s0low <- 0.2     # first-year survival at low population density
s1low <- 0.5     # survival of older sparrows at low population density
b <- 0.01        # strength of density dependence
f <- 6/2         # fecundity: females produced per female

for (i in 1:4){
  tspan <- tspans[i]
  for (j in 1:sims){
    N <- matrix(0,tspan,1)
    N[1] <- 10
    for(t in 2:tspan){
      s0 <- s0low*exp(-b*N[t-1])
      s1 <- s1low*exp(-b*N[t-1])
      N[t] <- rpois(1,N[t-1]*s0*f)+ rbinom(1,N[t-1],s1)
    }
    if(min(N) < threshold) ext_ind[j,i] = 1
    if(min(N) < threshold) ext_tm[j,i] = which.min(N)
    if(min(N) > threshold) QSD_n[j,i] = N[tspan]
  }
  sumtable[i,1] <- mean(ext_ind[,i])
  ext_tmtemp <- subset(ext_tm,ext_tm[,i] > 0)
  sumtable[i,2] <- mean(ext_tmtemp[,i])
  QSD_ntemp <- subset(QSD_n,QSD_n[,i] > 0)
  sumtable[i,3] <- mean(QSD_ntemp[,i])  # note: variance and range also possible
}

sumtable

################################################################################
# Part 3 of code: examining the impact of alternative vital rate values        #
# on stochastic population outcomes.                                           #
################################################################################
# Specify the number of simulations, time steps for each simulation, and the
# pseudo-extinction threshold
sims <- 10000
tspan <- 1000
threshold <- 1

# Storage for indicators on each simulation determining whether or not the
# population ever dropped below the pseudo-extinction threshold.
ext_ind <- matrix(0,sims,1)

# Storage for quasi-stationary distribution of abundance for simulations that
# do not go extinct
QSD_n <- matrix(0,sims,1)

# define parameters: note that s0low is defined within the time loop below
s1low <- 0.5     # survival of older sparrows at low population density
b <- 0.01        # strength of density dependence
f <- 6/2         # fecundity: females produced per female

for (j in 1:sims){
  N <- matrix(0,tspan,1)
  N[1] <- 10
  for(t in 2:tspan){
    s0low <- 0.2              # first-year survival at low population density
    #s0low <- rbeta(1,18.5,55.5)     # env stochasticity for first-year survival
    s0 <- s0low*exp(-b*N[t-1])
    s1 <- s1low*exp(-b*N[t-1])
    N[t] <- rpois(1,N[t-1]*s0*f)+ rbinom(1,N[t-1],s1)
  }
  if(min(N) < threshold) ext_ind[j,1] = 1
  if(min(N) > threshold) QSD_n[j,1] = N[tspan]
}
ext_prob <- mean(ext_ind[,1])
QSD_ntemp <- subset(QSD_n,QSD_n[,1] > 0)
QSD_mean <- mean(QSD_ntemp[,1])
ext_prob
QSD_mean


# Aside: this code can be used to calculate shape parameters of a
# Beta distribution with specified mean and variance
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams(0.25,0.0025)

