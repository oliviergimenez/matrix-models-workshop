# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# Exercise 8: The Ricker model of density dependence

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Part 1 of code: the Ricker model                                             #
################################################################################
tspan <- 100
N <- matrix(0,tspan,1)
N[1] <- 1
K <- 100
r <- 0.25
for (t in 2:tspan){
  N[t] <- N[t-1]*exp(r*(1-N[t-1]/K))
}

par(mar = c(5, 6, 4, 2),mgp = c(3.5, 1, 0))
plot(1:tspan, N, xlab = list("Time",cex=1.5),ylab = list("Abundance",cex=1.5),
  col = "black",type = "l",lwd = 2,cex.axis = 1.1)

################################################################################
# Part 2 of code: Bifurcation diagram for the Ricker model                     #
################################################################################
plot(-1, -1, xlim = c(0.1,4), ylim = c(0,600), xlab = "r", ylab = "N")
r <- seq(0.1, 4, by = 0.01)
for (i in 1:length(r)) {
    N <- vector()
    N[1] <- 50
    for (t in 2:tspan) {
      N[t] <- N[t-1]*exp(r[i]*(1-N[t-1]/K))
    }
    uval <- unique(N[50:tspan])
    points(rep(r[i], length(uval)), uval, cex = 0.1, pch = 19)
}

################################################################################
# Part 3 of code: DD Metapopulation model and bifurcation diagram              #
################################################################################
tspan <- 100
K <- 100
r <- 0.8
p <- 0.8
q <- 0.6
n <- matrix(0,2,tspan)
n[,1] <- rep(1,2)         # initial vector of abundance;
                          # row 1 for source, row 2 for sink
for (t in 2:tspan){
  A <- matrix(c(
    exp(r*(1-n[1,t-1]/K))*p, 0,
    exp(r*(1-n[1,t-1]/K))*(1-p), q),nrow=2, byrow=T)
  n[,t] <- A%*%n[,t-1]
}

par(mar = c(5, 6, 4, 2),mgp = c(3.5, 1, 0))
plot(1:tspan, n[1,], xlab = list("Time",cex=1.5),ylab = list("Abundance",cex=1.5),
  col = "black",type = "l",lwd = 2,cex.axis = 1.1)
lines(1:tspan, n[2,], col = "red", lwd = 2)
legend(x = 50, y = 25, legend = c("Source", "Sink"), lty = c(1, 1),
  lwd = c(2, 2), col = c("black", "red"), bty = "n", cex = 1.2)

### Bifurcation diagram for the sink population
plot(-1, -1, xlim = c(0.1,5), ylim = c(0,300), xlab = "r",ylab="Sink Abundance")
r <- seq(0.1, 5, by = 0.01)
for (i in 1:length(r)) {
  n <- matrix(0,2,tspan)
  n[,1] <- c(40,30)
  for (t in 2:tspan) {
    A <- matrix(c(
      exp(r[i]*(1-n[1,t-1]/K))*p, 0,
      exp(r[i]*(1-n[1,t-1]/K))*(1-p), q),nrow=2, byrow=T)
    n[,t] <- A%*%n[,t-1]
    }
    uval <- unique(n[2,50:tspan])
    points(rep(r[i], length(uval)), uval, cex = 0.1, pch = 19)
}

