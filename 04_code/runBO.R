rm(list = ls())

if(!library(pacman, logical.return = T)) install.packages("pacman")
pacman::p_load(tidyverse, hetGP, laGP, plgp, HDInterval, matrixcalc)
par(mfrow = c(1,1))


##########################################################
##########################################################
##########################################################
# Add scale and nugget 
# the main function for illustration and its curve
simpleFunction <- function(x) -6*dnorm(x,4) - 7.5*dnorm(x,7)+2

curve(simpleFunction, xlim = c(0,10), ylim = c(-3,3))
# using standard R numerical algorithm to optimise the function
op<- optimise(simpleFunction, interval = c(0,10)); op
abline(v = op$minimum, h = op$objective)

# using one-step-ahead updates 
set.seed(1)

# prior
X_pr = seq(0,10, length = 100)|> sort()
Y_pr = simpleFunction(X_pr)
D_pr = plgp::distance(X_pr)
S_pr = exp(-D_pr) + sqrt(.Machine$double.eps)
plot(X_pr, Y_pr, type = "l", ylim = c(-3,3), lwd = 3, col = "gray30")

# data
n = 1
X_dat <- sample(X_pr[X_pr<5],n)
Y_dat = simpleFunction(X_dat)
D_dat = plgp::distance(X_dat)

#---

S_dat = exp(-D_dat) + diag(.Machine$double.eps,n)

# posterior
D_pr_dat = plgp::distance(X_pr, X_dat)
S_pr_dat = exp(-D_pr_dat)

mu = S_pr_dat%*%solve(S_dat)%*%as.matrix(Y_dat)
covmat = S_pr - S_pr_dat%*%solve(S_dat)%*%t(S_pr_dat)
lower = mu - qnorm(.975)*sqrt(diag(covmat))
upper = mu + qnorm(.975)*sqrt(diag(covmat))


lines(X_pr, mu, col = "red")
points(X_dat, Y_dat, pch = 19, col = "grey40")
lines(X_pr, lower, col = "red", lty = 2)
lines(X_pr, upper, col = "red", lty = 2)
###

f_min = min(Y_dat)
EI = (f_min - mu)*pnorm((f_min-mu)/diag(covmat)) + 
  diag(covmat)*dnorm((f_min-mu)/diag(covmat))
#plot(X_pr, EI, type = "l") 
ind = which.max(EI)
X_dat = c(X_dat, X_pr[ind])

repeat{
  Y_dat = simpleFunction(X_dat)
#  X_pr1 = runif(1000, 0,10)|> sort()

  D_dat = plgp::distance(X_dat)
  S_dat = exp(-D_dat) + diag(1e-8, length(X_dat))
  
  # posterior
  D_pr_dat = plgp::distance(X_pr, X_dat)
  S_pr_dat = exp(-D_pr_dat)
  
  mu = S_pr_dat%*%solve(S_dat)%*%as.matrix(Y_dat)
  covmat = S_pr - S_pr_dat%*%solve(S_dat)%*%t(S_pr_dat)
  
  
  Sys.sleep(0.7)
  #par(mfrow = c(1,1))
  plot(X_pr, Y_pr, type = "l", ylim = c(-3,3), col = "gray30", lwd = 3, 
       main = paste0(op$minimum,"; ",X_dat[which.min(Y_dat)]) )
  abline(v = op$minimum)
  
  points(X_dat, Y_dat, pch = 16, col = "blue")
  points(X_dat[length(X_dat)], Y_dat[length(Y_dat)], pch = 16, col = "red")
  
  lower = mu - qnorm(.975)*sqrt(diag(covmat))
  upper = mu + qnorm(.975)*sqrt(diag(covmat))
  lines(X_pr, mu, col = "red", lty = 1)
  lines(X_pr, lower, col = "blue", lty = 3)
  lines(X_pr, upper, col = "blue", lty = 3)
  abline(v = X_dat[which.min(Y_dat)], lty = 2, col = "deeppink")
  
  # plot(X_pr, EI, type = "l")
  if(length(X_dat)==10){
      break
    
    
  }
  
  
  f_min = min(Y_dat)
  EI = (f_min - mu)*pnorm((f_min-mu)/diag(covmat)) + 
    diag(covmat)*dnorm((f_min-mu)/diag(covmat))
  #plot(X_pr, EI, type = "l") 
  ind = which(max(EI) - EI == 0)[1]
  X_dat = c(X_dat, X_pr[ind])
 # print(c(X_dat[length(X_dat)],Y_dat[length(Y_dat)]))

}
cbind(X_dat, Y_dat)[order(Y_dat),]

############
#######################
##############################################

##########################################################
##########################################################
##########################################################







