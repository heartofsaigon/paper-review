
rm(list = ls())

if(!library(pacman, logical.return = T)) install.packages("pacman")
pacman::p_load(tidyverse, hetGP, laGP, plgp, HDInterval, matrixcalc)
par(mfrow = c(1,1))

##########################################################
##########################################################
##########################################################
# Add scale and nugget 
# the main function for illustration and its curve
simpleFunction <- function(x) -6*dnorm(x,4) - 7.5*dnorm(x,7)+1.3

curve(simpleFunction, xlim = c(0,10), ylim = c(-3,3))
# using standard R numerical algorithm to optimise the function
op<- optimise(simpleFunction, interval = c(0,10)); op
abline(v = op$minimum, h = op$objective)

# using one-step-ahead updates 
set.seed(1)

# prior - before observing any data point
X_pr = seq(0, 10, length = 100)|> sort()
D_pr = plgp::distance(X_pr) # (x1 - x2)^2
S_pr = exp(-D_pr) + diag(.Machine$double.eps, 100) # exponentiated sqaured distance + jitter

lower =  -qnorm(.975)*sqrt(diag(S_pr))
upper =   qnorm(.975)*sqrt(diag(S_pr))

curve(simpleFunction, xlim = c(0,10), ylim = c(-3,3), lwd = 3)
lines(X_pr, lower, col = "red", lty = 2)
lines(X_pr, upper, col = "red", lty = 2)


# data
X_dat <- runif(1,0, 10)|> sort() # randomly choose 1 point
Y_dat = simpleFunction(X_dat) # evaluate new point

D_dat = plgp::distance(X_dat) 
S_dat = exp(-D_dat) + diag(.Machine$double.eps,length(X_dat))
# posterior
D_pr_dat = plgp::distance(X_pr, X_dat)  # distance between prior and new data point
S_pr_dat = exp(-D_pr_dat) #  cross covariance matrix. 
covmat = S_pr - S_pr_dat%*%solve(S_dat)%*%t(S_pr_dat) # predictive covariance. no require evaluation at new point 
mu = S_pr_dat%*%solve(S_dat)%*%as.matrix(Y_dat) # the mean requires evaluating new point

lower = mu - qnorm(.975)*sqrt(diag(covmat))
upper = mu + qnorm(.975)*sqrt(diag(covmat))


curve(simpleFunction, xlim = c(0,10), ylim = c(-3,3), lwd = 3)
points(X_dat, Y_dat, pch = 19, col = "grey40")
lines(X_pr, mu, col = "red")
lines(X_pr, lower, col = "red", lty = 2)
lines(X_pr, upper, col = "red", lty = 2)
###

f_min = min(Y_dat)
EI = (f_min - mu)*pnorm((f_min-mu)/diag(covmat)) + 
  diag(covmat)*dnorm((f_min-mu)/diag(covmat))
#plot(X_pr, EI, type = "l") 
ind = which.max(EI)
X_dat = c(X_dat, X_pr[ind])

stopit = 0
repeat{
  Y_dat = simpleFunction(X_dat)
  
  D_dat = plgp::distance(X_dat)
  S_dat = exp(-D_dat) + diag(.Machine$double.eps, length(X_dat))
  
  # posterior
  D_pr_dat = plgp::distance(X_pr, X_dat)
  S_pr_dat = exp(-D_pr_dat)
  
  mu = S_pr_dat%*%solve(S_dat)%*%as.matrix(Y_dat)
  covmat = S_pr - S_pr_dat%*%solve(S_dat)%*%t(S_pr_dat)
  
  
  Sys.sleep(0.7)
  #par(mfrow = c(1,1))
  curve(simpleFunction, xlim = c(0,10), ylim = c(-3,3), lwd = 3)
  
  points(as.matrix(X_dat), as.matrix(Y_dat), pch = 16, col = "red")
  points(X_dat[length(X_dat)], Y_dat[length(Y_dat)], pch = 16, col = "blue")
  
  
  # lower = mu - qnorm(.975)*sqrt(diag(covmat))
  # upper = mu + qnorm(.975)*sqrt(diag(covmat))
  lines(X_pr, mu, col = "red", lty = 1)
  
  interval<-
    rmvnorm(2000, as.matrix(mu), covmat)|>
    apply(2, \(i) quantile(i, c(0.025, 0.975)))
  lower = interval[1,]
  upper = interval[2,]
  
  lines(X_pr, lower, col = "red", lty = 3)
  lines(X_pr, upper, col = "red", lty = 3)
  
  
  
  
  f_min = min(Y_dat)
  EI = (f_min - mu)*pnorm((f_min-mu)/diag(covmat)) + 
    diag(covmat)*dnorm((f_min-mu)/diag(covmat))
  #plot(X_pr, EI, type = "l") 
  
  ind = which.max(EI)
  j=1
  while(X_pr[ind] %in% X_dat){
    j =j+1
    ind = order(max(EI)- EI)[2]
    if(j>5){
      stopit = 1
      break
    }
  }
  if(stopit ==1) break
  X_dat = c(X_dat, X_pr[ind])
  
}
cbind(X_dat, Y_dat)[order(Y_dat),]
############
#######################
##############################################

##########################################################
##########################################################
##########################################################