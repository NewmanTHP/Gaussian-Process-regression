library(ggplot2)
library(data.table)
library(MASS)
library(corrplot)
library(pracma)

## creating dataset ------------------------------------------------------------- 
y.original <- c(11, 12, 12.5, 13, 12.3, 11.2, 11.9, 12.4, 11.4, 13.1, 13.2, 12.9, 13.3)
y = y.original-mean(y.original) # de-meaning data
x <- c(1, 2, 3, 3.7, 4.2, 5.25, 7, 9, 9.9, 10.8, 11.6, 12, 13.1)
data <- as.data.table(cbind(x,y))

# creating test input data
test.inputs = seq(0,14,length.out=500)

## covariance functions --------------------------------------------------------

# squared exponential 
se_kernel <- function(x1, x2, zigma = 1, l = 1) {
  zigma^2 * exp(- (x1 - x2)^2 / (2 * l^2))
}


# Matern 
matern_kernel <- function(x1, x2, nu = 1.5, zigma = 1, l = 1) {
  if (!(nu %in% c(0.5, 1.5, 2.5))) {
    stop("p must be equal to 0.5, 1.5 or 2.5")
  }
  p <- nu - 0.5
  d <- abs(x1 - x2)
  if (p == 0) {
    zigma^2 * exp(- d / l)
  } else if (p == 1) {
    zigma^2 * (1 + sqrt(3)*d/l) * exp(- sqrt(3)*d/l)
  } else {
    zigma^2 * (1 + sqrt(5)*d/l + 5*d^2 / (3*l^2)) * exp(-sqrt(5)*d/l)
  }
}

## creating covariance matrices ------------------------------------------------

# hyperparameters:
l=0.8
nu = 1.5
zigma = 1.3

#covariance matrix K(X*, X*)
covariance.matrix.xstar.xstar = matrix(0,500,500)
for (i in 1:length(test.inputs)) {
  for (j in 1:length(test.inputs)) {
    
    covariance.matrix.xstar.xstar[i,j] <- se_kernel(l=l,test.inputs[i], test.inputs[j])
    
  }
}

#covariance matrix K(X,X) 
covariance.matrix.x.x = matrix(0,13,13)
for (i in 1:length(data$y)) {
  for (j in 1:length(data$y)) {
    
    covariance.matrix.x.x[i,j] <- se_kernel(l=l,data$x[i], data$x[j])
    
  }
}

#covariance matrix K(X,X*) 
covariance.matrix.x.xstar = matrix(0,13,500)
for (i in 1:length(data$y)) {
  for (j in 1:length(test.inputs)) {
    
    covariance.matrix.x.xstar[i,j] <- se_kernel(l=l,data$x[i], test.inputs[j])
    
  }
}

#covariance matrix K(X*,X) 
covariance.matrix.xstar.x = matrix(0,500,13)
for (i in 1:length(test.inputs)) {
  for (j in 1:length(data$y)) {
    
    covariance.matrix.xstar.x[i,j] <- se_kernel(l=l,test.inputs[i], data$x[j])
    
  }
}



################################################################################
#####                              PLOTS                                   #####
################################################################################

### Figure 1 :  Discrete time series data plot containing 13 observations
plot(x,y.original, xlim=c(0,14), ylim=c(10,14), pch=19, xlab = "time")


### Figure 2 :  Three randomly drawn functions from GP prior
# (Figure 2 in the associated paper is produced using 50 test inputs instead of 500)
matrix.of.random.values.from.GP.prior = matrix(0,500,1000)
set.seed(420)
for (i in 1:1000) {
  dat = mvrnorm(1, rep(0, times=500), covariance.matrix.xstar.xstar)
  for (j in 1:length(test.inputs)) {
    matrix.of.random.values.from.GP.prior[j,i] = dat[j]
  }
}

matrix.of.random.values.from.GP.prior = as.data.frame(matrix.of.random.values.from.GP.prior)
matrix.of.random.values.from.GP.prior$mean <- apply(matrix.of.random.values.from.GP.prior, 1, mean)
matrix.of.random.values.from.GP.prior$sd <- apply(matrix.of.random.values.from.GP.prior, 1, sd)
matrix.of.random.values.from.GP.prior$lowerbound <- matrix.of.random.values.from.GP.prior$mean - 2*matrix.of.random.values.from.GP.prior$sd
matrix.of.random.values.from.GP.prior$upperbound <- matrix.of.random.values.from.GP.prior$mean + 2*matrix.of.random.values.from.GP.prior$sd

lines(x=test.inputs, y=matrix.of.random.values.from.GP.prior$lowerbound, col=1)
lines(x=test.inputs, y=matrix.of.random.values.from.GP.prior$upperbound, col=1)

# plotting
set.seed(420)
# plotting area
plot(y=mvrnorm(1, rep(0, times=500), covariance.matrix.xstar.xstar), x=test.inputs, pch=19, cex = 0.1, col=4, xlim=c(0,14), ylim=c(-3,3), ylab="y", xlab = "time")
# shaded area
polygon(c(test.inputs, rev(test.inputs)), c(matrix.of.random.values.from.GP.prior$lowerbound, rev(matrix.of.random.values.from.GP.prior$upperbound)),
        col = "#E8E2E2", aplha=0.4)
# prior samples
j=1
set.seed(420)
for (i in 1:3){
  if(i==1){
    points(y=mvrnorm(1, rep(0, times=500), covariance.matrix.xstar.xstar), x=test.inputs, pch=19, cex = 0.4, col=4)
    j = j+1
  }else{
    lines(y=mvrnorm(1, rep(0, times=500), covariance.matrix.xstar.xstar), x=test.inputs,col=j, lwd=2)
    j = j+1
  }
}



### Figure 3 : Three randomly drawn functions from GP posterior

## A. noise free GP
sigma = 1e-6

# Cholesky factorisation 
L = t(chol(covariance.matrix.x.x + diag(x=sigma,nrow = length(data$y),ncol = length(data$y))))
alpha = mldivide(t(L), mldivide(L,data$y))
v = mldivide(L,covariance.matrix.x.xstar)
mean.pred = t(covariance.matrix.x.xstar)%*%alpha
cov.pred = covariance.matrix.xstar.xstar - t(v)%*%v

matrix.of.random.values.from.GP.posterior = matrix(0,500,100)
set.seed(420)
for (i in 1:100) {
  dat = mvrnorm(1, mu=mean.pred, Sigma = cov.pred)
  for (j in 1:length(test.inputs)) {
    matrix.of.random.values.from.GP.posterior[j,i] = dat[j]
  }
}

matrix.of.random.values.from.GP.posterior = as.data.frame(matrix.of.random.values.from.GP.posterior)
matrix.of.random.values.from.GP.posterior$mean <- apply(matrix.of.random.values.from.GP.posterior, 1, mean)
matrix.of.random.values.from.GP.posterior$sd <- apply(matrix.of.random.values.from.GP.posterior, 1, sd)
matrix.of.random.values.from.GP.posterior$lowerbound <- matrix.of.random.values.from.GP.posterior$mean - 2*matrix.of.random.values.from.GP.posterior$sd
matrix.of.random.values.from.GP.posterior$upperbound <- matrix.of.random.values.from.GP.posterior$mean + 2*matrix.of.random.values.from.GP.posterior$sd


set.seed(420)
# plotting area
plot(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, pch=19, cex = 0.1, ylab="y", xlab = "time", ylim = c(-2,2))
# shaded area
polygon(c(test.inputs, rev(test.inputs)), c(matrix.of.random.values.from.GP.posterior$lowerbound, rev(matrix.of.random.values.from.GP.posterior$upperbound)),
        col = "#E8E2E2", aplha=0.4)
# posterior samples
j=1
set.seed(420)
for (i in 1:3){
  if(i==1){
    points(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, pch=19, cex = 0.4, col=4)
    j = j+1
  }else{
    lines(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, lwd=2, col=j)
    j = j+1
  }
}
# data observed
points(x,y, pch=19, cex=0.8)


## B. with noise GP
sigma = 1e-2

# Cholesky factorisation 
L = t(chol(covariance.matrix.x.x + diag(x=sigma,nrow = length(data$y),ncol = length(data$y))))
alpha = mldivide(t(L), mldivide(L,data$y))
v = mldivide(L,covariance.matrix.x.xstar)
mean.pred = t(covariance.matrix.x.xstar)%*%alpha
cov.pred = covariance.matrix.xstar.xstar - t(v)%*%v

matrix.of.random.values.from.GP.posterior = matrix(0,500,100)
set.seed(420)
for (i in 1:100) {
  dat = mvrnorm(1, mu=mean.pred, Sigma = cov.pred)
  for (j in 1:length(test.inputs)) {
    matrix.of.random.values.from.GP.posterior[j,i] = dat[j]
  }
}

matrix.of.random.values.from.GP.posterior = as.data.frame(matrix.of.random.values.from.GP.posterior)
matrix.of.random.values.from.GP.posterior$mean <- apply(matrix.of.random.values.from.GP.posterior, 1, mean)
matrix.of.random.values.from.GP.posterior$sd <- apply(matrix.of.random.values.from.GP.posterior, 1, sd)
matrix.of.random.values.from.GP.posterior$lowerbound <- matrix.of.random.values.from.GP.posterior$mean - 2*matrix.of.random.values.from.GP.posterior$sd
matrix.of.random.values.from.GP.posterior$upperbound <- matrix.of.random.values.from.GP.posterior$mean + 2*matrix.of.random.values.from.GP.posterior$sd

set.seed(420)
# plotting area
plot(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, pch=19, cex = 0.1, ylab="y", xlab = "time", ylim = c(-2,2))
# shaded area
polygon(c(test.inputs, rev(test.inputs)), c(matrix.of.random.values.from.GP.posterior$lowerbound, rev(matrix.of.random.values.from.GP.posterior$upperbound)),
        col = "#E8E2E2", aplha=0.4)
# posterior samples
j=1
set.seed(420)
for (i in 1:3){
  if(i==1){
    points(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, pch=19, cex = 0.4, col=4)
    j = j+1
  }else{
    lines(y=mvrnorm(1, mu=mean.pred, Sigma = cov.pred), x=test.inputs, lwd=2, col=j)
    j = j+1
  }
}
# data observed
points(x,y, pch=19, cex=0.8)

