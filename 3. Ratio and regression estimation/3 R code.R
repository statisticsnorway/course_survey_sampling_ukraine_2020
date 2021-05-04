# -----------------------------------------------
# R-code for Solution ark 3, Ratio and regression estimation. Autumn 2020
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------

install.packages("survey")
library(survey)
# ----------
# Exercise 2
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'agsrs.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

yi <- sample.data$acres92
x1i <- sample.data$acres87
x2i <- sample.data$farms87

tx1 <- 964470625
tx2 <- 2087759
x1bar <- mean(x1i)
x1bar

x2bar <- mean(x2i)
x2bar

n <- nrow(sample.data)
N <- 3078

# a).
par(mfrow=c(1,2))
plot(x1i,yi,ylab='acres92',xlab='acres87')
plot(x2i,yi,ylab='acres92',xlab='farms87')

# Sample correlations
cor(x1i,yi)
cor(x2i,yi)

# b).
# With survey package
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
ratio.est1 <- svyratio(~acres92,~acres87,design=des)
ratio.est2 <- svyratio(~acres92,~farms87,design=des)

# Ratios
Bhat1 <- coef(ratio.est1)
Bhat1
Bhat2 <- coef(ratio.est2)
Bhat2

# Ratio estimates
Bhat1*tx1
Bhat2*tx2


# Or
predict(ratio.est1, total=tx1)$total
predict(ratio.est2, total=tx2)$total

# Variance of residuals
e1i <- yi-coef(ratio.est1)*x1i
s2.e1 <- var(e1i)

e2i <- yi-coef(ratio.est2)*x2i
s2.e2 <- var(e2i)

# SEs of the ratio estimates
N*sqrt((1-n/N)*s2.e1/n)
N*sqrt((1-n/N)*s2.e2/n)

# SE with survey package
# These are higher than those above as
# the ratio tx/(N*xbar) is used in estimation of se(Bhat)
SE(ratio.est1)*tx1
SE(ratio.est2)*tx2


# 
sqrt(ratio.est1$var*tx1*tx1/(tx1^2/(x1bar*N)^2))
sqrt(ratio.est2$var*tx2*tx2/(tx2^2/(x2bar*N)^2))

# c).
# Regression estimation
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
svytotal(~acres92, calibrate(des, ~farms87, pop=c(N,tx2)))

# Alternatively
reg.est <- svyglm(acres92~farms87, design=des,data=sample.data)
coef.reg <- coef(reg.est)
B0hat <- coef.reg[1]
B1hat <- coef.reg[2]

# Regression estimate
(B0hat+B1hat*tx2/N)*N

ei <- yi-B0hat-B1hat*x2i
s2.e <- var(ei)

# SE of the regression estimate
N*sqrt((1-n/N)*s2.e/n)

# Standard error of simple sample mean
mean(yi)*N
N*sqrt(1-n/N)*sd(yi)/sqrt(n)

# ----------
# Exercise 3
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'counties.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

N <- 3141
n <- nrow(sample.data)

xi <- sample.data$totpop

xbar <- mean(xi)

tx <- 255077536


# a).
par(mfrow=c(1,2))
with(sample.data,hist(physician,xlab='number of physicians',main=''))
with(sample.data,hist(log(physician+1),xlab='log of number of physicians',main=''))

# b).
# SRS estimate:N*ybar
yi <- sample.data$physician
mean(yi)*N

# SE of ybar
N*sqrt((1-n/N)*var(yi)/n)



# Minimum sample size required for a normal distribution
library(e1071)
28+25*skewness(yi)^2

# c).
sample.data$County[yi==max(yi)]
# Number of physicians versus total population 
par(mfrow=c(1,2))
plot(xi,yi,xlab='total population',ylab='number of physicians')
# When Cook county is excluded
select <- yi<max(yi)
plot(xi[select],yi[select],xlab='total population',ylab='number of physicians')

# d).
# Ratio estimation
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
ratio.est <- svyratio(~physician,~totpop,design=des)

# Ratios
Bhat <- coef(ratio.est)
Bhat


# Ratio estimates
Bhat*tx



# Variance of residuals
ei <- yi-coef(ratio.est)*xi
s2.e <- var(ei)


# SE of the ratio estimate
N*sqrt((1-n/N)*s2.e/n)

# Using an alternative variance estimator (Lohr, 2019, p.126)
N*sqrt((1-n/N)*s2.e/n)*tx/xbar/N


# This is equivalent to
sqrt(ratio.est$var*tx*tx)

# Regression estimation
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
svytotal(~physician, calibrate(des, ~totpop, pop=c(N,tx)))

# Alternatively
reg.est <- svyglm(physician~totpop, design=des,data=sample.data)
coef.reg <- coef(reg.est)
B0hat <- coef.reg[1]
B1hat <- coef.reg[2]

# Regression estimate
(B0hat+B1hat*tx/N)*N

ei <- yi-B0hat-B1hat*xi
s2.e <- var(ei)

# SE of the regression estimate
N*sqrt((1-n/N)*s2.e/n)

#f)
xi <- sample.data$landarea

xbar <- mean(xi)

tx <- 3536278

# f.a).
par(mfrow=c(1,1))
with(sample.data,hist(farmpop,xlab='farm population',main=''))

# f.b).
# SRS estimate:N*ybar
yi <- sample.data$farmpop
mean(yi)*N

# SE of ybar
N*sqrt((1-n/N)*var(yi)/n)

# f.c).
# Farm population versus land area
plot(xi,yi,xlab='land area',ylab='farm population')
cor(xi,yi)


# f.d).
# Ratio estimation
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
ratio.est <- svyratio(~farmpop,~landarea,design=des)

# Ratios
Bhat <- coef(ratio.est)
Bhat


# Ratio estimates
Bhat*tx


# Variance of residuals
ei <- yi-coef(ratio.est)*xi
s2.e <- var(ei)


# SE of the ratio estimate
N*sqrt((1-n/N)*s2.e/n)

# Using an alternative variance estimator (Lohr, 2019, p.126)
N*sqrt((1-n/N)*s2.e/n)*tx/xbar/N


# This is equivalent to
sqrt(ratio.est$var*tx*tx)


# Regression estimation
des <- svydesign(data=sample.data,prob=~1,ids=1:n,fpc=rep(n/N,n))
svytotal(~farmpop, calibrate(des, ~landarea, pop=c(N,tx)))

# Alternatively
reg.est <- svyglm(farmpop~landarea, design=des,data=sample.data)
coef.reg <- coef(reg.est)
B0hat <- coef.reg[1]
B1hat <- coef.reg[2]

# Regression estimate
(B0hat+B1hat*tx/N)*N

ei <- yi-B0hat-B1hat*xi
s2.e <- var(ei)

# SE of the regression estimate
N*sqrt((1-n/N)*s2.e/n)

# ----------
# Exercise 5
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
pop.data <- read.csv(paste(libname,'pop_industry.csv',sep=''),sep=",")
head(pop.data)

#Define known values from problem statement
sample.size <- 25
N <- nrow(pop.data)
#1)

#a.
#Get indexes of those units which have more than 50 employees
idx_more50 <- seq(1,N)[pop.data$emplnr > 50]
length(idx_more50)
#get total value for the companies which have more than 50 employees
t_more50 <- sum(pop.data$turnover[idx_more50])
t_more50

##Since we need to repeat (b) - (e) 10 times, let's create vectors which will remember results from each iteration and run a loop
results.t_R <- c()
results.t_ex <- c()
results.SE_R <- c()
results.SE_ex <- c()
for(i in seq(1,10)){
  #b.
  #Get indexes of those units which have less or equal 50 employees
  idx_leq50 <- seq(1,N)[pop.data$emplnr <= 50]
  #Take an SRS from those which have less or equal 50 employees. We set seed to get the same results
  set.seed(2020 + i)
  idx_s <- sample(idx_leq50, sample.size - length(idx_more50))
  idx_s
  
  #Ratio estimator for those who has less or equal 50 employees
  t_R_leq50 <- sum(pop.data$emplnr[idx_leq50])* mean(pop.data$turnover[idx_s]) / mean(pop.data$emplnr[idx_s])
  
  #c.
  #Ratio estimator for the population
  t_R <- t_more50 + t_R_leq50
  results.t_R <- c(results.t_R, t_R)
  
  #Standard error for ratio estimator
  n_s <- length(idx_s)
  n_more50 <- length(idx_more50)
  
  B_est <- mean(pop.data$turnover[idx_s]) / mean(pop.data$emplnr[idx_s])
  s_e2 <- 1/(n_s - 1) * sum((pop.data$turnover[idx_s] - B_est * pop.data$emplnr[idx_s])^2)
  s_e <- sqrt(s_e2)
  
  SE_R <- (N - n_more50)*sqrt(1-n_s/(N - n_more50))* 
    mean(pop.data$emplnr[idx_leq50]) / mean(pop.data$emplnr[idx_s]) * s_e / sqrt(n_s)
  results.SE_R <- c(results.SE_R, SE_R)
  
  #d.
  #Expnasion estimator
  t_ex <- t_more50 + (N - n_more50) * mean(pop.data$turnover[idx_s])
  results.t_ex <- c(results.t_ex, t_ex)
  
  #SE for expansion estimator
  s2 <- 1/(n_s - 1) * sum((pop.data$turnover[idx_s] - mean(pop.data$turnover[idx_s]))^2) 
  s <- sqrt(s2)
  
  SE_ex <- (N - n_more50)*sqrt(1-n_s/(N - n_more50)) * s * 1/sqrt(n_s)
  results.SE_ex <- c(results.SE_ex, SE_ex)
  
}

#Ratio estimator for 10 iterations 
results.t_R

#Confidense inerval for Ratio estimator
CI_R <- results.t_R + matrix(results.SE_R, ncol = 1) %*% matrix(c(-1.96, 1.96), ncol = 2)
CI_R

#Expantion estimator for 10 iterations
results.t_ex

#Confidence interval for expansion estimator
CI_ex <- results.t_ex + matrix(results.SE_ex, ncol = 1) %*% matrix(c(-1.96, 1.96), ncol = 2)
CI_ex

#Combine the table to fill in
cbind(results.t_R, CI_R, results.t_ex, CI_ex)

#Calculate empiric standard errors
sqrt(var(results.t_R))
sqrt(var(results.t_ex))