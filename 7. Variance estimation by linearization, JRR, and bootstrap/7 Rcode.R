# -----------------------------------------------
# Solution ark #7. Variance estimation
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------

# Load packages required
library(survey)
library(sampling)

# ----------
# Exercise 1
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'srs30.csv',sep=''),sep=',',header=FALSE)
head(sample.data)
dim(sample.data)

N <- 100
n <- nrow(sample.data)

yi <- sample.data[,1]


# Sample mean
ybar <- mean(yi)
ybar

# 1).
# Sample mean at jth jackknife replicate
ybar_j <- NULL
for(j in 1:30){
  temp <- mean(yi[-j])
  ybar_j <- c(ybar_j,temp)
}

# Mean of jackknife estimates
mean(ybar_j)

# Sample variance
s2 <- var(yi)
s2

# Jackknife variance
varhatJK <- (n-1)*sum((ybar_j-mean(ybar_j))^2)/n
varhatJK

# varhatJK is equiv. to s2/30
s2/30

# With fpc
varhatJKfpc <- varhatJK*(1-n/N)
varhatJKfpc

# 2).
# Bootstrap variance
Bvect <- seq(25,1000,by=25)
length(Bvect)

varhatbootMC <- NULL
varhatboot <- NULL

for(l in 1:length(Bvect)){
  
  B <- Bvect[l]
  
  thetabstarMC <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,n,replace=TRUE)
    temp.thetabstarMC <- mean(yi[s]) 
    thetabstarMC <- c(thetabstarMC,temp.thetabstarMC)
  }
  
  varhatbootMC <- c(varhatbootMC,var(thetabstarMC))  
  
  # The rescaling bootstrap the following, the naive bootstrap
  # is eqiv. to the rescaling boot, when n-1 units selected WR
  m <- n-1
  
  thetabstar <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,m,replace=TRUE)
    wiboot <- NULL
    for(j in 1:n){
      ri <- sum(s %in% j)
      temp <- n/(n-1)*ri*N/n
      wiboot <- c(wiboot,temp)
    }
    temp.thetabstar <- sum(yi*wiboot)/sum(wiboot) 
    thetabstar <- c(thetabstar,temp.thetabstar)
  }
  
  varhatboot <- c(varhatboot,var(thetabstar)) 
}


varhat.lim <- range(c(varhatbootMC,varhatboot,varhatJK,varhatJKfpc))

par(mfrow=c(1,1))
plot(Bvect,varhatboot,type='l',lty=1, xlim=range(Bvect), ylim=varhat.lim,
     ylab='Bootstrap estimate of variance',xlab='Number of replicates',col='blue')
lines(Bvect,varhatbootMC,type='l',lty=1,col='red')
abline(h=varhatJK)
text(x=500,y=0.52,labels=substitute(paste(hat(V)[JK],(bar(y)),' = 0.5326054',sep='')))
abline(h=varhatJKfpc,lty=1)
text(x=500,y=0.385,labels=expression(paste((1-f),hat(V)[JK],(bar(y)),' = 0.3728238',sep='')))

legend(775,0.75,legend=c("bootstrap, n","bootstrap, n-1"), col=c("red","blue"),lty=c(1,1),ncol=1)


# ----------
# Exercise 2
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(),"/")
sample.data <- read.csv(paste(libname,'agsrs.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

yi <- sample.data$acres92
xi <- sample.data$acres87

ybar <- mean(yi)
xbar <- mean(xi)


n <- nrow(sample.data)
N <- 3078


# Ratio estimator
Bhat <- ybar/xbar
Bhat

# Variance estimation with linearization
ei <- yi-Bhat*xi

# SEs of the ratio estimates#Since Xbar is known, 
# it can be used here instead of sample mean
varhatLIN <- (1-n/N)*var(ei)/n/xbar/xbar
varhatLIN 


# Jackknife
Bhat_j <- NULL
xi_j <- NULL
yi_j <- NULL
for(j in 1:n){
  xi_j <- xi[-j]
  yi_j <- yi[-j]
  
  temp <- mean(yi_j)/mean(xi_j)
  Bhat_j <- c(Bhat_j,temp)
}

mean(Bhat_j)

# Variance estimation with jackknife
varhatJK <- (n-1)*sum((Bhat_j-mean(Bhat_j))^2)/n
varhatJK 

# With fpc
varhatJKfpc <- varhatJK*(1-n/N)
varhatJKfpc

# Bootstrap
# Variance estimation with bootstrap
Bvect <- seq(25,1000,by=25)
length(Bvect)

varhatbootMC <- NULL
varhatboot <- NULL

for(l in 1:length(Bvect)){
  
  B <- Bvect[l]
  
  
  thetabstarMC <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,n,replace=TRUE)
    temp.thetabstarMC <- mean(yi[s])/mean(xi[s])
    thetabstarMC <- c(thetabstarMC,temp.thetabstarMC)
  }
  
  varhatbootMC <- c(varhatbootMC,var(thetabstarMC))  
  
  m <- n-1
  
  thetabstar <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,m,replace=TRUE)
    wiboot <- NULL
    for(j in 1:n){
      ri <- sum(s %in% j)
      temp <- n/(n-1)*ri*N/n
      wiboot <- c(wiboot,temp)
    }
    xbar.boot <- sum(xi*wiboot)/sum(wiboot) 
    ybar.boot <- sum(yi*wiboot)/sum(wiboot)
    temp.thetabstar <- ybar.boot/xbar.boot
    thetabstar <- c(thetabstar,temp.thetabstar)
  }
  
  varhatboot <- c(varhatboot,var(thetabstar)) 
}


varhat.lim <- range(c(varhatbootMC,varhatboot,varhatJK,varhatJKfpc,varhatLIN))

par(mfrow=c(1,1))
plot(Bvect,varhatboot,type='l',lty=1, xlim=range(Bvect), ylim=varhat.lim,
     ylab='Bootstrap estimate of variance',xlab='Number of replicates',col='blue')
lines(Bvect,varhatbootMC,type='l',lty=1,col='red')
abline(h=varhatJK)
text(x=500,y=3.68e-05,labels=substitute(paste(hat(V)[JK],(hat(B)),' = 3.707e-05',sep='')))
abline(h=varhatJKfpc,lty=1)
text(x=500,y=3.37e-05,labels=substitute(paste((1-f),' ',hat(V)[JK],(hat(B)),' = 3.346e-05',sep='')))
abline(h=varhatLIN)
text(x=500,y=3.28e-05,labels=substitute(paste(hat(V)[L],(hat(B)),' = 3.307e-05',sep='')))

legend(775,2.7e-05,legend=c("bootstrap, n-1","bootstrap, n"), col=c("blue","red"),lty=c(1,1),ncol=1)




# ----------
# Exercise 3
# ----------
N <- 1132
Xbar <- 10.3

yi <- c(125,119,83,85,99,117,69,133,154,168,61,80,114,147,122,106,82,88,97,99) 
xi <- c(12,11.4,7.9,9,10.5,7.9,7.3,10.2,11.7,11.3,5.7,8,10.3,12,9.2,8.5,7,10.7,9.3,8.2) 

n <- length(yi)
n

ybar <- mean(yi)
xbar <- mean(xi)

# Jackknife
# Regression estimation
B1hat <- sum((xi-xbar)*(yi-ybar))/sum((xi-xbar)^2)
B1hat

B0hat <- ybar-B1hat*xbar
B0hat

# Estimate from regression estimation
ybarreg <- B0hat+B1hat*Xbar
ybarreg

# Linearization
# SE of the regression estimator
ei <- (yi-ybar)-B1hat*(xi-xbar)

sereg <- sqrt((1-n/N)*var(ei)/n)
sereg

# With survey package
# Here, g-weights are used (see Lohr 2019, p.459): alternative variance estimator (Vhat.2)
des <- svydesign(data=data.frame(cbind(yi,xi),id=1:n,fpc=n/N), prob=~1,ids=~id,fpc=~fpc)
reg.est <- svyglm(yi~xi, design=des)
pop <- data.frame(xi=Xbar)


seregalt <- SE(svytotal(~yi, calibrate(des, ~xi, pop=c(1, Xbar))))
seregalt

# Jackknife
ybarreg_j <- NULL
xi_j <- NULL
yi_j <- NULL
for(j in 1:n){
  xi_j <- xi[-j]
  yi_j <- yi[-j]
  
  temp.B1hat <- sum((xi_j-mean(xi_j))*(yi_j-mean(yi_j)))/sum((xi_j-mean(xi_j))^2)
  temp.B0hat <- mean(yi_j)-temp.B1hat*mean(xi_j)
  temp <- temp.B0hat+temp.B1hat*Xbar
  ybarreg_j <- c(ybarreg_j,temp)
}

mean(ybarreg_j)

# Variance estimation with jackknife
seJK <- sqrt((n-1)*sum((ybarreg_j-mean(ybarreg_j))^2)/n)
seJK

# With fpc
seJKfpc <- sqrt(1-n/N)*seJK
seJKfpc

# Bootstrap
# Variance estimation with bootstrap
Bvect <- seq(25,1000,by=25)
length(Bvect)

sebootMC <- NULL
seboot <- NULL

for(l in 1:length(Bvect)){
  
  B <- Bvect[l]
  
  thetabstarMC <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,n,replace=TRUE)
    temp.B1hat <- sum((xi[s]-mean(xi[s]))*(yi[s]-mean(yi[s])))/sum((xi[s]-mean(xi[s]))^2)
    temp.B0hat <- mean(yi[s])-temp.B1hat*mean(xi[s])
    temp.thetabstarMC<- temp.B0hat+temp.B1hat*Xbar
    thetabstarMC <- c(thetabstarMC,temp.thetabstarMC)
  }
  
  sebootMC <- c(sebootMC,sd(thetabstarMC))  
  
  m <- n-1
  
  thetabstar <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,m,replace=TRUE)
    wiboot <- NULL
    for(j in 1:n){
      ri <- sum(s %in% j)
      temp <- n/(n-1)*ri*N/n
      wiboot <- c(wiboot,temp)
    }
    xbar.boot <- sum(xi*wiboot)/sum(wiboot) 
    ybar.boot <- sum(yi*wiboot)/sum(wiboot)
    temp.B1hat <- sum(wiboot*(xi-xbar.boot)*(yi-ybar.boot))/sum(wiboot*(xi-xbar.boot)^2)
    temp.B0hat <- ybar.boot-temp.B1hat*xbar.boot
    temp.thetabstar<- temp.B0hat+temp.B1hat*Xbar 
    thetabstar <- c(thetabstar,temp.thetabstar)
  }
  
  seboot <- c(seboot,sd(thetabstar)) 
}

se.lim <- range(c(sebootMC,seboot,seJK,seJKfpc,sereg,seregalt))

par(mfrow=c(1,1))
plot(Bvect,seboot,type='l',lty=1, xlim=range(Bvect), ylim=se.lim,
     ylab='Bootstrap estimate of the standard error',xlab='Number of replicates',col='blue')
lines(Bvect,sebootMC,type='l',lty=1,col='red')
abline(h=seJK)
text(x=100,y=5.42,labels=substitute(paste(se[JK],(bar(y)[reg]),' = 5.3846',sep='')))
abline(h=seJKfpc,lty=1)
text(x=600,y=5.30,labels=expression(paste((1-f),se[JK],(bar(y)[reg]),' = 5.3368',sep='')))
abline(h=seregalt)
text(x=600,y=5.03,labels=substitute(paste(se[Lalt],(bar(y)[reg]),' = 5.0654',sep='')))
abline(h=sereg)
text(x=600,y=4.00,labels=substitute(paste(se[L],(bar(y)[reg]),' = 3.9622',sep='')))

legend(775,4.15,legend=c("bootstrap, n-1","bootstrap, n"), col=c("blue","red"),lty=c(1,1),ncol=1)

# ----------
# Exercise 4
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'certify.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

# Treat the respondents as the sample
N <- 18609
n <- nrow(sample.data)

w0 <- N/n
w0

# 
sample.data$phd <- 1
sample.data$phd[sample.data$college!='P'] <- 2
sample.data$phd <- as.factor(sample.data$phd)

sample.data$work <- 3
sample.data$work[sample.data$workenv=='I'] <- 1
sample.data$work[sample.data$workenv=='A'] <- 2
sample.data$work <- as.factor(sample.data$work)


# Iterative raking with "sampling" package
# create a matirx of dummy variables based on phd and work
a <- model.matrix(~-1+sample.data$phd+sample.data$work)

# Compute g-weights, which are the coefficients that shall be used to adjust the initial weights
# Final weights = initial weights * g-weights
cal.g <- calib(a,d=rep(w0,n),total=c(10235,8374,6327,6885),method="raking")

# Raking weights
wirak <- cal.g*w0


# Estimated number of members opposing certification
Yhatrak <- sum(wirak[sample.data$certify==5])
Yhatrak

# Bootstrap
# Variance estimation with bootstrap
Bvect <- seq(25,1000,by=25)
length(Bvect)

m <- n-1

w0istar <- w0*n/m

varhatboot <- NULL

for(l in 1:length(Bvect)){
  
  B <- Bvect[l]
  
  thetabstar <- NULL
  for(i in 1:B){
    set.seed(170420+i)
    s <- sample(n,m,replace=TRUE)
    data.boot <- sample.data[s,]
    
    # Iterative raking with "sampling" package
    a.boot <- model.matrix(~-1+data.boot$phd+data.boot$work)
    cal.boot <- calib(a.boot,d=rep(w0istar,m),total=c(10235,8374,6327,6885),method="raking")
    # Raking weights
    wiboot <- cal.boot*w0istar
    
    # Estimated number of members opposing certification
    temp.thetabstar <- sum(wiboot[data.boot$certify==5])
    thetabstar <- c(thetabstar,temp.thetabstar)
  }
  
  varhatboot <- c(varhatboot,var(thetabstar)) 
}


varhat.lim <- range(varhatboot)

par(mfrow=c(1,1))
plot(Bvect,varhatboot,type='l',lty=1, xlim=range(Bvect), ylim=varhat.lim,
     ylab='Bootstrap estimate of variance',xlab='Number of replicates')

# Standard CI based on bootstrap estimate of the standard error
Yhatrak+qnorm(c(0.025,0.0975),0,1)*sqrt(varhatboot[length(Bvect)])


# Confidence interval based on percentiles from bootstrapping
quantile(thetabstar,c(0.025,0.0975))


