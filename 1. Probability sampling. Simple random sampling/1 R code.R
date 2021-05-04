# -----------------------------------------------
# R-code for Solution ark 1, Probability Sampling. Simple Random Sampling. Autumn 2020
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------

# ----------
# Exercise 3: Parts 3-5
# ----------
# Function to extract all subsets of a vector for a given size
# Alternatively: combn(N,n,simplify=FALSE)
enum.choose <- function(x, k) {
  if(k > length(x)) stop('k > length(x)')
  if(choose(length(x), k)==1){
    list(as.vector(combn(x, k)))
  } else {
    cbn <- combn(x, k)
    lapply(seq(ncol(cbn)), function(i) cbn[,i])
  }
}

# Population set
yi <- c(1,2,4,4,7,7,7,8)

# Population size
N <- 8
# Sample size
n <- 3


# For SRS without replacement # ----------
# Number of all possible samples
all.nr.samples.WOR <- choose(N,n)
all.nr.samples.WOR
# All possible samples
all.subsets.WOR <- enum.choose(1:N,n)
length(all.subsets.WOR)

# Obtain ybar_s WOR
ybar.WOR <- NULL
for(i in 1:all.nr.samples.WOR){
  s <- all.subsets.WOR[[i]]
  yi.s <- yi[s]
  temp <- mean(yi.s)
  ybar.WOR <- c(ybar.WOR,temp) 
}

# For SRS with replacement # ----------
# Number of all possible samples
all.nr.samples.WR <- N^n
all.nr.samples.WR
# All possible samples
all.subsets.WR <- expand.grid(rep(list(1:N), n))
dim(all.subsets.WR)

# Obtain ybar_s WR
ybar.WR <- NULL
for(i in 1:all.nr.samples.WR){
  s <- unlist(all.subsets.WR[i,])
  yi.s <- yi[s]
  temp <- mean(yi.s)
  ybar.WR <- c(ybar.WR,temp) 
}

# Plot histogram
par(mfrow=c(1,2))
hist(ybar.WOR,seq(min(ybar.WOR)-1,max(ybar.WOR)+1,1),freq=FALSE)
hist(ybar.WR,seq(min(ybar.WR)-1,max(ybar.WR)+1,1),freq=FALSE)

# Calculate expected value
mean(ybar.WOR)
mean(ybar.WR)

# True value
mean(yi)

# Calculate sampling variance
var(ybar.WOR)*(all.nr.samples.WOR-1)/(all.nr.samples.WOR)
var(ybar.WR)*(all.nr.samples.WR-1)/(all.nr.samples.WR)

# Calculate sampling variance analytically
S2 <- var(yi)
var.srswor <- (1-n/N)*S2/n
var.srswor

var.srswr <- S2*(N-1)/N/n
var.srswr

# ----------
# Exercise 4
# ----------
# Frequency table
pubi <- seq(0, 10, by=1)
freqi <- c(28,4,3,4,4,2,1,0,2,1,1)

# Publications in sample
yi <- rep(pubi,freqi)

# Popuation size
N <- 807
# Sample size
n <- 50

# 1).
# Plot histogram
h <- hist(yi,seq(min(yi)-1,max(yi)+1,1),freq=FALSE)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE, main = "Histogram of referred publications (%)",xlab="Referred publications",ylab="Percentage (%)")


# 2).
# Estimate the mean number of publications per faculy member 
mean(yi)


# s2
s2 <- var(yi)
s2

# SE
sqrt((1-n/N)*s2/n)


# 3).
# Skewness
library(e1071)
skewness(yi)

# Calculate an estimate of skewness manually
mu3 <- sum((yi-mean(yi))^3)/n
mu3
mu2 <- s2
mu2

gamma <- mu3/mu2^(3/2)
gamma

# nmin
ceiling(28+25*gamma^2)

# 4).
# Proportion of faculty members with no referred publications
phat <- sum(yi==0)/n
phat

se <- sqrt((1-n/N)*phat*(1-phat)/(n-1))
se

# 95% CI
phat+c(-1.96,1.96)*se

