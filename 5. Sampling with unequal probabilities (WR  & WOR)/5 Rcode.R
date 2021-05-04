# -----------------------------------------------
# Solution ark # 5.Sampling with unequal probabilities
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------
# ----------
# Exercise 2
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'azcounties.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

N <- nrow(sample.data)
ti <- sample.data$housing
Mi <- sample.data$population

n <- 1

cor(ti,Mi)

# 1).
# Inclusion probabilities
pi <- Mi/sum(Mi)
summary(pi)

# that for each possible sample
that <- ti/pi
that


# Variance of that
t <- sum(ti)
t
sum(pi*(that-t)^2)


# 2).
# Inclusion probabilities under SRS
pi.srs <- n/N

# that for each possible sample
that.srs <- ti/pi.srs
that.srs


# Variance of that
sum(pi.srs*(that.srs-t)^2)

# -----------
# Exercise 3
# -----------
N <- 27
id.nr <- c(14,23,9,14,16,6,14,19,21,11)
n <- length(id.nr)

Mi <- c(65,25,48,65,2,62,65,62,61,41)

yij <- list(c(3, 0, 0, 4),c(2, 1, 2, 0),c(0, 0, 1, 0),c(2, 0, 1, 0),c(2, 0),
            c(0, 2, 2, 5),c(1, 0, 0, 3),c(4, 1, 0, 0),c(2, 2, 3, 1),c(2, 5, 12, 3))

pi <- c(0.0805452, 0.0309789, 0.0594796, 0.0805452, 0.0024783, 0.0768278,
        0.0805452, 0.0768278, 0.0755886, 0.0508055) 

mi <- sapply(yij,length)

# Estimated totals for each academic unit
tihat <- sapply(yij,sum)*Mi/mi
tihat

# Estimated total
that <- sum(tihat/pi)/n
that

# Standard error of the estimated total
tihat.breve <- tihat/pi
sd(tihat.breve)
sd(tihat.breve)/sqrt(n)

# ----------
# Exercise 4
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), "/")
sample.data <- read.csv(paste(libname,'classpps.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

# Load data
incprob.data <- read.csv(paste(libname,'classppsjp.csv',sep=''),sep=',')
head(incprob.data)
dim(incprob.data)

# Merge two data
sample.data <- merge(sample.data,incprob.data[,names(incprob.data)!='clssize'],by='class',all=T)
sample.data

N <- 15
n <- 5
mi <- 4
M0 <- 647

# Find that.HT
tihat.HT.data <- with(sample.data,aggregate(hours*clssize,by=list(class),FUN=mean))
names(tihat.HT.data) <- c('class','tihat.HT')
tihat.HT <- tihat.HT.data$tihat.HT
tihat.HT


incprob.data <- merge(incprob.data,tihat.HT.data,by='class',all=T)
incprob.data

# 1).
# Estimate the HT variance of that.HT
# Variance estimates of tihat
si2.data <- with(sample.data,aggregate(hours,by=list(class),FUN=var))
names(si2.data) <- c('class','si2')            

incprob.data <- merge(incprob.data,si2.data,by='class',all=T)
incprob.data

#  
incprob.data$varhattihat <- with(incprob.data,clssize^2*(1-mi/clssize)*si2/mi)
incprob.data$varhattihat


# Variance estimate of sum(ti/pi)
# 
incprob.data$psuid <- c(3,1,4,2,5)
incprob.data <- incprob.data[order(incprob.data$psuid),]

temp.HT <- 0
for(i in 1:n){
  temp.pi <- incprob.data$SelectionProb[i]
  temp.ti <- incprob.data$tihat.HT[i]
  for(j in 1:n){
    temp.pj <- incprob.data$SelectionProb[j]
    temp.tj <- incprob.data$tihat.HT[j]
    if(i==j){temp.pij <- temp.pi
    tempi.HT <- (1-temp.pi)*temp.ti/temp.pi*temp.ti/temp.pi}
    if(i!=j){temp.pij <- incprob.data[i,j+4]
    tempi.HT <- (1-temp.pi*temp.pj/temp.pij)*temp.ti/temp.pi*temp.tj/temp.pj}
    temp.HT <- temp.HT + tempi.HT
  }
}
temp.HT

# The HT estimate of the variance of that.HT
varhatthat.HT <- temp.HT + with(incprob.data,sum(varhattihat*SamplingWeight))
varhatthat.HT

# The SYG estimate of the variance of that.HT
temp.SYG <- 0
for(i in 1:n){
  temp.pi <- incprob.data$SelectionProb[i]
  temp.ti <- incprob.data$tihat.HT[i]
  for(j in 1:n){
    temp.pj <- incprob.data$SelectionProb[j]
    temp.tj <- incprob.data$tihat.HT[j]
    if(i==j){tempi.SYG <- 0}
    if(i!=j){temp.pij <- incprob.data[i,j+4]
    tempi.SYG <- -(1-temp.pi*temp.pj/temp.pij)/2*(temp.ti/temp.pi-temp.tj/temp.pj)^2}
    temp.SYG <- temp.SYG + tempi.SYG
  }
}
temp.SYG

varhatthat.SYG <- temp.SYG + with(incprob.data,sum(varhattihat*SamplingWeight))
varhatthat.SYG

# 2).
# With-replacement variance estimator of that.HT
that.HT  <- with(incprob.data,sum(tihat.HT*SamplingWeight))
that.HT

varhatthat.WR <- with(incprob.data,sum((tihat.HT/SelectionProb-that.HT/n)^2)*n/(n-1))
varhatthat.WR

# Finite population correction adjusted WR variance estimate
varhatthat.WR*(1-n/N)


# 3).
# Standard errors
sqrt(varhatthat.HT)
sqrt(varhatthat.SYG)
sqrt(varhatthat.WR)
sqrt(varhatthat.WR*(1-n/N))


