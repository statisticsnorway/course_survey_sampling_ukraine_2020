# -----------------------------------------------
# Solution ark # 4. Cluster sampling (one-stage, two-stage). Autumn 2020
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------


# ----------
# Exercise 3
# ----------
N <- 45
n <- 6

Mi <- c(52,19,37,39,8,14)
yij <- list(c(146,180,251,152,72,181,171,361,73,186),c(99,101,52,121),c(199,179,98,63,126,87,62),
            c(226,129,57,46,86,43,85,165),c(12,23),c(87,43,59))

yij


# Summary statistics for each cluster
# Number sample supermarkets in each cluster
mi <- sapply(yij,length)
mi

# Mean of cases sold in each cluster
yibar <- sapply(yij,mean)
round(yibar,2)

# SD of cases sold in each cluster
si <- sapply(yij,sd)
round(si,4) 

# Plot of the data
Iij <- as.factor(rep(1:6,mi))
plot(Iij,unlist(yij),xlab='City',ylab='Number of cases sold')

# Estimate of the total number of cases sold
tihat <- Mi*yibar
round(tihat,2)

that <- sum(tihat*N/n)
that

# Estimate of the average number of cases sold per supermarket
ybarr <- sum(tihat)/sum(Mi)
ybarr

# Variance estimate of that
st2 <- var(tihat)
st2

varyibar <- Mi^2*(1-mi/Mi)*si^2/mi
varyibar
sum(varyibar)

varhat.that <- N^2*(1-n/N)*st2/n + N/n*sum(varyibar)
varhat.that
sqrt(varhat.that)

# Variance estimate of ybarr
mbar <- sum(Mi)/n
mbar

s2.yibar <- sum(Mi*Mi*(yibar-ybarr)^2)/(n-1)
s2.yibar

varhat.ybarr <- (1-n/N)*s2.yibar/n/mbar/mbar + 1/n/N/mbar/mbar*sum(varyibar)
varhat.ybarr
sqrt(varhat.ybarr)

# ----------
# Exercise 5
# ----------
# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(), '/')
sample.data <- read.csv(paste(libname,'measles.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

n <- length(unique(sample.data$school))
N <- 46

# 1).
kij <- as.numeric(sample.data$form==1)
rij <- as.numeric(sample.data$returnf==1)

sum.data1 <- with(sample.data,aggregate(cbind(Mitotal,mi),by=list(school),FUN=mean))
names(sum.data1) <- c('school','Mi','mi')

sum.data2 <- with(sample.data,aggregate(cbind(kij,rij),by=list(school),FUN=sum))
names(sum.data2) <- c('school','ki','ri')

sum.data2$yibar <- with(sum.data2,ri/ki)

sum.data <- merge(sum.data1,sum.data2,by='school',all=T)
sum.data

# 2).
with(sum.data,Mi/ki)

# 3).
tihat <- with(sum.data,yibar*Mi)
sum(tihat)

sum.Mi <- sum(sum.data$Mi)
sum.Mi

# Estimate of the percentage of parents returned a consent form
ybarr <- sum(tihat)/sum.Mi
ybarr

# Calculation of standard error of ybarr
mbar <- sum.Mi/n
mbar

# Residuals from ratio estimation
ei <- with(sum.data,Mi*(yibar-ybarr))
# SD of the residuals
se2 <- var(ei)
se2


# SD of forms returned
si2 <- with(sum.data,ki*yibar*(1-yibar)/(ki-1))
si2

# 
varyibar <- with(sum.data,Mi^2*(1-ki/Mi)*si2/ki)
varyibar
sum(varyibar)


# Variance estimate of ybarr
varhat.ybarr <- (1-n/N)*se2/n/mbar/mbar + 1/n/N/mbar/mbar*sum(varyibar)
varhat.ybarr

se.ybarr <- sqrt(varhat.ybarr)
se.ybarr

# 95% CI
ybarr+qnorm(c(0.025,0.975),0,1)*se.ybarr


# 4).
sum.ki <- sum(sum.data$ki)
sum.ri <- sum(sum.data$ri)

phatsrs <- sum.ri/sum.ki
phatsrs

# Variance under SRS
varhat.phatsrs <- phatsrs*(1-phatsrs)/(sum.ki-1)
varhat.phatsrs

# Design effect
deffhat <- varhat.ybarr/varhat.phatsrs
deffhat

