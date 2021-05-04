# -----------------------------------------------
# Solution ark #6. Survey nonresponse
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------
# ----------
# Exercise 2
# ----------
y <- c(600, 520, 620, 500, 380, 460, 450, 250, 400, 780)
s <- seq(1,10)

#(1)
set.seed(2020)
simp <- sample(s, 10, replace = TRUE)
simp

yimp <- y[simp]
yimp

ycomp <- c(y, yimp)
ycomp

mean(ycomp)
mean(y)

var(ycomp)
v <- var(ycomp)/length(ycomp)
se <- sqrt(v)
se

CIimp <- mean(ycomp) + c(-1.96, 1.96) * se
CIimp

#(2)
mean(y)

vresp <- var(y)/length(y)
seresp <- sqrt(vresp)
seresp

CIresp = mean(y) + c(-1.96, 1.96) * seresp
CIresp

# ----------
# Exercise 4
# ----------
# Load package "sampling"
library(sampling)

# Load data
# You should change libname depending on where the data is located on your PC
libname <- paste0(getwd(),"/")
sample.data <- read.csv(paste(libname,'certify.csv',sep=''),sep=',')
head(sample.data)
dim(sample.data)

N <- 18609
n <- N

nR <- nrow(sample.data)
nR

ph.edu <- c(0.55,0.38,0.07)
sum(ph.edu)


ph.workpl <- c(0.29,0.34,0.11,0.26)
sum(ph.workpl)



#(1).
nh.edu <- round(n*ph.edu,0)
nh.edu

nh.workpl <- round(n*ph.workpl)
nh.workpl

# Number of respondent by the groups
with(sample.data,table(college))
nhR.edu <- c(sum(sample.data$college=="P"),sum(sample.data$college=="M"),
             nR-sum(sample.data$college %in% c("P","M")))
nhR.edu


with(sample.data,table(workenv))
nhR.workpl <- c(sum(sample.data$workenv=="I"),sum(sample.data$workenv=="A"),
                sum(sample.data$workenv=="G"),
                nR-sum(sample.data$workenv %in% c("I","A","G")))
nhR.workpl

# Response rates by the groups
rrh.edu <- round(nhR.edu/nh.edu*100,0)
rrh.edu

rrh.workpl <- round(nhR.workpl/nh.workpl*100,0)
rrh.workpl

# (2).
sample.data$phd <- 1
sample.data$phd[sample.data$college!='P'] <- 2
sample.data$phd <- as.factor(sample.data$phd)

sample.data$work <- 3
sample.data$work[sample.data$workenv=='I'] <- 1
sample.data$work[sample.data$workenv=='A'] <- 2
sample.data$work <- as.factor(sample.data$work)

T0 <- with(sample.data,table(phd,work))

w0 <- N/nR
w0

T1 <- T0*w0
round(T1,1)

round(apply(T1,1,sum),1)
round(apply(T1,2,sum),1)

T2 <- T1
T2[1,] <- T1[1,]*round(n*ph.edu[1],0)/apply(T1,1,sum)[1]
T2[2,] <- T1[2,]*(n-round(n*ph.edu[1],0))/apply(T1,1,sum)[2]
T2

round(T2,1)
round(apply(T2,1,sum),1)
round(apply(T2,2,sum),1)

T3 <- T2
T3[,1] <- T2[,1]*round(n*ph.workpl[1],0)/apply(T2,2,sum)[1]
T3[,2] <- T2[,2]*round(n*ph.workpl[2],0)/apply(T2,2,sum)[2]
T3[,3] <- T2[,3]*(n-round(n*ph.workpl[1],0)-round(n*ph.workpl[2],0))/apply(T2,2,sum)[3]
T3

round(T3,1)
round(apply(T3,1,sum),1)
round(apply(T3,2,sum),1)

# Iterative raking with "sampling" package
# create a matirx of dummy variables based on phd and work
a <- model.matrix(~-1+sample.data$phd+sample.data$work)
head(a)

# Compute g-weights, which are the coefficients that shall be used to adjust the initial weights
# Final weights = initial weights * g-weights
cal.g <- calib(a,d=rep(w0,nR),total=c(10235,8374,6327,6885),method="raking")

# Raking weights
wtilde <- cal.g*w0

# Raking weights by phd X work
with(sample.data,aggregate(wtilde,by=list(phd,work),FUN=mean))

# Estimated cell totals after raking adjustment: phd X work
with(sample.data,aggregate(wtilde,by=list(phd,work),FUN=sum))


# Estimated marginal totals after raking adjustment: phd
with(sample.data,aggregate(wtilde,by=list(phd),FUN=sum))


# Estimated marginal totals after raking adjustment: work
with(sample.data,aggregate(wtilde,by=list(work),FUN=sum))



# (3).
# Estimation of the percantages for the question regarding the issue of certification
# Without weights
round(table(sample.data$certify)/nR*100,1)

# With raking weights
round(with(sample.data,aggregate(wtilde,by=list(certify),FUN=sum))$x/sum(wtilde)*100,1)



