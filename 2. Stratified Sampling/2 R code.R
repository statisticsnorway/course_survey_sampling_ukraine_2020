# -----------------------------------------------
# R-code for Solution ark 2, Stratified Sampling. Autumn 2020
# -----------------------------------------------
# Melike Oguz-Alper, Tatsiana Pekarskaya
# -----------------------------------------------

# ----------
# Exercise 3
# ----------
# 1).
Nh <- c(102,310,217,178)

freq_pub <- c(0:8)


y1i <- rep(freq_pub,c(1,2,0,1,0,2,0,1,0))
y2i <- rep(freq_pub,c(10,2,0,1,2,1,1,0,2))
y3i <- rep(freq_pub,c(9,0,1,0,2,0,1,0,0))
y4i <- rep(freq_pub,c(8,2,0,1,0,0,0,0,0))

nh <- c(length(y1i),length(y2i),length(y3i),length(y4i))
nh
sum(nh)

sh2 <- c(var(y1i),var(y2i),var(y3i),var(y4i))
sh2

yhbar <- c(mean(y1i),mean(y2i),mean(y3i),mean(y4i))
yhbar

that <- sum(Nh*yhbar)
that

varstr <- sum(Nh^2*(1-nh/Nh)*sh2/nh)
varstr

sestr <- sqrt(varstr)
sestr

# 3).
hatp_str <- 1/sum(Nh) * sum(Nh * c(1,10,9,8)/ nh)
hatp_str

hatp_h <- c(1,10,9,8) / nh
hatV_str <- sum(1/sum(Nh)^2 * Nh^2 * (1 - nh/Nh) * hatp_h * (1 - hatp_h) / (nh - 1))
hatV_str

sqrt(hatV_str)