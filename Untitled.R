install.packages("evir")
install.packages("evd")
install.packages("gsl")
install.packages("evmix")
library(evir)
library(evd)
library(gsl)
library(evmix)

appl <- read.csv("AAPL.csv")
appl_adjts <- ts(data=appl["Adj.Close"])

# calculate negative daily log-returns
n <- length(appl_adjts)  
Xt <- -log(appl_adjts[-1]/appl_adjts[-n])  

# export negative daily log-returns as .csv file
Xt_df <- data.frame(Xt)
write.csv(Xt_df, "C:\\Users\\unowen\\Desktop\\APPLndl.csv")

####### part (b) Block maximum based method ######

### part (b1) ###
Mblocks <- c()
max_n <- 22

# note that the maximum i-th block is 11 since any block
# greater than 11 will have Xt values that is out of bounds
block_max = 0
for (i in 1:11) {
  # find max Xt in block i
  for (j in 1:max_n) {
    if (Xt[max_n*(i-1)+j] > block_max) {
      block_max = Xt[max_n*(i-1)+j]
    }
  }
  Mblocks <- c(Mblocks, block_max)
  block_max = 0
}

# create histogram of the n-block sequence
hist(Mblocks, 
     main = "Histogram for 22-block sequence",
     xlab = "Maximum negative log-return from block",
     cex.main = 0.8, col="firebrick")

### part (b2) ###
# fit sequence into GEV distribution
# so GEV cumulative function has parameters loc, scale, shape
Mblock_GEV <- fgev(Mblocks)
Mblock_par <- Mblock_GEV$estimate

### part (b3) ###
# QQ-plot analysis for goodness of fitting
# we see that the points fit the regression line
# well so it looks that the GED is a good fit for the
# data

# get quantiles of the inverse CDF of GEV distribution
qqx <- qgev(ppoints(length(Mblocks)), 
            loc = Mblock_par[1], 
            scale = Mblock_par[2], 
            shape = Mblock_par[3])

# get sample quantiles of Mblocks probabilities
qqy <- quantile(Mblocks, probs = ppoints(length(Mblocks)))

# qqplot with line y=x
plot(qqx, qqy, 
     main = "QQ-plot",
     xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles",
     cex.main = 0.8, pch = 16)
abline(0, 1, col = "red")


####### part (c) Threshold exceedance based method ######

### part (c1) ###
# Create sample mean excess plot and get appropriate threshold u
# we will take u = 0.01 since it looks that the points start to 
# oscillate for threshold > -0.01
meplot(Xt, main = "Sample Mean Excess Plot",
       cex.main = 0.8, pch = 16)
abline(v=-0.01, col="red")
GPD_mu <- -0.01

### part (c2) ###
# fit sequence into GPD and we get the
# shape and scale parameters to define the GPD
Mblock_GPD <- fgpd(Mblocks, u = -0.01)
GPD_shape <- Mblock_GPD$xi
GPD_scale <- Mblock_GPD$sigmau

### part (c3) ###
# survival function of F(u) = P(X > u)
F_u <- 1 - (length(Xt[Xt <= GPD_mu])/251)

# Estimate for VaR_0.99, since alpha = 0.99 > F(u) we get
GPD_var99 <- GPD_mu + (GPD_scale / GPD_shape)*
                       ( ((1 - 0.99)/F_u)^(-GPD_shape) - 1 )

# Estimate for CVaR_0.99, since shape < 1 and alpha = 0.99 > F(u) we get
GPD_cvar99 <- (GPD_var99/(1-GPD_shape)) + ((GPD_scale - GPD_shape*GPD_mu)/(1-GPD_shape))


####### part (d) Hill's method ######

### part (d1) ###
# Produce a Hill plot
hillplot(Xt)


### part (d2) ###
# Develop a Hill's estimator by choosing a k
# it looks that when k = 25, the Hill plot is flat
# around this point
k <- 25
abline(v=k, col="red")

Xt_des <- sort(Xt, decreasing = TRUE)
Xt_hill <- Xt_des[1:k]
Xt_hill_est <- ((1/k)*(sum(log(Xt_hill))) - log(Xt_hill[k]))^(-1)

### part (d3) ###
# Using the hill estimator, calculate VaR_0.99 and CVaR_0.99

# we use the survival function from the estimated distribution to
# calculate VaR_0.99 
hill_var99 <- Xt_hill[k] * ( (251*0.01/k)^(-1 / Xt_hill_est) ) 

# using VaR_0.99, we calculate CVaR_0.99 as:
hill_cvar99 <- hill_var99 + (( ((Xt_hill[k])^(Xt_hill_est)) * k)/(0.01*251)) *
                            (-(hill_var99)^(1-Xt_hill_est))/(1-Xt_hill_est)



