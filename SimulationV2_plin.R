###### Note: Must run detrending to get the y_hats to adjust to y
rm(list = ls())
sim_ctys <- 100 # Number of county level yields per year to simulate
sim_year <- 100 # Number of years worth of data to simulate
Ms        <- 1:100 # Index of MC simulations
#J        <- 100 # Number of approximate quantiles to obtain for each Monte Carlo simulation
K        <- 1000 # Number of Simulations for the 
Q        <- 4 # 1 + number of quantiles used in determining the B-spline knots (i.e. using 1/Q to (Q-1)/Q quantiles as the knots)
degree   <- 3 # Degree used in Bspline bases

source("BSplineFunctions.R") # Get detrended yield and price and functions to run B-spline quantile regression
library(tidyverse)
library(splines)
library(readxl)
library(MASS)
library(sn)
library(Peacock.test)
library(fasano.franceschini.test)
library(xtable)

set.seed(10062018) # For reproducability
taus <- c(0.1,0.25,0.5,0.75,0.9)# Fixed values of tau used for the quantiles
tau_p <- runif(K,0.01,0.99) # Fixed values of tau_p used in the singular simulation
tau_y <- runif(K,0.01,0.99) # Fixed values of tau_y used in the singular simulation

# Simulate detrended stocks, s: Based on method of moments estimation of s
s <- rbeta(sim_year,7,44)
s_kn <- quantile(s,seq(1/Q,(Q-1)/Q,1/Q)) # knots for Bspline bases
Bs <- bs(s,degree=degree,knots=s_kn,Boundary.knots=c(0,1),intercept = TRUE) # Bspline basis used in quantile regression for p
s_rep <- rep(s,sim_ctys) # replicate s for quantile regression for y
Bs_rep <- bs(s_rep,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE) # Bspline basis used in quantile regression for y
s_tru_quants <- qbeta(c(0.1,seq(1/Q,(Q-1)/Q,1/Q),0.9),7,44) # true quantiles used for 1 simulation of joint density (1st Monte Carlo sim)
Bs_tru_quants <- bs(rep(s_tru_quants,each=101),degree=degree,knots=s_kn,Boundary.knots=c(0,1),intercept = TRUE) # Bspline basis used in joint density simulation
Bs_est <- bs(seq(0,1,0.01),degree=degree,knots=s_kn,Boundary.knots=c(0,1),intercept = TRUE) # Bspline basis used in producing plots for simulation evaluation
lambdas <- 0.01 # Set smoothing parameter lambda (previously calculated using cross validation techniques)
n_lambdas <- length(lambdas) # Only used in 

# Obtain Difference matrix in the penalty term
Ds <- matrix(0,nrow=ncol(Bs)-2,ncol=ncol(Bs))

for (j in 1:nrow(Ds)){
  for (k in 1:ncol(Ds)){
    if (k-j == 0 || k-j == 2){
      Ds[j,k] = 1
    } else if(k-j == 1){
      Ds[j,k] = -2
    }
  }
}

Dps <- matrix(0,nrow=2*(ncol(Bs))-2,ncol=2*(ncol(Bs)))

for (j in 1:nrow(Dps)){
  for (k in 1:ncol(Dps)){
    if (k-j == 0 || k-j == 2){
      Dps[j,k] = 1
    } else if(k-j == 1){
      Dps[j,k] = -2
    }
  }
}

# Matrices used to save the simulation data
# p_app_lin <- array(0,dim=c(nrow(Bs_est),length(taus),M))
# p_app_non <- array(0,dim=c(nrow(Bs_est),length(taus),M))
# y_app_lin_p_lin <- array(0,dim=c(nrow(Bs_est)*(Q+1),length(taus),M))
# y_app_lin_p_non <- array(0,dim=c(nrow(Bs_est)*(Q+1),length(taus),M))
# y_app_non_p_lin <- array(0,dim=c(nrow(Bs_est)*(Q+1),length(taus),M))
# y_app_non_p_non <- array(0,dim=c(nrow(Bs_est)*(Q+1),length(taus),M))


# Set functions for conditional mean and standard deviation of p and s
# Based on linear regression and residual histograms of p
# a_p and a_y are to ensure mean and standard deviation of residuals are 0 and 1, respectively
m_p_lin <- function(s){0.2 - 0.4*s}
s_p_lin <- function(s){0.5 - 0.5*s}
# m_p_non <- function(s){-0.2+0.4*exp(-2*s)}
# s_p_non <- function(s){0.5*exp(-2*s)}
# 
a_p <- 3 # Skewness parameter because detrended prices are empirically right skewed
d_p <- a_p / sqrt(1 + a_p^2)  #To ensure eps_p terms have theoretical zero mean and unit variance

m_y_non <- function(p,s){-25 + 14.45*exp(p) + 22.18*s}
#m_y_non <- function(p,s){-5.31 + 4.67*p - 9.01*s + 136.04*p^2 + 52.56*s^2}
s_y <- 33
a_y <- -3 # Skewness parameter because detrended prices are empirically left skewed
d_y <- a_y / sqrt(1 + a_y^2) #Used to ensure eps_y terms have theoretical zero mean and unit variance

p_eps <- qsn(tau_p,
             xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
             omega=sqrt((1-2*d_p^2/pi)^-1),
             alpha=a_p)
p_tru <- cbind(p_eps*s_p_lin(s_tru_quants[1]) + m_p_lin(s_tru_quants[1]),
               p_eps*s_p_lin(s_tru_quants[2]) + m_p_lin(s_tru_quants[2]),
               p_eps*s_p_lin(s_tru_quants[3]) + m_p_lin(s_tru_quants[3]),
               p_eps*s_p_lin(s_tru_quants[4]) + m_p_lin(s_tru_quants[4]),
               p_eps*s_p_lin(s_tru_quants[5]) + m_p_lin(s_tru_quants[5]))

y_eps <- qsn(tau_y,
             xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
             omega=sqrt((1-2*d_y^2/pi)^-1),
             alpha=a_y)
y_tru <- cbind(y_eps*s_y + m_y_non(p_tru[,1],s_tru_quants[1]),
               y_eps*s_y + m_y_non(p_tru[,2],s_tru_quants[2]),
               y_eps*s_y + m_y_non(p_tru[,3],s_tru_quants[3]),
               y_eps*s_y + m_y_non(p_tru[,4],s_tru_quants[4]),
               y_eps*s_y + m_y_non(p_tru[,5],s_tru_quants[5]))

p_vals_s1 <- NULL
p_vals_s2 <- NULL
p_vals_s3 <- NULL
p_vals_s4 <- NULL
p_vals_s5 <- NULL

for (m in Ms){
  set.seed(160793+m) # For reproducibility
  eps_p <- rsn(sim_year,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)
  eps_y <- rsn(sim_year*sim_ctys,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)
  
  p_tilde <- m_p_lin(s) + s_p_lin(s)*eps_p
  y_tilde <- m_y_non(p_tilde,s) + s_y*eps_y
  
  lambda <- lambdas[1]
  
  # For the 100th Monte Carlo simulation, we simulate an entire joint density function of y and p given a particular s
  # Our main goal is to approximate that conditional join density function
  # if (m == 1){
    p_sim <- NULL
    y_sim <- NULL
    # start_time <- Sys.time()
    for (k in 1:K){
      # Simulate, for a given s,"K" values from the joint density function of y and p using quantile regression
      if(k %% 5 == 0){
        print(paste0("Running Simulation: ",m,".",k))
      }
      lambda <- optim.lambda(unique(p_tilde),Bs,D(ncol(Bs)),tau_p[k])$estimate
      p_beta_cond[k,] <- quant.reg(unique(p_tilde),Bs,D(ncol(Bs)),lambda,tau_p[k])$b
      p_tilde_sim <- quant.reg(p_tilde,Bs,Ds,lambda,tau_p[k])$b
      # Simulate y given p and s
      p_kn <- quantile(p_tilde_rep,seq(1/Q,(Q-1)/Q,1/Q))
      Bp_rep <- bs(p_tilde_rep,degree=degree,knots=p_kn,Boundary.knots = c(-2,2),intercept = TRUE)
      Bs_rep <- bs(s_rep,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
      Bps_rep <- cbind(Bp_rep,Bs_rep)
      y_tilde_sim <- quant.reg(y_tilde,Bps_rep,Dps,lambda=1e-8,tau_y[k])$b
      p_sim <- cbind(p_sim,p_tilde_sim)
      y_sim <- cbind(y_sim,y_tilde_sim)
    }
    # end_time <- Sys.time()
  # }
    # p_vals_s1 <- c(p_vals_s1,fasano.franceschini.test(cbind(p_sim[1,],y_sim[1,]),cbind(p_tru[,1],y_tru[,1]),nPermute=1000)$p.value)
    # p_vals_s2 <- c(p_vals_s2,fasano.franceschini.test(cbind(p_sim[2,],y_sim[2,]),cbind(p_tru[,2],y_tru[,2]),nPermute=1000)$p.value)
    # p_vals_s3 <- c(p_vals_s3,fasano.franceschini.test(cbind(p_sim[3,],y_sim[3,]),cbind(p_tru[,3],y_tru[,3]),nPermute=1000)$p.value)
    # p_vals_s4 <- c(p_vals_s4,fasano.franceschini.test(cbind(p_sim[4,],y_sim[4,]),cbind(p_tru[,4],y_tru[,4]),nPermute=1000)$p.value)
    # p_vals_s5 <- c(p_vals_s5,fasano.franceschini.test(cbind(p_sim[5,],y_sim[5,]),cbind(p_tru[,5],y_tru[,5]),nPermute=1000)$p.value)
}
# saveRDS(p_vals_s1,"p_vals_s1_plin.rds")
# saveRDS(p_vals_s2,"p_vals_s2_plin.rds")
# saveRDS(p_vals_s3,"p_vals_s3_plin.rds")
# saveRDS(p_vals_s4,"p_vals_s4_plin.rds")
# saveRDS(p_vals_s5,"p_vals_s5_plin.rds")
saveRDS(p_sim,"p_beta_plin.rds")
saveRDS(y_sim,"y_beta_plin.rds")
# saveRDS(p_tru,"p_tru_plin.rds")
# saveRDS(y_tru,"y_tru_plin.rds")

# p_tru <- rsn(K,
#              xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
#              omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)*s_p_non(s_tru_quants[2]) + m_p_non(s_tru_quants[2])
# y_tru <- rsn(K,
#              xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
#              omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)*s_y + m_y_non(p_tru,s_tru_quants[2])
# peacock2(cbind(p_sim[1,],y_sim[1,]),cbind(p_tru,y_tru))
# fasano.franceschini.test(cbind(p_sim[1,],y_sim[1,]),cbind(p_tru,y_tru))
# 
# p_tru <- rsn(10000,
#              xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
#              omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)*s_p_non(s_tru_quants[3]) + m_p_non(s_tru_quants[3])
# y_tru <- rsn(10000,
#              xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
#              omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)*s_y + m_y_non(p_tru,s_tru_quants[3])
# peacock2(cbind(p_sim[2,],y_sim[2,]),cbind(p_tru,y_tru))
# fasano.franceschini.test(cbind(p_sim[2,],y_sim[2,]),cbind(p_tru,y_tru))
# 
# p_tru <- rsn(K,
#              xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
#              omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)*s_p_non(s_tru_quants[4]) + m_p_non(s_tru_quants[4])
# y_tru <- rsn(K,
#              xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
#              omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)*s_y + m_y_non(p_tru,s_tru_quants[4])
# peacock2(cbind(p_sim[3,],y_sim[3,]),cbind(p_tru,y_tru))
# # fasano.franceschini.test(cbind(p_sim[3,],y_sim[3,]),cbind(p_tru,y_tru))$p.value
# xtable(rbind(summary(p_vals_s1),summary(p_vals_s2),summary(p_vals_s3)),digits=4)


