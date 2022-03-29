###### Note: Must run detrending to get the y_hats to adjust to y
rm(list = ls())
sim_ctys <- 500 # Number of county level yields per year to simulate
sim_year <- 100 # Number of years worth of data to simulate
M        <- 100 # Number of Monte Carlo Simulations
J        <- 100 # Number of approximate quantiles to obtain for each Monte Carlo simulation
K        <- 5000 # Number of Simulations for the 
Q        <- 4 # 1 + number of quantiles used in determining the B-spline knots (i.e. using 1/Q to (Q-1)/Q quantiles as the knots)
degree   <- 3 # Degree used in Bspline bases

source("Corn.R") # Get detrended yield and price and functions to run B-spline quantile regression
library(tidyverse)
library(splines)
library(readxl)
library(MASS)
library(sn)

set.seed(10062018) # For reproducability
taus <- seq(1/J,(J-1)/J,1/J)# Fixed values of tau used for the quantiles
tau_p <- runif(K) # Fixed values of tau_p used in the singular simulation
tau_y <- runif(K) # Fixed values of tau_y used in the singular simulation

# Simulate detrended stocks, s: Based on method of moments estimation of s
s <- rbeta(sim_year,1.5,6)
s_kn <- quantile(s,seq(1/Q,(Q-1)/Q,1/Q)) # knots for Bspline bases
Bs <- bs(s,degree=degree,knots=s_kn,Boundary.knots=c(0,1),intercept = TRUE) # Bspline basis used in quantile regression for p
s_rep <- rep(s,sim_ctys) # replicate s for quantile regression for y
Bs_rep <- bs(s_rep,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE) # Bspline basis used in quantile regression for y
s_tru_quants <- qbeta(c(0.1,seq(1/Q,(Q-1)/Q,1/Q),0.9),1.5,6) # true quantiles used for 1 simulation of joint density (1st Monte Carlo sim)
Bs_tru_quants <- bs(rep(s_tru_quants,101),degree=degree,knots=s_kn,Boundary.knots=c(0,1),intercept = TRUE) # Bspline basis used in joint density simulation
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
p_app_lin <- array(0,dim=c(nrow(Bs_est),J-1,M))
p_app_non <- array(0,dim=c(nrow(Bs_est),J-1,M))
y_app_lin_p_lin <- array(0,dim=c(nrow(Bs_est)*(Q+1),J-1,M))
y_app_lin_p_non <- array(0,dim=c(nrow(Bs_est)*(Q+1),J-1,M))
y_app_non_p_lin <- array(0,dim=c(nrow(Bs_est)*(Q+1),J-1,M))
y_app_non_p_non <- array(0,dim=c(nrow(Bs_est)*(Q+1),J-1,M))

for (m in 1:M){
  set.seed(160793+m) # For reproducibility
  
  # Set functions for conditional mean and standard deviation of p and s
  # Based on linear regression and residual histograms of p
  # a_p and a_y are to ensure mean and standard deviation of residuals are 0 and 1, respectively
  m_p_lin <- function(s){-0.5 + 0.5*s}
  s_p_lin <- function(s){0.3 + 0.1*s}
  m_p_non <- function(s){-0.5+0.5*exp(-s)}
  s_p_non <- function(s){0.15+0.25*exp(-3*s)}
  # 
  a_p <- 3 # Skewness parameter because detrended prices are empirically right skewed
  d_p <- a_p / sqrt(1 + a_p^2)  #To ensure eps_p terms have theoretical zero mean and unit variance
  eps_p <- rsn(sim_year,-d_p*sqrt(2/pi),(1-2*d_p^2/pi)^-1,a_p)
  
  m_y_lin <- function(p,s){-1.11 + 14.45*p + 22.18*s}
  m_y_non <- function(p,s){-5.31 + 4.67*p - 9.01*s + 136.04*p^2 + 52.56*s^2}
  s_y <- 31
  a_y <- -3 # Skewness parameter because detrended prices are empirically left skewed
  d_y <- a_y / sqrt(1 + a_y^2) #To ensure eps_y terms have theoretical zero mean and unit variance
  eps_y <- rsn(sim_year*sim_ctys,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)
  
  # Run simulations for all four possible combinations of linear and non-linear p and y
  for (l in 0:3){
    # if l (integer divide) 2 is 0, then p is linear, else p is non-linear
    if (l %/% 2 == 0){
      p_tilde <- m_p_lin(s) + s_p_lin(s)*eps_p
    } else {
      p_tilde <- m_p_non(s) + s_p_non(s)*eps_p
    }
    # if l mod 2 is 0, then y is linear, else y is non-linear
    if (l %% 2 == 0){
      y_tilde <- m_y_lin(p_tilde,s) + s_y*eps_y
    } else {
      y_tilde <- m_y_non(p_tilde,s) + s_y*eps_y
    }
    for (v in 1:n_lambdas){
      # Set smoothing parameter lambda
      lambda <- lambdas[v]
      p_sim <- NULL
      y_sim <- NULL
      for (tau in taus){
        print(paste0("Getting ",tau," Quantile from Simulation ",m,".",l," with lambda = ",lambdas[v]))
        # Obtain tau^th conditional quantiles for p for given s
        p_tm <- quant.reg(p_tilde,Bs,Ds,lambda,tau)
        p_est <- Bs_est %*% p_tm$b
        # Obtain tau^th conditional quantiles for y for given p and s
        p_tilde_rep <- rep(p_tilde,sim_ctys)
        p_kn <- quantile(p_tilde_rep,seq(1/Q,(Q-1)/Q,1/Q))
        Bp_rep <- bs(p_tilde_rep,degree=degree,knots=p_kn,Boundary.knots=c(-10,10),intercept = TRUE) # Bspline basis used in quantile regression for y
        Bp_est <- bs(rep(seq(-10,10,0.2),5),degree=degree,knots=p_kn,Boundary.knots=c(-10,10),intercept = TRUE) # Bspline basis used in producing plots for simulation evaluation
        Bps_rep <- cbind(Bp_rep,Bs_rep)
        Bps_est <- cbind(Bp_est,Bs_tru_quants)
        y_tm <- quant.reg(y_tilde,Bps_rep,Dps,lambda,tau)
        y_est <- Bps_est %*% y_tm$b
        p_sim <- cbind(p_sim,p_est)
        y_sim <- cbind(y_sim,y_est)
      }
      
      # Save the approximations (Note, p_app_lin and p_app_non are conducted twice, so no need to redo)
      if (l == 0){
        p_app_lin[,,m] <- p_sim
        y_app_lin_p_lin[,,m] <- y_sim
      } else if (l == 1) {
        y_app_non_p_lin[,,m] <- y_sim
      } else if (l == 2) {
        p_app_non[,,m] <- p_sim
        y_app_lin_p_non[,,m] <- y_sim
      } else {
        y_app_non_p_non[,,m] <- y_sim
      }
    }
    # For the 1st Monte Carlo simulation, we simulate an entire joint density function of y and p given a particular s
    # Our main goal is to approximate that conditional join density function
    if (m == 1){
      p_sim <- NULL
      y_sim <- NULL
      for (k in 1:K){
        # Simulate, for a given s,"K" values from the joint density function of y and p using quantile regression
        print(paste0("Running Simulation: ",l,".",k))
        Bs_quants <- bs(s_tru_quants[2:4],degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
        # Simulate p given s
        if (l %% 2 == 0){
          p_tilde_sim <- Bs_quants %*% quant.reg(p_tilde,Bs,Ds,lambda,tau_p[k])$b
        } else if (l == 1){
          p_tilde_sim <- p_app_sim_lin[,k]
        } else {
          p_tilde_sim <- p_app_sim_non[,k]
        }
        # Simulate y given p and s
        p_tilde_rep <- rep(p_tilde,sim_ctys)
        p_kn <- quantile(p_tilde_rep,seq(1/Q,(Q-1)/Q,1/Q))
        Bp <- bs(as.vector(p_tilde_sim),degree=degree,knots=p_kn,Boundary.knots = c(-10,10),intercept = TRUE) # Bspline basis used to get simulated y for a given tau_y
        Bp_rep <- bs(p_tilde_rep,degree=degree,knots=p_kn,Boundary.knots = c(-10,10),intercept = TRUE)
        Bs_rep <- bs(s_rep,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
        Bps_quants <- cbind(Bp,Bs_quants)
        Bps_rep <- cbind(Bp_rep,Bs_rep)
        y_tilde_sim <- Bps_quants %*% quant.reg(y_tilde,Bps_rep,Dps,lambda,tau_y[k])$b
        p_sim <- cbind(p_sim,p_tilde_sim)
        y_sim <- cbind(y_sim,y_tilde_sim)
      }
      if (l == 0){
        p_app_sim_lin <- p_sim
        y_app_sim_lin_p_lin <- y_sim
      } else if (l == 1) {
        y_app_sim_non_p_lin <- y_sim
      } else if (l == 2) {
        p_app_sim_non <- p_sim
        y_app_sim_lin_p_non <- y_sim
      } else {
        y_app_sim_non_p_non <- y_sim
      }
    }
  }
}

