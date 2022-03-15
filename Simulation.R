###### Note: Must run detrending to get the y_hats to adjust to y
rm(list = ls())
sim_ctys <- 500 # Number of county level yields per year to simulate
sim_year <- 100 # Number of years worth of data to simulate
M        <- 100 # Number of Monte Carlo Simulations
J        <- 100 # Number of approximate quantiles to obtain for each Monte Carlo simulation
K        <- 5000 # Number of Simulations for the 
Q        <- 4 # 1 + number of quantiles used in determining the B-spline knots (i.e. using 1/Q to (Q-1)/Q quantiles as the knots)

source("Corn.R") # Get detrended yield and price and functions to run B-spline quantile regression
library(tidyverse)
library(splines)
library(readxl)
library(MASS)
library(sn)

set.seed(10062018) # For reproducability
tau <- seq(1/J,1,1/J) - 1/(2*J) # Fixed values of tau used for the quantiles

# Simulate detrended stocks, s: Based on method of moments estimation of s
s <- rbeta(sim_year,1.5,6)
s_tru_quants <- qbeta(seq(1/Q,(Q-1)/Q,1/Q),1.5,6)
s_rep <- rep(s,sim_ctys)
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

Dps <- matrix(0,nrow=2*ncol(Bs)-2,ncol=2*ncol(Bs))

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
p_tru_lin <- array(0,dim=c(Q-1,J,M))
p_app_lin <- array(0,dim=c(Q-1,J,M))
p_tru_non <- array(0,dim=c(Q-1,J,M))
p_app_non <- array(0,dim=c(Q-1,J,M))
y_tru_lin_p_lin <- array(0,dim=c((Q-1)^2,J,M))
y_app_lin_p_lin <- array(0,dim=c((Q-1)^2,J,M))
y_tru_lin_p_non <- array(0,dim=c((Q-1)^2,J,M))
y_app_lin_p_non <- array(0,dim=c((Q-1)^2,J,M))
y_tru_non_p_lin <- array(0,dim=c((Q-1)^2,J,M))
y_app_non_p_lin <- array(0,dim=c((Q-1)^2,J,M))
y_tru_non_p_non <- array(0,dim=c((Q-1)^2,J,M))
y_app_non_p_non <- array(0,dim=c((Q-1)^2,J,M))

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
  a_p <- 3
  d_p <- a_p / sqrt(1 + a_p^2)
  eps_p <- rsn(sim_year,-d_p*sqrt(2/pi),(1-2*d_p^2/pi)^-1,a_p)
  
  m_y_lin <- function(p,s){-1.11 + 14.45*p + 22.18*s}
  m_y_non <- function(p,s){-5.31 + 4.67*p - 9.01*s + 136.04*p^2 + 52.56*s^2}
  s_y <- 31
  a_y <- -3
  d_y <- a_y / sqrt(1 + a_y^2)
  eps_y <- rsn(sim_year*sim_ctys,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)
  
  # Run simulations for all four possible combinations of linear and non-linear p and y
  for (l in 0:3){
    # if l (integer divide) 2 is 0, then p is linear, else p is non-linear
    if (l %/% 2 == 0){
      p_tilde <- m_p_lin(s) + s_p_lin(s)*eps_p
      # True quantiles for which to compare our simulation results
      p_tru_lin[,,m] <- matrix(m_p_lin(s_tru_quants) + s_p_lin(s_tru_quants)*
                                 qsn(rep(tau,each=Q-1),
                                            -d_p*sqrt(2/pi),
                                            (1-2*d_p^2/pi)^-1,
                                            a_p),
                               nrow=Q-1)
      # True quantiles used for calculation of B-spline basis for quantile regression of y
      p_tru_quants <- m_p_lin(rep(s_tru_quants,each=Q-1)) + 
        s_p_lin(rep(s_tru_quants,each=Q-1))*qsn(seq(1/Q,(Q-1)/Q,1/Q),
                                                -d_p*sqrt(2/pi),
                                                (1-2*d_p^2/pi)^-1,
                                                a_p)
    } else {
      p_tilde <- m_p_non(s) + s_p_non(s)*eps_p
      # True quantiles for which to compare our simulation results
      p_tru_non[,,m] <- matrix(m_p_non(s_tru_quants) + s_p_non(s_tru_quants)*
                                 qsn(rep(tau,each=Q-1),
                                                 -d_p*sqrt(2/pi),
                                                 (1-2*d_p^2/pi)^-1,
                                                 a_p),
                               nrow=Q-1)
      
      # True quantiles used for calculation of B-spline basis for quantile regression of y
      p_tru_quants <- m_p_non(rep(s_tru_quants,each=Q-1)) + 
        s_p_non(rep(s_tru_quants,each=Q-1))*qsn(seq(1/Q,(Q-1)/Q,1/Q),
                                                -d_p*sqrt(2/pi),
                                                (1-2*d_p^2/pi)^-1,
                                                a_p)
    }
    # if l mod 2 is 0, then y is linear, else y is non-linear
    if (l %% 2 == 0){
      y_tilde <- m_y_lin(p_tilde,s) + s_y*eps_y
      if (l %/% 2 == 0){
        # True quantiles for which to compare our simulation results
        y_tru_lin_p_lin[,,m] <- matrix(m_y_lin(p_tru_quants,rep(s_tru_quants,each=Q-1)) + s_y*qsn(rep(tau,each=(Q-1)^2),
                                                 -d_y*sqrt(2/pi),
                                                 (1-2*d_y^2/pi)^-1,
                                                 a_y),
                               nrow=(Q-1)^2)
      } else {
        # True quantiles for which to compare our simulation results
        y_tru_lin_p_non[,,m] <-matrix(m_y_lin(p_tru_quants,rep(s_tru_quants,each=Q-1)) + s_y*qsn(rep(tau,each=(Q-1)^2),
                                                                                                        -d_y*sqrt(2/pi),
                                                                                                        (1-2*d_y^2/pi)^-1,
                                                                                                        a_y),
                                      nrow=(Q-1)^2)
      }
    } else {
      y_tilde <- m_y_non(p_tilde,s) + s_y*eps_y
      if (l %/% 2 == 0){
        # True quantiles for which to compare our simulation results
        y_tru_non_p_lin[,,m] <- matrix(m_y_non(p_tru_quants,rep(s_tru_quants,each=Q-1)) + s_y*qsn(rep(tau,each=(Q-1)^2),
                                                                                                         -d_y*sqrt(2/pi),
                                                                                                         (1-2*d_y^2/pi)^-1,
                                                                                                         a_y),
                                       nrow=(Q-1)^2)
      } else {
        # True quantiles for which to compare our simulation results
        y_tru_non_p_non[,,m] <- matrix(m_y_non(p_tru_quants,rep(s_tru_quants,each=Q-1)) + s_y*qsn(rep(tau,each=(Q-1)^2),
                                                                                                         -d_y*sqrt(2/pi),
                                                                                                         (1-2*d_y^2/pi)^-1,
                                                                                                         a_y),
                                       nrow=(Q-1)^2)
      }
    }
    for (v in 1:n_lambdas){
      # Set smoothing parameter lambda
      lambda <- lambdas[v]
      p_sim <- NULL
      y_sim <- NULL
      for (j in 1:J){
        print(paste0("Getting Quantile ",j," from Simulation ",m,".",l," with lambda = ",lambdas[v]))
        
        # Obtain ((j-1)/J + 0.005)th conditional quantiles for p for given s in Bs_quants
        # Bs is used in obtaining quantile regression coefficients
        # Bs_quants is used with regression coefficients to obtain approximate quantiles
        s_kn <- quantile(s,seq(1/Q,(Q-1)/Q,1/Q))
        Bs <- bs(s,degree=3,knots=s_kn,intercept = TRUE)
        Bs_quants <- bs(s_tru_quants,degree=3,knots=s_kn,Boundary.knots = c(min(s),max(s)),intercept = TRUE)
        tau_p <- tau[j]
        p_tm <- quant.reg(p_tilde,Bs,Ds,lambda,tau_p)
        p_tilde_sim <- Bs_quants %*% p_tm$b
        
        # Obtain ((j-1)/J + 0.005)th conditional quantiles for p for given s in Bs_quants
        # Bp_rep is used in obtaining quantile regression coefficients
        # Bp_quants is used with regression coefficients to obtain approximate quantiles
        # Need to extend Bp and Bs because there are `sim_ctys` county level yields in a given simulated year
        tau_y <- tau[j]
        p_tilde_rep <- rep(p_tilde,sim_ctys)
        p_kn <- quantile(p_tilde_rep,seq(1/Q,(Q-1)/Q,1/Q))
        Bp_quants <- bs(p_tru_quants,degree=3,knots=p_kn,Boundary.knots = c(min(p_tilde),max(p_tilde)),intercept = TRUE)
        Bp_rep <- bs(p_tilde_rep,degree=3,knots=p_kn,intercept = TRUE)
        Bs_quants <- bs(rep(s_tru_quants,each=Q-1),degree=3,knots=s_kn,Boundary.knots = c(min(s),max(s)),intercept = TRUE)
        Bs_rep <- bs(s_rep,degree=3,knots=s_kn,intercept = TRUE)
        Bps_quants <- cbind(Bp_quants,Bs_quants)
        Bps_rep <- cbind(Bp_rep,Bs_rep)
        y_tm <- quant.reg(y_tilde,Bps_rep,Dps,lambda,tau_y)
        y_tilde_sim <- Bps_quants %*% y_tm$b
        p_sim <- cbind(p_sim,p_tilde_sim)
        y_sim <- cbind(y_sim,y_tilde_sim)
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
      if (l == 0){
        p_tru_sim_lin <- (m_p_lin(s_tru_quants) + s_p_lin(s_tru_quants)*rsn((Q-1)*K,-d_p*sqrt(2/pi),(1-2*d_p^2/pi)^-1,a_p)) %>% matrix(nrow=Q-1)
        y_tru_sim_lin_p_lin <- (m_y_lin(as.vector(p_tru_sim_lin),rep(s_tru_quants,K)) + s_y*rsn((Q-1)*K,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)) %>% matrix(nrow=(Q-1))
      } else if (l == 1) {
        y_tru_sim_non_p_lin <- (m_y_non(as.vector(p_tru_sim_lin),rep(s_tru_quants,K)) + s_y*rsn((Q-1)*K,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)) %>% matrix(nrow=(Q-1))
      } else if (l == 2) {
        p_tru_sim_non <- (m_p_non(s_tru_quants) + s_p_non(s_tru_quants)*rsn((Q-1)*K,-d_p*sqrt(2/pi),(1-2*d_p^2/pi)^-1,a_p)) %>% matrix(nrow=Q-1)
        y_tru_sim_lin_p_non <- (m_y_lin(as.vector(p_tru_sim_non),rep(s_tru_quants,K)) + s_y*rsn((Q-1)*K,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)) %>% matrix(nrow=(Q-1))
      } else {
        y_tru_sim_non_p_non <- (m_y_non(as.vector(p_tru_sim_non),rep(s_tru_quants,K)) + s_y*rsn((Q-1)*K,-d_y*sqrt(2/pi),(1-2*d_y^2/pi)^-1,a_y)) %>% matrix(nrow=(Q-1))
      }
      p_sim <- NULL
      y_sim <- NULL
      for (k in 1:K){
        print(paste0("Running Simulation: ",l,".",k))
        s_kn <- quantile(s,seq(1/Q,(Q-1)/Q,1/Q))
        Bs <- bs(s,degree=3,knots=s_kn,intercept = TRUE)
        Bs_quants <- bs(s_tru_quants,degree=3,knots=s_kn,Boundary.knots = c(min(s),max(s)),intercept = TRUE)
        if (l %% 2 == 0){
          tau_p <- runif(1)
          p_tilde_sim <- Bs_quants %*% quant.reg(p_tilde,Bs,Ds,lambda,tau_p)$b
        } else if (l == 1){
          p_tilde_sim <- p_app_sim_lin[,k]
        } else {
          p_tilde_sim <- p_app_sim_non[,k]
        }
        tau_y <- runif(1)
        p_tilde_rep <- rep(p_tilde,sim_ctys)
        p_kn <- quantile((p_tilde_rep),seq(1/Q,(Q-1)/Q,1/Q))
        Bp <- bs((as.vector(p_tilde_sim)),degree=3,knots=p_kn,Boundary.knots = c(min(p_tilde),max(p_tilde)),intercept = TRUE)
        Bp_rep <- bs((p_tilde_rep),degree=3,knots=p_kn,intercept = TRUE)
        Bs_quants <- bs(s_tru_quants,degree=3,knots=s_kn,Boundary.knots = c(min(s),max(s)),intercept = TRUE)
        Bs_rep <- bs(s_rep,degree=3,knots=s_kn,intercept = TRUE)
        Bps_quants <- cbind(Bp,Bs_quants)
        Bps_rep <- cbind(Bp_rep,Bs_rep)
        y_tilde_sim <- Bps_quants %*% quant.reg(y_tilde,Bps_rep,Dps,lambda,tau_y)$b
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

