# Running our dataset technique on the corn crop dataset
rm(list=ls())
library(readxl)
library(splines)
library(MASS)
library(tidyverse)
library(usmap)
library(tools)
library(latex2exp)
source("BSplineFunctions.R")

corn_price      <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_futures_price_Dec")[-1,]
corn_stocks     <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_ending_stock")
corn_natl_yield <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_production_nationwide")
corn_cty_yield  <- read_excel("data/futures_price_yield_200206.xlsx",sheet = "corn_yield_bu_acre")
corn_IV         <- read_excel("data/Implied_Volatility.xlsx",sheet = "Corn_IV")
corn_non_yield  <- read_excel("data/Corn_NonIrrigatedYield.xlsx")
GDPDEF <- read.csv("GDPDEF.csv")

corn_price      <- corn_price  %>% filter(year >= 1990 & year <= 2018)
corn_stocks     <- corn_stocks %>% filter(year >= 1989 & year <= 2017)
corn_natl_yield <- corn_natl_yield %>% filter(Year >= 1989 & Year <= 2017)
corn_cty_yield  <- corn_cty_yield  %>% filter(Year >= 1980 & Year <= 2018) %>%
  filter(! (State %in% c("NEBRASKA","KANSAS")))
corn_cty_yield  <- rbind(corn_cty_yield,corn_non_yield) %>% 
  arrange(Year,State)

D <- function(d){
  Ds <- matrix(0,nrow=d-2,ncol=d)
  for (j in 1:nrow(Ds)){
    for (k in 1:ncol(Ds)){
      if (k-j == 0 || k-j == 2){
        Ds[j,k] = 1
      } else if(k-j == 1){
        Ds[j,k] = -2
      }
    }
  }
  return(Ds)
}

kron <- function(A,B){
  A[,rep(seq(ncol(A)), each = ncol(B))] * B[,rep(seq(ncol(B)),ncol(A))]
}

set.seed(10062018) # For reproducability
Q <- 4
K <- 1000
degree <- 3
states <- toTitleCase(tolower(unique(corn_cty_yield$State)))
tau_p <- runif(K) # Fixed values of tau_p used in the singular simulation
tau_y <- runif(K) # Fixed values of tau_y used in the singular simulation

t_tilde <- corn_cty_yield$Year - min(corn_cty_yield$Year) + 1
Bt <- bs(t_tilde,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
lambda <- 0.5

Bc <- model.matrix(~0+State,corn_cty_yield)
Bct <- kron(Bc,Bt)
Bc2020 <- kron(diag(12),matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE),each=12),nrow=12))

prem_calc <- function(p,y,phi,p_harv,y_bar){
  mean(pmax(phi*p_harv*y_bar - p*y,0))
}

y <- corn_cty_yield$`Yield (Bu/Acre)`
y_beta <- mean.reg(y,Bct,D(ncol(Bct)),lambda)
y_hat <- Bct %*% y_beta
corn_cty_yield$y_tilde <- y - as.vector(y_hat)
y_adj_2020 <- Bc2020 %*% y_beta
y_bar <- sapply(1990:2018,function(t){
  tmp <- corn_cty_yield %>%
    filter(Year >= t-10 & Year <= t-1) %>%
    group_by(State) %>%
    summarise(bar = mean(y_tilde))
  return(tmp$bar)
}) + matrix(rep(y_adj_2020,29),nrow=12)

s <- corn_stocks$`Ending Stocks (Million Bushels)`
s_hat <- loess(`Production Measured in Million Bushels`~Year,
               corn_natl_yield)
s_tilde <- s/s_hat$fitted
s_kn <- quantile(s_tilde,seq(1/Q,(Q-1)/Q,1/Q))
Bs <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
Gs <- matrix(0,nrow=ncol(Bs),ncol=ncol(Bs))
for(t in 1:29){
  Gs <- Gs + 1/29*(Bs[t,] %*% t(Bs[t,]))
}

corn_cty_yield <- corn_cty_yield %>% filter(Year >= 1990)
p_fut <- as.numeric(corn_price$`Feb_price (dollar)`)*GDPDEF[17:45,3]
p_harv <- as.numeric(corn_price$`Dec_price (dollar)`)*GDPDEF[17:45,3]
corn_price$p <- log(p_harv)
n_ctys_per_year <- corn_cty_yield %>%
  group_by(Year) %>%
  summarise(n = n())

t_tilde <- corn_cty_yield$Year - min(corn_cty_yield$Year) + 1
Bt <- bs(t_tilde,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
Bt_uni <- bs(unique(t_tilde),knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
# lambda <- 0.5
Gt <- matrix(0,nrow=ncol(Bt_uni),ncol=ncol(Bt_uni))
for(t in 1:29){
  Gt <- Gt + 1/29*(Bt_uni[t,] %*% t(Bt_uni[t,]))
}

p_beta <- mean.reg(corn_price$p,Bt_uni,D(ncol(Bt_uni)),lambda)
p_hat <- Bt_uni %*% p_beta
p_tilde <- rep(corn_price$p - as.vector(p_hat),n_ctys_per_year$n)
var_p_tilde <- var(unique(p_tilde))
var_p_hat <- sapply(1:29,function(t){
  (var_p_tilde/29*t(Bt_uni[t,]) %*% MASS::ginv(Gt) %*% Bt_uni[t,])[1,1]
})

Bc <- model.matrix(~0+State,corn_cty_yield)
Bct <- kron(Bc,Bt)
B2020 <- matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)))
Bc2020 <- kron(diag(12),matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE),each=12),nrow=12))
y_beta <- mean.reg(corn_cty_yield$`Yield (Bu/Acre)`,Bct,D(ncol(Bct)),lambda)
y_hat <- Bct %*% y_beta
y_tilde <- corn_cty_yield$`Yield (Bu/Acre)` - as.vector(y_hat)
y_hat_2020 <- Bc2020 %*% y_beta
var_y_tilde <- corn_cty_yield %>% 
  mutate(y_tilde = y_tilde) %>% 
  group_by(State) %>% 
  summarise(var_y_tilde = var(y_tilde),
            n = n())
var_y_hat <- var_y_tilde$var_y_tilde/var_y_tilde$n * 
  (t(B2020) %*% MASS::ginv(Gt) %*% B2020)[1,1]

s_tilde <- rep(s_tilde,n_ctys_per_year$n)
s_kn <- quantile(unique(s_tilde),seq(1/Q,(Q-1)/Q,1/Q))
Bs_rep <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
p_kn <- quantile(unique(p_tilde),seq(1/Q,(Q-1)/Q,1/Q))
Bp_rep <- bs(p_tilde,degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)

IV <- corn_IV$`Volatility Factor`[1:29]*sqrt(9/12)/100
IV_hat <- Bs %*% mean.reg(IV,Bs,D(ncol(Bs)),lambda)
var_IV_tilde <- var(IV - IV_hat)
var_IV_hat <- sapply(1:29,function(t){
  (var_IV_tilde/29*t(Bs[t,]) %*% MASS::ginv(Gs) %*% Bs[t,])[1,1]
})

p_fut_hat <- Bs %*% mean.reg(p_fut,Bs,D(ncol(Bs)),lambda)
var_p_fut_tilde <- var(p_fut - p_fut_hat)
var_p_fut_hat <- sapply(1:29,function(t){
  (var_p_fut_tilde/29*t(Bs[t,]) %*% MASS::ginv(Gs) %*% Bs[t,])[1,1]
})

p_beta_cond <- array(NA,dim=c(K,7))
p_beta_uncond <- array(NA,dim=c(K,1))
y_beta_cond <- array(NA,dim=c(K,14*12))
y_beta_uncond <- array(NA,dim=c(K,7*12))

for (k in 1:K){
  print(paste0("Simulating Corn Value ",k))
  # Simulate p given s for conditional density
  p_beta_cond[k,] <- quant.reg(unique(p_tilde),Bs,D(ncol(Bs)),lambda=.45,tau_p[k])$b
  # Simulate y given p and s for conditional density
  i = 0
  for (state in unique(corn_cty_yield$State)){
    i <- i + 1
    Bps_rep <- cbind(Bp_rep,Bs_rep)[corn_cty_yield$State==state,]
    y_beta_cond[k,((i-1)*14+1):(i*14)] <- quant.reg(y_tilde[corn_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
  }
  # Simulate p for unconditional density
  p_beta_uncond[k,] <- quantile(unique(p_tilde),tau_p[k])
  # Simulate y given p for unconditional density
  i = 0
  for (state in unique(corn_cty_yield$State)){
    i <- i + 1
    Bps_rep <- Bp_rep[corn_cty_yield$State==state,]
    y_beta_uncond[k,((i-1)*7+1):(i*7)] <- quant.reg(y_tilde[corn_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
  }
}

saveRDS(p_beta_uncond,"p_beta_uncond_Corn.RDS")
saveRDS(y_beta_uncond,"y_beta_uncond_Corn_NonIrr.RDS")
saveRDS(p_beta_cond,"p_beta_cond_Corn.RDS")
saveRDS(y_beta_cond,"y_beta_cond_Corn_NonIrr.RDS")

B = 50
y_beta_cond_jk <- array(NA,dim=c(K,14*12,B))
y_beta_uncond_jk <- array(NA,dim=c(K,7*12,B))

dat <- data.frame(Year = corn_cty_yield$Year,
                  State = corn_cty_yield$State,
                  County = corn_cty_yield$County,
                  s_tilde = s_tilde,
                  p_tilde = p_tilde,
                  y_tilde = y_tilde) %>%
  arrange(State,County) %>%
  group_by(State) %>%
  mutate(grp = rep(1:B,length.out=n()))

for (b in 1:B){
  tmp <- dat[dat$grp != b,]
  
  s_tilde <- tmp$s_tilde
  p_tilde <- tmp$p_tilde
  y_tilde <- tmp$y_tilde
  Bs_rep <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
  Bp_rep <- bs(p_tilde,degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  
  for (k in 1:K){
    print(paste0("Simulating Corn Value ",k," for Jackknife ",b))
    # Simulate y given p and s for conditional density
    i = 0
    for (state in unique(corn_cty_yield$State)){
      i <- i + 1
      Bps_rep <- cbind(Bp_rep,Bs_rep)[tmp$State==state,]
      y_beta_cond_jk[k,((i-1)*14+1):(i*14),b] <- quant.reg(y_tilde[tmp$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    }
    
    # Simulate y given p for uconditional density
    i = 0
    for (state in unique(corn_cty_yield$State)){
      i <- i + 1
      Bps_rep <- Bp_rep[tmp$State==state,]
      y_beta_uncond_jk[k,((i-1)*7+1):(i*7),b] <- quant.reg(y_tilde[tmp$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    }
  }
}

saveRDS(y_beta_uncond_jk,"y_beta_uncond_jk_Corn_NonIrr.RDS")
saveRDS(y_beta_cond_jk,"y_beta_cond_jk_Corn_NonIrr.RDS")
