# Running our dataset technique on the soybean crop dataset
rm(list=ls())
library(readxl)
library(splines)
library(MASS)
library(tidyverse)
library(usmap)
library(tools)
library(latex2exp)
source("BSplineFunctions.R")

soybean_price      <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "soybean_futures_price_Nov")[-1,]
soybean_stocks     <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "soybean_ending_stock")
soybean_natl_yield <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "soybean_production_nationwide")
soybean_cty_yield  <- read_excel("data/futures_price_yield_200206.xlsx",sheet = "soybean_yield_bu_acre")
soybean_IV         <- read_excel("data/Implied_Volatility.xlsx",sheet = "Soybean_IV")
soybean_non_yield  <- read_excel("data/Soybean_NonIrrigatedYield.xlsx")
GDPDEF <- read.csv("GDPDEF.csv")

soybean_price      <- soybean_price  %>% filter(year >= 1990 & year <= 2018)
soybean_stocks     <- soybean_stocks %>% filter(year >= 1989 & year <= 2017)
soybean_natl_yield <- soybean_natl_yield %>% filter(Year >= 1989 & Year <= 2017)
soybean_cty_yield  <- soybean_cty_yield  %>% filter(Year >= 1980 & Year <= 2018) %>%
  filter(! (State %in% c("NEBRASKA","KANSAS")))
soybean_cty_yield  <- rbind(soybean_cty_yield,soybean_non_yield)%>% 
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
states <- toTitleCase(tolower(unique(soybean_cty_yield$State)))
tau_p <- runif(K) # Fixed values of tau_p used in the singular simulation
tau_y <- runif(K) # Fixed values of tau_y used in the singular simulation

t_tilde <- soybean_cty_yield$Year - min(soybean_cty_yield$Year) + 1
Bt <- bs(t_tilde,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
lambda <- 0.5

Bc <- model.matrix(~0+State,soybean_cty_yield)
Bct <- kron(Bc,Bt)
Bc2020 <- kron(diag(12),matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE),each=12),nrow=12))

prem_calc <- function(p,y,phi,p_harv,y_bar){
  mean(pmax(phi*p_harv*y_bar - p*y,0))
}

y <- soybean_cty_yield$`Yield (Bu/Acre)`
y_beta <- mean.reg(y,Bct,D(ncol(Bct)),lambda)
y_hat <- Bct %*% y_beta
soybean_cty_yield$y_tilde <- y - as.vector(y_hat)
y_adj_2020 <- Bc2020 %*% y_beta
y_bar <- sapply(1990:2018,function(t){
  tmp <- soybean_cty_yield %>%
    filter(Year >= t-10 & Year <= t-1) %>%
    group_by(State) %>%
    summarise(bar = mean(y_tilde))
  return(tmp$bar)
}) + matrix(rep(y_adj_2020,29),nrow=12)

s <- soybean_stocks$`Ending Stocks \r\n(Million Bushels)`
s_hat <- loess(`Production Measured in Million Bushels`~Year,
               soybean_natl_yield)
s_tilde <- s/s_hat$fitted
s_kn <- quantile(s_tilde,seq(1/Q,(Q-1)/Q,1/Q))
Bs <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
Gs <- matrix(0,nrow=ncol(Bs),ncol=ncol(Bs))
for(t in 1:29){
  Gs <- Gs + 1/29*(Bs[t,] %*% t(Bs[t,]))
}

soybean_cty_yield <- soybean_cty_yield %>% filter(Year >= 1990)
p_fut <- as.numeric(soybean_price$`Feb_price (dollar)`)*GDPDEF[17:45,3]
p_harv <- as.numeric(soybean_price$`Nov_price (dollar)`)*GDPDEF[17:45,3]
soybean_price$p <- log(p_harv)
n_ctys_per_year <- soybean_cty_yield %>%
  group_by(Year) %>%
  summarise(n = n())

t_tilde <- soybean_cty_yield$Year - min(soybean_cty_yield$Year) + 1
Bt <- bs(t_tilde,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
Bt_uni <- bs(unique(t_tilde),knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
# lambda <- 0.5
Gt <- matrix(0,nrow=ncol(Bt_uni),ncol=ncol(Bt_uni))
for(t in 1:29){
  Gt <- Gt + 1/29*(Bt_uni[t,] %*% t(Bt_uni[t,]))
}

p_beta <- mean.reg(soybean_price$p,Bt_uni,D(ncol(Bt_uni)),lambda)
p_hat <- Bt_uni %*% p_beta
p_tilde <- rep(soybean_price$p - as.vector(p_hat),n_ctys_per_year$n)
var_p_tilde <- var(unique(p_tilde))
var_p_hat <- sapply(1:29,function(t){
  (var_p_tilde/29*t(Bt_uni[t,]) %*% MASS::ginv(Gt) %*% Bt_uni[t,])[1,1]
})

Bc <- model.matrix(~0+State,soybean_cty_yield)
Bct <- kron(Bc,Bt)
B2020 <- matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)))
Bc2020 <- kron(diag(12),matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE),each=12),nrow=12))
y_beta <- mean.reg(soybean_cty_yield$`Yield (Bu/Acre)`,Bct,D(ncol(Bct)),lambda)
y_hat <- Bct %*% y_beta
y_tilde <- soybean_cty_yield$`Yield (Bu/Acre)` - as.vector(y_hat)
y_hat_2020 <- Bc2020 %*% y_beta
var_y_tilde <- soybean_cty_yield %>% 
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

IV <- soybean_IV$`Volatility Factor`[1:29]*sqrt(9/12)/100
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

# p_beta_cond <- array(NA,dim=c(K,7))
# p_beta_uncond <- array(NA,dim=c(K,1))
# y_beta_cond <- array(NA,dim=c(K,14*12))
# y_beta_uncond <- array(NA,dim=c(K,7*12))
p_beta_cond <- readRDS("p_beta_cond_Soybean.RDS")
y_beta_cond <- readRDS("y_beta_cond_Soybean.RDS")
y_beta_cond_jk <- readRDS("y_beta_cond_jk_Soybean.RDS")
p_beta_uncond <- readRDS("p_beta_uncond_Soybean.RDS")
y_beta_uncond <- readRDS("y_beta_uncond_Soybean.RDS")
y_beta_uncond_jk <- readRDS("y_beta_uncond_jk_Soybean.RDS")

for (k in 1:K){
  print(paste0("Simulating Soybean Value ",k))
  # Simulate p given s for conditional density
  y_beta_sim <- NULL
  # lambda <- optim.lambda(unique(p_tilde),Bs,D(ncol(Bs)),tau_p[k])$estimate
  # p_beta_cond[k,] <- quant.reg(unique(p_tilde),Bs,D(ncol(Bs)),lambda,tau_p[k])$b
  # # Simulate y given p and s for conditional density
  # Bp <- bs(as.vector(p_tilde_sim),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  # Bps <- cbind(Bp,Bs)
  i = 0
  for (state in unique(soybean_cty_yield$State)){
  # for (state in c("NEBRASKA","KANSAS")){
    i <- i + 1
    # if (state == "KANSAS"){
    #   i <- 4
    # } else {
    #   i <- 8
    # }
    Bps_rep <- cbind(Bp_rep,Bs_rep)[soybean_cty_yield$State==state,]
    # tmp <- Bps %*% quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    # y_tilde_sim <- c(y_tilde_sim,as.vector(tmp))
    y_beta_cond[k,((i-1)*14+1):(i*14)] <- quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
  }
  
  # Bps <- kron(Bp,Bs)
  # Bps_rep <- kron(Bp_rep,Bs_rep)
  # y_tilde_sim <- Bps %*% quant.reg(y_tilde,Bps_rep,D(49),lambda,tau_y[k])$b
  # y_sim_cond_kron <- cbind(y_sim_cond_kron,y_tilde_sim)
  
  # Simulate p given s for conditional density
  y_beta_sim <- NULL
  p_beta_uncond[k,] <- quantile(unique(p_tilde),tau_p[k])
  # Simulate y given p and s for conditional density
  # Bp <- bs(rep(p_tilde_sim,29),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  # Bps <- Bp
  i = 0
  for (state in unique(soybean_cty_yield$State)){
  # for (state in c("NEBRASKA","KANSAS")){
    i <- i + 1
    # if (state == "KANSAS"){
    #   i <- 4
    # } else {
    #   i <- 8
    # }
    Bps_rep <- Bp_rep[soybean_cty_yield$State==state,]
    # tmp <- Bps %*% quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    # y_tilde_sim <- c(y_tilde_sim,as.vector(tmp))
    y_beta_uncond[k,((i-1)*7+1):(i*7)] <- quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
  }
  
  # y_tilde_sim <- Bp %*% quant.reg(y_tilde,Bp_rep,D(7),lambda,tau_y[k])$b
  # y_sim_uncond_kron <- cbind(y_sim_uncond_kron,y_tilde_sim)
}

saveRDS(p_beta_uncond,"p_beta_uncond_Soybean_NonIrr.RDS")
saveRDS(y_beta_uncond,"y_beta_uncond_Soybean_NonIrr.RDS")
saveRDS(p_beta_cond,"p_beta_cond_Soybean_NonIrr.RDS")
saveRDS(y_beta_cond,"y_beta_cond_Soybean_NonIrr.RDS")

# sd_p_uncond <- apply(p_sim_uncond,1,sd)
# sd_p_cond <- apply(p_sim_cond,1,sd)
# 
# cor_uncond <- sapply(1:12,function(s){
#   sapply(1:29,function(t){
#     cor(exp(p_hat[t] + p_sim_uncond[t,]),y_sim_uncond[29*(s-1)+t,])
#   })
# }) %>% matrix(nrow=29)
# 
# cor_cond <- sapply(1:12,function(s){
#   sapply(1:29,function(t){
#     cor(exp(p_hat[t] + p_sim_cond[t,]),y_sim_cond[29*(s-1)+t,])
#   })
# }) %>% matrix(nrow=29)
# 
# dt <- data.frame(State = rep(states,each=29*2),
#                  s = rep(unique(s_tilde),each=2,times=12),
#                  Year = rep(1990:2018,each=2,times=12),
#                  Method = rep(c("Proposed","Current"),times=12*29),
#                  phi = as.character(rep(c(0.7),times=12*29*2)))
# 
# dt$prem <- apply(dt,1,function(x){
#   s <- which(states==x[1])
#   t <- as.numeric(x[3]) - 1989
#   p <- exp(p_hat[t] + p_sim_cond[t,])
#   y <- y_adj_2020[s,1] + y_sim_cond[29*(s-1)+t,]
#   if (x[4] == "Current"){
#     set.seed(10062018)
#     y <- y_adj_2020[s,1] + y_sim_uncond[29*(s-1)+t,]
#     cor <- cor_uncond[t,s]
#     eps <- mvrnorm(K,mu=c(0,0),Sigma=matrix(c(1,cor,cor,1),nrow=2))
#     p <- quantile(p,pnorm(eps[,1]))
#     y <- quantile(y,pnorm(eps[,2]))
#   }
#   prem_calc(p,y,as.numeric(x[5]),p_harv[t],y_bar[s,t])
# })
# 
# prems_prop <- (dt %>%
#     filter(phi == 0.7) %>%
#     group_by(State,Year) %>%
#     summarise(PremDiff=prem[Method=="Proposed"]) %>%
#     ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
# 
# prems_curr <- (dt %>%
#                  filter(phi == 0.7) %>%
#                  group_by(State,Year) %>%
#                  summarise(PremDiff=prem[Method=="Current"]) %>%
#                  ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
# 
# saveRDS(sd_p_uncond,"Soybean_sd_p_uncond.RDS")
# saveRDS(sd_p_cond,"Soybean_sd_p_cond.RDS")
# saveRDS(cor_uncond,"Soybean_cor_uncond.RDS")
# saveRDS(cor_cond,"Soybean_cor_cond.RDS")
# saveRDS(prems_prop,"Soybean_prems_prop.RDS")
# saveRDS(prems_curr,"Soybean_prems_curr.RDS")

B = 50

dat <- data.frame(Year = soybean_cty_yield$Year,
                  State = soybean_cty_yield$State,
                  County = soybean_cty_yield$County,
                  s_tilde = s_tilde,
                  p_tilde = p_tilde,
                  y_tilde = y_tilde) %>%
  arrange(State,County) %>%
  group_by(State) %>%
  mutate(grp = rep(1:B,length.out=n()))

# sd_p_cond_jk <- matrix(NA,nrow=29,ncol=B)
# cor_cond_jk <- array(NA,dim=c(29,12,B))
# sd_p_uncond_jk <- matrix(NA,nrow=29,ncol=B)
# cor_uncond_jk <- array(NA,dim=c(29,12,B))
# 
# prems_prop_jk <- array(NA,dim=c(29,12,B))
# prems_curr_jk <- array(NA,dim=c(29,12,B))

# y_beta_cond_jk <- array(NA,dim=c(K,14*12,B))
# y_beta_uncond_jk <- array(NA,dim=c(K,7*12,B))

for (b in 1:B){
  tmp <- dat[dat$grp != b,]
  
  s_tilde <- tmp$s_tilde
  p_tilde <- tmp$p_tilde
  y_tilde <- tmp$y_tilde
  # s_kn <- quantile(unique(s_tilde_tmp),seq(1/Q,(Q-1)/Q,1/Q))
  # Bs <- bs(unique(s_tilde),degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
  Bs_rep <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
  # p_kn <- quantile(unique(p_tilde_tmp),seq(1/Q,(Q-1)/Q,1/Q))
  Bp_rep <- bs(p_tilde,degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  
  for (k in 1:K){
    print(paste0("Simulating Soybean Value ",k," for Jackknife ",b))
    # Simulate y given p and s for conditional density
    y_beta_sim <- NULL
    # Bp <- bs(as.vector(p_tilde_sim),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
    # Bps <- cbind(Bp,Bs)
    i = 0
    for (state in unique(soybean_cty_yield$State)){
    # for (state in c("NEBRASKA","KANSAS")){
      i <- i + 1
      # if (state == "KANSAS"){
      #   i <- 4
      # } else {
      #   i <- 8
      # }
      Bps_rep <- cbind(Bp_rep,Bs_rep)[tmp$State==state,]
      # tmp <- Bps %*% quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
      # y_tilde_sim <- c(y_tilde_sim,as.vector(tmp))
      y_beta_cond_jk[k,((i-1)*14+1):(i*14),b] <- quant.reg(y_tilde[tmp$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    }
    
    # Bps <- kron(Bp,Bs)
    # Bps_rep <- kron(Bp_rep,Bs_rep)
    # y_tilde_sim <- Bps %*% quant.reg(y_tilde,Bps_rep,D(49),lambda,tau_y[k])$b
    # y_sim_cond_kron <- cbind(y_sim_cond_kron,y_tilde_sim)
    
    # Simulate p given s for conditional density
    y_beta_sim <- NULL
    # p_beta_sim <- rep(quantile(unique(p_tilde),tau_p[k]),length(unique(s_tilde)))
    # Simulate y given p and s for conditional density
    # Bp <- bs(as.vector(p_tilde_sim),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
    # Bps <- Bp
    i = 0
    for (state in unique(soybean_cty_yield$State)){
    # for (state in c("NEBRASKA","KANSAS")){
      i <- i + 1
      # if (state == "KANSAS"){
      #   i <- 4
      # } else {
      #   i <- 8
      # }
      Bps_rep <- Bp_rep[tmp$State==state,]
      # tmp <- Bps %*% quant.reg(y_tilde[soybean_cty_yield$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
      # y_tilde_sim <- c(y_tilde_sim,as.vector(tmp))
      y_beta_uncond_jk[k,((i-1)*7+1):(i*7),b] <- quant.reg(y_tilde[tmp$State==state],Bps_rep,D(ncol(Bps_rep)),lambda=1e-8,tau_y[k])$b
    }
    
    # y_tilde_sim <- Bp %*% quant.reg(y_tilde,Bp_rep,D(7),lambda,tau_y[k])$b
    # y_sim_uncond_kron <- cbind(y_sim_uncond_kron,y_tilde_sim)
  }
  
  
  
  # sd_p_cond_jk[,b] <-  apply(p_sim_cond,1,sd)
  # sd_p_uncond_jk[,b] <-  apply(p_sim_uncond,1,sd)
  # 
  # cor_uncond <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     cor(exp(p_hat[t] + p_sim_uncond[t,]),y_sim_uncond[29*(s-1)+t,])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # cor_uncond_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     cor(exp(p_hat[t] + p_sim_uncond[t,]),y_sim_uncond[29*(s-1)+t,])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # cor_cond_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     cor(exp(p_hat[t] + p_sim_cond[t,]),y_sim_cond[29*(s-1)+t,])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # dt <- data.frame(State = rep(states,each=29*2),
  #                  s = rep(unique(s_tilde),each=2,times=12),
  #                  Year = rep(1990:2018,each=2,times=12),
  #                  Method = rep(c("Proposed","Current"),times=12*29),
  #                  phi = as.character(rep(c(0.7),times=12*29*2)))
  # 
  # dt$prem <- apply(dt,1,function(x){
  #   s <- which(states==x[1])
  #   t <- as.numeric(x[3]) - 1989
  #   p <- exp(p_hat[t] + p_sim_cond[t,])
  #   y <- y_adj_2020[s,1] + y_sim_cond[29*(s-1)+t,]
  #   if (x[4] == "Current"){
  #     set.seed(10062018)
  #     y <- y_adj_2020[s,1] + y_sim_uncond[29*(s-1)+t,]
  #     cor <- cor_uncond[t,s]
  #     eps <- mvrnorm(K,mu=c(0,0),Sigma=matrix(c(1,cor,cor,1),nrow=2))
  #     p <- quantile(p,pnorm(eps[,1]))
  #     y <- quantile(y,pnorm(eps[,2]))
  #   }
  #   prem_calc(p,y,as.numeric(x[5]),p_harv[t],y_bar[s,t])
  # })
  # 
  # prems_prop_jk[,,b] <- (dt %>%
  #                     filter(phi == 0.7) %>%
  #                     group_by(State,Year) %>%
  #                     summarise(PremDiff=prem[Method=="Proposed"]) %>%
  #                     ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
  # 
  # prems_curr_jk[,,b] <- (dt %>%
  #                          filter(phi == 0.7) %>%
  #                          group_by(State,Year) %>%
  #                          summarise(PremDiff=prem[Method=="Current"]) %>%
  #                          ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
}

saveRDS(y_beta_uncond_jk,"y_beta_uncond_jk_Soybean_NonIrr.RDS")
saveRDS(y_beta_cond_jk,"y_beta_cond_jk_Soybean_NonIrr.RDS")