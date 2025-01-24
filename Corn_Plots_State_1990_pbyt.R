# Running our dataset technique on the corn crop dataset
rm(list=ls())
library(readxl)
library(splines)
library(MASS)
library(tidyverse)
library(usmap)
library(tools)
library(latex2exp)
library(truncnorm)
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
B2020 <- matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)))
Bc2020 <- kron(diag(12),matrix(rep(bs(41,knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE),each=12),nrow=12))

prem_calc <- function(p,y,phi,p_harv,y_bar){
  mean(pmax(phi*p_harv*y_bar - p*y,0))
}


Bt_uni <- bs(unique(t_tilde),knots = c(10,20,30,40),Boundary.knots=c(0,50),intercept=TRUE)
# lambda <- 0.5
Gt <- matrix(0,nrow=ncol(Bt_uni),ncol=ncol(Bt_uni))
for(t in 11:39){
  Gt <- Gt + 1/29*(Bt_uni[t,] %*% t(Bt_uni[t,]))
}

y <- corn_cty_yield$`Yield (Bu/Acre)`
y_beta <- mean.reg(y,Bct,D(ncol(Bct)),lambda)
y_hat <- Bct %*% y_beta
corn_cty_yield$y_tilde <- y - as.vector(y_hat)
y_hat_2020 <- Bc2020 %*% y_beta
y_bar <- sapply(1990:2018,function(t){
  tmp <- corn_cty_yield %>%
    filter(Year >= t-10 & Year <= t-1) %>%
    group_by(State) %>%
    summarise(bar = mean(y_tilde))
  return(tmp$bar)
}) + matrix(rep(y_hat_2020,29),nrow=12)
var_y_tilde <- corn_cty_yield %>% 
  mutate(y_tilde = y_tilde) %>% 
  group_by(State) %>% 
  summarise(var_y_tilde = var(y_tilde),
            n = n())
var_y_hat <- var_y_tilde$var_y_tilde/var_y_tilde$n * 
  (t(B2020) %*% MASS::ginv(Gt) %*% B2020)[1,1]

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
Bt_uni <- Bt_uni[11:39,]
p_fut <- as.numeric(corn_price$`Feb_price (dollar)`)*GDPDEF[17:45,3]
p_harv <- as.numeric(corn_price$`Dec_price (dollar)`)*GDPDEF[17:45,3]
corn_price$p <- log(p_harv)
n_ctys_per_year <- corn_cty_yield %>%
  group_by(Year) %>%
  summarise(n = n())

p_beta <- mean.reg(corn_price$p,Bt_uni,D(ncol(Bt_uni)),lambda)
p_hat <- Bt_uni %*% p_beta
p_tilde <- rep(corn_price$p - as.vector(p_hat),n_ctys_per_year$n)
var_p_tilde <- var(unique(p_tilde))
var_p_hat <- sapply(1:29,function(t){
  (var_p_tilde/29*t(Bt_uni[t,]) %*% MASS::ginv(Gt) %*% Bt_uni[t,])[1,1]
})

s_tilde <- rep(s_tilde,n_ctys_per_year$n)
s_kn <- quantile(unique(s_tilde),seq(1/Q,(Q-1)/Q,1/Q))
Bs_rep <- bs(s_tilde,degree=degree,knots=s_kn,Boundary.knots = c(0,1),intercept = TRUE)
p_kn <- quantile(unique(p_tilde),seq(1/Q,(Q-1)/Q,1/Q))
Bp_rep <- bs(p_tilde,degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)

IV <- corn_IV$`Volatility Factor`[1:29]*sqrt(9/12)/100
IV_hat <- Bs %*% mean.reg(IV,Bs,D(ncol(Bs)),lambda)
var_IV_tilde <- var(IV - IV_hat)
var_IV_hat <- sapply(1:29,function(t){
  (var_IV_tilde/29*(t(Bs[t,]) %*% MASS::ginv(Gs) %*% Bs[t,]))[1,1]
})

p_fut_hat <- Bs %*% mean.reg(p_fut,Bs,D(ncol(Bs)),lambda)
var_p_fut_tilde <- var(p_fut - p_fut_hat)
var_p_fut_hat <- sapply(1:29,function(t){
  (var_p_fut_tilde/29*(t(Bs[t,]) %*% MASS::ginv(Gs) %*% Bs[t,]))[1,1]
})

p_beta_cond <- readRDS("p_beta_cond_Corn.RDS")
y_beta_cond <- readRDS("y_beta_cond_Corn_NonIrr.RDS")
y_beta_cond_jk <- readRDS("y_beta_cond_jk_Corn_NonIrr.RDS")
p_beta_uncond <- readRDS("p_beta_uncond_Corn.RDS")
y_beta_uncond <- readRDS("y_beta_uncond_Corn_NonIrr.RDS")
y_beta_uncond_jk <- readRDS("y_beta_uncond_jk_Corn_NonIrr.RDS")

states <- toTitleCase(tolower(unique(corn_cty_yield$State)))

p_fut_sim <- rnorm(29*K,p_fut_hat,sqrt(var_p_fut_hat)) %>%
  matrix(ncol=K)
p_hat_sim <- rnorm(29*K,p_hat,sqrt(var_p_hat)) %>% 
  matrix(ncol=K)
y_hat_sim <- rnorm(29*12*K,y_hat_2020[,1],sqrt(var_y_hat)) %>% 
  matrix(ncol=K)
y_bar_sim <- rnorm(29*12*K,as.vector(t(y_bar)),sqrt(var_y_hat)) %>% 
  matrix(ncol=K)
IV_sim <- rtruncnorm(29*K,IV_hat,sqrt(var_IV_hat),a=0) %>% 
  matrix(ncol=K)

p_sim_cond <- Bs %*% t(p_beta_cond)
p_sim_uncond <- matrix(p_beta_uncond,nrow=29,ncol=K,byrow=TRUE)
p_sim_IV <- (rnorm(29*K,0,IV_sim) %>% matrix(ncol=K)) + log(p_fut_sim) - p_hat_sim
y_sim_cond <- array(NA,dim=c(29*12,K))
y_sim_uncond <- array(NA,dim=c(29*12,K))
y_sim_cond_IV <- array(NA,dim=c(29*12,K))
y_sim_uncond_IV <- array(NA,dim=c(29*12,K))

for (k in 1:K){
  Bp <- bs(as.vector(p_sim_cond[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  Bp_uncond <- bs(as.vector(p_sim_uncond[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  Bp_IV <- bs(as.vector(p_sim_IV[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
  Bps <- cbind(Bp,Bs)
  Bps_IV <- cbind(Bp_IV,Bs)
  for (i in 1:length(states)){
    y_sim_cond[(29*(i-1)+1):(29*i),k] <- Bps %*% y_beta_cond[k,((i-1)*14+1):(i*14)]
    y_sim_uncond[(29*(i-1)+1):(29*i),k] <- Bp_uncond %*% y_beta_uncond[k,((i-1)*7+1):(i*7)]
    y_sim_cond_IV[(29*(i-1)+1):(29*i),k] <- Bps_IV %*% y_beta_cond[k,((i-1)*14+1):(i*14)]
    y_sim_uncond_IV[(29*(i-1)+1):(29*i),k] <- Bp_IV %*% y_beta_uncond[k,((i-1)*7+1):(i*7)]
  }
}

sd_p_cond <- apply(exp(p_hat_sim + p_sim_cond),1,sd)
sd_p_uncond <- apply(exp(p_hat_sim + p_sim_uncond),1,sd)

cor_uncond <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    cor(exp(p_hat_sim[t,] + p_sim_uncond[t,]),
        y_hat_sim[29*(s-1)+t,] + y_sim_uncond[29*(s-1)+t,])
  })
}) %>% matrix(nrow=29)

cor_cond <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    cor(exp(p_hat_sim[t,] + p_sim_cond[t,]),
        y_hat_sim[29*(s-1)+t,] + y_sim_cond[29*(s-1)+t,])
  })
}) %>% matrix(nrow=29)

prems_curr_phi0.7 <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    prem_calc(exp(p_hat_sim[t,] + p_sim_IV[t,]),
              y_hat_sim[29*(s-1)+t,] + y_sim_uncond_IV[29*(s-1)+t,],
              0.7,
              p_fut[t],
              y_bar[s,t])
  })
})

prems_curr_phi0.85 <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    prem_calc(exp(p_hat_sim[t,] + p_sim_IV[t,]),
              y_hat_sim[29*(s-1)+t,] + y_sim_uncond_IV[29*(s-1)+t,],
              0.85,
              p_fut[t],
              y_bar[s,t])
  })
})

# prems_prop_phi0.7 <- sapply(1:12,function(s){
#   sapply(1:29,function(t){
#     prem_calc(exp(p_hat_sim[t,] + p_sim_IV[t,]),
#               y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
#               0.7,
#               p_fut[t],
#               y_bar[s,t])
#   })
# })
# 
# prems_prop_phi0.85 <- sapply(1:12,function(s){
#   sapply(1:29,function(t){
#     prem_calc(exp(p_hat_sim[t,] + p_sim_IV[t,]),
#               y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
#               0.85,
#               p_fut[t],
#               y_bar[s,t])
#   })
# })

prems_prop_phi0.7 <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    prem_calc(exp(p_hat_sim[t,] + p_sim_IV[t,]),
              y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
              0.7,
              p_fut[t],
              y_bar[s,t])
  })
})

prems_prop_phi0.85 <- sapply(1:12,function(s){
  sapply(1:29,function(t){
    prem_calc(p_sim_cond[t,],
              y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
              0.85,
              p_fut[t],
              y_bar[s,t])
  })
})

B = 50
cor_cond_jk <- array(NA,dim=c(29,12,B))
cor_uncond_jk <- array(NA,dim=c(29,12,B))
cor_cond_IV_jk <- array(NA,dim=c(29,12,B))
cor_uncond_IV_jk <- array(NA,dim=c(29,12,B))

# prems_curr_phi0.7_jk <- array(NA,dim=c(29,12,B))
# prems_curr_phi0.9_jk <- array(NA,dim=c(29,12,B))
# prems_prop_phi0.7_jk <- array(NA,dim=c(29,12,B))
# prems_prop_phi0.9_jk <- array(NA,dim=c(29,12,B))
# 
for (b in 1:B){
  print(b)
  y_sim_cond <- array(NA,dim=c(29*12,K))
  y_sim_uncond <- array(NA,dim=c(29*12,K))
  for (k in 1:K){
    Bp <- bs(as.vector(p_sim_cond[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
    Bp_uncond <- bs(as.vector(p_sim_uncond[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
    # Bp_IV <- bs(as.vector(p_sim_IV[,k]),degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
    Bps <- cbind(Bp,Bs)
    # Bps_IV <- cbind(Bp_IV,Bs)
    for (i in 1:length(states)){
      y_sim_cond[(29*(i-1)+1):(29*i),k] <- Bps %*% y_beta_cond_jk[k,((i-1)*14+1):(i*14),b]
      y_sim_uncond[(29*(i-1)+1):(29*i),k] <- Bp_uncond %*% y_beta_uncond_jk[k,((i-1)*7+1):(i*7),b]
    #   y_sim_cond_IV[(29*(i-1)+1):(29*i),k] <- Bps_IV %*% y_beta_cond_jk[k,((i-1)*14+1):(i*14),b]
    #   y_sim_uncond_IV[(29*(i-1)+1):(29*i),k] <- Bp_IV %*% y_beta_uncond_jk[k,((i-1)*7+1):(i*7),b]
    }
  }
  cor_uncond_jk[,,b] <- sapply(1:12,function(s){
    sapply(1:29,function(t){
      cor(exp(p_hat_sim[t,] + p_sim_uncond[t,]),
          y_hat_sim[29*(s-1)+t,] + y_sim_uncond[29*(s-1)+t,])
    })
  }) %>% matrix(nrow=29)

  cor_cond_jk[,,b] <- sapply(1:12,function(s){
    sapply(1:29,function(t){
      cor(exp(p_hat_sim[t,] + p_sim_cond[t,]),
          y_hat_sim[29*(s-1)+t,] + y_sim_cond[29*(s-1)+t,])
    })
  }) %>% matrix(nrow=29)

  # prems_curr_phi0.7_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     prem_calc(p_fut_sim[t,] * exp(p_hat_sim[t,] + p_sim_IV[t,]),
  #               y_hat_sim[29*(s-1)+t,] + y_sim_uncond_IV[29*(s-1)+t,],
  #               0.7,
  #               p_fut_sim[t,],
  #               y_bar[s,t])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # prems_curr_phi0.9_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     prem_calc(p_fut_sim[t,] * exp(p_hat_sim[t,] + p_sim_IV[t,]),
  #               y_hat_sim[29*(s-1)+t,] + y_sim_uncond_IV[29*(s-1)+t,],
  #               0.9,
  #               p_fut_sim[t,],
  #               y_bar[s,t])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # prems_prop_phi0.7_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     prem_calc(p_fut_sim[t,] * exp(p_hat_sim[t,] + p_sim_IV[t,]),
  #               y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
  #               0.7,
  #               p_fut_sim[t,],
  #               y_bar[s,t])
  #   })
  # }) %>% matrix(nrow=29)
  # 
  # prems_prop_phi0.9_jk[,,b] <- sapply(1:12,function(s){
  #   sapply(1:29,function(t){
  #     prem_calc(p_fut_sim[t,] * exp(p_hat_sim[t,] + p_sim_IV[t,]),
  #               y_hat_sim[29*(s-1)+t,] + y_sim_cond_IV[29*(s-1)+t,],
  #               0.9,
  #               p_fut_sim[t,],
  #               y_bar[s,t])
  #   })
  # }) %>% matrix(nrow=29)
}

# dat <- data.frame(State = rep(states,each=29*2),
#                   s = rep(unique(s_tilde),each=2,times=12),
#                   Year = rep(1990:2018,each=2,times=12),
#                   Method = rep(c("Proposed","Current"),times=12*29),
#                   phi = as.character(rep(c(0.7),times=12*29*2)))
# 
# dat$prem <- apply(dat,1,function(x){
#   s <- which(states==x[1])
#   t <- as.numeric(x[3]) - 1989
#   if (x[4] == "Current"){
#     p_tilde <- qnorm(tau_p,0,IV[t])
#     y_tilde <- rep(0,K)
#     for (k in 1:K){
#       Bp <- bs(p_tilde[k],degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
#       y_tilde[k] <- Bp %*% y_beta_uncond[k,((s-1)*7+1):(s*7)]
#     }
#     p <- p_harv[t] * exp(p_tilde)
#     # p <- 4 * exp(p_tilde)
#     y <- y_adj_2020[s,1] + y_tilde
#   } else {
#     IV_tmp <- qnorm(tau_IV,IV_hat[t],sqrt(IV_sigma2))
#     p_harv_tmp <- qnorm(tau_phat,p_harv_hat[t],sqrt(p_harv_sigma2))
#     p_tilde <- qnorm(tau_p,0,IV_tmp)
#     # p_tilde <- qnorm(tau_p,0,IV[t])
#     y_tilde <- rep(0,K)
#     for (k in 1:K){
#       Bp <- bs(p_tilde[k],degree=degree,knots=p_kn,Boundary.knots = c(-1,1),intercept = TRUE)
#       Bps <- c(Bp[1,],Bs[t,])
#       y_tilde[k] <- t(Bps) %*% y_beta_cond[k,((s-1)*14+1):(s*14)]
#     }
#     p <- p_harv_tmp * exp(p_tilde)
#     # p <- 4 * exp(p_tilde)
#     y <- y_adj_2020[s,1] + y_tilde
#   }
#   prem_calc(p,y,as.numeric(x[5]),p_harv[t],y_bar[s,t])
# })

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
# prems_prop <- (dat %>%
#                  filter(phi == 0.7) %>%
#                  group_by(State,Year) %>%
#                  summarise(PremDiff=prem[Method=="Proposed"]) %>%
#                  ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
# 
# prems_curr <- (dat %>%
#                  filter(phi == 0.7) %>%
#                  group_by(State,Year) %>%
#                  summarise(PremDiff=prem[Method=="Current"]) %>%
#                  ungroup %>% select(PremDiff))$PremDiff %>% matrix(nrow=29)
#
# corn_cty_yield %>% 
#   group_by(Year) %>% 
#   summarise(AvgYield = mean(`Yield (Bu/Acre)`)) %>% 
#   rename(year=Year) %>%
#   left_join(corn_stocks) %>%
#   rename(Year=year) %>%
#   ggplot() +
#   geom_line(aes(x=year,y=AvgYield))+
#   geom_line(aes(x=year,y=s_tilde))+
#   scale_y_continuous()
#   theme_bw()
# ggsave("Corn_Yield_Plot.png")

p_sd_p <- data.frame(s=unique(s_tilde),
                p_cond=sd_p_cond) %>%
  ggplot()+
  # geom_line(aes(x=s,y=p_cond))+
  geom_point(aes(x=s,y=p_cond)) +
  geom_smooth(aes(x=s,y=p_cond),colour="black",se=FALSE) +
  xlab(expression(tilde(s)[t]))+
  ylab(expression(paste(hat(sd),"(",p["2020,t"],"|",tilde(s)[t],")")))+
  # geom_hline(aes(yintercept=sd(p_sim_uncond)),colour="red",linetype=2)+
  theme_bw()

p_IV <- data.frame(s=unique(s_tilde),
           # p_cond=apply(p_sim_cond,1,sd),
           IV = IV,
           hat = IV_hat) %>%
  #pivot_longer(cols=-s,names_to="Type",values_to="Value") %>%
  ggplot()+
  # geom_line(aes(x=s,y=p_cond))+
  geom_point(aes(x=s,y=IV)) +
  geom_smooth(aes(x=s,y=hat),colour="black",se=FALSE) +
  xlab(expression(tilde(s)[t]))+
  ylab(expression(paste(hat(sd),"(",p["2020,t"],"|",tilde(s)[t],")")))+
  # geom_hline(aes(yintercept=sd(p_sim_uncond)),colour="red",linetype=2)+
  theme_bw()
ggsave("Figures/Corn_IV_state.pdf",p_IV,width=7,height=4)

# p_sd_p_IV <- data.frame(s=unique(s_tilde),
#                      p_cond=sd_p_IV) %>%
#   ggplot()+
#   # geom_line(aes(x=s,y=p_cond))+
#   geom_point(aes(x=s,y=p_cond)) +
#   geom_smooth(aes(x=s,y=p_cond)) +
#   xlab(expression(tilde(s)[t]))+
#   ylab(expression(paste("sd(",p["2020,t"],"|",tilde(s)[t],")")))+
#   # geom_hline(aes(yintercept=sd(p_sim_uncond)),colour="red",linetype=2)+
#   theme_bw()
# p_sd <- grid.arrange(p_sd_p,p_sd_p_IV,ncol=2)
ggsave("Figures/Corn_sd_state_NonIrr.pdf",p_sd_p,width=7,height=4)
  
group.colors <- c(Conditional = "#000000",
                  Unconditional = "red")#Unconditional = "#4aa02c")

RMA_cor <- c(-.4,-.35,-.4,
             -.25,0,-.3,
             -.35,-.35,0,
             -.25,0,-.2)

p_cor_yp <- data.frame(state=rep(states,each=29),
           s=rep(unique(s_tilde),times=12),
           Conditional=as.vector(cor_cond),
           Unconditional=as.vector(cor_uncond),
           low=as.vector(cor_cond - 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B)),
           upp=as.vector(cor_cond + 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B))) %>%
  gather("Type","Value",-c(s,state,low,upp)) %>%
  # gather("Type","Value",-c(s,state)) %>%
  ggplot()+
  #geom_line(aes(x=s,y=Value,linetype=Type,colour=Type))+
  # geom_point(aes(x=s,y=Value,colour=Type))+
  geom_smooth(aes(x=s,y=Value,colour=Type,linetype=Type),se=FALSE)+
  #geom_point(aes(x=s,y=Conditional))+
  scale_colour_manual(values=group.colors)+
  geom_smooth(aes(x=s,y=low),colour="black",se=FALSE,linetype="dotted")+
  geom_smooth(aes(x=s,y=upp),colour="black",se=FALSE,linetype="dotted")+
  #geom_hline(yintercept=0,colour="red")+
  facet_wrap(.~state,nrow=4,scales="free")+
  xlab(expression(tilde(s)[t]))+
  ylab(expression(paste(hat(cor),"(",y["2020,jt"],",",p["2020,t"],"|",tilde(s)[t],")")))+
  theme_bw() +
  theme(legend.position = "none")

# p_cor_yp_IV <- data.frame(state=rep(states,each=29),
#                     s=rep(unique(s_tilde),times=12),
#                     Conditional=as.vector(cor_cond_IV),
#                     Unconditional=as.vector(cor_uncond_IV)) %>%
#   #low=as.vector(cor_cond - 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B)),
#   #upp=as.vector(cor_cond + 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B))) %>%
#   # gather("Type","Value",-c(s,state,low,upp)) %>%
#   gather("Type","Value",-c(s,state)) %>%
#   ggplot()+
#   #geom_line(aes(x=s,y=Value,linetype=Type,colour=Type))+
#   geom_point(aes(x=s,y=Value,colour=Type))+
#   geom_smooth(aes(x=s,y=Value,colour=Type))+
#   #geom_point(aes(x=s,y=Conditional))+
#   scale_colour_manual(values=group.colors)+
#   # geom_point(aes(x=s,y=low),colour="blue")+
#   # geom_point(aes(x=s,y=upp),colour="blue")+
#   #geom_hline(yintercept=0,colour="red")+
#   facet_wrap(.~state,nrow=4,scales="free")+
#   xlab(expression(tilde(s)[t]))+
#   ylab(expression(paste(hat(cor),"(",y["2020,jt"],",",p["2020,t"],"|",tilde(s)[t],")")))+
#   theme_bw() +
#   theme(legend.position = "none")

# p_cor <- grid.arrange(p_cor_yp,p_cor_yp_IV,ncol=2)
ggsave("Figures/Corn_cor_state_NonIrr.pdf",p_cor_yp,width=8,height=6)

p_cor_yp <- data.frame(state=rep(states,each=29),
                       s=rep(unique(s_tilde),times=12),
                       Conditional=as.vector(cor_cond),
                       # Unconditional=as.vector(cor_uncond),
                       Unconditional = rep(RMA_cor,each=29),
                       low=as.vector(cor_cond - 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B)),
                       upp=as.vector(cor_cond + 1.96*sqrt(apply(cor_cond_jk,c(1,2),var)*(B-1)^2/B))) %>%
  gather("Type","Value",-c(s,state,low,upp)) %>%
  # gather("Type","Value",-c(s,state)) %>%
  ggplot()+
  #geom_line(aes(x=s,y=Value,linetype=Type,colour=Type))+
  # geom_point(aes(x=s,y=Value,colour=Type))+
  geom_smooth(aes(x=s,y=Value,colour=Type,linetype=Type),se=FALSE)+
  #geom_point(aes(x=s,y=Conditional))+
  scale_colour_manual(values=group.colors)+
  geom_smooth(aes(x=s,y=low),colour="black",se=FALSE,linetype="dotted")+
  geom_smooth(aes(x=s,y=upp),colour="black",se=FALSE,linetype="dotted")+
  #geom_hline(yintercept=0,colour="red")+
  facet_wrap(.~state,nrow=4,scales="free")+
  xlab(expression(tilde(s)[t]))+
  ylab(expression(paste(hat(cor),"(",y["2020,jt"],",",p["2020,t"],"|",tilde(s)[t],")")))+
  theme_bw() +
  theme(legend.position = "none")
ggsave("Figures/Corn_cor_state_NonIrr_Poster.pdf",p_cor_yp,width=8,height=7)

# dt <- corn_price %>%
#   rename(Year=year) %>%
#   full_join(corn_cty_yield) %>%
#   mutate(y = y_tilde + as.vector(Bc %*% y_adj_2020)) %>%
#   mutate(state = toTitleCase(tolower(State))) %>%
#   dplyr::select(Year,state,p,y)
# 
# uncond_pred <- NULL
# cond_pred <- NULL
#
# for (s in 1:12){
#   for (t in 1:29){
#     state <- states[s]
#     dens_uncond <- kde2d(exp(p_hat$fitted[t]+p_sim_uncond[t,]),y_sim_uncond[29*(s-1)+t,] + y_adj_2020[s,1],lims=c(2.4,8.5,22,300))
#     dens_cond <- kde2d(exp(p_hat$fitted[t]+p_sim_cond[t,]),y_sim_cond[29*(s-1)+t,] + y_adj_2020[s,1],lims=c(2.4,8.5,22,300))
#     uncond_curr <- interp(x=rep(dens_uncond$x,times=25),
#                           y=rep(dens_uncond$y,each=25),
#                           z=as.vector(dens_uncond$z),
#                           xo=dt$p[dt$Year==1989+t & dt$state==state] %>% unique(),
#                           yo=dt$y[dt$Year==1989+t & dt$state==state])$z %>% as.vector() %>% log %>% sum
#     uncond_pred <- c(uncond_pred,uncond_curr)
#     cond_curr <- interp(x=rep(dens_cond$x,times=25),
#                         y=rep(dens_cond$y,each=25),
#                         z=as.vector(dens_cond$z),
#                         xo=dt$p[dt$Year==1989+t & dt$state==state] %>% unique(),
#                         yo=dt$y[dt$Year==1989+t & dt$state==state],
#                         linear=TRUE)$z %>% as.vector() %>% log %>% sum
#     cond_pred <- c(cond_pred,cond_curr)
#   }
# }

# Posterior predictive of unconditional distribution -154162.7
# Posterior predictive of conditional distribution -146503.8

appender <- function(string) {
  TeX(ifelse(string==1,"RMA",paste("$\\tilde{s}_t = $", string)))
}

dt <- data.frame(state=rep(states,4),
                 s_tilde=c(round(rep(unique(s_tilde)[order(unique(s_tilde))[round(29*c(.2,.5,.8))]],each=12),2),rep(1,12)),
                 values=c(as.vector(t(cor_cond[order(unique(s_tilde))[round(29*c(.2,.5,.8))],])),RMA_cor))
p2 <- plot_usmap(data=dt,include=states)+
  scale_fill_gradient(
    limits = c(min(dt$values),max(c(dt$values,0))),
    low="red",high="white",
    name = expression(paste(hat(cor)))) +
  facet_wrap(.~s_tilde,labeller=as_labeller(appender,
                                            default = label_parsed),
             nrow=2) +
  theme(legend.position = c(1,0.3))
ggsave("Figures/Corn_cor_statemap_1990_NonIrr.pdf",p2,scale=2,width=3,height=2)
ggsave("Figures/Corn_cor_statemap_NonIrr_Poster.pdf",p2,width=11,height=8)

# t <- which.min(s_tilde)
# p5 <- data.frame(State = rep(c("Michigan","Minnesota","Wisconsin"),each=K),
#            Conditional = c(y_adj_2020[5,1] + y_sim_cond[29*(5-1)+t,],
#                            y_adj_2020[6,1] + y_sim_cond[29*(6-1)+t,],
#                            y_adj_2020[12,1] + y_sim_cond[29*(12-1)+t,]),
#            Unconditional = c(y_adj_2020[5,1] + y_sim_uncond[29*(5-1)+t,],
#                              y_adj_2020[6,1] + y_sim_uncond[29*(6-1)+t,],
#                              y_adj_2020[12,1] + y_sim_uncond[29*(12-1)+t,])) %>%
#   gather(key="Type",value="Value",-State) %>%
#   ggplot(aes(x=Value,fill=Type)) +
#   geom_density(alpha=0.5) +
#   facet_wrap(.~State,nrow=1,scales="free") +
#   xlab(expression(y["2020,jt"])) +
#   ylab(expression(paste("g(",y["2020,jt"],")"))) +
#   theme_bw()
# ggsave("Figures/Corn_state_y_dist.pdf",p5,width=8,height=3)

# p5 <- dat %>% # data.frame(state=rep(states,each=29),
#                  s=rep(unique(s_tilde),times=12),
#                  # Value=as.vector(prems_prop),
#                  Proposed=as.vector(prems_prop),
#                  Current=as.vector(prems_curr)) %>% #,
#                  # low=as.vector(prems_prop - 1.96*sqrt(apply(prems_prop_jk,c(1,2),var)*(B-1)^2/B)),
#                  # upp=as.vector(prems_prop + 1.96*sqrt(apply(prems_prop_jk,c(1,2),var)*(B-1)^2/B))) %>%
  #gather("Type","Value",-c(s,state,low,upp)) %>%
#  pivot_longer(cols=c(Proposed,Current),names_to="Type",values_to="Value") %>%
#   ggplot() +
#   #geom_line(aes(x=s,y=prem,colour=Method)) +
#   #geom_line(aes(x=s,y=Value,linetype=Type,colour=Type))+
#   #scale_colour_manual(values=group.colors)+
#   geom_line(aes(x=s,y=prem))+
#   # geom_line(aes(x=s,y=low),linetype=3)+
#   # geom_line(aes(x=s,y=upp),linetype=3)+
#   #annotate(geom="rect",xmin=0.06,xmax=0.09,ymin=-Inf,ymax=Inf,alpha=0.2) +
#   #annotate(geom="rect",xmin=0.135,xmax=0.15,ymin=-Inf,ymax=Inf,alpha=0.2) +
#   geom_hline(aes(yintercept=0),colour="red",linetype=2) +
#   #ylim(c(0,20)) +
#   facet_wrap(.~State,nrow=4,scales="free")+
#   theme_bw()+
#   xlab(expression(tilde(s)[t]))+
#   ylab("Premium ($ per acre)")
#   #ylab("Premium Difference($ per acre)")
# ggsave("Figures/Corn_prem_state_1990_phi0.7.pdf",p5,width=8,height=5)


prem_data_0.7 <- data.frame(State=rep(unique(corn_cty_yield$State),each=29),
           s=rep(unique(s_tilde),times=12),
           Year=rep(1990:2018,times=12),
           p_fut=rep(p_fut,times=12),
           p_harv=rep(p_harv,times=12),
           y_bar=as.vector(t(y_bar)),
           y_hat=rep(y_hat_2020,each=29),
           prems_3chan = as.vector(t(prems_prop_phi0.7)),
           prems_2chan = as.vector(t(prems_curr_phi0.7))) %>%
  right_join(corn_cty_yield[,c("State","County","Year","y_tilde")],by=c("State","Year")) %>% 
  mutate(y = y_hat + y_tilde) %>%
  mutate(LR_num = pmax(0.7*p_fut*y_bar - p_harv*y,0))

prem_0.7 <- prem_data_0.7 %>%
  group_by(State,Year) %>%
  summarise(sum_LR_num = sum(LR_num),
            sum_prems_2 = sum(prems_2chan),
            sum_prems_3 = sum(prems_3chan)) %>%
  ungroup %>%
  mutate(LR_3_cede = ifelse(sum_LR_num > sum_prems_3,sum_LR_num/sum_prems_3,1),
         LR_3_retain = ifelse(sum_LR_num > sum_prems_3,1,sum_LR_num/sum_prems_3),
         LR_2_cede = ifelse(sum_LR_num > sum_prems_2,sum_LR_num/sum_prems_2,1),
         LR_2_retain = ifelse(sum_LR_num > sum_prems_2,1,sum_LR_num/sum_prems_2)) %>%
  mutate(D = ifelse(is.nan((LR_2_cede/LR_2_retain)/(LR_3_cede/LR_3_retain)),
                    1,
                    (LR_2_cede/LR_2_retain)/(LR_3_cede/LR_3_retain))) %>% 
  group_by(State) %>%
  summarise(D_star = length(which(D > 1))) %>%
  ungroup %>%
  mutate(p_val = pbinom(D_star,29,0.5,lower.tail=FALSE))

prem_data_0.85 <- data.frame(State=rep(unique(corn_cty_yield$State),each=29),
                        s=rep(unique(s_tilde),times=12),
                        Year=rep(1990:2018,times=12),
                        p_fut=rep(p_fut,times=12),
                        p_harv=rep(p_harv,times=12),
                        y_bar=as.vector(t(y_bar)),
                        y_hat=rep(y_hat_2020,each=29),
                        prems_3chan = as.vector(t(prems_prop_phi0.85)),
                        prems_2chan = as.vector(t(prems_curr_phi0.85))) %>%
  right_join(corn_cty_yield[,c("State","County","Year","y_tilde")],by=c("State","Year")) %>% 
  mutate(y = y_hat + y_tilde) %>%
  mutate(LR_num = pmax(0.85*p_fut*y_bar - p_harv*y,0))

prem_0.85 <- prem_data_0.85 %>%
  group_by(State,Year) %>%
  summarise(sum_LR_num = sum(LR_num),
            sum_prems_2 = sum(prems_2chan),
            sum_prems_3 = sum(prems_3chan)) %>%
  ungroup %>%
  mutate(LR_3_cede = ifelse(sum_LR_num > sum_prems_3,sum_LR_num/sum_prems_3,1),
         LR_3_retain = ifelse(sum_LR_num > sum_prems_3,1,sum_LR_num/sum_prems_3),
         LR_2_cede = ifelse(sum_LR_num > sum_prems_2,sum_LR_num/sum_prems_2,1),
         LR_2_retain = ifelse(sum_LR_num > sum_prems_2,1,sum_LR_num/sum_prems_2)) %>%
  mutate(D = ifelse(is.nan((LR_2_cede/LR_2_retain)/(LR_3_cede/LR_3_retain)),
                    1,
                    (LR_2_cede/LR_2_retain)/(LR_3_cede/LR_3_retain))) %>% 
  group_by(State) %>%
  summarise(D_star = length(which(D > 1))) %>%
  ungroup %>%
  mutate(p_val = pbinom(D_star,29,0.5,lower.tail=FALSE))

  left_join(prem_0.7,prem_0.85,by="State") %>%
    select(-contains("D_star")) %>%
    mutate(State = states,
           p_val.x = ifelse(p_val.x == 1.00,0.9999,p_val.x),
           p_val.y = ifelse(p_val.y == 1.00,0.9999,p_val.y)) %>%
    rename('70% Coverage' = p_val.x,
           '85% Coverage' = p_val.y) %>%
    select(State,'70% Coverage','85% Coverage') %>%
  xtable::xtable(digits=4) %>%
  print(include.rownames = FALSE)
  
  prem_data_0.7 <- data.frame(State=rep(unique(corn_cty_yield$State),each=29),
                              s=rep(unique(s_tilde),times=12),
                              Year=rep(1990:2018,times=12),
                              p_fut=rep(p_fut,times=12),
                              p_harv=rep(p_harv,times=12),
                              y_bar=as.vector(t(y_bar)),
                              y_hat=rep(y_hat_2020,each=29),
                              prems_3chan = as.vector(t(prems_prop_phi0.7)),
                              prems_2chan = as.vector(t(prems_curr_phi0.7))) %>%
    right_join(corn_cty_yield[,c("State","County","Year","y_tilde")],by=c("State","Year")) %>% 
    mutate(y = y_hat + y_tilde) %>%
    mutate(LR_num = pmax(0.7*p_fut*y_bar - p_harv*y,0))
  
  prem_0.7 <- prem_data_0.7 %>%
    group_by(State,Year) %>%
    summarise(sum_LR_num = sum(LR_num),
              sum_prems_2 = sum(prems_2chan),
              sum_prems_3 = sum(prems_3chan)) %>%
    ungroup %>%
    mutate(LR_3_cede = ifelse(sum_LR_num > sum_prems_3,sum_LR_num/sum_prems_3,1),
           LR_3_retain = ifelse(sum_LR_num > sum_prems_3,1,sum_LR_num/sum_prems_3),
           LR_2_cede = ifelse(sum_LR_num > sum_prems_2,sum_LR_num/sum_prems_2,1),
           LR_2_retain = ifelse(sum_LR_num > sum_prems_2,1,sum_LR_num/sum_prems_2)) %>%
    mutate(D = ifelse(is.nan((LR_3_cede/LR_3_retain)/(LR_2_cede/LR_2_retain)),
                      1,
                      (LR_3_cede/LR_3_retain)/(LR_2_cede/LR_2_retain))) %>% 
    group_by(State) %>%
    summarise(D_star = length(which(D > 1))) %>%
    ungroup %>%
    mutate(p_val = pbinom(D_star,29,0.5,lower.tail=FALSE))
  
  prem_data_0.85 <- data.frame(State=rep(unique(corn_cty_yield$State),each=29),
                               s=rep(unique(s_tilde),times=12),
                               Year=rep(1990:2018,times=12),
                               p_fut=rep(p_fut,times=12),
                               p_harv=rep(p_harv,times=12),
                               y_bar=as.vector(t(y_bar)),
                               y_hat=rep(y_hat_2020,each=29),
                               prems_3chan = as.vector(t(prems_prop_phi0.85)),
                               prems_2chan = as.vector(t(prems_curr_phi0.85))) %>%
    right_join(corn_cty_yield[,c("State","County","Year","y_tilde")],by=c("State","Year")) %>% 
    mutate(y = y_hat + y_tilde) %>%
    mutate(LR_num = pmax(0.85*p_fut*y_bar - p_harv*y,0))
  
  prem_0.85 <- prem_data_0.85 %>%
    group_by(State,Year) %>%
    summarise(sum_LR_num = sum(LR_num),
              sum_prems_2 = sum(prems_2chan),
              sum_prems_3 = sum(prems_3chan)) %>%
    ungroup %>%
    mutate(LR_3_cede = ifelse(sum_LR_num > sum_prems_3,sum_LR_num/sum_prems_3,1),
           LR_3_retain = ifelse(sum_LR_num > sum_prems_3,1,sum_LR_num/sum_prems_3),
           LR_2_cede = ifelse(sum_LR_num > sum_prems_2,sum_LR_num/sum_prems_2,1),
           LR_2_retain = ifelse(sum_LR_num > sum_prems_2,1,sum_LR_num/sum_prems_2)) %>%
    mutate(D = ifelse(is.nan((LR_3_cede/LR_3_retain)/(LR_2_cede/LR_2_retain)),
                      1,
                      (LR_3_cede/LR_3_retain)/(LR_2_cede/LR_2_retain))) %>% 
    group_by(State) %>%
    summarise(D_star = length(which(D > 1))) %>%
    ungroup %>%
    mutate(p_val = pbinom(D_star,29,0.5,lower.tail=FALSE))
  
  left_join(prem_0.7,prem_0.85,by="State") %>%
    select(-contains("D_star")) %>%
    mutate(State = states,
           p_val.x = ifelse(p_val.x == 1.00,0.9999,p_val.x),
           p_val.y = ifelse(p_val.y == 1.00,0.9999,p_val.y)) %>%
    rename('70% Coverage' = p_val.x,
           '85% Coverage' = p_val.y) %>%
    select(State,'70% Coverage','85% Coverage') %>%
    xtable::xtable(digits=4) %>%
    print(include.rownames = FALSE)