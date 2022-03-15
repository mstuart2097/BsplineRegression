library(readxl)
library(splines)
library(tidyverse)
library(MASS)
source("BSplineFunctions.R")

corn_price      <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_futures_price_Dec")[-1,]
corn_stocks     <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_ending_stock")
corn_natl_yield <- read_excel("data/futures_price_yield_200308.xlsx",sheet = "corn_production_nationwide")
corn_cty_yield  <- read_excel("data/futures_price_yield_200206.xlsx",sheet = "corn_yield_bu_acre")

corn_price      <- corn_price  %>% filter(year >= 1973 & year <= 2018)
corn_stocks     <- corn_stocks %>% filter(year >= 1973 & year <= 2018)
corn_natl_yield <- corn_natl_yield %>% filter(Year >= 1973 & Year <= 2018)
corn_cty_yield  <- corn_cty_yield  %>% filter(Year >= 1973 & Year <= 2018)

t_tilde <- corn_cty_yield$Year - 1972

n_ctys_by_year <- sapply(corn_cty_yield$Year 
                         %>% unique,
                         function(t){length(which(corn_cty_yield$Year == t))})

p_tilde <- rep(log(as.numeric(corn_price$`Dec_price (original)`)
               /as.numeric(corn_price$`Feb_price (original)`)),
           times = n_ctys_by_year)

p_pri <- rep(as.numeric(corn_price$`Feb_price (original)`),
             times = n_ctys_by_year)

s_hat <- loess(`Production Measured in Million Bushels`~Year,
               corn_natl_yield)

plot(corn_natl_yield$Year,corn_natl_yield$`Production Measured in Million Bushels`)
lines(predict(s_hat),x=corn_natl_yield$Year)

s <- rep(corn_stocks$`Ending Stocks (Million Bushels)`
         /s_hat$fitted,
         times = n_ctys_by_year) %>% as.vector

Bs <- bs(t_tilde,knots = c(10,20,30,40))
lambda <- 3
Ds <- matrix(0,nrow=ncol(Bs),ncol=ncol(Bs))

for (j in 1:nrow(Ds)){
  for (k in 1:ncol(Ds)){
    if (k-j == 0 || k-j == 2){
      Ds[j,k] = 1
    } else if(k-j == 1){
      Ds[j,k] = -2
    }
  }
}

y <- as.vector(corn_cty_yield$`Yield (Bu/Acre)`)

y_hat <- Bs %*% mean.reg(y,Bs,Ds,lambda)
y_tilde <- as.vector(corn_cty_yield$`Yield (Bu/Acre)` - y_hat)

