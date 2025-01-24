rm(list=ls())
library(tidyverse)
library(grid)
library(splines)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(MASS)
library(ggpubr)
library(sn)
library(quantreg)
library(directlabels)
library(latex2exp)

dt <- data.frame(x = seq(-3,3,0.01))
dt$`0.01` <- dt$x*(0.01 - (dt$x < 0))
dt$`0.15` <- dt$x*(0.15 - (dt$x < 0))
dt$`0.29` <- dt$x*(0.29 - (dt$x < 0))
dt$`0.43` <- dt$x*(0.43 - (dt$x < 0))
dt$`0.57` <- dt$x*(0.57 - (dt$x < 0))
dt$`0.71` <- dt$x*(0.71 - (dt$x < 0))
dt$`0.85` <- dt$x*(0.85 - (dt$x < 0))
dt$`0.99` <- dt$x*(0.99 - (dt$x < 0))
dt %>% 
  gather("tau","value",-x) %>%
  ggplot() +
  geom_line(aes(x=x,y=value,linetype=tau,color=tau)) +
  xlab(TeX("u")) +
  ylab(TeX("$\\rho_{\\tau}(u)$"))+
  #ylab(expression(paste(rho[tau](u)))) + 
  labs(color=TeX("$\\tau$"),linetype=TeX("$\\tau$")) +
  theme_bw()
ggsave("check.pdf", width = 5, height = 4)

set.seed(8675309)
x <- runif(100)
y <- rnorm(100,2*x + sin(2*pi*x),x)
#y <- rnorm(1000,x,x)
Bx <- bs(x,knots=quantile(x,c(0.25,0.5,0.75)),degree=3,
         Boundary.knots=c(0,1),intercept = TRUE)
dt <- as.data.frame(cbind(y,x,Bx))
names(dt) <- c("y","x","Bx_1","Bx_2","Bx_3","Bx_4","Bx_5","Bx_6","Bx_7")
dt$`5th Quantile` <- rq(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt,tau=0.05)$fitted.values
dt$`25th Quantile` <- rq(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt,tau=0.25)$fitted.values
dt$`50th Quantile` <- rq(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt,tau=0.5)$fitted.values
dt$`75th Quantile` <- rq(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt,tau=0.75)$fitted.values
dt$`95th Quantile` <- rq(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt,tau=0.95)$fitted.values
dt$Mean <- lm(y ~ 0+Bx_1+Bx_2+Bx_3+Bx_4+Bx_5+Bx_6+Bx_7,dt)$fitted.values
dt %>% dplyr::select(!contains("Bx")) %>%
  gather("label","value",-c(y,x)) %>%
  mutate(label=factor(label,levels=c("5th Quantile","25th Quantile","50th Quantile","75th Quantile","95th Quantile","Mean"))) %>%
  mutate(regtype=ifelse(label=="Mean","Mean","Quantile")) %>%
  mutate(line_lab = if_else(x==max(x), as.character(label), NA_character_)) %>%
  filter(x < 0.95) %>%
  ggplot()+
  geom_point(aes(x=x,y=y),size=0.5) +
  geom_line(aes(x=x,y=value,group=label,color=regtype))+
  scale_color_manual(values=c('Red','Black'))+
  xlim(c(0,1.25))+
  geom_dl(aes(x=x,y=value,label = label), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  labs(color="")+
  theme_bw()+
  theme(legend.position="none")
ggsave("quantreg.pdf", width = 5, height = 4)

rm(list=ls())
m_p_lin <- function(s){0.2 - 0.4*s}
s_p_lin <- function(s){0.5 - 0.5*s}
m_p_non <- function(s){-0.2+0.4*exp(-2*s)}
s_p_non <- function(s){0.5*exp(-2*s)}
# 
a_p <- 3 # Skewness parameter because detrended prices are empirically right skewed
d_p <- a_p / sqrt(1 + a_p^2)  #To ensure eps_p terms have theoretical zero mean and unit variance

m_y_non <- function(p,s){-25 + 14.45*exp(p) + 22.18*s}
#m_y_non <- function(p,s){-5.31 + 4.67*p - 9.01*s + 136.04*p^2 + 52.56*s^2}
s_y <- 33
a_y <- -3 # Skewness parameter because detrended prices are empirically left skewed
d_y <- a_y / sqrt(1 + a_y^2) #Used to ensure eps_y terms have theoretical zero mean and unit variance

#Used to ensure eps_y terms have theoretical zero mean and unit variance
Q <- 4
s_tru_quants <- qbeta(c(0.1,seq(1/Q,(Q-1)/Q,1/Q),0.9),7,44)

K <- 5000
set.seed(10062018) # For reproducability
taus <- c(0.1,0.25,0.5,0.75,0.9)# Fixed values of tau used for the quantiles
tau_p <- runif(K) # Fixed values of tau_p used in the singular simulation
tau_y <- runif(K) # Fixed values of tau_y used in the singular simulation
s_sim <- rbeta(100,7,44)

# p_app_lin <- readRDS("p_app_lin.RDS")
# p_app_non <- readRDS("p_app_non.RDS")
# y_app_lin_p_lin <- readRDS("y_app_lin_p_lin.RDS")
# y_app_lin_p_non <- readRDS("y_app_lin_p_non.RDS")
# y_app_non_p_lin <- readRDS("y_app_non_p_lin.RDS")
# y_app_non_p_non <- readRDS("y_app_non_p_non.RDS")

appender <- function(string) {
  TeX(paste("$\\tau_p = $", string)) 
}
# 
# apprx = p_app_lin
# dt <- data.frame(s=seq(0,0.5,0.01),
#                  med=apply(apprx,c(1,2),quantile,probs=0.5),
#                  low=apply(apprx,c(1,2),quantile,probs=0.025),
#                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>%
#   mutate(`tru.1`=m_p_lin(s)+s_p_lin(s)*qsn(0.1,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.2`=m_p_lin(s)+s_p_lin(s)*qsn(0.25,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.3`=m_p_lin(s)+s_p_lin(s)*qsn(0.5,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.4`=m_p_lin(s)+s_p_lin(s)*qsn(0.75,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.5`=m_p_lin(s)+s_p_lin(s)*qsn(0.9,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)) %>%
#   gather(key="var",value="value",-s)
# dt <- dt %>%
#   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
#   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
#   dplyr::select(-var)
# dt <- dt %>%
#   mutate(tau=replace(tau,tau==1,0.1),
#          tau=replace(tau,tau==2,0.25),
#          tau=replace(tau,tau==3,0.5),
#          tau=replace(tau,tau==4,0.75),
#          tau=replace(tau,tau==5,0.9))
# dt <- dt %>%
#   spread(key="Type",value="value")
# p <- ggplot(dt %>% filter(s >= 0.05 & s <= 0.2)) +
#   geom_line(aes(x=s,y=med),color="red") +
#   geom_line(aes(x=s,y=low),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=upp),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=tru)) +
#   facet_wrap(tau~.,scales="free_y",nrow=2,labeller=as_labeller(appender,
#                                                                default = label_parsed)) +
#   xlab(expression(tilde(s))) +
#   ylab(expression(paste(q[tau[p]],"(",tilde(s),")"))) +
#   scale_x_continuous(breaks=c(0.075,0.175))+
#   theme_bw()
# ggsave("p_lin.pdf",p, width = 6, height = 4)
# 
# apprx = p_app_non
# dt <- data.frame(s=seq(0,0.5,0.01),
#                  med=apply(apprx,c(1,2),quantile,probs=0.5),
#                  low=apply(apprx,c(1,2),quantile,probs=0.025),
#                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>%
#   mutate(`tru.1`=m_p_non(s)+s_p_non(s)*qsn(0.1,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.2`=m_p_non(s)+s_p_non(s)*qsn(0.25,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.3`=m_p_non(s)+s_p_non(s)*qsn(0.5,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.4`=m_p_non(s)+s_p_non(s)*qsn(0.75,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p),
#          `tru.5`=m_p_non(s)+s_p_non(s)*qsn(0.9,xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),omega=sqrt((1-2*d_p^2/pi)^-1),alpha=a_p)) %>%
#   gather(key="var",value="value",-s)
# dt <- dt %>%
#   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
#   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
#   dplyr::select(-var)
# dt <- dt %>%
#   mutate(tau=replace(tau,tau==1,0.1),
#          tau=replace(tau,tau==2,0.25),
#          tau=replace(tau,tau==3,0.5),
#          tau=replace(tau,tau==4,0.75),
#          tau=replace(tau,tau==5,0.9))
# dt <- dt %>%
#   spread(key="Type",value="value")
# p <- ggplot(dt %>% filter(s >= 0.05 & s <= 0.2)) +
#   geom_line(aes(x=s,y=med),color="red") +
#   geom_line(aes(x=s,y=low),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=upp),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=tru)) +
#   facet_wrap(tau~.,scales="free_y",nrow=2,labeller=as_labeller(appender, 
#                                                                default = label_parsed)) + 
#   xlab(expression(tilde(s))) +
#   ylab(expression(paste(q[tau[p]],"(",tilde(s),")"))) +
#   scale_x_continuous(breaks=c(0.075,0.175))+
#   theme_bw()
# ggsave("p_non.pdf",p, width = 6, height = 4)
# 
# appender2 <- function(string){
#   if(string == "0.1" || string == "0.25" || string == "0.5" || string == "0.75" || string == "0.9")  {
#     TeX(paste("$\\tau_y = $",string))
#   } else {
#     TeX(paste("$\\tilde{s} = $",string))
#   }
# }
# # 
# # apprx = y_app_lin_p_lin
# # dt <- data.frame(p=rep(seq(-3,3,0.06),5),
# #                  s=rep(round(s_tru_quants,3),each=101),
# #                  med=apply(apprx,c(1,2),quantile,probs=0.5),
# #                  low=apply(apprx,c(1,2),quantile,probs=0.025),
# #                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>% 
# #   mutate(`tru.1`=m_y_lin(p,s)+s_y*qsn(0.1,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.2`=m_y_lin(p,s)+s_y*qsn(0.25,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.3`=m_y_lin(p,s)+s_y*qsn(0.5,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.4`=m_y_lin(p,s)+s_y*qsn(0.75,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.5`=m_y_lin(p,s)+s_y*qsn(0.9,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)) %>%
# #   gather(key="var",value="value",-c("p","s"))
# # dt <- dt %>%
# #   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
# #   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
# #   dplyr::select(-var)
# # dt <- dt %>%
# #   mutate(tau=replace(tau,tau==1,0.1),
# #          tau=replace(tau,tau==2,0.25),
# #          tau=replace(tau,tau==3,0.5),
# #          tau=replace(tau,tau==4,0.75),
# #          tau=replace(tau,tau==5,0.9))
# # dt <- dt %>%
# #   spread(key="Type",value="value")
# # 
# # p <- ggplot(dt %>% filter(p >= -0.5 & p <= 1)) +
# #   geom_line(aes(x=p,y=med),color="red") +
# #   geom_line(aes(x=p,y=low),color="red",linetype="dashed") +
# #   geom_line(aes(x=p,y=upp),color="red",linetype="dashed") +
# #   geom_line(aes(x=p,y=tru)) +
# #   facet_grid(tau~s,scales="free_y",labeller=as_labeller(appender2, 
# #                                                         default = label_parsed)) + 
# #   xlab(expression(tilde(p))) +
# #   ylab(expression(paste(q[tau[y]],"(",tilde(p),",",tilde(s),")"))) +
# #   theme_bw()
# # ggsave("y_lin_p_lin.pdf",p, width = 6, height = 4)
# #
# # p <- ggplot(dt %>% filter(p >= -0.5 & p <= 1 & s == median(s))) +
# #   geom_line(aes(x=p,y=med-tru),color="red") +
# #   geom_line(aes(x=p,y=low-tru),color="red",linetype="dashed") +
# #   geom_line(aes(x=p,y=upp-tru),color="red",linetype="dashed") +
# #   geom_hline(yintercept = 0)+
# #   facet_wrap(tau~.,scales="free_y",nrow=2,labeller=as_labeller(appender2, 
# #                                                                default = label_parsed))+
# #   xlab(expression(tilde(p))) +
# #   ylab(expression(paste(q[tau[y]],"(",tilde(p),",",tilde(s),")"))) +
# #   theme_bw()
# # ggsave("y_lin_p_lin_zoom.pdf",p, width = 6, height = 4)
# # 
# # apprx = y_app_lin_p_non
# # dt <- data.frame(p=rep(seq(-0.5,1,0.02),5),
# #                  s=rep(round(s_tru_quants,3),each=76),
# #                  med=apply(apprx,c(1,2),quantile,probs=0.5),
# #                  low=apply(apprx,c(1,2),quantile,probs=0.025),
# #                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>% 
# #   mutate(`tru.1`=m_y_lin(p,s)+s_y*qsn(0.1,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.2`=m_y_lin(p,s)+s_y*qsn(0.25,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.3`=m_y_lin(p,s)+s_y*qsn(0.5,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.4`=m_y_lin(p,s)+s_y*qsn(0.75,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
# #          `tru.5`=m_y_lin(p,s)+s_y*qsn(0.9,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)) %>%
# #   gather(key="var",value="value",-c("p","s"))
# # dt <- dt %>%
# #   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
# #   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
# #   dplyr::select(-var)
# # dt <- dt %>%
# #   mutate(tau=replace(tau,tau==1,0.1),
# #          tau=replace(tau,tau==2,0.25),
# #          tau=replace(tau,tau==3,0.5),
# #          tau=replace(tau,tau==4,0.75),
# #          tau=replace(tau,tau==5,0.9))
# # dt <- dt %>%
# #   spread(key="Type",value="value")
# # p <- ggplot(dt) +
# #   geom_line(aes(x=p,y=med),color="red") +
# #   geom_line(aes(x=p,y=low),color="red",linetype="dashed") +
# #   geom_line(aes(x=p,y=upp),color="red",linetype="dashed") +
# #   geom_line(aes(x=p,y=tru)) +
# #   facet_grid(tau~s,scales="free_y",labeller=as_labeller(appender2, 
# #                                                         default = label_parsed)) + 
# #   xlab(expression(tilde(p))) +
# #   ylab(expression(paste(q[tau[y]],"(",tilde(p),",",tilde(s),")"))) +
# #   theme_bw()
# # ggsave("y_lin_p_non.pdf",p, width = 6, height = 4)
# 
# apprx = y_app_non_p_lin
# dt <- data.frame(p=rep(seq(-0.5,1,0.02),5),
#                  s=rep(round(s_tru_quants,3),each=76),
#                  med=apply(apprx,c(1,2),quantile,probs=0.5),
#                  low=apply(apprx,c(1,2),quantile,probs=0.025),
#                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>%
#   mutate(`tru.1`=m_y_non(p,s)+s_y*qsn(0.1,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.2`=m_y_non(p,s)+s_y*qsn(0.25,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.3`=m_y_non(p,s)+s_y*qsn(0.5,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.4`=m_y_non(p,s)+s_y*qsn(0.75,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.5`=m_y_non(p,s)+s_y*qsn(0.9,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)) %>%
#   gather(key="var",value="value",-c("p","s"))
# dt <- dt %>%
#   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
#   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
#   dplyr::select(-var)
# dt <- dt %>%
#   mutate(tau=replace(tau,tau==1,0.1),
#          tau=replace(tau,tau==2,0.25),
#          tau=replace(tau,tau==3,0.5),
#          tau=replace(tau,tau==4,0.75),
#          tau=replace(tau,tau==5,0.9))
# dt <- dt %>%
#   spread(key="Type",value="value")
# p <- ggplot(dt) +
#   geom_line(aes(x=p,y=med),color="red") +
#   geom_line(aes(x=p,y=low),color="red",linetype="dashed") +
#   geom_line(aes(x=p,y=upp),color="red",linetype="dashed") +
#   geom_line(aes(x=p,y=tru)) +
#   facet_grid(tau~s,scales="free_y",labeller=as_labeller(appender2,
#                                                         default = label_parsed)) +
#   xlab(expression(tilde(p))) +
#   ylab(expression(paste(q[tau[y]],"(",tilde(p),",",tilde(s),")"))) +
#   scale_x_continuous(breaks=c(-0.25,0.25,0.75)) +
#   theme_bw()
# ggsave("y_non_p_lin.pdf",p, width = 8, height = 6)
# 
# apprx = y_app_non_p_non
# dt <- data.frame(p=rep(seq(-0.5,1,0.02),5),
#                  s=rep(round(s_tru_quants,3),each=76),
#                  med=apply(apprx,c(1,2),quantile,probs=0.5),
#                  low=apply(apprx,c(1,2),quantile,probs=0.025),
#                  upp=apply(apprx,c(1,2),quantile,probs=0.975)) %>% 
#   mutate(`tru.1`=m_y_non(p,s)+s_y*qsn(0.1,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.2`=m_y_non(p,s)+s_y*qsn(0.25,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.3`=m_y_non(p,s)+s_y*qsn(0.5,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.4`=m_y_non(p,s)+s_y*qsn(0.75,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y),
#          `tru.5`=m_y_non(p,s)+s_y*qsn(0.9,xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),omega=sqrt((1-2*d_y^2/pi)^-1),alpha=a_y)) %>%
#   gather(key="var",value="value",-c("p","s"))
# dt <- dt %>%
#   mutate(Type=sapply(strsplit(dt$var,"\\."),"[[",1)) %>%
#   mutate(tau=(as.numeric(sapply(strsplit(dt$var,"\\."),"[[",2)))) %>%
#   dplyr::select(-var)
# dt <- dt %>%
#   mutate(tau=replace(tau,tau==1,0.1),
#          tau=replace(tau,tau==2,0.25),
#          tau=replace(tau,tau==3,0.5),
#          tau=replace(tau,tau==4,0.75),
#          tau=replace(tau,tau==5,0.9))
# dt <- dt %>%
#   spread(key="Type",value="value")
# p <- ggplot(dt) +
#   geom_line(aes(x=p,y=med),color="red") +
#   geom_line(aes(x=p,y=low),color="red",linetype="dashed") +
#   geom_line(aes(x=p,y=upp),color="red",linetype="dashed") +
#   geom_line(aes(x=p,y=tru)) +
#   facet_grid(tau~s,scales="free_y",labeller=as_labeller(appender2, 
#                                                         default = label_parsed)) + 
#   xlab(expression(tilde(p))) +
#   ylab(expression(paste(q[tau[y]],"(",tilde(p),",",tilde(s),")"))) +
#   scale_x_continuous(breaks=c(-0.25,0.25,0.75)) +
#   theme_bw()
# ggsave("y_non_p_non.pdf",p, width = 8, height = 6)

m = 1

p_app_sim_lin <- readRDS(paste0("saves/p_app_sim_lin",m,".RDS"))
p_app_sim_non <- readRDS(paste0("saves/p_app_sim_non",m,".RDS"))
y_app_sim_non_p_lin <- readRDS(paste0("saves/y_app_sim_non_p_lin",m,".RDS"))
y_app_sim_non_p_non <- readRDS(paste0("saves/y_app_sim_non_p_non",m,".RDS"))

# #Sys.setenv(PATH=paste0("C:\\Users\\matsu\\AppData\\Local\\Programs\\orca;", Sys.getenv("PATH")))
# 
appender <- function(string) {
  TeX(paste("$\\tilde{s} = $", string))
}

p_seq <- seq(-0.9,0.9,0.05)
y_seq <- seq(-70,70,5)
p_len <- length(p_seq)
y_len <- length(y_seq)
dt <- data.frame(p = rep(p_seq,y_len),
                 y = rep(y_seq,each=p_len))
dt$val1 <- dt$val2 <- dt$val3 <- 0
# m <- list(
#   l = 0,
#   r = 0,
#   b = 0,
#   t = 0,
#   pad = 1
# )

for (i in c(86,68,98)){
  s <- s_sim[i]
  p_app <- p_app_sim_lin[i,]
  y_app <- y_app_sim_non_p_lin[i,]
  app <- kde2d(p_app,y_app,n=50)
  g_app <- apply(dt,1,function(k){
    x1 <- app$x[max(which(app$x <= k[1]))]
    x2 <- app$x[max(which(app$x <= k[1])) + 1]
    y1 <- app$y[max(which(app$y <= k[2]))]
    y2 <- app$y[max(which(app$y <= k[2])) + 1]
    z11 <- app$z[max(which(app$x <= k[1])),max(which(app$y <= k[2]))]
    z12 <- app$z[max(which(app$x <= k[1])),max(which(app$y <= k[2])) + 1]
    z21 <- app$z[max(which(app$x <= k[1])) + 1,max(which(app$y <= k[2]))]
    z22 <- app$z[max(which(app$x <= k[1])) + 1,max(which(app$y <= k[2])) + 1]
    p1 <- (x2 - k[1])/(x2 - x1)*z11 + (k[1] - x1)/(x2 - x1)*z21
    p2 <- (x2 - k[1])/(x2 - x1)*z12 + (k[1] - x1)/(x2 - x1)*z22
    g <- (y2 - k[2])/(y2 - y1)*p1 + (k[2] - y1)/(y2 - y1)*p1
    return(g)
  })

  dt$g <- g_tru <- dsn((dt$p - m_p_lin(s))/s_p_lin(s),
               xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
               omega=sqrt((1-2*d_p^2/pi)^-1),
               alpha=a_p)/s_p_lin(s)*
    dsn((dt$y - m_y_non(dt$p,s))/s_y,
        xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
        omega=sqrt((1-2*d_y^2/pi)^-1),
        alpha=a_y)/s_y

  if (i == 86){
    dt$val1 <- (g_app - g_tru)^2
  } else if (i == 68){
    dt$val2 <- (g_app - g_tru)^2
  } else {
    dt$val3 <- (g_app - g_tru)^2
  }
}

# ggplot(dt %>%
#          filter(p >= -0.35 & p <= 0.35 & y >= -70 & y <= 40) %>%
#          gather(type,value,-c(p,y)) %>%
#          filter(type != "g") %>%
#          mutate(type=replace(type,type=="val1","0.093")) %>%
#          mutate(type=replace(type,type=="val2","0.173")) %>%
#          mutate(type=replace(type,type=="val3","0.281"))) +
#   geom_tile(aes(x = p,
#                 y = y,fill=value)) +
#   scale_fill_gradient2(
#                        limits = c(0,4),
#                        low="blue",mid="white",high="red") +
#   facet_grid(type~.)+
#   xlab(expression(tilde(p)[t])) + ylab(expression(tilde(y)[t])) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   labs(fill=expression(paste("(",hat(g) - g,")/",g)))
# ggsave("p_lin_y_non_heat.pdf",width=6,height=8)

ggplot(dt %>%
         filter(p >= -0.35 & p <= 0.35 & y >= -70 & y <= 40) %>%
         gather(type,value,-c(p,y)) %>%
         filter(type != "g") %>%
         mutate(type=replace(type,type=="val1","0.093")) %>%
         mutate(type=replace(type,type=="val2","0.173")) %>%
         mutate(type=replace(type,type=="val3","0.281"))) +
  geom_tile(aes(x = p,
                y = y,fill=value)) +
  scale_fill_gradient2(low="white",high="red") +
  facet_grid(type~.,labeller=as_labeller(appender,
                                         default = label_parsed))+
  xlab(expression(tilde(p)[t])) + ylab(expression(tilde(y)[paste("jt")])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(fill=expression(paste("(",hat(g) - g,")^2")))
ggsave("p_lin_y_non_heat2.pdf",width=6,height=8)

for (i in c(86,68,98)){
  s <- s_sim[i]
  p_app <- p_app_sim_non[i,]
  y_app <- y_app_sim_non_p_non[i,]
  app <- kde2d(p_app,y_app,n=50)
  g_app <- apply(dt,1,function(k){
    x1 <- app$x[max(which(app$x <= k[1]))]
    x2 <- app$x[max(which(app$x <= k[1])) + 1]
    y1 <- app$y[max(which(app$y <= k[2]))]
    y2 <- app$y[max(which(app$y <= k[2])) + 1]
    z11 <- app$z[max(which(app$x <= k[1])),max(which(app$y <= k[2]))]
    z12 <- app$z[max(which(app$x <= k[1])),max(which(app$y <= k[2])) + 1]
    z21 <- app$z[max(which(app$x <= k[1])) + 1,max(which(app$y <= k[2]))]
    z22 <- app$z[max(which(app$x <= k[1])) + 1,max(which(app$y <= k[2])) + 1]
    p1 <- (x2 - k[1])/(x2 - x1)*z11 + (k[1] - x1)/(x2 - x1)*z21
    p2 <- (x2 - k[1])/(x2 - x1)*z12 + (k[1] - x1)/(x2 - x1)*z22
    g <- (y2 - k[2])/(y2 - y1)*p1 + (k[2] - y1)/(y2 - y1)*p1
    return(g)
  })

  dt$g <- g_tru <- dsn((dt$p - m_p_non(s))/s_p_non(s),
               xi=-sqrt((1-2*d_p^2/pi)^-1)*d_p*sqrt(2/pi),
               omega=sqrt((1-2*d_p^2/pi)^-1),
               alpha=a_p)/s_p_non(s)*
    dsn((dt$y - m_y_non(dt$p,s))/s_y,
        xi=-sqrt((1-2*d_y^2/pi)^-1)*d_y*sqrt(2/pi),
        omega=sqrt((1-2*d_y^2/pi)^-1),
        alpha=a_y)/s_y

  if (i == 1){
    dt$val1 <- (g_app - g_tru)^2
  } else if (i == 2){
    dt$val2 <- (g_app - g_tru)^2
  } else {
    dt$val3 <- (g_app - g_tru)^2
  }
}

# ggplot(dt %>%
#          filter(p >= -0.35 & p <= 0.35 & y >= -70 & y <= 40) %>%
#          gather(type,value,-c(p,y)) %>%
#          filter(type != "g") %>%
#          mutate(type=replace(type,type=="val1","0.093")) %>%
#          mutate(type=replace(type,type=="val2","0.173")) %>%
#          mutate(type=replace(type,type=="val3","0.281"))) +
#   geom_tile(aes(x = p,
#                 y = y,fill=value)) +
#   scale_fill_gradient2(
#                        limits = c(0,4),
#                        low="blue",mid="white",high="red") +
#   facet_grid(type~.)+
#   xlab(expression(tilde(p)[t])) + ylab(expression(tilde(y)[t])) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   labs(fill=expression(paste("(",hat(g) - g,")/",g)))
# ggsave("p_non_y_non_heat.pdf",width=6,height=8)

ggplot(dt %>%
         filter(p >= -0.35 & p <= 0.35 & y >= -70 & y <= 40) %>%
         gather(type,value,-c(p,y)) %>%
         filter(type != "g") %>%
         mutate(type=replace(type,type=="val1","0.093")) %>%
         mutate(type=replace(type,type=="val2","0.173")) %>%
         mutate(type=replace(type,type=="val3","0.281"))) +
  geom_tile(aes(x = p,
                y = y,fill=value)) +
  scale_fill_gradient2(low="white",high="red") +
  facet_grid(type~.,labeller=as_labeller(appender,
                                         default = label_parsed))+
  xlab(expression(tilde(p)[t])) + ylab(expression(tilde(y)[paste("jt")])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(fill=expression(paste("(",hat(g) - g,")^2")))
ggsave("p_non_y_non_heat2.pdf",width=6,height=8)

# y_lin_p_lin_app <- NULL
# y_lin_p_non_app <- NULL
# y_non_p_lin_app <- NULL
# y_non_p_non_app <- NULL
# 
# for (m in 1:100){
#   print(paste0("m = ",m))
#   p_app_sim_lin <- try(readRDS(paste0("saves/p_app_sim_lin",m,".RDS")),silent=TRUE)
#   p_app_sim_non <- try(readRDS(paste0("saves/p_app_sim_non",m,".RDS")),silent=TRUE)
#   y_app_sim_lin_p_lin <- try(readRDS(paste0("saves/y_app_sim_lin_p_lin",m,".RDS")),silent=TRUE)
#   y_app_sim_lin_p_non <- try(readRDS(paste0("saves/y_app_sim_lin_p_non",m,".RDS")),silent=TRUE)
#   y_app_sim_non_p_lin <- try(readRDS(paste0("saves/y_app_sim_non_p_lin",m,".RDS")),silent=TRUE)
#   y_app_sim_non_p_non <- try(readRDS(paste0("saves/y_app_sim_non_p_non",m,".RDS")),silent=TRUE)
#   y_lin_p_lin_tmp <- NULL
#   y_lin_p_non_tmp <- NULL
#   y_non_p_lin_tmp <- NULL
#   y_non_p_non_tmp <- NULL
#   if (class(p_app_sim_lin) != "try-error" &&
#       class(p_app_sim_non) != "try-error" &&
#       class(y_app_sim_lin_p_lin) != "try-error" &&
#       class(y_app_sim_non_p_lin) != "try-error" &&
#       class(y_app_sim_lin_p_non) != "try-error" && 
#       class(y_app_sim_non_p_non) != "try-error"){
#     for (t in 1:100){
#       print(paste0("t = ",t))
#       y_lin_p_lin_tmp <- c(y_lin_p_lin_tmp,cor(p_app_sim_lin[t,],y_app_sim_lin_p_lin[t,]))
#       y_lin_p_non_tmp <- c(y_lin_p_non_tmp,cor(p_app_sim_non[t,],y_app_sim_lin_p_non[t,]))
#       y_non_p_lin_tmp <- c(y_non_p_lin_tmp,cor(p_app_sim_lin[t,],y_app_sim_non_p_lin[t,]))
#       y_non_p_non_tmp <- c(y_non_p_non_tmp,cor(p_app_sim_non[t,],y_app_sim_non_p_non[t,]))
#     }
#     y_lin_p_lin_app <- cbind(y_lin_p_lin_app,y_lin_p_lin_tmp)
#     y_lin_p_non_app <- cbind(y_lin_p_non_app,y_lin_p_non_tmp)
#     y_non_p_lin_app <- cbind(y_non_p_lin_app,y_non_p_lin_tmp)
#     y_non_p_non_app <- cbind(y_non_p_non_app,y_non_p_non_tmp)
#   }
# }
# 
# y_lin_p_lin_tru <- 14.45*s_p_lin(s)/s_y
# y_lin_p_non_tru <- 14.45*s_p_non(s)/s_y
# s_p <- (1 - 2*d_p^2/pi)^(-1/2)
# m_p <- -s_p*d_p*sqrt(2/pi)
# dM <- function(t){
#   2*exp(m_p*t+s_p^2*t^2/2)*
#     (pnorm(s_p*d_p*t)*(m_p+s_p^2*t) + dnorm(s_p*d_p*t)*s_p*d_p)
# }
# y_non_p_lin_tru <- 14.45*exp(m_p_lin(s))*dM(s_p_lin(s))/s_y
# y_non_p_non_tru <- 14.45*exp(m_p_non(s))*dM(s_p_non(s))/s_y
# 
# appender <- function(string) {
#   if (string == "y_lin"){
#     TeX(paste("Linear ","$\\tilde{y}$"))
#   } else {
#     TeX(paste("Non-linear ","$\\tilde{y}$"))
#   }
# }
# 
# y_lin_p_lin <- data.frame(s = s,
#                           p = "p_lin",
#                           y = "y_lin",
#                           med = apply(y_lin_p_lin_app,1,quantile,probs=0.5),
#                           low = apply(y_lin_p_lin_app,1,quantile,probs=0.025),
#                           upp = apply(y_lin_p_lin_app,1,quantile,probs=0.975),
#                           tru = y_lin_p_lin_tru)
# 
# y_lin_p_non <- data.frame(s = s,
#                           p = "p_non",
#                           y = "y_lin",
#                           med = apply(y_lin_p_non_app,1,quantile,probs=0.5),
#                           low = apply(y_lin_p_non_app,1,quantile,probs=0.025),
#                           upp = apply(y_lin_p_non_app,1,quantile,probs=0.975),
#                           tru = y_lin_p_non_tru)
# 
# y_non_p_lin <- data.frame(s = s,
#                           p = "p_lin",
#                           y = "y_non",
#                           med = apply(y_non_p_lin_app,1,quantile,probs=0.5),
#                           low = apply(y_non_p_lin_app,1,quantile,probs=0.025),
#                           upp = apply(y_non_p_lin_app,1,quantile,probs=0.975),
#                           tru = y_non_p_lin_tru)
# 
# y_non_p_non <- data.frame(s = s,
#                           p = "p_non",
#                           y = "y_non",
#                           med = apply(y_non_p_non_app,1,quantile,probs=0.5),
#                           low = apply(y_non_p_non_app,1,quantile,probs=0.025),
#                           upp = apply(y_non_p_non_app,1,quantile,probs=0.975),
#                           tru = y_non_p_non_tru)
# 
# as.data.frame(rbind(y_lin_p_lin,y_lin_p_non,y_non_p_lin,y_non_p_non)) %>%
#   ggplot() +
#   geom_line(aes(x=s,y=med),color="red") +
#   geom_line(aes(x=s,y=low),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=upp),color="red",linetype="dashed") +
#   geom_line(aes(x=s,y=tru)) +
#   facet_grid(y~p,scales="free_y",labeller=as_labeller(appender, 
#                                                       default = label_parsed)) +
#   xlab(expression(tilde(s))) +
#   ylab(expression(paste("cor(",tilde(y),",",tilde(p),"|",tilde(s),")"))) +
#   theme_bw()
# ggsave("cor.pdf",width=6,height=8)

