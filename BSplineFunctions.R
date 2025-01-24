
# Function to find moments of Y and P using B-spline
# data is the imported data
# B is the B-spline matrix (you can use function bs in packages "splines" to calculate this)
# D is the m-th difference matrix used in the penalty term
# lambda is the smoothing parameter

# var.loss <- function(data,B,D,b,lambda1,lambda2){
#   Delta <- tail(b,length(b)-2) - 
#     2*head(b[tail(b,length(b) - 1)],length(b)-2) +
#     head(b,length(b)-2)
#   w <- (Delta < 0)*1
#   resid <- data - B %*% b
#   loss <- sum(resid^2) + lambda1/2 * (b %*% t(D) %*% D %*% b) + lambda2 * w*Delta^2
#   return(loss)
# }
# 

# Functions to find
# Calculate Euclidean Norm
sigmoid <- function(x){exp(x)/(1+exp(x))}

# Performs the quantile regression using the above loss function
mean.reg <- function(data,B,D,lambda){
  BtB <- t(B) %*% B
  DtD <- t(D) %*% D
  return(solve(BtB + lambda/2 * DtD) %*% t(B) %*% data)
}

# Calculated the loss function for the quantile regression for a particular quantile (tau)
quant.loss=function(data,B,D,b,lambda,tau){
  resid=data-B%*%b
  U=tau-as.numeric(resid<0)
  loss=sum(resid*U)+lambda/2*t(as.vector(b))%*%t(D)%*%D%*%as.vector(b)
  return(loss)
}

# Performs the quantile regression using the above loss function
quant.reg <- function(data,B,D,lambda,tau){
  n <- length(data)
  BtB <- t(B) %*% B
  DtD <- t(D) %*% D
  b <- solve(BtB + lambda/2 * DtD) %*% t(B) %*% data
  sigma2 <- mean((data - B%*%b)^2)
  b.init <- qnorm(tau,b,sqrt(sigma2))
  beta <- nlm(f=quant.loss,p=b.init,data=data,B=B,D=D,lambda=lambda,tau=tau)
  b <- beta$estimate
  resid <- data-B%*%b
  U <- tau-as.numeric(resid<0)
  loss <- sum(resid*U)
  return(list(b=b,fit=loss))
}

# Extracts the GACV value for a given lambda
GACV.fn <- function(data,B,D,lambda,tau){
  n <- length(data)
  BtB <- t(B) %*% B
  DtD <- t(D) %*% D
  b <- solve(BtB + lambda/2 * DtD) %*% t(B) %*% data
  sigma2 <- mean((data - B%*%b)^2)
  b.init <- qnorm(tau,b,sqrt(sigma2))
  beta <- nlm(f=quant.loss,p=b.init,data=data,B=B,D=D,lambda=lambda,tau=tau)
  b <- beta$estimate
  resid <- data-B%*%b
  U <- tau-as.numeric(resid<0)
  loss <- sum(resid*U)
  trace_hat <- sum(diag(B %*% solve(BtB + lambda * DtD) %*% t(B)))
  GACV <- loss/(n - trace_hat)
  #return(list(b=b,fit=loss,vcov=solve(beta$hessian)))
  return(GACV)
}

optim.lambda <- function(data,B,D,tau,start=0.1){
  lambda <- nlm(f=GACV.fn,p=start,data=data,B=B,D=D,tau=tau,iterlim = 1000)
  j <- 1
  while(lambda$estimate < 0){
    lambda <- nlm(f=GACV.fn,p=start+0.1*j,data=data,B=B,D=D,tau=tau,iterlim = 1000)
    j <- j + 1
  }
  # if(lambda$code == 3){
  #   lambda <- nlm(f=GACV.fn,p=lambda$estimate,data=data,B=B,D=D,tau=tau,iterlim = 1000)
  # }
  return(lambda)
}

# Simulation based on quantile regression
# Input is the Bspline matrix, the matrix of regression coefficients,
# and the vector of quantiles 
# sim.quant.reg <- function(B,beta,tau){
#   quants <- B %*% beta
#   d <- nrow(B)
#   n <- length(tau)
#   u <- runif(d)
#   res <- NULL
#   for (i in 1:d){
#     r <- u[i]
#     tau.low <- ifelse(r < tau[1],1,
#                       ifelse(r > tau[n],n - 1,
#                              which.max(tau[between(tau,0,r)])))
#     tau.upp <- tau.low + 1
#     b.low <- ifelse(tau.low==1,
#                     beta[,1],
#                     beta[,1] + apply(exp(beta[,2:tau.low]),1,cumsum))
#     b.upp <- beta[,1] + apply(exp(beta[,2:tau.upp]),1,cumsum)
#     quants.low <- B[i,] %*% b.low
#     quants.upp <- B[i,] %*% b.upp
#     slope <- (quants.upp - quants.low)/(tau[tau.upp] - tau[tau.low])
#     intercept <- quants.low - slope*tau[tau.low]
#     res <- c(res,intercept+slope*r)
#   }
#   return(res)
# }



