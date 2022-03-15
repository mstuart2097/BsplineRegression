# Calculated the loss function for the quantile regression for a particular quantile (tau)
quant.loss=function(data,B,D,b,lambda,tau){
  # B is the B-spline matrix with degree `p` and given knots
  # b is the quantile regression coefficients
  # D is the difference matrix used in the penalty term
  # lambda is the smoothing parameter of the penalty term
  # tau is the quantile used in the simulation
  resid=data-B%*%b
  U=tau-as.numeric(resid<0)
  loss=sum(resid*U)+lambda/2*t(as.vector(b))%*%t(D)%*%D%*%as.vector(b)
  return(loss)
}

# Performs the quantile regression using the above loss function
quant.reg <- function(data,B,D,lambda,tau){
  BtB <- t(B) %*% B
  DtD <- t(D) %*% D
  b.init <- solve(BtB + lambda/2 * DtD) %*% t(B) %*% data
  beta <- nlm(f=quant.loss,p=b.init,data=data,B=B,D=D,lambda=lambda,tau=tau)
  b <- beta$estimate
  resid=data-B%*%b
  U=tau-as.numeric(resid<0)
  loss=sum(resid*U)
  return(list(b=b,fit=loss))
}