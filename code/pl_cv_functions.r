# Functions to run particle learning for local level DLM

dlpred.ll <- function(y, x, suff.x, theta, lambda=1, rb=F)
{
  log.pred = ifelse(rb, dnorm(y,suff.x[1],sqrt(suff.x[2]+theta*(1+lambda)),log=T), dnorm(y, x[1], sqrt(theta*(1+lambda)),log=T))
  return(log.pred)
}

revo.ll <- function(y, x, suff.x, theta, lambda=1, rb=F)
{
  omega = lambda / (1+lambda)
  mu = ifelse(rb, omega*(y + suff.x[1]/lambda), omega*(y + x[1]/lambda))
  tau = ifelse(rb, sqrt(omega*(theta + suff.x[2]/(lambda*(1+lambda)))), sqrt(omega*theta))
  return(c(rnorm(1,mu,tau), x[1]))
}

rprior.ll <- function(a0=1,b0=1,lambda=1)
{
  theta = 1 / rgamma(1, a0, b0)
  x = c(rnorm(1, 0, sqrt(theta)), 0)
  suff.theta = c(a0+.5, b0 + .5*x[1]^2)
  suff.x = c(0, theta)
  theta = 1 / rgamma(1, suff.theta[1], suff.theta[2])
  return(list(x=x,theta=theta,suff.x=suff.x,suff.theta=suff.theta))
}

rmove.ll <- function(suff.theta)
{
  1 / rgamma(1,suff.theta[1],suff.theta[2])
}

smap.theta.ll <- function(suff.theta, y, x, lambda=1)
{
  suff.theta[1] = suff.theta[1] + 1
  suff.theta[2] = suff.theta[2] + .5*((y-x[1])^2 + (1/lambda)*((x[1]-x[2])^2))
  return(suff.theta)
}

smap.state.ll <- function(suff.x, y, theta, lambda=1)
{
  A = (suff.x[2]+theta*lambda)/(suff.x[2]+theta*(1+lambda))
  suff.x[1] = A*y + (1-A)*suff.x[1]
  suff.x[2] = A*theta
  return(suff.x)
}