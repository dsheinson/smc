###########################################################################
# Functions associated with DLM formulation of regression with AR(p) errors
###########################################################################

dllik <- function(y, x, F = 1, beta = 0, U = 0, sigma2 = 1, V = 1)
{
  d = length(beta)
  p = length(x)
  q = length(y)
  U = as.matrix(U); F = as.matrix(F); V = as.matrix(V)
  stopifnot(dim(U)[1] == q & dim(F)[1] == q & dim(U)[2] == d & dim(F)[2] == p & dim(V)[1] == q & dim(V)[1] == dim(V)[2])
  
  return(mvdnorm(y,U%*%beta + F%*%x,sigma2*V,log=TRUE))
}

revo <- function(x, phi = 1, sigma2 = 1, W = 1)
{
  p = length(x)
  W = as.matrix(W)
  stopifnot(length(phi) == p & dim(W)[1] == p & dim(W)[2] == p)
  
  G = makeG(phi)
  return(rep(G%*%x + t(chol(W))%*%rnorm(p,0,sqrt(sigma2)),1))
}

###################
# Utility functions
###################

mvdnorm <- function(x, mu, Sigma, log = TRUE)
{
  d = length(x)
  Sigma = as.matrix(Sigma)
  stopifnot(length(mu) == d & dim(Sigma)[1] == d & dim(Sigma)[1] == d)
  
  log.lik = -(d/2)*log(2*pi) - .5*log(det(Sigma)) - .5*t(x - mu)%*%solve(Sigma)%*%(x - mu)
  if(log) return(log.lik) else return(exp(log.lik))
}

makeG <- function(phi) {
  m = length(phi)
  G <- diag(0,m)
  G[,1] <- phi
  G[-m,-1] <- diag(1,m-1)
  return(G)
}