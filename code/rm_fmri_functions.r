source("dlm_mcmc_functions.r")
source("reg_mcmc_functions.r")

# Functions to pass into resample-move particle filter for dlm regression model
# data y is a scalar
# state x is a (p+1)-length vector with evolving state and time index
# fixed parameters theta is a (d + p + 2)-length vector with regression coefficients beta, AR parameters phi, measurement noise variance sigma2m and AR white noise variance sigma2s

dllik.dlm <- function(y, x, theta, U, u = NULL, dr = TRUE)
{
  if(!is.null(u)) stopifnot(length(u) == dim(U)[1])
  p = length(x) - 1
  t = x[p+1]
  d = dim(U)[2]
  beta = theta[1:d]
  if(dr) F = c(u[t],rep(0,p-1)) else F = c(1,rep(0,p-1))
  mu = sum(U[t,]*beta) + sum(F*x[1:p])
  return(dnorm(y,mu,sqrt(theta[d+p+1]),log=TRUE))
}

revo.dlm <- function(x, theta)
{
  p = length(x) - 1
  d = length(theta) - 2 - p
  phi = theta[(d+1):(d+p)]
  G = makeG(phi)
  mu = rep(G %*% x[1:p],1)
  mu[1] = rnorm(1, mu[1], sqrt(theta[d+p+2]))
  return(c(mu, x[p+1]+1))
}

rprior.dlm <- function(b0 = c(0,0), B0 = 1e6*diag(2), phi0 = 0, Phi0 = 1e6, am0 = 1e-6, bm0 = 1e-6, as0 = 1e-6, bs0 = 1e-6, m0=0)
{
  p = length(phi0)
  d = length(b0)
  beta = rep(b0 + t(chol(B0))%*%rnorm(d),1)
  phi = rep(phi0 + t(chol(Phi0))%*%rnorm(p),1)
  while(!is.stationary(phi)) phi = rep(phi0 + t(chol(Phi0))%*%rnorm(p),1)
  sigma2m = 1 / rgamma(1,am0,bm0)
  sigma2s = 1 / rgamma(1,as0,bs0)
  C0 = makeC0(phi)
  x0 = rep(m0 + t(chol(C0)) %*% rnorm(p, 0, sqrt(sigma2s)), 1)
  return(list(x=c(x0,0),theta=c(beta,phi,sigma2m,sigma2s)))
}

rmove.dlm <- function(y, x, theta, psi, prior, n.iter = 1, smooth = TRUE, store.smooth = TRUE)
{ 
  # Get dimensions and adjust parameters for MCMC functions
  p = dim(x)[1] - 1
  d = length(theta) - 2 - p
  t = dim(x)[2] - 1
  x = matrix(x[1:p,],nr=p)
  theta.list = list(beta = theta[1:d], phi = theta[(d+1):(d+p)], sigma2m = theta[d+p+1], sigma2s = theta[d+p+2])
  y = matrix(y, nr = 1)
  nt = dim(y)[2]
  psi$F = array(psi$F[,,1:nt], c(1,p,nt))
  psi$U = array(psi$U[,,1:nt], c(1,d,nt))
  
  for(k in 1:n.iter)
  {
    # Sample smoothed states by FFBS?
    if(smooth) x.smooth = sample.states(y, x, theta.list, psi, prior) else x.smooth = x
    
    # Sample theta from full conditional
    theta.list$beta = sample.beta(y, x.smooth, theta.list, psi, prior)
    theta.list$phi = sample.phi(y, x.smooth, theta.list, psi, prior)$phi
    theta.list$sigma2m = sample.sigma2m(y, x.smooth, theta.list, psi, prior)
    theta.list$sigma2s = sample.sigma2s(y, x.smooth, theta.list, psi, prior)
    theta = c(theta.list$beta, theta.list$phi, theta.list$sigma2m, theta.list$sigma2s)
    
    # Return smoothed states?
    if(store.smooth) x = rbind(x.smooth,0:t) else x = rbind(x,0:t)
  }
  return(list(state=x,theta=theta))
}