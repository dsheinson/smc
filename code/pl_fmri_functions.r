source("dlm_mcmc_functions.r")

# Functions to run particle learning for AR(p) DLM

# theta = (beta, phi, sigma2s, sigma2m) (d + p + 2)-length vector
# y is a scalar, univariate observation
# x = (x1,x2,...,xp, t) (last element to to track time, since U_t and F_t change over time)
# U is an nt by d matrix of known covariates
# F is an nt-length vector of known state coefficient matrices

dlpred.ar <- function(y, x, suff.x, theta, U, F)
{
  d = dim(U)[2]
  p = length(x) - 1
  Gx = sum(x[1:p]*theta[(d+1):(d+p)])
  U = t(U[x[p+1]+1,])
  F = F[x[p+1]+1]
  mu = (U%*%theta[1:d])[1,1] + F*Gx
  tau = sqrt((F^2)*theta[d+p+1] + theta[d+p+2])
  log.pred = dnorm(y, mu, tau, log=T)
  return(log.pred)
}

revo.ar <- function(y, x, suff.x, theta, U, F)
{
  d = dim(U)[2]
  p = length(x) - 1
  Gx = sum(x[1:p]*theta[(d+1):(d+p)])
  U = t(U[x[p+1]+1,])
  F = F[x[p+1]+1]
  omega = 1/((F^2)/theta[d+p+2] + 1/theta[d+p+1])
  mu = omega*((y - (U%*%theta[1:d])[1,1])*F*(1/theta[d+p+2]) + Gx/theta[d+p+1])
  tau = sqrt(omega)
  x1 = rnorm(1,mu,tau)
  if(p > 1) x1 = c(x1, x[1:(p-1)])
  return(c(x1,x[p+1]+1))
}

rprior.ar <- function(prior)
{
  sigma2m = 1 / rgamma(1, prior$am0, prior$bm0)  
  d = length(prior$b0)
  cl = chol(prior$B0)
  beta = as.numeric(prior$b0 + t(cl)%*%rnorm(d, 0, sqrt(sigma2m)))
  sigma2s = 1 / rgamma(1, prior$as0, prior$bs0)
  p = length(prior$phi0[[1]])
  cl = chol(prior$Phi0[[1]])
  phi = as.numeric(prior$phi0[[1]] + t(cl)%*%rnorm(p, 0, sqrt(sigma2s)))
  cl = chol(prior$C0)
  x0 = as.numeric(prior$m0 + t(cl)%*%rnorm(p, 0, 1))
  if(p > 1) x.pre = makex0(matrix(x0), phi) else x.pre = x0
  suff.x = c(prior$m0, as.numeric(prior$C0), 0)
  suff.theta = c(rep(0, d + d^2 + 1 + p), as.numeric(x.pre%*%t(x.pre)), 0, prior$am0+.5, prior$as0+.5)
  return(list(x=c(x.pre,0), theta = c(beta,phi,sigma2s,sigma2m),suff.x=suff.x,suff.theta=suff.theta))
}

rmove.ar <- function(suff.theta, prior)
{
  d = length(prior$b0)
  p = length(prior$phi0[[1]])
  B0.prec = solve(prior$B0)
  Bt = solve(matrix(suff.theta[(d+1):(d+d^2)],nr=d) + B0.prec)
  bt = as.numeric(Bt%*%(suff.theta[1:d]+B0.prec%*%prior$b0))
  bmt = .5*(suff.theta[d + d^2 + 1] + t(prior$b0)%*%B0.prec%*%prior$b0 - t(bt)%*%solve(Bt)%*%bt) + prior$bm0
  sigma2m = 1 / rgamma(1, suff.theta[d + d^2 + p + p^2 + 3], as.numeric(bmt))
  clbeta = chol(Bt)
  beta = bt + t(clbeta)%*%rnorm(d, 0, sqrt(sigma2m))
  Phi0.prec = solve(prior$Phi0[[1]])
  Phit = solve(matrix(suff.theta[(d + d^2 + 2 + p):(d + d^2 + 1 + p + p^2)],nr=p) + Phi0.prec)
  phit = as.numeric(Phit%*%(suff.theta[(d + d^2 + 2):(d + d^2 + 1 + p)]+Phi0.prec%*%prior$phi0[[1]]))
  bst = .5*(suff.theta[d + d^2 + p + p^2 + 2] + t(prior$phi0[[1]])%*%Phi0.prec%*%prior$phi0[[1]] - t(phit)%*%solve(Phit)%*%phit) + prior$bs0
  sigma2s = 1 / rgamma(1, suff.theta[d + d^2 + p + p^2 + 4], as.numeric(bst))
  clphi = chol(Phit)
  phi = phit + as.numeric(t(clphi)%*%rnorm(p, 0, sqrt(sigma2s)))
  return(c(beta,phi,sigma2s,sigma2m))
}

smap.theta.ar <- function(suff.theta, y, x.new, x.curr, U, F)
{
  d = dim(U)[2]
  p = length(x.new) - 1
  U = t(U[x.new[p+1],])
  F = F[x.new[p+1]]
  suff.theta[1:d] = suff.theta[1:d] + t(U)%*%(y - F*x.new[1])
  suff.theta[(d+1):(d+d^2)] = suff.theta[(d+1):(d+d^2)] + as.numeric(t(U)%*%U)
  suff.theta[d + d^2 + 1] = suff.theta[d + d^2 + 1] + (y - F*x.new[1])^2
  suff.theta[(d + d^2 + 2):(d + d^2 + p + 1)] = suff.theta[(d + d^2 + 2):(d + d^2 + 1 + p)] + x.new[1]*x.curr[1:p]
  suff.theta[(d + d^2 + p + 2):(d + d^2 + p + p^2 + 1)] = suff.theta[(d + d^2 + p + 2):(d + d^2 + p + p^2 + 1)] + as.numeric(x.curr[1:p]%*%t(x.curr[1:p]))
  suff.theta[(d + d^2 + p + p^2 + 2):(d + d^2 + p + p^2 + 4)] = suff.theta[(d + d^2 + p + p^2 + 2):(d + d^2 + p + p^2 + 4)] + c(sum(x.new[1:p]*x.new[1:p]), .5, .5)
  return(suff.theta)
}

smap.state.ar <- function(suff.x, y, theta, U, F)
{
  return(suff.x)
}

# Utility functions

rprior.convert <- function(cov, d, sd.fac = 1)
{
  p = length(cov$center) - d - 2
  b0 = cov$center[1:d]
  B0 = (sd.fac^2)*diag(d)%*%cov$cov[(1:d),(1:d)]
  phi0 = cov$center[(d+1):(d+p)]
  Phi0 = (sd.fac^2)*cov$cov[(d+1):(d+p),(d+1):(d+p)]
  b = cov$center[(d+p+1):(d+p+2)] / ((sd.fac^2)*diag(cov$cov[(d+p+1):(d+p+2),(d+p+1):(d+p+2)]))
  a = b*cov$center[(d+p+1):(d+p+2)]
  return(list(b0=b0, B0=B0, phi0=list(phi0), Phi0=list(Phi0), as0 = a[1], bs0 = b[1], am0 = a[2], bm0 = b[2]))
}