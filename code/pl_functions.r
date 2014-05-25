source("dlm_mcmc_functions.r")

# Functions to run particle learning for local level DLM

dlpred.ll <- function(y, x, suff.x, theta, lambda=1, rb=F)
{
  log.pred = ifelse(rb,dnorm(y,suff.x[1],sqrt(suff.x[2]+theta*(1+lambda)),log=T),dnorm(y, x, sqrt(theta*(1+lambda)), log=T))
  return(log.pred)
}

revo.ll <- function(y, x, suff.x, theta, lambda=1, rb=F)
{
  omega = lambda / (1+lambda)
  mu = omega*(y + (rb*suff.x[1] + (1-rb)*x)/lambda)
  tau = sqrt(omega*(theta + rb*(suff.x[2]/(lambda*(1+lambda)))))
  return(rnorm(1,mu,tau))
}

rprior.ll <- function(a0=1,b0=1)
{
  theta = 1 / rgamma(1, a0, b0)
  x = rnorm(1, 0, sqrt(theta))
  suff.theta = c(a0+.5,b0+x^2)
  suff.x = c(0, theta)
  return(list(x=x,theta=theta,suff.x=suff.x,suff.theta=suff.theta))
}

rmove.ll <- function(theta, suff.theta)
{
  1 / rgamma(1,suff.theta[1],suff.theta[2])
}

smap.theta.ll <- function(suff.theta, y, x.new, x.curr, lambda=1)
{
  suff.theta[1] = suff.theta[1] + 1
  suff.theta[2] = suff.theta[2] + .5*((y-x.new)^2 + (1/lambda)*(x.new-x.curr)^2)
  return(suff.theta)
}

smap.state.ll <- function(suff.x, y, theta, lambda=1)
{
  A = (suff.x[2]+theta*lambda)/(suff.x[2]+theta*(1+lambda))
  suff.x[1] = A*y + (1-A)*suff.x[1]
  suff.x[2] = A*theta
  return(suff.x)
}

# Functions to run particle learning for AR(1) DLM
# theta = (beta, phi, sigma2s, sigma2m)
# x = (x, t) (to track time, since U_t and F_t change over time)
# U is an nt by d matrix of covariates
# F is an nt length vector of known state coefficients

dlpred.ar <- function(y, x, suff.x, theta, U, F, rb=T)
{
  d = dim(U)[2]
  mu = t(U[x[2],])%*%theta[1:d] + F[x[2]]*theta[d+1]*(rb*suff.x[1] + (1-rb)*x[1])
  tau = sqrt((F[x[2]]^2)*(rb*((theta[d+1]^2)*suff.x[2]) + theta[d+2]) + theta[d+3])
  log.pred = dnorm(y, mu, tau, log=T)
  return(log.pred)
}

revo.ar <- function(y, x, suff.x, theta, U, F, rb=T)
{
  d = dim(U)[2]
  omega = 1/((F[x[2]]^2)*(1/theta[d+3]) + 1/theta[d+2])
  mu = omega*((y - t(U[x[2],])%*%theta[1:d])*F[x[2]]*(1/theta[d+3]) + theta[d+1]*(rb*suff.x[1] + (1-rb)*x[1])*(1/theta[d+2]))
  tau = sqrt(rb*suff.x[2]*((theta[d+3]*theta[d+1])/((F[x[2]]^2)*theta[d+2] + theta[d+3]))^2 + omega)
  return(rnorm(1,mu,tau))
}

rprior.ar.joint <- function(prior)
{
  sigma2m = 1 / rgamma(1, prior$am0, prior$bm0)  
  d = length(prior$b0)
  cl = chol(prior$B0)
  beta = as.numeric(prior$b0 + t(cl)%*%rnorm(d, 0, sqrt(sigma2m))%*%cl)
  sigma2s = 1 / rgamma(1, prior$as0, prior$bs0)
  phi = rnorm(1, prior$phi0[[1]], sqrt(sigma2s*prior$Phi0[[1]][1,1]))
  while(!is.stationary(phi)) phi = rnorm(1, prior$phi0[[1]], sqrt(sigma2s*prior$Phi0[[1]][1,1]))
  x0 = c(rnorm(1, 0, sqrt(sigma2s / (1 - phi^2))), 0)
  suff.x = c(0, sigma2s / (1 - phi^2), 0)
  suff.theta = c(0, 0, 0, 0, x0^2, 0, prior$am0, prior$as0 + .5, x0)
  return(list(x=x0, theta = c(beta,phi,sigma2s,sigma2m),suff.x=suff.x,suff.theta=suff.theta))
}

rprior.ar.marg <- function(prior)
{
  d = length(prior$b0)
  cl = chol(prior$B0)
  beta = as.numeric(prior$b0 + t(cl)%*%rnorm(d, 0, 1)%*%cl)
  phi = rnorm(1, prior$phi0[[1]], sqrt(prior$Phi0[[1]][1,1]))
  while(!is.stationary(phi)) phi = rnorm(1, prior$phi0[[1]], sqrt(prior$Phi0[[1]][1,1]))
  sigma2s = 1 / rgamma(1, prior$as0, prior$bs0)
  sigma2m = 1 / rgamma(1, prior$am0, prior$bm0)
  x0 = c(rnorm(1, 0, sqrt(sigma2s / (1 - phi^2))), 0)
  suff.x = c(0, sigma2s / (1 - phi^2), 0)
  suff.theta = c(numeric(d), 0, 0, 0, x0^2, 0, prior$am0, prior$as0 + .5, x0)
  return(list(x=x0, theta = c(beta,phi,sigma2s,sigma2m),suff.x=suff.x,suff.theta=suff.theta))
}

rmove.ar.joint <- function(theta, suff.theta, prior)
{
  d = length(prior$b0)
  B0.prec = solve(prior$B0)
  Bt = solve(suff.theta[2] + B0.prec)
  bt = as.numeric(Bt%*%(suff.theta[1]+B0.prec%*%prior$b0))
  bmt = .5*(suff.theta[3] + t(prior$b0)%*%B0.pred%*%prior$b0 + t(bt)%*%solve(Bt)%*%bt) + prior$bm0
  theta[d+3] = 1 / rgamma(1, suff.theta[7], bmt)
  clbeta = chol(Bt)
  theta[1:d] = bt + t(clbeta)%*%rnorm(d, 0, sqrt(theta[d+3]))%*%clbeta
  Phi0.prec = solve(prior$Phi0)
  Phit = solve(suff.theta[5] + Phi0.prec)
  phit = as.numeric(Phit%*%(suff.theta[4]+Phi0.prec%*%prior$phi0))
  bst = .5*(suff.theta[6] + t(prior$phi0)%*%Phi0.pred%*%prior$phi0 + t(phit)%*%solve(Phit)%*%phit) + prior$bs0
  sigma2s = 1 / rgamma(1, suff.theta[8], bst)
  clphi = chol(Phit)
  phi = phit + as.numeric(t(clphi)%*%rnorm(1, 0, sqrt(sigma2s))%*%clphi)
  while(!is.stationary(phi)) phi = phit + as.numeric(t(clphi)%*%rnorm(1, 0, sqrt(sigma2s))%*%clphi)
  logMH = Psi.C0(suff.theta[9], 0, phi, sigma2s) - Psi.C0(suff.theta[9], 0, theta[d+1], theta[d+2])
  u = log(runif(1))
  if(u < logMH) theta[(d+1):(d+2)] = c(phi,sigma2s)
  return(theta)
}

rmove.ar.marg <- function(theta, suff.theta, prior)
{
  d = length(prior$b0)
  B0.prec = solve(prior$B0)
  Bt = solve(suff.theta[2]/theta[d+3] + B0.prec)
  bt = as.numeric(Bt%*%(suff.theta[1]/theta[d+3]+B0.prec%*%prior$b0))
  clbeta = chol(Bt)
  theta[1:d] = bt + t(clbeta)%*%rnorm(d, 0, 1)%*%clbeta
  SSy = t(theta[1:d])%*%suff.theta[2]%*%theta[1:d] - 2*t(theta[1:d])%*%suff.theta[1] + suff.theta[3]
  bmt = .5*as.numeric(SSy) + prior$bm0
  theta[d+3] = 1 / rgamma(1, suff.theta[7], bmt)
  Phi0.prec = solve(prior$Phi0)
  Phit = solve(suff.theta[5]/theta[d+2] + Phi0.prec)
  phit = as.numeric(Phit%*%(suff.theta[4]/theta[d+2]+Phi0.prec%*%prior$phi0))
  clphi = chol(Phit)
  phi = phit + as.numeric(t(clphi)%*%rnorm(1, 0, 1)%*%clphi)
  while(!is.stationary(phi)) phi = phit + as.numeric(t(clphi)%*%rnorm(1, 0, 1)%*%clphi)
  logMH = Psi.C0(suff.theta[9], 0, phi, theta[d+2]) - Psi.C0(suff.theta[9], 0, theta[d+1], theta[d+2])
  u = log(runif(1))
  if(u < logMH) theta[d+1] = phi
  SSx = suff.theta[6] + suff.theta[5]*theta[d+1]^2 - 2*theta[d+1]*suff.theta[4]
  bst = .5*(SSx + (1-theta[d+1]^2)*suff.theta[9]^2) + prior$bs0
  theta[d+2] = 1 / rgamma(1, suff.theta[8], bst)
  return(theta)
}

smap.theta.ar <- function(suff.theta, y, x.new, x.curr, U, F)
{
  suff.theta[1] = suff.theta[1] + U[x.new[2],]%*%(y - F[x.new[2]]*x.new)
  suff.theta[2] = suff.theta[2] + U[x.new[2],]%*%t(U[x.new[2],])
  suff.theta[3] = suff.theta[3] + (y - F[x.new[2]]*x.new)^2
  suff.theta[4:8] = suff.theta[4:8] + c(x.new*x.curr, x.curr^2, x.new^2,.5,.5)
  return(suff.theta)
}

smap.state.ar <- function(suff.x, y, theta, U, F)
{
  d = dim(U)[2]
  A = ((F[suff.x[3]]^2)*((theta[d+1]^2)*suff.x[2]+theta[d+2]))/((F[suff.x[3]]^2)*((theta[d+1]^2)*suff.x[2]+theta[d+2]) + theta[d+3])
  suff.x[1] = A*(y-t(U[suff.x[3],])%*%theta[1:d])/F[suff.x[3]] + (1-A)*theta[d+1]*suff.x[1]
  suff.x[2] = ((theta[d+1]^2)*suff.x[2] + theta[d+2])*(1-A)
  suff.x[3] = suff.x[3] + 1
  return(suff.x)
}