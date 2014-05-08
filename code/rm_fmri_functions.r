source("dlm_mcmc_functions.r")
source("reg_mcmc_functions.r")
source("dlm_ar_functions.r")

# y is a scalar
# x is a 2-length vector of state and time index
# theta is a length-5 vector (beta0, beta1, phi, sigma2s, sigma2m)
# u is an nt-length vector of the predicted BOLD response

dllik.M101 <- function(y, x, theta, u) dnorm(y, theta[1] + (theta[2] + x[1])*u[x[2]], sqrt(theta[5]),log=TRUE)
dllik.M011 <- function(y, x, theta, u) dnorm(y, theta[1] + theta[2]*u[x[2]] + x[1], sqrt(theta[5]),log=TRUE)

revo <- function(x, theta) c(rnorm(1, theta[3]*x[1], sqrt(theta[4])), x[2]+1)

rprior <- function(b0 = c(900,5), B0 = c(5,5), phi0 = 0.8, Phi0 = 0.2, as0 = 1, bs0 = 1, am0 = 1, bm0 = 1)
{
  beta = rnorm(2, b0, sqrt(B0))
  phi = rnorm(1, phi0, sqrt(Phi0))
  while(!is.stationary(phi)) phi = rnorm(1, phi0, sqrt(Phi0))
  sigma2s = 1 / rgamma(1, as0, bs0)
  sigma2m = 1/ rgamma(1, am0, bm0)
  x0 = c(rnorm(1, 0, sqrt(sigma2s / (1 - phi^2))), 0)
  return(list(x=x0, theta = c(beta,phi,sigma2s,sigma2m)))
}

rmcmc <- function(y, x, theta, u, mod, n.iter = 1, smooth = TRUE, store.smooth = TRUE)
{ 
  y = matrix(y, nr = 1)
  nt = dim(y)[2]
  psi = list();
  psi$V = 1
  psi$F = array(u, c(1,1,nt))
  if(mod == "M101") psi$F[1,1,] = u[1:nt]
  psi$U = array(rbind(1,u[1:nt]), c(1,2,nt))
  xt = x[2,]
  x = matrix(x[1,], nr = 1)
  theta.list = list(beta = theta[1:2], phi = theta[3], sigma2s = theta[4], sigma2m = theta[5])
  prior = list(b0 = c(900,5), B0 = 5*diag(2), phi0 = list(0.8), Phi0 = list(matrix(0.2)), as0 = 1, bs0 = 1, am0 = 1, bm0 = 1, m0 = 0)
  
  for(k in 1:n.iter)
  {
    # Sample smoothed states by FFBS?
    if(smooth) x.smooth = sample.states(y, x, theta.list, psi, prior, FALSE) else x.smooth = x
    
    # Sample theta from full conditional
    theta.list$beta = sample.beta(y, x.smooth, theta.list, psi, prior)
    theta.list$phi = sample.phi(y, x.smooth, theta.list, psi, prior, FALSE)$phi
    theta.list$sigma2m = sample.sigma2m(y, x.smooth, theta.list, psi, prior)
    theta.list$sigma2s = sample.sigma2s(y, x.smooth, theta.list, psi, prior, FALSE)
    theta = c(theta.list$beta, theta.list$phi, theta.list$sigma2m, theta.list$sigma2s)
    
    # Return smoothed states?
    if(store.smooth) x = rbind(x.smooth,xt) else x = rbind(x,xt)
  }
  return(list(state=x,theta=theta))
}