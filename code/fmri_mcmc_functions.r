source("dlm_ar_functions.r")

fmri_mcmc <- function(y, psi, prior, initial, mcmc.details, steps, progress=TRUE, print.iter=FALSE) {
  # Initial values of x and theta
  if(missing(initial)) {
    # Set up initial values    
  } else {
    x = initial$x
    theta = initial$theta
  }
  
  # Check dimensions
  check = check.dim(y, x, theta, psi, prior)
  y = check$y
  x = check$x
  theta = check$theta
  psi = check$psi
  prior = check$prior
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  d = dim(psi$U)[2]
  
  # MCMC details
  if(missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0
  } else {
    n.thin = mcmc.details$n.thin
    n.sims = mcmc.details$n.sims
    n.burn = mcmc.details$n.burn
  }
  if (missing(steps)) {
    steps=c('beta','sigma2m','phi','sigma2s','x')
  } 
  
  # save structures
  n.iter <- (n.sims - n.burn)%/%n.thin
  keep.beta <- matrix(NA, n.iter, d)
  keep.sigma2m <- rep(NA, n.iter)
  keep.phi <- matrix(NA, n.iter, p)
  keep.sigma2s <- rep(NA, n.iter)
  keep.x  <- array(NA, c(n.iter, p, nt + 1))
  accept.phi <- 0
  
  # Run mcmc
  if(progress & !print.iter) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if(print.iter & (i %% n.thin == 0))
    {
      print(i)
    } else if(progress) {
      setTxtProgressBar(pb,i)
    }
    
    if('beta' %in% steps) theta$beta = sample.beta(y, x, theta, psi, prior)
    if('sigma2m' %in% steps) theta$sigma2m = sample.sigma2m(y, x, theta, psi, prior)
    if('phi' %in% steps){
      samp.phi = sample.phi(y, x, theta, psi, prior)
      theta$phi = samp.phi$phi
      accept.phi = accept.phi + samp.phi$accept
    }
    if('sigma2s' %in% steps) theta$sigma2s = sample.sigma2m(y, x, theta, psi, prior)
    if('x' %in% steps) x = sample.states(y, x, theta, psi, prior)
    
    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.beta[ii,] = theta$beta
      keep.sigma2m[ii] = theta$sigma2m
      keep.phi[ii,] = theta$phi
      keep.sigma2s[ii] = theta$sigma2s
      keep.x[ii,,] = x
    }
  }
  
  return(list(beta = keep.beta, sigma2m = keep.sigma2m, phi = keep.phi, sigma2s = keep.sigma2s, x = keep.x, accept.phi=accept.phi, mcmc.details=list(n.sims=n.sims,n.thin=n.thin,n.burn=n.burn)))
}

# Functions to sample from full conditional distributions
# Inputs:
# y is a q by nt matrix of data
# x is a p by nt + 1 matrix of states
# theta is d + p + 2 vector of fixed parameters
# psi is a list with components V q by q covariance matrix, U a q by d by nt array of beta covariates, and F a q by p by nt array of state covariates
# prior is a list with components b0, B0 (mean and covariance on normal prior for beta), phi0, Phi0 (mean and covariance on truncated normal prior for phi), am0, bm0 (shape and rate for inverse-gamma prior on sigma2m), as0, bs0 (shape and rate for inverse-gamma prior on sigma2s), and m0 (mean on normal prior for initial state; covariance matrix is constructed using makeC0() to ensure initial state comes from a stationary distribution)

sample.beta <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  d = dim(psi$U)[2]
  Vinv = solve(psi$V)
  
  # Subtract states from data
  e = makee(y, psi$F, x)
  
  # Calculate UVU^-1 and UVe
  UVU <- matrix(0,nr=d,nc=d)
  UVe <- rep(0,d)
  for(k in 1:nt)
  {
    Uk = matrix(psi$U[,,k],nr=q,nc=d)
    UVU = UVU + t(Uk)%*%Vinv%*%Uk
    UVe = UVe + t(Uk)%*%Vinv%*%e[,k]
  }
    
  # If sigma2m = 0, betas are completely determined given the states
  if(theta$sigma2m == 0)
  {    
    return(solve(UVU)%*%UVe)
  } else { # sample beta from full conditional (normal distribution)
    B0.prec = solve(prior$B0)
    Bn = solve((1/theta$sigma2m)*UVU + B0.prec)
    bn = Bn%*%((1/theta$sigma2m)*UVe + B0.prec%*%prior$b0)
    return(rep(bn + t(chol(Bn))%*%rnorm(d,0,1),1))
  }
}

sample.sigma2m <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  Vinv = solve(psi$V)

  # Substract states and inputs from data
  E = makeE(y, psi$U, theta$beta, psi$F, x)
  
  # Calculate SSy
  SSy = 0
  for(i in 1:nt) SSy = SSy + t(E[,i])%*%Vinv%*%E[,i]
  
  # Sample from inverse gamma
  amn = nt*q/2 + prior$am0
  bmn =  SSy/2 + prior$bm0
  return(1/rgamma(1,amn,bmn))
}

sample.phi <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(x)[1]

  # Make X
  X = makeXtilde(x, theta$phi)
  x1 = x[1,2:(nt+1)]
  XX = t(X)%*%X
  Xx1 = t(X)%*%x1
  
  # Calculate phin and Phin and sample phi
  Phi0.prec = solve(prior$Phi0)
  Phin = solve((1/theta$sigma2s)*XX + Phi0.prec)
  phin = Phin%*%((1/theta$sigma2s)*Xx1 + Phi0.prec%*%prior$phi0)
  phi.p = rep(phin + t(chol(Phin))%*%rnorm(p),1) 
  while(!is.stationary(phi.p)) phi.p = rep(t(chol(Phin))%*%rnorm(p,phin,1),1)
  
  # Perform MH step
  accept = FALSE
  logMH <- Psi(x[,1], prior$m0, phi.p, theta$sigma2s) - Psi(x[,1], prior$m0, theta$phi, theta$sigma2s)
  if (log(runif(1)) < logMH)
  {
    theta$phi <- phi.p
    accept = TRUE
  }
  return(list(phi=theta$phi,accept=accept))
}

sample.sigma2s <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(X)[1]
  q = dim(y)[1]
  C0 = makeC0(phi)
  
  # Calculate SSx
  X = makeXtilde(x, theta$phi)
  e = x[1,2:(nt+1)] - X%*%theta$phi
  SSx = t(e)%*%e
  
  # Sample from inverse gamma
  asn = (p/2)*(nt + 1) + prior$as0
  bsn = (SSx + t(x[,1]-psi$m0)%*%solve(C0)%*%(x[,1]-psi$m0))/2 + prior$bs0
  return(1/rgamma(1,asn,bsn))
}

sample.states <- function(y, x, theta, psi, prior)
{
  # Extract dlm components
  p = length(theta$phi)
  G = makeG(theta$phi)
  V = theta$sigma2m*psi$V
  sf = c(1,rep(0,p-1))
  W = theta$sigma2s*sf%*%t(sf)
  C0 = makeC0(theta$phi)
  return(ffbs(y, psi$U, theta$beta, psi$F, G, V, W, prior$m0, C0))
}

###################
# Utility Functions
###################

save.iteration <- function(current.iter, n.burn, n.thin) {
  if (current.iter<n.burn | ((current.iter-n.burn)%%n.thin)!=0) {
    return(0)
  } else {
    return((current.iter-n.burn)/n.thin)
  }
}

# check.dim() checks that y, x, and components of theta, psi, and prior are of correct dimensions
check.dim <- function(y, x, theta, psi, prior)
{
  # Data and states
  if(is.null(dim(y))) y = matrix(y, nr=1)
  y = as.matrix(y)
  q = dim(y)[1]
  nt = dim(y)[2]
  if(is.null(dim(x))) x = matrix(x, nr=1)
  x = as.matrix(x)
  p = dim(x)[1]
  stopifnot(dim(x)[2] == nt + 1)
  
  # Unknown parameters
  if(is.null(theta$beta)) beta = 0 else beta = rep(theta$beta, 1)
  d = length(beta)
  if(is.null(theta$phi)) phi = rep(0,p) else phi = rep(theta$phi, 1)
  stopifnot(length(phi) == p)
  if(is.null(theta$sigma2m)) sigma2m = 0 else sigma2m = theta$sigma2m[1]
  if(is.null(theta$sigma2s)) sigma2s = 0 else sigma2s = theta$sigma2s[1]
  theta = list(beta = beta, phi = phi, sigma2m = sigma2m, sigma2s = sigma2s)
  
  # Known parameters
  if(is.null(psi$V)) V = as.matrix(diag(q)) else V = as.matrix(psi$V)
  stopifnot(dim(V)[1] == q & dim(V)[2] == q)
  if(is.null(psi$U)) U = array(0,c(q,d,nt))
  if(is.matrix(U))
  {
    stopifnot(dim(U)[1] == q & dim(U)[2] == d)
    U = array(U,dim=c(q,d,nt))
  } else {
    stopifnot(length(dim(U) == 3) & dim(U)[1] == q & dim(U)[2] == d & dim(U)[3] == nt)
  }
  if(is.null(F)) F = array(0,c(q,p,nt))
  if(is.matrix(F))
  {
    stopifnot(dim(F)[1] == q & dim(F)[2] == p)
    F = array(F,dim=c(q,p,nt))
  } else {
    stopifnot(length(dim(F) == 3) & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  }
  psi = list(V = V, U = U, F = F)
    
  # Prior hyperparameters
  checked.prior = list()
  checked.prior$m0 = rep(prior$m0,1)
  stopifnot(length(checked.prior$m0) == p)
  if(!is.null(prior$b0) & !is.null(prior$B0))
  {
    checked.prior$b0 = rep(prior$b0,1)
    checked.prior$B0 = as.matrix(prior$B0)
    stopifnot(length(checked.prior$b0) == d & dim(checked.prior$B0)[1] == d & dim(checked.prior$B0)[2] == d)
  }
  if(!is.null(prior$phi0) & !is.null(prior$Phi0))
  {
    checked.prior$phi0 = rep(prior$phi0,1)
    checked.prior$Phi0 = as.matrix(prior$Phi0)
    stopifnot(length(checked.prior$phi0) == p & dim(checked.prior$Phi0)[1] == p & dim(checked.prior$Phi0)[2] == p)
  }
  if(!is.null(prior$am0) & !is.null(prior$bm0))
  {
    checked.prior$am0 = rep(prior$am0,1)[1]
    checked.prior$bm0 = rep(prior$bm0,1)[1]
  }
  if(!is.null(prior$as0) & !is.null(prior$bs0))
  {
    checked.prior$as0 = rep(prior$as0,1)[1]
    checked.prior$bs0 = rep(prior$bs0,1)[1]
  }
  prior = checked.prior
     
  return(list(y=y,x=x,theta=theta,psi=psi,prior=prior))
}

ffbs <- function(y, U, beta, F, G, V, W, m0, C0)
{
  nt = dim(y)[2]
  p = length(m0)
  q = dim(y)[1]
  d = length(beta)
  stopifnot((dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(C0)[1] == dim(C0)[2]))
  stopifnot((dim(C0)[1] == p) & (dim(G)[1] == p) & (dim(F)[1] == q) & (dim(F)[2] == p) & (dim(F)[3] == nt) & (dim(V)[1] == q) & (dim(W)[1] == p))
  stopifnot(dim(U)[1] == q & dim(U)[2] == d & dim(U)[3] == nt)
            
  # Initialize output arrays
  m = matrix(NA, nr = p, nc = nt + 1)
  C = array(NA, dim = c(p, p, nt + 1))
  f = matrix(NA, nr = q, nc = nt)
  Q = array(NA, dim = c(q, q, nt))
  A = matrix(NA, nr = p, nc = nt)
  R = array(NA, dim = c(p, p, nt))
  m[1,] = m0; C[,,1] = C0
  
  # Forward-filtering
  for(i in 1:nt)
  {
    Ui = matrix(U[,,i],nr=q,nc=d); Fi = matrix(F[,,i],nr=q,nc=p)
    A[,i] = G%*%m[,i]; R[,,i] = G%*%C[,,i]%*%t(G) + W
    f[,i] = Ui%*%beta + Fi%*%A[,i]; Q[,,i] = Fi%*%R[,,i]%*%t(Fi) + V
    e = y[,i] - f[,i]; Qinv = solve(Q[,,i])
    RFQinv = R[,,i]%*%t(Fi)%*%Qinv
    m[,i+1] = A[,i] + RFQinv%*%e
    C[,,i+1] = R[,,i] - RFQinv%*%Fi%*%R[,,i]
  }
  
  # Backward-sampling
  x = matrix(NA, nr=p, nc=nt+1)
  x[,nt+1] = t(chol(C[,,nt+1]))%*%rnorm(p,m[,nt+1],1)
  for(i in nt:1)
  {
    CGRinv = C[,,i]%*%t(G)%*%solve(R[,,i])
    h = m[,i] + CGRinv%*%(x[,i+1] - A[,i])
    H = C[,,i] - CGRinv%*%G%*%C[,,i]
    x[,i] = t(chol(H))%*%rnorm(p, h, 1)
  }
  return(x)
}

makee <- function(y, F, x)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  stopifnot(dim(x)[2] == nt+1 & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  
  e = matrix(NA, nr=q, nc=nt)
  for(i in 1:nt) e[,i] = y[,i] - matrix(F[,,i],nr=q,nc=p)%*%x[,i+1]
  return(e)
}

makeE <- function(y, U, beta, F, x)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  d = length(beta)
  stopifnot(dim(U)[1] == q & dim(U)[2] == d & dim(U)[3] == nt)
  stopifnot(dim(x)[2] == nt+1 & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  
  E = matrix(NA, nr=q, nc=nt)
  for(i in 1:nt) E[,i] = y[,i] - matrix(U[,,i],nr=q,nc=d)%*%beta - matrix(F[,,i],nr=q,nc=p)%*%x[,i+1]
  return(E)
}

makeXtilde <- function(x, phi)
{
  nt = dim(x)[2] - 1
  p = dim(x)[1]
  stopifnot(p == length(phi))
  Xtilde = matrix(NA, nr = nt, nc = p)
  
  if(p == 1)
  {
    Xtilde[,p] = x[,1:nt] 
  } else {
    for(i in 1:nt) if(i == 1) Xtilde[i,] = makex0(x, phi) else Xtilde[i,] = c(x[1,i], Xtilde[i-1,1:(p-1)])
  }
  return(Xtilde)
}

makex0 <- function(x, phi)
{
  nt = dim(x)[2] - 1
  p = dim(x)[1]
  stopifnot(p == length(phi) & p > 1)
  
  x0 = rep(NA,p)
  x0[1] = x[1,1]
  xj = x[,1]
  for(i in 1:(p-1))
  {
    xi = xj[p] / phi[p]
    if(i < p - 1)
    {
      xj.new = rep(NA,p-i-1)
      for(j in p:(i+2)) xj.new[j] = xj[j-1] - phi[j-1]*xi
      xj = xj.new
    }
    x0[i+1] = xi
  }
  return(x0)
}

is.stationary <- function(phi)
{
  return(all(Mod(polyroot(c(1,-phi)))>1))
}

Psi <- function(x0, m0, phi, sigma2s, log=TRUE) 
{
  p = length(x0)
  stopifnot((p == length(m0)) & (p == length(phi)))
  
  C0 <- makeC0(phi)
  C0.prec <- solve(C0)
  logPsi <- -.5*(log(det(C0))-(1/sigma2s)*t(x0-m0)%*%C0.prec%*%(x0-m0))/2
  return(ifelse(log, logPsi, exp(logPsi)))
}