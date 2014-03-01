# y is a q by nt matrix of data
# x is a p by nt + 1 matrix of states
# theta is d + p + 2 vector of fixed parameters
# psi is a list with components V q by q covariance matrix, U a q by d by nt array of beta covariates, and F a q by p by nt array of state covariates
# prior is a list with components b0, B0 (mean and covariance on normal prior for beta), phi0, Phi0 (mean and covariance on truncated normal prior for phi), am0, bm0 (shape and rate for inverse-gamma prior on sigma2m), as0, bs0 (shape and rate for inverse-gamma prior on sigma2s), m0, and C0 (mean and covariance on normal prior for initial state)

sample.beta <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  Vinv = solve(psi$V)
  
  # Subtract states from data
  e = makee(y, psi$F, x)
  
  # Calculate UVU^-1 and UVe
  UVU <- matrix(0,nr=d,nc=d)
  UVe <- rep(0,d)
  for(k in 1:nt)
  {
    UVU = UVU + t(psi$U[,,k])%*%Vinv%*%psi$U[,,k]
    UVe = UVe + t(psi$U[,,k])%*%Vinv%*%e[,k]
  }
  UVUinv = solve(UVU)
    
  # If sigma2m = 0, betas are completely determined given the states
  if(theta$sigma2m == 0)
  {    
    return(UVUinv%*%UVe)
  } else { # sample beta from full conditional (normal distribution)
    B0.prec = solve(prior$B0)
    Bn = solve((1/theta$sigma2m)*UVUinv + B0.prec)
    bn = Bn%*%((1/theta$sigma2m)*UVe + B0.prec%*%prior$b0)
    return(rep(t(chol(Bn))%*%rnorm(q,bn,1),1))
  }
}

sample.sigma2m <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  p = dim(X)[1]
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
  XXinv = solve(t(X)%*%X)
  Xx1 = t(X)%*%x1
  
  # Calculate phin and Phin and sample phi
  Phi0.prec = solve(prior$Phi0)
  Phin = solve((1/theta$sigma2s)*XXinv + Phi0.prec)
  phin = Phin%*%((1/theta$sigma2s)*Xx1 + Phi0.prec%*%prior$phi0)
  phi.p = rep(t(chol(Phin))%*%rnorm(p,phin,1),1) 
  while(!is.stationary(phi.p)) phi.p = rep(t(chol(Phi_n))%*%rnorm(p,phi_n,1),1)
  
  # Perform MH step
  logMH <- Psi(x[,1], psi$m0, phi.p, theta$sigma2s) - Psi(x[,1], psi$m0, theta$phi, theta$sigma2s)
  if (log(runif(1)) < logMH) theta$phi <- phi.p
  return(theta$phi)
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

###################
# Utility Functions
###################

# check.dim() checks that y, x, and components of theta, psi, and prior are of correct dimensions
check.dim <- function(y, x, theta, psi, prior)
{
  # Data and states
  y = as.matrix(y)
  q = dim(y)[1]
  nt = dim(y)[2]
  x = as.matrix(x)
  p = dim(x)[1]
  stopifnot(dim(x)[2] == nt + 1)
  
  # Unknown parameters
  if(is.null(theta$beta)) beta = 0 else beta = rep(theta$beta, 1)
  d = length(beta)
  if(is.null(theta$phi) phi = rep(0,p) else phi = rep(theta$phi, 1)
  stopifnot(length(phi == p))
  if(is.null(theta$sigma2m)) sigma2m = 0 else sigma2m = theta$sigma2m[1]
  if(is.null(theta$sigma2s)) sigma2s = 0 else sigma2s = theta$sigma2s[1]
  theta = list(beta = beta, phi = phi, sigma2m = sigma2m, sigma2s)
  
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
  checked.prior$C0 = as.matrix(prior$C0)
  stopifnot(dim(checked.prior$C0)[1] == p & dim(checked.prior$C0)[2] == p)
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
  if(!is.null(prior$am0) & !is.null(bm0))
  {
    checked.prior$am0 = rep(prior$am0,1)[1]
    checked.prior$bm0 = rep(prior$bm0,1)[1]
  }
  if(!is.null(prior$as0) & !is.null(bs0))
  {
    checked.prior$as0 = rep(prior$as0,1)[1]
    checked.prior$bs0 = rep(prior$bs0,1)[1]
  }
  prior = checked.prior
     
  return(list(y=y,x=x,theta=theta,psi=psi,prior=prior))
}

makee <- function(y, F, x)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  stopifnot(dim(x)[2] == nt+1 & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  
  e = matrix(NA, nr=q, nc=nt)
  for(i in 1:nt) e[,i] = y[,i] - F[,,i]%*%x[,i+1]
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
  for(i in 1:nt) E[,i] = y[,i] - U[,,i]%*%beta - F[,,i]%*%x[,i+1]
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
    Xtilde[,p] == x[,1:nt] 
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

makeG <- function(phi)
{
  m = length(phi)
  G <- diag(0,m)
  G[,1] <- phi
  G[-m,-1] <- diag(1,m-1)
  return(G)
}

is.stationary <- function(phi)
{
  return(all(Mod(polyroot(c(1,-phi)))>1))
}

makeC0 <- function(phi)
{
  p <- length(phi)
  G <- makeG(phi)
  f <- c(1,rep(0,p-1))
  vv <- solve(diag(p^2)-G%x%G)%*%as.numeric(f%*%t(f))
  nr <- sqrt(length(vv))
  return(matrix(vv,nr,nr))
}

Psi <- function(x0, m0, phi, sigma2s, log=TRUE) 
{
  p = length(x0)
  stopifnot((p == length(m0)) & (p == length(phi)))
  
  C0 <- makeC0(phi)
  C0.prec <- solve(C0)
  logPsi <- -(log(det(C0))-(1/sigma2s)*t(x0-m0)%*%C0.prec%*%(x0-m0))/2
  return(ifelse(log, logPsi, exp(logPsi)))
}