source("dlm_ar_functions.r")

########## Gibbs sampling functions ####################
# dlm.mcmc runs an MCMC algorithm for a DLM where the state follows an AR(p) process
# Inputs:
# y is an nt-length vector or q by nt matrix of data
# psi is a list with components V a q by q covariance matrix, U a q by d by nt array of beta covariates, and F a q by p by nt array of state covariates
# prior is a list with components b0, B0 (mean and covariance on normal prior for beta), phi0, Phi0 (lists with elements phi0[[i]] and Phi0[[i]] giving the mean and covariance on truncated normal prior for phi[[i]]), am0, bm0 (shape and rate for inverse-gamma prior on sigma2m), as0, bs0 (vectors with elements as0[i] and bs0[i] giving the shape and rate for inverse-gamma prior on sigma2s[i]), and m0 (mean on normal prior for initial state; covariance matrix is constructed using makeC0() to ensure initial state comes from a stationary distribution)
# initial (optional) is a list with components x a p by nt + 1 state matrix and theta a list with components beta (d-length vector), sigma2m (scalar), phi (P-length list of phi vectors), and sigma2s (P-length vector of white noise variance of AR components)
# same is a boolean indicating whether all AR components in the model are assumed to have the same model parameters
# mcmc.details (optional) a list with scalar components n.sims (total number of MCMC iterations), n.burn (burn in iterations), and n.thin (save every n.thin-th iteration after burn-in)
# steps (optional) a character vector list which parameters to sample in the MCMC
# progress a boolean indicating if a progress bar should be displayed
# print.iter a boolean indicating if the iteration number should be printed (only if progress = FALSE)

dlm.ar.mcmc <- function(y, psi, prior, initial, same = FALSE, mcmc.details, steps, progress=TRUE, print.iter=FALSE) {
  # Get dimensions of data and parameters
  if(!is.matrix(y)) y = matrix(y, nr = 1)
  if(!is.list(prior$phi0)) prior$phi0 = list(prior$phi0)
  nt = dim(y)[2]
  q = dim(y)[1]
  p = sapply(prior$phi0, length)
  P = length(p)
  d = length(prior$b0)
  
  # Set up initial values
  if(missing(initial)) {
    x = matrix(rnorm(sum(p)*(nt+1)),nr=sum(p),nc=nt+1)
    phi.temp = list(); length(phi.temp) = P
    for(i in 1:P) phi.temp[[i]] = rep(0,p[i])
    theta = list(beta = rep(0,d), sigma2m = 1, phi = phi.temp, sigma2s = rep(1,P))
    rm(phi.temp)
  } else {
    x = initial$x
    theta = initial$theta
    check = check.dim(y, x, theta, psi, prior)
    x = check$x
    theta = check$theta
    psi = check$psi
    prior = check$prior
    rm(check)
  }
  
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
  keep.phi = list(); length(keep.phi) = P
  for(i in 1:P) keep.phi[[i]] <- matrix(NA, n.iter, p[i])
  keep.sigma2s <- matrix(NA, n.iter, P)
  keep.x  <- array(NA, c(n.iter, sum(p), nt + 1))
  accept.phi <- list(); length(accept.phi) = P
  
  # Run mcmc
  if(progress) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if(print.iter & !progress & (i %% n.thin == 0))
    {
      print(i)
    } else if(progress) {
      setTxtProgressBar(pb,i)
    }
    
    if('beta' %in% steps) theta$beta = sample.beta(y, x, theta, psi, prior)
    if('sigma2m' %in% steps) theta$sigma2m = sample.sigma2m(y, x, theta, psi, prior)
    if('phi' %in% steps){
      samp.phi = sample.phi(y, x, theta, psi, prior)
      for(j in 1:P)
      {
        theta$phi[[j]] = samp.phi$phi[[j]]
        if(i > n.burn) accept.phi[[j]] = c(accept.phi[[j]], samp.phi$accept[[j]])
      }
    }
    if('sigma2s' %in% steps) theta$sigma2s = sample.sigma2s(y, x, theta, psi, prior)
    if('x' %in% steps) x = sample.states(y, x, theta, psi, prior)
    
    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.beta[ii,] = theta$beta
      keep.sigma2m[ii] = theta$sigma2m
      for(j in 1:P) keep.phi[[j]][ii,] = theta$phi[[j]]
      keep.sigma2s[ii,] = theta$sigma2s
      keep.x[ii,,] = x
    }
  }
  
  return(list(beta = keep.beta, sigma2m = keep.sigma2m, phi = keep.phi, sigma2s = keep.sigma2s, x = keep.x, accept.phi=accept.phi, mcmc.details=list(n.sims=n.sims,n.thin=n.thin,n.burn=n.burn), initial=initial))
}

# Functions to sample from full conditional distributions
# Inputs:
# y is a q by nt matrix of data
# x a p by nt + 1 state matrix 
# theta a list with components beta (d-length vector), sigma2m (scalar), phi (p-length vector), and sigma2s (scalar)
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
  e = makeE1(y, psi$F, x)
  
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
  E = makeE2(y, psi$U, theta$beta, psi$F, x)
  
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
  p = sapply(prior$phi0, length)
  P = length(p)
  accept = list(); length(accept) = P
  
  for(j in 1:P)
  {
    # Make X
    cj = ifelse(j == 1, 1, sum(p[1:(j-1)])+1)
    xj = x[cj:(cj+p[j]-1),]
    if(!is.matrix(xj)) xj = matrix(xj, nr = 1)
    X = makeXtilde(xj, theta$phi[[j]])
    x1 = x[cj,2:(nt+1)]
    XX = t(X)%*%X
    Xx1 = t(X)%*%x1
  
    # Calculate phin and Phin and sample phi
    Phi0.prec = solve(prior$Phi0[[j]])
    Phin = solve((1/theta$sigma2s[j])*XX + Phi0.prec)
    phin = Phin%*%((1/theta$sigma2s[j])*Xx1 + Phi0.prec%*%prior$phi0[[j]])
    phi.p = rep(phin + t(chol(Phin))%*%rnorm(p[j]),1)
  
    # Reject phi.p until stationary
    while(!is.stationary(phi.p))
    {
      phi.p = rep(phin + t(chol(Phin))%*%rnorm(p),1)
      accept[[j]] = c(accept[[j]], FALSE)
    }
    accept[[j]] = c(accept[[j]], TRUE)
  
    # Perform MH step
    logMH <- Psi.C0(x[cj:(cj+p[j]-1),1], prior$m0[cj:(cj+p[j]-1)], phi.p, theta$sigma2s[j]) - Psi.C0(x[cj:(cj+p[j]-1),1], prior$m0[cj:(cj+p[j]-1)], theta$phi[[j]], theta$sigma2s[j])
    if (log(runif(1)) < logMH)
    {
      theta$phi[[j]] <- phi.p
      accept[[j]] = c(accept[[j]], TRUE)
    } else {
      accept[[j]] = c(accept[[j]], FALSE)
    }
  }
  return(list(phi=theta$phi,accept=accept))
}

sample.sigma2s <- function(y, x, theta, psi, prior)
{
  nt = dim(y)[2]
  q = dim(y)[1]
  p = sapply(prior$phi0, length)
  P = length(p)

  for(j in 1:P)
  {
    # Calculate SSx
    C0 = makeC0(theta$phi[[j]])
    cj = ifelse(j == 1, 1, sum(p[1:(j-1)])+1)  
    xj = x[cj:(cj+p[j]-1),]
    if(!is.matrix(xj)) xj = matrix(xj, nr = 1)
    X = makeXtilde(xj, theta$phi[[j]])
    x1 = x[cj,2:(nt+1)]
    e = x1 - X%*%theta$phi[[j]]
    SSx = t(e)%*%e
  
    # Sample from inverse gamma
    asn = (p[j]/2)*(nt + 1) + prior$as0[j]
    bsn = (SSx + t(x[cj:(cj+p[j]-1),1]-prior$m0[cj:(cj+p[j]-1)])%*%solve(C0)%*%(x[cj:(cj+p[j]-1),1]-prior$m0[cj:(cj+p[j]-1)]))/2 + prior$bs0[j]
    theta$sigma2s[j] = 1/rgamma(1,asn,bsn)
  }
  return(theta$sigma2s)
}

sample.states <- function(y, x, theta, psi, prior)
{
  # Extract dlm components
  p = sapply(prior$phi0, length)
  P = length(p)
  G = makeG(theta$phi[[1]])
  sf = c(1,rep(0,p[1]-1))
  W = theta$sigma2s[1]*sf%*%t(sf)
  C0 = theta$sigma2s[1]*makeC0(theta$phi[[1]])
  if(P > 1)
  {
    for(j in 2:P)
    {
      G = bdiag(list(G,makeG(theta$phi[[j]])))
      sf = c(1,rep(0,p[j]-1))
      Wj = theta$sigma2s[j]*sf%*%t(sf)
      W = bdiag(list(W,Wj))
      C0 = bdiag(list(C0,theta$sigma2s[j]*makeC0(theta$phi[[j]])))
    }
  }
  V = theta$sigma2m*psi$V

  # Forward filtering backward sampling
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
  # Data assumed to be a q by nt matrix
  q = dim(y)[1]
  nt = dim(y)[2]
  
  # States assumed to be either a vector or a matrix
  if(!is.matrix(x)) x = matrix(x, nr=1)
  stopifnot(is.matrix(x))
  p.sum = dim(x)[1]
  stopifnot(dim(x)[2] == nt + 1)
  prior$m0 = rep(prior$m0,1)
  stopifnot(length(prior$m0) == p.sum)
  
  # Unknown parameters and priors
  # beta
  prior$b0 = rep(prior$b0,1)
  prior$B0 = as.matrix(prior$B0)
  stopifnot(length(prior$b0) == dim(prior$B0)[1] & dim(prior$B0)[1] == dim(prior$B0)[2])
  d = length(prior$b0)
  if(is.null(theta$beta)) theta$beta = rep(0,d) else theta$beta = rep(theta$beta, 1)
  stopifnot(length(theta$beta) == d)
  # sigma2m
  prior$am0 = rep(prior$am0,1)[1]
  prior$bm0 = rep(prior$bm0,1)[1]
  if(is.null(theta$sigma2m)) theta$sigma2m = 0 else theta$sigma2m = theta$sigma2m[1]
  # phi
  stopifnot(is.list(prior$phi0) & is.list(prior$Phi0))
  P = length(prior$phi0)
  stopifnot(P == length(prior$Phi0))
  p = rep(NA, P)
  for(j in 1:length(prior$phi0))
  {
    prior$phi0[[j]] = rep(prior$phi0[[j]],1)
    prior$Phi0[[j]] = as.matrix(prior$Phi0[[j]])
    p[j] = length(prior$phi0[[j]])
    stopifnot(dim(prior$Phi0[[j]])[1] == p[j] & dim(prior$Phi0[[j]])[2] == p[j])
  }
  stopifnot(sum(p) == p.sum)
  if(is.null(theta$phi))
  {
    theta$phi = list(); length(theta$phi) = P
    for(j in 1:P) theta$phi[[j]] = rep(0,p[j])
  } else if(!is.list(theta$phi)) {
    theta$phi = list(theta$phi)
    stopifnot(P == 1 & p[1] == length(theta$phi[[1]]))
  }
  # sigma2s
  stopifnot(length(prior$as0) == P & length(prior$bs0) == P)
  if(is.null(theta$sigma2s)) theta$sigma2s = rep(0, P) else theta$sigma2s = rep(theta$sigma2s, 1)
  stopifnot(length(theta$sigma2s) == P)
  
  # Known parameters
  if(is.null(psi$V)) psi$V = as.matrix(diag(q)) else psi$V = as.matrix(psi$V)
  stopifnot(dim(psi$V)[1] == q & dim(psi$V)[2] == q)
  if(is.null(psi$U)) psi$U = array(0,c(q,d,nt))
  if(is.matrix(psi$U)) 
  {
    stopifnot(dim(psi$U)[1] == q & dim(psi$U)[2] == d)
    psi$U = array(psi$U,dim=c(q,d,nt))
  } else {
    stopifnot(length(dim(psi$U) == 3) & dim(psi$U)[1] == q & dim(psi$U)[2] == d & dim(psi$U)[3] == nt)
  }
  if(is.null(psi$F)) psi$F = array(0,c(q,p.sum,nt))
  if(is.matrix(psi$F))
  {
    stopifnot(dim(psi$F)[1] == q & dim(psi$F)[2] == p.sum)
    psi$F = array(psi$F,dim=c(q,p,nt))
  } else {
    stopifnot(length(dim(psi$F) == 3) & dim(psi$F)[1] == q & dim(psi$F)[2] == p.sum & dim(psi$F)[3] == nt)
  } 
  return(list(x=x,theta=theta,psi=psi,prior=prior))
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
  m[,1] = m0; C[,,1] = C0

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
  oo = 1:p
  chol.C0 = try(chol(C[,,nt+1]),silent=TRUE)
  if(class(chol.C0) == "try-error")
  {
    chol.C0 = chol(C[,,nt+1],pivot=TRUE)
    oo = order(attr(chol.C0, "pivot"))
  }
  x[,nt+1] = m[,nt+1] + t(as.matrix(chol.C0)[,oo])%*%rnorm(p,0,1)
  for(i in nt:1)
  {
    CGRinv = C[,,i]%*%t(G)%*%solve(R[,,i])
    h = m[,i] + CGRinv%*%(x[,i+1] - A[,i])
    H = C[,,i] - CGRinv%*%G%*%C[,,i]
    
    oo = 1:p
    chol.H = try(chol(H),silent=TRUE)
    if(class(chol.H) == "try-error")
    {
      chol.H = chol(H,pivot=TRUE)
      oo = order(attr(chol.H, "pivot"))
    }
    x[,i] = h + t(as.matrix(chol.H)[,oo])%*%rnorm(p, 0, 1)
  }
  return(x)
}

makeE1 <- function(y, F, x)
{
  nt = dim(y)[2]
  p = dim(x)[1]
  q = dim(y)[1]
  stopifnot(dim(x)[2] == nt+1 & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  
  e = matrix(NA, nr=q, nc=nt)
  for(i in 1:nt) e[,i] = y[,i] - matrix(F[,,i],nr=q,nc=p)%*%x[,i+1]
  return(e)
}

makeE2 <- function(y, U, beta, F, x)
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

Psi.C0 <- function(x0, m0, phi, sigma2s, log=TRUE) 
{
  p = length(x0)
  stopifnot((p == length(m0)) & (p == length(phi)))
  
  C0 <- makeC0(phi)
  C0.prec <- solve(C0)
  logPsi <- -.5*log(det(C0))-.5*(1/sigma2s)*t(x0-m0)%*%C0.prec%*%(x0-m0)
  return(ifelse(log, logPsi, exp(logPsi)))
}

mvdnorm <- function(x, mu, Sigma, log = TRUE)
{
  d = length(x)
  Sigma = as.matrix(Sigma)
  stopifnot(length(mu) == d & dim(Sigma)[1] == d & dim(Sigma)[1] == d)
  
  log.lik = -(d/2)*log(2*pi) - .5*log(det(Sigma)) - .5*t(x - mu)%*%solve(Sigma)%*%(x - mu)
  if(log) return(log.lik) else return(exp(log.lik))
}