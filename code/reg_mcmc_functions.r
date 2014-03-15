source("dlm_ar_functions.r")

########## Gibbs sampling functions ####################
# reg.ar.mcmc runs an MCMC algorithm for a linear regression model where the error term follows an AR(p) process
# y is an n-length vector
# X is an n by q design matrix
# prior is a list with components b0, B0 (mean and covariance for normal prior on beta), phi0, Phi0 (mean and covariance for (truncated) normal prior on phi), v0, and d0 (shape and rate for IG prior on sigma2)
# initial (optional) is a list with components beta (q-length vector), phi (p-length vector), and sigma2 (scalar)
# mcmc.details (optional) a list with scalar components n.sims (total number of MCMC iterations), n.burn (burn in iterations), and n.thin (save every n.thin-th iteration after burn-in)
# steps (optional) a character vector list which parameters to sample in the MCMC
# progress a boolean indicating if a progress bar should be displayed
# print.iter a boolean indicating if the iteration number should be printed (only if progress = FALSE)

reg.ar.mcmc <- function(y, X, prior, initial, mcmc.details, steps, progress=TRUE, print.iter=FALSE) {
  n <- length(y)
  p <- length(prior$phi0)
  q <- ncol(X)
  
  # Deal with missing arguments
  if (missing(initial)) {
    psi <- list(beta=rep(0,q), phi=rep(0,p), sigma2s=1)
  } else {
    psi <- initial
  }
  if (missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0  
  } else {
    n.thin = mcmc.details$n.thin
    n.sims = mcmc.details$n.sims
    n.burn = mcmc.details$n.burn
  }
  if (missing(steps)) {
    steps=c('beta','phi','sigma2s')
  } 
  
  # save structures
  n.iter <- (n.burn-n.sims)%/%n.thin
  keep.beta   <- matrix(NA, n.iter, q)
  keep.phi    <- matrix(NA, n.iter, p)
  keep.sigma2s <- rep(   NA, n.iter)
  accept.phi <- c()
  
  # Run mcmc
  if(progress) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if(print.iter & !progress & (i %% n.thin == 0))
    {
      print(i)
    } else if(progress) {
      setTxtProgressBar(pb,i)
    }

    if ('beta'  %in%steps) psi$beta   <- sample.ar.beta(  y, X, psi, prior)
    if ('phi'   %in%steps)
    {
      samp.phi <- sample.ar.phi(   y, X, psi, prior)
      psi$phi <- samp.phi$phi
      if(i > n.burn) accept.phi <- c(accept.phi, samp.phi$accept)
    }
    if ('sigma2'%in%steps) psi$sigma2s <- sample.ar.sigma2s(y, X, psi, prior)
    
    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.beta  [ii,] <- psi$beta
      keep.phi   [ii,] <- psi$phi
      keep.sigma2[ii ] <- psi$sigma2s
    }
  }
  
  return(list(beta=keep.beta, phi=keep.phi, sigma2=keep.sigma2s, accept.phi=accept.phi, mcmc.details=list(n.sims=n.sims,n.thin=n.thin,n.burn=n.burn), initial=initial))
}

# Functions to sample from full conditional distributions
# Inputs:
# y is an n-length vector
# X is an n by q design matrix
# psi has components beta, phi, sigma2 (current values)
# prior is a list with components b0, B0 (mean and covariance for beta), phi0, Phi0 (mean and covariance for phi), v0, and d0 (shape and rate for sigma2)

sample.ar.beta <- function(y, X, psi, prior)
{  
  tildey <- make.tildey(y,psi$phi)
  tildeX <- make.tildeX(X,psi$phi)
  
  B0.prec = solve(prior$B0)
  Bn <- solve(B0.prec+t(tildeX)%*%tildeX/psi$sigma2)
  bn <- Bn%*%(B0.prec%*%prior$b0 + t(tildeX)%*%tildey/psi$sigma2)
  
  return(rep(bn+t(chol(Bn.inv))%*%rnorm(length(psi$beta))),1)
}

sample.ar.phi <- function(y, X, psi, prior) 
{  
  # Calculate phin and Phin and sample phi
  E <- makeE(y,X,psi$beta,length(psi$phi))
  e <- makee(y,X,psi$beta,length(psi$phi))
  tau <- 1/psi$sigma2
  Phi0.prec = solve(prior$Phi0)
  Phin <- solve(Phi0.prec + tau*(t(E)%*%E))
  phin <- Phin%*%(Phi0.prec%*%prior$phi0 + tau*(t(E)%*%e))
  phi.p <- rep(phin + t(chol(Phin))%*%rnorm(length(psi$phi)),1)
  
  # Reject phi.p until stationary
  accept <- c()
  while(!is.stationary(phi.p))
  {
    phi.p <- mvrnorm(1, phin, solve(Phin))
    accept = c(accept,FALSE)
  }
  accept = c(accept, TRUE)
  
  # Perform MH step
  logMH <- Psi(y,X,psi, phi.p) - Psi(y,X,psi, psi$phi)
  if (log(runif(1)) < logMH) psi$phi <- phi.p
  return(list(phi=psi$phi,accept=accept))
}

sample.ar.sigma2s <- function(y, X, psi, prior)
{
  return(1/rgamma(1,prior$v0+length(y)/2, prior$d0+d1(y,X,psi)/2))
}


save.iteration <- function(current.iter, n.burn, n.thin) {
  if (current.iter<n.burn | ((current.iter-n.burn)%%n.thin)!=0) {
    return(0)
  } else {
    return((current.iter-n.burn)/n.thin)
  }
}

################# Utility functions #####################

makeOmega <- function(phi, theta=NULL) {
  m <- length(phi)
  G <- makeG(phi)
  if (is.null(theta)) {
    f <- rep(1,m)
  } else {
    stopifnot(length(phi)==length(theta))
    f <- makef(theta)
  }
  vv <- solve(diag(m^2)-G%x%G)%*%as.numeric(f%*%t(f))
  nr <- sqrt(length(vv))
  return(matrix(vv,nr,nr))
}

makeSigmap <- function(phi) return( makeOmega(phi) )

makeQ <- function(phi) return(chol(makeSigmap(phi))) 

makee <- function(y,X,beta, p, n=length(y)) return(as.numeric(y-X%*%beta)[seq(p+1,n)])

makeE <- function(y,X,beta, p, n=length(y)) { 
  e <- makee(y,X,beta,0) # get entire error vector 
  E <- matrix(NA, n-p, p)
  for (i in (p+1):n) {
    E[i-p,] <- e[seq(i-1,i-p)] # lag 1 up to lag p errors
  }
  return(E)
}

makey1 <- function(y,p) return(y[1:p]) 

makeX1 <- function(X, p) return(X[1:p,])

Psi <- function(y,X,psi, phi=NULL, log=T) {
  if (!is.null(phi)) psi$phi <- phi
  
  Sigmap <- makeSigmap(psi$phi)
  Sigmap.inv <- solve(Sigmap)
  p <- length(psi$phi)
  y1 <- makey1(y,p)
  X1 <- makeX1(X,p)
  yXb <- y1 - X1%*%psi$beta
  
  logPsi <- -( log(det(Sigmap))+t(yXb) %*% Sigmap.inv %*%yXb/psi$sigma2 )/2
  
  return(ifelse(log, logPsi, exp(logPsi)))
}

make.tildey1 <- function(y, phi, p=length(phi)) {
  Q <- makeQ(phi)
  y1 <- makey1(y,p)
  return(solve(Q)%*%y1)
}

make.tildey2 <- function(y, phi, p=length(phi)) {
  n <- length(y)
  tildey2 <- rep(NA,n-p)
  for (i in 1:(n-p))
    tildey2[i] <- sum(y[p+i],-phi*y[seq(p+i-1,i)])
  return(tildey2)
}

make.tildey <- function(y, phi, p=length(phi)) return(c(make.tildey1(y,phi,p), make.tildey2(y,phi,p)))


make.tildeX1 <- function(X, phi, p=length(phi)) {
  Q <- makeQ(phi)
  X1 <- makeX1(X,p)
  return(solve(Q)%*%X1)
}

make.tildeX2 <- function(X, phi, p=length(phi)) {
  n <- nrow(X)
  q <- ncol(X)
  tildeX2 <- matrix(NA, n-p, q)
  for (i in 1:(n-p)) {
    for (j in 1:q) {
      tildeX2[i,j] <- sum(X[p+i,j],-phi*X[seq(p+i-1,i),j])
    }
  }
  return(tildeX2)
}


make.tildeX <- function(X, phi, p=length(phi)) return(rbind(make.tildeX1(X,phi,p), make.tildeX2(X,phi,p)))

d1 <- function(y,X,psi) {
  tildey <- make.tildey(y, psi$phi)
  tildeX <- make.tildeX(X, psi$phi)
  
  return(sum((tildey-tildeX%*%psi$beta)^2))
}
