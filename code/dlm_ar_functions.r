# Functions to help construct parameters of a DLM formulation of regression with AR(p) errors and measurement noise

makeG <- function(phi) {
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

bdiag <- function(mylist)
{
  if(length(mylist) > 0)
  {
    ps = sapply(mylist,dim)
    nrow = sum(ps[1,])
    ncol = sum(ps[2,])
    C0 = matrix(0,nr=nrow,nc=ncol)
    row.index = 1; col.index = 1
    for(i in 1:length(mylist))
    {
      if(all(ps[,i] != 0))
      {
        C0[row.index:(row.index + ps[1,i] - 1),col.index:(col.index + ps[2,i] - 1)] = mylist[[i]]
        row.index = row.index + ps[1,i]
        col.index = col.index + ps[2,i]
      }
    }
    return(C0)
  } else { return(matrix(nr=0,nc=0)) }
}

# Functions to build a DLM model using R package 'dlm'
# par = c(logit(phi1), ..., logit(phip), log(sigma2s), log(sigma2m))

build_M101 <- function(par,V,U,nd=NULL) # assume all error components are AR(1)
{
  require(dlm)
  d = dim(U)[2]
  lp = length(par)
  p = (lp - 1) / 2
  stopifnot(p %% 1 == 0)
   
  FF = matrix(rep(0,d+p),nr=1)
  JFF = matrix(1:(d+p),nr=1)
  if(p == 1) U.add = matrix(U[1,d,]) else U.add = t(U[1,(d-p+1):d,])
  X = cbind(t(U[1,,]),U.add)
  dd = dim(X)[2]
  V=exp(par[lp])*V
  GG = diag(dd)
  for(j in 1:p) GG[dd-p+j,dd-p+j] = unlogit(par[j],-1,1)
  W = 0*diag(dd)
  for(j in 1:p) W[dd-p+j,dd-p+j] = exp(par[p+j])
  m0 = rep(0,dd)
  C0 = 1e6*diag(dd)
  for(j in 1:p) C0[dd-p+j,dd-p+j] = exp(par[p+j])/(1-unlogit(par[j],-1,1)^2)
  mod = dlm(FF=FF,V=V,GG=GG,W=W,m0=m0,C0=C0,JFF=JFF,X=X)
  return(mod)
}

build_M101_same <- function(par,V,U,nd) # assume all error components are AR(1)
{
  require(dlm)
  d = dim(U)[2]
  stopifnot(length(par) == 3 & nd <= d)
  
  FF = matrix(rep(0,d+nd),nr=1)
  JFF = matrix(1:(d+nd),nr=1)
  if(nd == 1) U.add = matrix(U[1,d,]) else U.add = t(U[1,(d-nd+1):d,])
  X = cbind(t(U[1,,]),U.add)
  dd = dim(X)[2]
  V=exp(par[3])*V
  GG = diag(dd)
  for(j in 1:nd) GG[dd-nd+j,dd-nd+j] = unlogit(par[1],-1,1)
  W = 0*diag(dd)
  for(j in 1:nd) W[dd-nd+j,dd-nd+j] = exp(par[2])
  m0 = rep(0,dd)
  C0 = 1e6*diag(dd)
  for(j in 1:nd) C0[dd-nd+j,dd-nd+j] = exp(par[2])/(1-unlogit(par[1],-1,1)^2)
  mod = dlm(FF=FF,V=V,GG=GG,W=W,m0=m0,C0=C0,JFF=JFF,X=X)
  return(mod)
}

build_M011 <- function(par,V,U,nd=NULL)
{
  require(dlm)
  d = dim(U)[2]
  
  FF = matrix(rep(0,d+1),nr=1)
  JFF = matrix(1:(d+1),nr=1)
  X = cbind(t(U[1,,]),1)
  V=exp(par[3])*V
  GG = diag(dim(X)[2]); GG[dim(X)[2],dim(X)[2]] = unlogit(par[1],-1,1)
  W = 0*diag(dim(X)[2]); W[dim(X)[2],dim(X)[2]] = exp(par[2])
  m0 = rep(0,dim(X)[2])
  C0 = 1e6*diag(dim(X)[2]); C0[dim(X)[2],dim(X)[2]] = exp(par[2])/(1-unlogit(par[1],-1,1)^2)
  mod = dlm(FF=FF,V=V,GG=GG,W=W,m0=m0,C0=C0,JFF=JFF,X=X)
  return(mod)
}

###################
# Utility functions
###################

logit <- function(x,a,b)
{
  u = (x - a) / (b - a)
  return(log(u / (1 - u)))
}

unlogit <- function(u,a,b)
{
  x = exp(u) / (1 + exp(u))
  return(x*(b-a) + a)
}

d.unlogit <- function(x,a,b)
{
  ((b-a)*exp(x)*(1+2*exp(x))) / (1+exp(x))^2
}