dr.sim <- function(nt, F, G, V, W, sigma)
{
  F =  as.matrix(F); G = as.matrix(G); V = as.matrix(V); W = as.matrix(W)
  stopifnot((dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(F)[1] == dim(V)[1]) & (dim(F)[2] == dim(G)[1]) & (dim(G)[1] == dim(W)[1]))
  p = dim(W)[1]; q = dim(V)[1]
  chol.W = chol(W); chol.V = chol(V)
  x0 = t(chol.W)%*%rnorm(p,0,sigma)
  x = matrix(NA, nr = nt + 1, nc = p)
  y = matrix(NA, nr = nt, nc = q)
  x[1,] = x0
  for(i in 1:nt)
  {
    x[i+1,] = G%*%x[i,] + t(chol.W)%*%rnorm(p,0,sigma)
    y[i,] = F%*%x[i+1,] + t(chol.V)%*%rnorm(q,0,sigma)
  }
  return(list(y = y, x = x))
}

dr.post <- function(y, F, G, V, W, a0, b0, m0, C0)
{
  y = as.matrix(y); m0 = as.matrix(m0); C0 = as.matrix(C0)
  F =  as.matrix(F); G = as.matrix(G); V = as.matrix(V); W = as.matrix(W)
  stopifnot((dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(F)[1] == dim(V)[1]) & (dim(F)[2] == dim(G)[1]) & (dim(G)[1] == dim(W)[1]))
  p = dim(W)[1]; q = dim(V)[1]
  nt = dim(y)[1]
  m = matrix(NA, nr = nt + 1, nc = p)
  C = array(NA, dim = c(p, p, nt + 1))
  a <- b <- rep(NA, nt + 1)
  f = matrix(NA, nr = nt, nc = q)
  Q = array(NA, dim = c(q, q, nt))
  m[1,] = m0; C[,,1] = C0; a[1] = a0; b[1] = b0
  for(i in 1:nt)
  {
    A = G%*%m[i,]; R = G%*%C[,,i]%*%t(G) + W
    f[i,] = F%*%A; Q[,,i] = F%*%R%*%t(F) + V  
    e = y[i,] - f[i,]; Qinv = solve(Q[,,i])
    RFQinv = R%*%t(F)%*%Qinv
    m[i+1,] = A + RFQinv%*%e
    C[,,i+1] = R - RFQinv%*%F%*%R
    a[i+1] = a[i] + 1/2
    b[i+1] = b[i] + (1/2)%*%e%*%Qinv%*%e
  }
  return(list(m = m, C = C, a = a, b = b, f = f, Q = Q))
}

# Assume all arguments are vectors of scalars
dr.prob <- function(y, f, Q, a, b)
{
  nt = length(y)
  log.marglik = 0
  for(i in 1:nt) log.marglik = log.marglik + dt.ns(y[i],f[i],sqrt(Q[i]*b[i]/a[i]), df = 2*a[i], log = TRUE)
  return(log.marglik)
}

dr.pf.prob <- function(out)
{
  nt = dim(out$increment)[2]
  log.marglik <- 0
  for(i in 1:nt) log.marglik = log.marglik + log(sum(exp(out$increment[,i])*out$weight[,i]))
  return(log.marglik)
}

###################
# Utility functions
###################

dt.ns <- function(x, mu, sigma, ...)
{
  dt((x - mu)/sigma, ...) / sigma
}