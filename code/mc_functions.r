# Function to simulate nt + 1 states (x_t) and nt observations (y_t) from a dlm with common state and observations variance factor sigma^2
# y_t = Fx_t + v_t, v_t iid N(0,sigma^2*V)
# x_t = Gx_t-1 + w_t, w_t iid N(0,sigma^2*W)
dlm.sim <- function(nt, F, G, V, W, sigma)
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

# Function to explicitly calculate sufficient statistics for the filtered distribution of the states, sigma^2, and one-step ahead predictions using Kalman-filter-like equations
# Prior distributions: sigma^2 ~ IG(a0,b0) and x_0|sigma^2 ~ N(m0,sigma^2*C0)
dlm.post <- function(y, F, G, V, W, a0, b0, m0, C0)
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
    a[i+1] = a[i] + q/2
    b[i+1] = b[i] + (1/2)%*%e%*%Qinv%*%e
  }
  return(list(m = m, C = C, a = a, b = b, f = f, Q = Q, A = A, R = R))
}

# Function to calculate the log marginal likelihood of the data for a DLM where both the states and observations are univariate and the states and observations have common variance factor sigma^2
# Assume all arguments are vectors; f and Q should have length equal to the data y while a and b should have length one greater than y (because they include prior values a0 and b0)
dlm.lmarglik <- function(y, f, Q, a, b)
{
  stopifnot(length(y) == length(f) & length(y) == length(Q) & length(y) == length(a) - 1 & length(a) == length(b))
  nt = length(y)
  log.marglik = 0
  for(i in 1:nt) log.marglik = log.marglik + dt((y[i]-f[i])/sqrt(Q[i]*b[i]/a[i]), df = 2*a[i], log=TRUE) - log(sqrt(Q[i]*b[i]/a[i]))
  return(log.marglik)
}

# Function to approximate the log marginal likelihood using particle filtering given a list returned by rm_pf()
pf.lmarglik <- function(out)
{
  nt = dim(out$increment)[2]
  log.marglik <- 0
  for(i in 1:nt) log.marglik = log.marglik + log(sum(exp(out$increment[,i])*out$weight[,i]))
  return(log.marglik)
}

# Function to calculate posterior model probabilities given a set of log marginal likelihoods and prior model probabilities
postModProbs <- function(lmarglik, priorModProbs)
{
  stopifnot(length(lmarglik) == length(priorModProbs))
  postModProbs.log = lmarglik + log(priorModProbs) - log(sum(exp(lmarglik)*priorModProbs))
  return(exp(postModProbs.log))
}