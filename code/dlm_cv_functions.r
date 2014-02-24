# Function to explicitly calculate sufficient statistics for the filtered distribution of the states, sigma^2, and one-step ahead predictions using Kalman-filter-like equations
# Prior distributions: sigma^2 ~ IG(a0,b0) and x_0|sigma^2 ~ N(m0,sigma^2*C0)
cv.post <- function(y, F, G, V, W, a0, b0, m0, C0)
{
  y = as.matrix(y)
  m0 = rep(m0,1)
  C0 = as.matrix(C0); G = as.matrix(G); V = as.matrix(V); W = as.matrix(W)
  stopifnot((dim(y)[1] == dim(V)[1]) & (dim(C0)[1] == dim(C0)[2]) & (dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(G)[1] == dim(W)[1]) & (dim(G)[1] == dim(C0)[1]))
  p = dim(W)[1]; q = dim(V)[1]
  nt = dim(y)[2]

  if(is.null(dim(F))) F = array(F,c(1,q,nt))
  if(is.matrix(F)) F = array(rep(F,nt),dim=c(dim(F)[1],dim(F)[2],nt)) else stopifnot(length(dim(F) == 3) & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  stopifnot((dim(F)[1] == dim(V)[1]) & (dim(F)[2] == dim(G)[1]) & (dim(F)[3] == nt))
  
  m = matrix(NA, nr = p, nc = nt + 1)
  C = array(NA, dim = c(p, p, nt + 1))
  a <- b <- rep(NA, nt + 1)
  f = matrix(NA, nr = q, nc = nt)
  Q = array(NA, dim = c(q, q, nt))
  A = matrix(NA, nr = p, nc = nt)
  R = array(NA, dim = c(p, p, nt))
  m[,1] = m0; C[,,1] = C0; a[1] = a0; b[1] = b0
  for(i in 1:nt)
  {
    A[,i] = G%*%m[,i]; R[,,i] = G%*%C[,,i]%*%t(G) + W
    f[,i] = F[,,i]%*%A[,i]; Q[,,i] = F[,,i]%*%R[,,i]%*%t(F[,,i]) + V  
    e = y[,i] - f[,i]; Qinv = solve(Q[,,i])
    RFQinv = R[,,i]%*%t(F[,,i])%*%Qinv
    m[,i+1] = A[,i] + RFQinv%*%e
    C[,,i+1] = R[,,i] - RFQinv%*%F[,,i]%*%R[,,i]
    a[i+1] = a[i] + q/2
    b[i+1] = b[i] + (1/2)%*%e%*%Qinv%*%e
  }
  return(list(m = m, C = C, a = a, b = b, f = f, Q = Q, A = A, R = R))
}

# Function to calculate the log marginal likelihood of the data for a DLM where the states and observations have common variance factor sigma^2
cv.lmarglik <- function(y, f, Q, a, b)
{
  y = as.matrix(y); f = as.matrix(f); a = rep(a,1); b = rep(b,1)
  stopifnot((length(dim(Q)) == 3) & (dim(y)[1] == dim(f)[1]) & (dim(y)[2] == dim(f)[2]) & (dim(y)[1] == dim(Q)[1]) & (dim(y)[2] == dim(Q)[3]) & (dim(Q)[1] == dim(Q)[2]) & (dim(y)[2] == length(a) - 1) & (length(a) == length(b)))
  nt = dim(y)[2]
  log.marglik = 0
  for(i in 1:nt) log.marglik = log.marglik + dmt(y[,i], f[,i], Q[,,i]*b[i]/a[i], 2*a[i])
  return(log.marglik)
}

#### Utility functions ####

# Calculate the density of the multivariate student-t distribution
dmt <- function(x, mu, Sigma, df, log = TRUE)
{
  Sigma = as.matrix(Sigma)
  p = dim(Sigma)[1]
  stopifnot((length(df) == 1) & (length(x) == p) & (length(mu) == p) & (dim(Sigma)[2] == p))

  log.lik = log(gamma((df+p)/2)) - log(gamma((df)/2)) - (p/2)*(log(df) + log(pi)) - .5*log(det(Sigma)) - ((df+p)/2)*log(1 + (1/df)*(t(x-mu)%*%solve(Sigma)%*%(x-mu)))
  if(log) return(log.lik) else return(exp(log.lik))
}