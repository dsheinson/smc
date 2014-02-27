# Function to simulate nt + 1 states (x_t) and nt observations (y_t) from a dlm
# y_t = U_t*beta + F_t*x_t + v_t, v_t iid N(0,V)
# x_t = Gx_t-1 + w_t, w_t iid N(0,W)
dlm.sim <- function(nt, F, G, V, W, m0, C0, beta, U = 0)
{
  # Initialize parameters and check dimensions
  m0 = rep(m0,1)
  G = as.matrix(G); V = as.matrix(V); W = as.matrix(W); C0 = as.matrix(C0)
  stopifnot((dim(C0)[1] == dim(C0)[2]) & (dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(G)[1] == dim(W)[1]) & (dim(G)[1] == dim(C0)[1]) & (dim(G)[1] == length(m0)))
  p = dim(W)[1]; q = dim(V)[1]
  
  if(is.null(dim(F))) F = array(F,c(q,p,nt))
  if(is.matrix(F)) F = array(rep(F,nt),dim=c(dim(F)[1],dim(F)[2],nt)) else stopifnot(length(dim(F) == 3) & dim(F)[1] == q & dim(F)[2] == p & dim(F)[3] == nt)
  
  if(missing(beta)) beta = rep(0,1)
  d = length(beta)
  if(is.null(dim(U))) U = array(U,c(q,d,nt))
  if(is.matrix(U))
  {
    stopifnot(dim(U)[1] == q & dim(U)[2] == d)
    U = array(U,dim=c(q,d,nt))
  } else {
    stopifnot(length(dim(U) == 3) & dim(U)[1] == q & dim(U)[2] == d & dim(U)[3] == nt)
  }
  
  # Simulate data
  try(chol.C0 = chol(C0),silent=TRUE)
  if(class(chol.C0) == "try-error") chol.C0 = chol(C0,pivot=TRUE)
  try(chol.V = chol(V),silent=TRUE)
  if(class(chol.V) == "try-error") chol.V = chol(V,pivot=TRUE)
  try(chol.W = chol(W),silent=TRUE)
  if(class(chol.W) == "try-error") chol.W = chol(W,pivot=TRUE)
  x0 = chol.C0%*%rnorm(p,m0,1)
  x = matrix(NA, nr = p, nc = nt + 1)
  y = matrix(NA, nr = q, nc = nt)
  x[,1] = x0
  for(i in 1:nt)
  {
    x[,i+1] = G%*%x[,i] + t(chol.W)%*%rnorm(p,0,1)
    y[,i] = U[,,i]%*%beta + F[,,i]%*%x[,i+1] + t(chol.V)%*%rnorm(q,0,1)
  }
  return(list(y = y, x = x, true.params = list(F=F,G=G,V=V,W=W,m0=m0,C0=C0,beta=beta,U=U)))
}