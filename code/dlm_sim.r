# Function to simulate nt + 1 states (x_t) and nt observations (y_t) from a dlm
# y_t = U_t*beta + F_t*x_t + v_t, v_t iid N(0,V)
# x_t = Gx_t-1 + w_t, w_t iid N(0,W)
dlm.sim <- function(nt, F, G, V, W, x0, beta, U = 0)
{
  # Initialize parameters and check dimensions
  x0 = rep(x0,1)
  G = as.matrix(G); V = as.matrix(V); W = as.matrix(W)
  stopifnot((dim(G)[1] == dim(G)[2]) & (dim(V)[1] == dim(V)[2]) & (dim(W)[1] == dim(W)[2]) & (dim(G)[1] == dim(W)[1]) & (dim(G)[1] == length(x0)))
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
  
  # Compute choleski decompositions of CO, V, W and check if semi-positive definite
  oo.V <- 1:q; oo.W <- 1:p
  chol.V = try(chol(V),silent=TRUE)
  if(class(chol.V) == "try-error")
  {
    chol.V = chol(V,pivot=TRUE)
    oo.V = order(attr(chol.V, "pivot"))
  }
  chol.W = try(chol(W),silent=TRUE)
  if(class(chol.W) == "try-error")
  {
    chol.W = chol(W,pivot=TRUE)
    oo.W = order(attr(chol.W, "pivot"))
  }

  # Simulate data
  x = matrix(NA, nr = p, nc = nt + 1)
  y = matrix(NA, nr = q, nc = nt)
  x[,1] = x0
  for(i in 1:nt)
  {
    x[,i+1] = G%*%x[,i] + t(chol.W[,oo.W])%*%rnorm(p,0,1)
    y[,i] = matrix(U[,,i],nr=q,nc=d)%*%beta + matrix(F[,,i],nr=q,nc=p)%*%x[,i+1] + t(chol.V[,oo.V])%*%rnorm(q,0,1)
  }
  return(list(y = y, x = x, true.params = list(F=F,G=G,V=V,W=W,beta=beta,U=U)))
}