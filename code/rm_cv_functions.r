dllik.cv <- function(y, x, theta, f = 1, v = 1) dnorm(y,f*x,sqrt(theta*v),log=TRUE)

revo.cv <- function(x, theta, g = 1, w = 1) rnorm(1,g*x,sqrt(theta*w))

rprior.cv <- function(a0=1,b0=1,m0=0,C0=1)
{
  theta = 1 / rgamma(1,a0,b0)
  x = rnorm(1,m0,sqrt(C0*theta))
  return(list(x=x,theta=theta))
}

# assume y and x are univariate 1 by nt and 1 by nt + 1 matrices, respectively
rm_mcmc <- function(y, x, theta, a0, b0, mydlm, L, n.iter, move.states = TRUE)
{
  # Set k = nt - L where L is the lag from the current time point
  nt = length(y)
  if(L >= nt)
  {
    k = 1
    L = nt - 1
  } else {
   k = nt - L
  }

  for(t in 1:n.iter)
  {
    # Sample smoothed states by FFBS
    if(move.states | n.iter > 1)
    {
      samp.x = sample.states(y, theta, mydlm, L)
      if(k > 1) samp.x = c(x[1:(k-1)], samp.x)
    } else {
      samp.x = x
    }
    
    # Sample theta from full conditional
    theta = sample.theta(y, samp.x, a0, b0, mydlm)
  }
  return(list(state=samp.x,theta=theta))
}

###################
# Utility functions
###################

sample.theta <- function(y, x, a0, b0, mydlm)
{
  K = length(y)
  F = mydlm$F; G = mydlm$G; V = mydlm$V; W = mydlm$W; m0 = mydlm$m0; C0 = mydlm$C0 
  a = a0 + K + 0.5
  b = (1/(2*V))*sum((y - x[2:(K+1)])^2) + (1/(2*W))*sum((x[2:(K+1)] - x[1:K])^2) + (1/(2*C0))*((x[1]-m0)^2) + b0
  return(1/rgamma(1,a,b))
}

sample.states <- function(y, theta, mydlm, L)
{
  # Initialize DLM
  nt = length(y)
  F = mydlm$F
  G = mydlm$G
  V = theta*mydlm$V
  W = theta*mydlm$W
  
  # Forward-filtering
  m = rep(NA, nt + 1)
  C = rep(NA, nt + 1)
  f = rep(NA, nt)
  Q = rep(NA, nt)
  A = rep(NA, nt)
  R = rep(NA, nt)
  m[1] = mydlm$m0; C[1] = theta*mydlm$C0
  for(i in 1:nt)
  {
    A[i] = G*m[i]; R[i] = G*C[i]*G + W
    f[i] = F*A[i]; Q[i] = F*R[i]*F + V  
    e = y[i] - f[i]; Qinv = 1 / Q[i]
    RFQinv = R[i]*F*Qinv
    m[i+1] = A[i] + RFQinv*e
    C[i+1] = R[i] - RFQinv*F*R[i]
  }
  
  # Backward-sampling
  x = rep(NA, L + 2)
  x[L + 2] = rnorm(1, m[nt + 1], sqrt(C[nt + 1]))
  for(i in (L+1):1)
  {
    h = m[nt - L - 1 + i] + C[nt - L - 1 + i]*G*(1/R[nt - L - 1 + i])*(x[i+1] - A[nt - L - 1 + i])
    H = C[nt - L - 1 + i] - C[nt - L - 1 + i]*G*(1/R[nt - L - 1 + i])*G*C[nt - L - 1 + i]
    x[i] = rnorm(1, h, sqrt(H))
  }
  return(x)
}