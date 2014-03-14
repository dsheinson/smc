dllik <- function(y, x, theta, f = 1, v = 1) dnorm(y,f*x,sqrt(theta*v),log=TRUE)

revo <- function(x, theta, g = 1, w = 1) rnorm(1,g*x,sqrt(theta*w))

rprior <- function(a=1,b=1,m0=0,C0=1)
{
  mytheta = 1 / rgamma(1,a,b)
  mystate = rnorm(1,m0,sqrt(C0*mytheta))
  return(list(x=mystate,theta=mytheta))
}

rm_mcmc <- function(y, x, theta, a0, b0, mydlm, n.iter)
{
  for(t in 1:n.iter)
  {
    # Sample theta from full conditional
    theta = sample.theta(y, x, theta, a0, b0, mydlm)
    
    # Sample states by FFBS
    x = sample.states(y, theta, mydlm)
  }
  return(list(state=x,theta=theta))
}

###################
# Utility functions
###################

sample.theta <- function(y, x, theta, a0, b0, mydlm)
{
  K = length(y)
  F = mydlm$F; G = mydlm$G; V = mydlm$V; W = mydlm$W; m0 = mydlm$m0; C0 = mydlm$C0 
  a = a0 + K + 0.5
  b = (1/(2*V))*sum((y - x[2:(K+1)])^2) + (1/(2*W))*sum((x[2:(K+1)] - x[1:K])^2) + (1/(2*C0))*((x[1]-m0)^2) + b0
  return(1/rgamma(1,a,b))
}

sample.states <- function(y, theta, mydlm)
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
  x = rep(NA, nt + 1)
  x[nt + 1] = rnorm(1, m[nt + 1], sqrt(C[nt + 1]))
  for(i in nt:1)
  {
    h = m[i] + C[i]*G*(1/R[i])*(x[i+1] - A[i])
    H = C[i] - C[i]*G*(1/R[i])*G*C[i]
    x[i] = rnorm(1, h, sqrt(H))
  }
  return(x)
}