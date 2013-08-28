dllik <- function(y, x, theta) dnorm(y,x,sqrt(theta),log=TRUE)

revo <- function(x, theta) rnorm(1,x,sqrt(theta))

rprior <- function(j)
{
  mytheta = 1 / rgamma(1,1,.25)
  mystate = rnorm(1,0,sqrt(mytheta))
  return(list(x=mystate,theta=mytheta))
}

require(dlm)
rm_mcmc <- function(y, x, theta, n.iter)
{
  for(t in 1:n.iter)
  {
    # Sample theta from full conditional
    theta = sample.theta(y, x, theta)

    # Sample states by FFBS
    mydlm = dlm(list(m0=0, C0=theta, FF=1, V=theta, GG=1, W=theta))
    x = dlmBSample(dlmFilter(y, mydlm))
  }
  return(list(state=x,theta=theta))
}

sample.theta <- function(y, x, theta)
{
  K = length(y)
  a = K + 1.5
  b = .5*(sum((y - x[2:(K+1)])^2) + sum((x[2:(K+1)] - x[1:K])^2) + x[1]^2) + 4
  return(1/rgamma(1,a,b))
}

sample.theta.mh <- function(y, x, theta)
{
  K = length(y)
  theta.proposal = 0
  while(theta.proposal == 0) theta.proposal = rgamma(1,theta^2,theta)
  ll.curr <- ll.proposal <- 0
  lp.curr <- dgamma(theta,1,.25,log=T)+dnorm(x[1],0,sqrt(theta),log=TRUE)
  lp.proposal <- dgamma(theta.proposal,1,.25,log=T)+dnorm(x[1],0,sqrt(theta.proposal),log=TRUE)
  for(k in 1:K)
  {
    ll.curr = ll.curr + dllik(y[k],x[k+1],theta)
    ll.proposal = ll.proposal + dllik(y[k],x[k+1],theta.proposal)
    lp.curr = lp.curr + dnorm(x[k+1],x[k],sqrt(theta),log=TRUE)
    lp.proposal = lp.proposal + dnorm(x[k+1],x[k],sqrt(theta.proposal),log=TRUE)
  }
  numer = ll.proposal + lp.proposal + dnorm(theta,theta.proposal,1,log=T)
  denom = ll.curr + lp.curr + dnorm(theta.proposal,theta,1,log=T)
  logMH = numer - denom
  if(log(runif(1)) < logMH) theta = theta.proposal
  return(theta)
}