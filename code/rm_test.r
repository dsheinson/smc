# Simulate data from a dlm with common state/observation variance
v = 1; F = 2; G = .95; W = 1; V = 1; x0 = rnorm(1,0,sqrt(v))
n = 200
x = rep(NA,n+1)
y = rep(NA,n)
x[1] = x0
for(i in 1:n)
{
  x[i+1] = G*x[i] + W*rnorm(1,0,v)
  y[i] = F*x[i] + V*rnorm(1,0,v)
}

# Plot data
gmin = min(x,y); gmax = max(x,y)
plot(0:n,x,ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
lines(1:n,y,lty=2)
legend("topright",legend=expression(x,y),lty=c(1,2))

## Compute marginal filtered distributions of states and unknown variance
m = C = a = b = rep(NA,n+1)
m[1] = 0; C[1] = 1; a[1] = 1; b[1] = .25
for(i in 1:n)
{
  A = G*m[i]; R = G*C[i]*G + W
  f = F*A; Q = F*R*F + V  
  e = y[i] - f; Qinv = 1/Q
  RFQinv = R*F*(1/Q)
  m[i+1] = A + RFQinv*e
  C[i+1] = R - RFQinv*F*R
  a[i+1] = a[i] + 1/2
  b[i+1] = b[i] + (1/2)*e*Qinv*e
}

## Perform MCMC using JAGS
require(rjags)
d1 = list(y=y, n=length(y))
mod1 = jags.model("dlmRW-model.txt", data=d1, n.chains=3, n.adapt=1e3)
samps1 = coda.samples(mod1, c("tau.e","x"), n.iter = 1e4)
#gelman.diag(samps1) # psrf should be around 1
#summary(samps1)

## Perform SMC
# kernel density particle filter
dllik <- function(y, x, theta) dnorm(y,2*x,sqrt(exp(theta)),log=TRUE)
revo <- function(x, theta) rnorm(1,.95*x,sqrt(exp(theta)))
rprior <- function(j)
{
  mytheta = rgamma(1,1,.25)
  mystate = rnorm(1,0,sqrt(mytheta))
  return(list(x=mystate,theta=log(mytheta)))
}
pstate <- function(x, theta) .95*x
source("kd_pf.r")
out.kd = kd_pf(y, dllik, pstate, revo, rprior, 100, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

# resample-move
# function to evaluate the logarithm of the likelihood of a new observation given the current state
dllik <- function(y, x, theta)
{
  if(!is.matrix(x)) x = matrix(x, 1)
  dnorm(y,2*x[1,dim(x)[2]],sqrt(theta),log=TRUE)
}

# function to sample from state evolution equation
revo <- function(x, theta)
{
  if(!is.matrix(x)) x = matrix(x, 1)
  rnorm(1,.95*x[1,dim(x)[2]],sqrt(theta))
}

# function to sample prior state and parameter
rprior <- function(j)
{
  mytheta = rgamma(1,1,.25)
  mystate = rnorm(1,0,sqrt(mytheta))
  return(list(x=mystate,theta=mytheta))
}

# function to move parameter by Metropolis-Hastings
rmove <- function(y, x, theta)
{
  n.iter = 1
  if(!is.matrix(x)) x = matrix(x, 1)
  if(!is.matrix(y)) y = matrix(y, 1)
  K = dim(y)[2]
  for(t in 1:n.iter)
  {
    theta.proposal = 0
    while(theta.proposal == 0) theta.proposal = rgamma(1,theta^2,theta)
#    x.proposal = rnorm(1,x[1,1],sqrt(theta.proposal))
    x.proposal = x[1,1]
    ll.curr <- ll.proposal <- 0
    lp.curr <- dgamma(theta,1,.25,log=T)+dnorm(x[1,1],0,sqrt(theta),log=TRUE)
    lp.proposal <- dgamma(theta.proposal,1,.25,log=T)+dnorm(x.proposal,0,sqrt(theta.proposal),log=TRUE)
    for(k in 1:K)
    {
#      x.proposal <- c(x.proposal, rnorm(1,x.proposal[k],sqrt(theta.proposal)))
      x.proposal <- c(x.proposal, x[1,k+1])
      ll.curr = ll.curr + dllik(y[1,k],x[1,k+1],theta)
      ll.proposal = ll.proposal + dllik(y[1,k],x.proposal[k+1],theta.proposal)
      lp.curr = lp.curr + dnorm(x[1,k+1],x[1,k],sqrt(theta),log=TRUE)
      lp.proposal = lp.proposal + dnorm(x.proposal[k+1],x.proposal[k],sqrt(theta.proposal),log=TRUE)
    }
    numer = ll.proposal + lp.proposal + dnorm(theta,theta.proposal,1,log=T)
    denom = ll.curr + lp.curr + dnorm(theta.proposal,theta,1,log=T)
    logMH = numer - denom
    if(log(runif(1)) < logMH)
    { 
      theta = theta.proposal
      x[1,] = x.proposal
    }
  }
  return(list(state=x,theta=theta))
}

# run resample-move particle filter
source("rm_pf.r")
out.rm = rm_pf(y, dllik, revo, rprior, rmove, 100, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

#########################################

# Plot 95% credible intervals of filtered states over time
gmin = min(x,y); gmax = max(x,y)
plot(0:n,x,ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
lines(1:n,y,lty=2)
legend("topright",legend=expression(x,y),lty=c(1,2))
# Kalman filter
lk = qt(0.025,2*a)*sqrt(C*(b/a)) + m
uk = qt(0.975,2*a)*sqrt(C*(b/a)) + m
lines(0:n,lk,col=2)
lines(0:n,uk,col=2)
# JAGS MCMC
lj = apply(samps1[[1]], 2, function(x) quantile(x,probs=0.025))
uj = apply(samps1[[1]], 2, function(x) quantile(x,probs=0.975))
lines(0:n,lj[-1],col=4)
lines(0:n,uj[-1],col=4)
# kernel density particle filter
require(Hmisc)
mystates = out.kd$state[1,,]
lkd = rep(NA,n+1)
ukd = rep(NA,n+1)
for(i in 1:dim(mystates)[2])
{
  lkd[i] = wtd.quantile(mystates[,i],out.kd$weight[,i],probs=0.025,normwt=TRUE)
  ukd[i] = wtd.quantile(mystates[,i],out.kd$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lkd,col=3)
lines(0:n,ukd,col=3)
# Resample-move SMC
require(Hmisc)
lr = rep(NA,n+1)
ur = rep(NA,n+1)
for(i in 1:length(out.rm$state))
{
  lr[i] = wtd.quantile(out.rm$state[[i]][1,,i],out.rm$weight[,i],probs=0.025,normwt=TRUE)
  ur[i] = wtd.quantile(out.rm$state[[i]][1,,i],out.rm$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr,col=6)
lines(0:n,ur,col=6)
legend("bottomleft",legend=c("KF","MCMC","KD","RM"),lty=c(1,1,1,1),col=c(2,4,3,6))

# Plot 95% credible intervals of unknown precision over time
# Analytic estimates
lk = qgamma(0.025,a,b)
uk = qgamma(0.975,a,b)
burn = 20
ymin = min(lk[-(1:burn)],v)
ymax = max(uk[-(1:burn)],v)
plot(0:n,lk,type="l",col=2,ylim=c(ymin,ymax),xlab=expression(t),ylab=expression(phi))
lines(0:n,uk,col=2)
abline(h=v)
# JAGS MCMC
abline(h=c(lj[1],uj[1]),col=4)
# Kernel density SMC
lkd = rep(NA,n+1)
ukd = rep(NA,n+1)
for(i in 1:(n+1))
{
  lkd[i] = wtd.quantile(1/exp(out.kd$theta[1,,i]),out.kd$weight[,i],probs=0.025,normwt=TRUE)
  ukd[i] = wtd.quantile(1/exp(out.kd$theta[1,,i]),out.kd$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lkd,col=3)
lines(0:n,ukd,col=3)
# resample-move
lr = rep(NA,n+1)
ur = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr[i] = wtd.quantile(1/out.rm$theta[1,,i],out.rm$weight[,i],probs=0.025,normwt=TRUE)
  ur[i] = wtd.quantile(1/out.rm$theta[1,,i],out.rm$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr,col=6)
lines(0:n,ur,col=6)
legend("topright",legend=c("KF","MCMC","KD","RM"),lty=c(1,1,1,1),col=c(2,4,3,6))