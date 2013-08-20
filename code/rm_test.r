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

## Perform SMC using kernel density particle filter
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

## Perform SMC using resample-move particle filter - move parameter only using Metropolis-Hastings
source("rm_test_functions.r")
rmove <- function(y, x, theta) return(list(state=x,theta=sample.theta(y, x, theta)))
source("rm_pf.r")
out.rm1 = rm_pf(y, dllik, revo, rprior, rmove, 100, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

## Perform SMC using resample-move particle filter - move states and parameter using Gibbs sampling
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
out.rm2 = rm_pf(y, dllik, revo, rprior, rmove, 100, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

#######
# Plots
#######

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
# Resample-move SMC - moving theta only
require(Hmisc)
lr1 = rep(NA,n+1)
ur1 = rep(NA,n+1)
for(i in 1:length(out.rm1$state))
{
  lr1[i] = wtd.quantile(out.rm1$state[[i]][1,,i],out.rm1$weight[,i],probs=0.025,normwt=TRUE)
  ur1[i] = wtd.quantile(out.rm1$state[[i]][1,,i],out.rm1$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr1,col=6)
lines(0:n,ur1,col=6)
# Resample-move SMC - moving states and theta
require(Hmisc)
lr2 = rep(NA,n+1)
ur2 = rep(NA,n+1)
for(i in 1:length(out.rm2$state))
{
  lr2[i] = wtd.quantile(out.rm2$state[[i]][1,,i],out.rm2$weight[,i],probs=0.025,normwt=TRUE)
  ur2[i] = wtd.quantile(out.rm2$state[[i]][1,,i],out.rm2$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr2,col=7)
lines(0:n,ur2,col=7)

legend("bottomleft",legend=c("KF","MCMC","KD","RM1","RM2"),lty=c(1,1,1,1,1),col=c(2,4,3,6,7))

# Plot 95% credible intervals of unknown precision over time
# Analytic estimates
x11()
lk = qgamma(0.025,a,b)
uk = qgamma(0.975,a,b)
ymin = 0.05
ymax = 1.95
#burn = 20
#ymin = min(lk[-(1:burn)],v)
#ymax = max(uk[-(1:burn)],v)
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
# resample-move - moving theta only
lr1 = rep(NA,n+1)
ur1 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr1[i] = wtd.quantile(1/out.rm1$theta[1,,i],out.rm1$weight[,i],probs=0.025,normwt=TRUE)
  ur1[i] = wtd.quantile(1/out.rm1$theta[1,,i],out.rm1$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr1,col=6)
lines(0:n,ur1,col=6)
# resample-move - moving theta only
lr2 = rep(NA,n+1)
ur2 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.025,normwt=TRUE)
  ur2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr2,col=7)
lines(0:n,ur2,col=7)

legend("topright",legend=c("KF","MCMC","KD","RM1","RM2"),lty=c(1,1,1,1,1),col=c(2,4,3,6,7))