# Simulate data from a dlm with common state/observation variance
v = 1; F = 1; G = 1; V = 1; W = 1
x0 = rnorm(1,0,sqrt(v))
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

## Compute true marginal filtered distributions of states and unknown variance
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

## Perform SMC using kernel density particle filter
dllik <- function(y, x, theta) dnorm(y,x,sqrt(exp(theta)),log=TRUE)
revo <- function(x, theta) rnorm(1,x,sqrt(exp(theta)))
rprior <- function(j)
{
  mytheta = 1 / rgamma(1,1,.25)
  mystate = rnorm(1,0,sqrt(mytheta))
  return(list(x=mystate,theta=log(mytheta)))
}
pstate <- function(x, theta) x
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

# Save data
save.image("../data/rm_test.rdata")
load("../data/rm_test.rdata")

#######
# Plots
#######

# Compute 95% credible intervals of filtered states over time
# True posterior
lk = qt(0.025,2*a)*sqrt(C*(b/a))
uk = qt(0.975,2*a)*sqrt(C*(b/a))

# kernel density particle filter
require(Hmisc)
mystates = apply(out.kd$state[1,,], 2, function(x) x - mean(x))
lkd = rep(NA,n+1)
ukd = rep(NA,n+1)
for(i in 1:dim(mystates)[2])
{
  lkd[i] = wtd.quantile(mystates[,i],out.kd$weight[,i],probs=0.025,normwt=TRUE)
  ukd[i] = wtd.quantile(mystates[,i],out.kd$weight[,i],probs=0.975,normwt=TRUE)
}

# Resample-move SMC - moving theta only
require(Hmisc)
lr1 = rep(NA,n+1)
ur1 = rep(NA,n+1)
for(i in 1:length(out.rm1$state))
{
  mystates = out.rm1$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr1[i] = wtd.quantile(mystates,out.rm1$weight[,i],probs=0.025,normwt=TRUE)
  ur1[i] = wtd.quantile(mystates,out.rm1$weight[,i],probs=0.975,normwt=TRUE)
}

# Resample-move SMC - moving states and theta
require(Hmisc)
lr2 = rep(NA,n+1)
ur2 = rep(NA,n+1)
for(i in 1:length(out.rm2$state))
{
  mystates = out.rm2$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.025,normwt=TRUE)
  ur2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.975,normwt=TRUE)
}

# Plot 95% credible intervals of filtered states over time
gmin = min(lk,lkd,lr1,lr2,lrj); gmax = max(uk,ukd,ur1,ur2,urj)
plot(0:n,rep(0,n+1),ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
lines(0:n,lk,col=2)
lines(0:n,uk,col=2)
lines(0:n,lkd,col=3)
lines(0:n,ukd,col=3)
lines(0:n,lr1,col=6)
lines(0:n,ur1,col=6)
lines(0:n,lr2,col=4)
lines(0:n,ur2,col=4)
legend("bottomright",legend=c("x","Post.","KD","RM1","RM2"),lty=rep(1,5),col=c(1,2,3,6,4))

# Compute p-values of test of particle filtered distribution versus true posterior
pval = matrix(NA,n+1,3)
for(i in 1:(n+1))
{
  make.t <- function(x) (x - m[i])/(sqrt(C[i]*(b[i]/a[i])))
  pval[i,1] = ks.test(make.t(out.kd$state[1,,i]),pt,df=2*a[i])$p.value
  pval[i,2] = ks.test(make.t(out.rm1$state[[i]][1,,i]),pt,df=2*a[i])$p.value
  pval[i,3] = ks.test(make.t(out.rm2$state[[i]][1,,i]),pt,df=2*a[i])$p.value
}

# Plot p-values over time
gmin = min(pval[-1,]); gmax = max(pval[-1,])
x11()
plot(0:n,pval[,1],ylim=c(gmin,gmax),type="l",col=3,xlab=expression(t),ylab="P-value")
lines(0:n,pval[,2],col=6)
lines(0:n,pval[,3],col=4)
legend("topleft",legend=c("KD","RM1","RM2"),lty=rep(1,3),col=c(3,6,4))

# Plot 95% credible intervals of unknown precision over time
# Analytic estimates
x11()
lk = qgamma(0.025,a,b)
uk = qgamma(0.975,a,b)
burn = 15
ymin = min(cbind(lk,v)[-(1:burn),])
ymax = max(cbind(uk,v)[-(1:burn),])
plot(0:n,lk,type="l",col=2,ylim=c(ymin,ymax),xlab=expression(t),ylab=expression(phi))
lines(0:n,uk,col=2)
abline(h=v)
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
# resample-move - moving theta and states
lr2 = rep(NA,n+1)
ur2 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.025,normwt=TRUE)
  ur2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr2,col=4)
lines(0:n,ur2,col=4)
legend("topright",legend=c("True precision","True posterior","KD","RM1","RM2"),lty=rep(1,5),col=c(1,2,3,6,4))