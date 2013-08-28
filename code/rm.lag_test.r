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

source("rm_pf.r")
source("rm_test_functions.r")

## Perform SMC using resample-move particle filter - move states and parameter using Gibbs sampling
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
out.rm25 = rm_pf(y, dllik, revo, rprior, rmove, 100, .25, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

## Perform SMC using resample-move particle filter - move states and parameter using Gibbs sampling
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
out.rm50 = rm_pf(y, dllik, revo, rprior, rmove, 100, .50, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

## Perform SMC using resample-move particle filter - move states and parameter using Gibbs sampling
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
out.rm75 = rm_pf(y, dllik, revo, rprior, rmove, 100, .75, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

## Perform SMC using resample-move particle filter - move states and parameter using Gibbs sampling
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
out.rm100 = rm_pf(y, dllik, revo, rprior, rmove, 100, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

# Save data
save.image("../data/rm.lag_test.rdata")

#######
# Plots
#######

# Compute 95% credible intervals of filtered states over time
# True posterior
lk = qt(0.025,2*a)*sqrt(C*(b/a))
uk = qt(0.975,2*a)*sqrt(C*(b/a))

# Resample-move SMC moving states and theta - lag 50
require(Hmisc)
lr25 = rep(NA,n+1)
ur25 = rep(NA,n+1)
for(i in 1:length(out.rm25$state))
{
  mystates = out.rm25$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr25[i] = wtd.quantile(mystates,out.rm25$weight[,i],probs=0.025,normwt=TRUE)
  ur25[i] = wtd.quantile(mystates,out.rm25$weight[,i],probs=0.975,normwt=TRUE)
}

# Resample-move SMC moving states and theta - lag 100
require(Hmisc)
lr50 = rep(NA,n+1)
ur50 = rep(NA,n+1)
for(i in 1:length(out.rm50$state))
{
  mystates = out.rm50$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr50[i] = wtd.quantile(mystates,out.rm50$weight[,i],probs=0.025,normwt=TRUE)
  ur50[i] = wtd.quantile(mystates,out.rm50$weight[,i],probs=0.975,normwt=TRUE)
}

# Resample-move SMC moving states and theta - lag 150
require(Hmisc)
lr75 = rep(NA,n+1)
ur75 = rep(NA,n+1)
for(i in 1:length(out.rm75$state))
{
  mystates = out.rm75$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr75[i] = wtd.quantile(mystates,out.rm75$weight[,i],probs=0.025,normwt=TRUE)
  ur75[i] = wtd.quantile(mystates,out.rm75$weight[,i],probs=0.975,normwt=TRUE)
}

# Resample-move SMC moving states and theta - no lag
require(Hmisc)
lr100 = rep(NA,n+1)
ur100 = rep(NA,n+1)
for(i in 1:length(out.rm100$state))
{
  mystates = out.rm100$state[[i]][1,,i]
  mystates = mystates - mean(mystates)
  lr100[i] = wtd.quantile(mystates,out.rm100$weight[,i],probs=0.025,normwt=TRUE)
  ur100[i] = wtd.quantile(mystates,out.rm100$weight[,i],probs=0.975,normwt=TRUE)
}

# Plot 95% credible intervals of filtered states over time
gmin = min(lk,lr25,lr50,lr75,lr100); gmax = max(uk,ur25,ur50,ur75,ur100)
plot(0:n,rep(0,n+1),ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
lines(0:n,lk,col=2)
lines(0:n,uk,col=2)
lines(0:n,lr25,col=3)
lines(0:n,ur25,col=3)
lines(0:n,lr50,col=4)
lines(0:n,ur50,col=4)
lines(0:n,lr75,col=6)
lines(0:n,ur75,col=6)
lines(0:n,lr100,col=7)
lines(0:n,ur100,col=7)
legend("bottomright",legend=c("x","Post.","Lag 50","Lag 100","Lag 150","Lag 200"),lty=rep(1,6),col=c(1,2,3,4,6,7))

# Compute p-values of test of particle filtered distribution versus true posterior
pval = matrix(NA,n+1,4)
for(i in 1:(n+1))
{
  make.t <- function(x) (x - m[i])/(sqrt(C[i]*(b[i]/a[i])))
  pval[i,1] = ks.test(make.t(out.rm25$state[[i]][1,,i]),pt,df=2*a[i])$p.value
  pval[i,2] = ks.test(make.t(out.rm50$state[[i]][1,,i]),pt,df=2*a[i])$p.value
  pval[i,3] = ks.test(make.t(out.rm75$state[[i]][1,,i]),pt,df=2*a[i])$p.value
  pval[i,4] = ks.test(make.t(out.rm100$state[[i]][1,,i]),pt,df=2*a[i])$p.value
}

# Plot p-values over time
gmin = min(pval[-1,]); gmax = max(pval[-1,])
x11()
plot(0:n,pval[,1],ylim=c(gmin,gmax),type="l",col=3,xlab=expression(t),ylab="P-value")
lines(0:n,pval[,2],col=4)
lines(0:n,pval[,3],col=6)
lines(0:n,pval[,4],col=7)

# Plot 95% credible intervals of unknown precision over time
# Analytic estimates
x11()
lk = qgamma(0.025,a,b)
uk = qgamma(0.975,a,b)
burn = 10
ymin = min(cbind(lk,v)[-(1:burn),])
ymax = max(cbind(uk,v)[-(1:burn),])
plot(0:n,lk,type="l",col=2,ylim=c(ymin,ymax),xlab=expression(t),ylab=expression(phi))
lines(0:n,uk,col=2)
abline(h=v)
# resample-move lag 25
lr25 = rep(NA,n+1)
ur25 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr25[i] = wtd.quantile(1/out.rm25$theta[1,,i],out.rm25$weight[,i],probs=0.025,normwt=TRUE)
  ur25[i] = wtd.quantile(1/out.rm25$theta[1,,i],out.rm25$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr25,col=3)
lines(0:n,ur25,col=3)
# resample-move lag 50
lr50 = rep(NA,n+1)
ur50 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr50[i] = wtd.quantile(1/out.rm50$theta[1,,i],out.rm50$weight[,i],probs=0.025,normwt=TRUE)
  ur50[i] = wtd.quantile(1/out.rm50$theta[1,,i],out.rm50$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr50,col=4)
lines(0:n,ur50,col=4)
# resample-move lag 75
lr75 = rep(NA,n+1)
ur75 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr75[i] = wtd.quantile(1/out.rm75$theta[1,,i],out.rm75$weight[,i],probs=0.025,normwt=TRUE)
  ur75[i] = wtd.quantile(1/out.rm75$theta[1,,i],out.rm75$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr75,col=6)
lines(0:n,ur75,col=6)
# resample-move lag 100
lr100 = rep(NA,n+1)
ur100 = rep(NA,n+1)
for(i in 1:(n+1))
{
  lr100[i] = wtd.quantile(1/out.rm100$theta[1,,i],out.rm100$weight[,i],probs=0.025,normwt=TRUE)
  ur100[i] = wtd.quantile(1/out.rm100$theta[1,,i],out.rm100$weight[,i],probs=0.975,normwt=TRUE)
}
lines(0:n,lr100,col=7)
lines(0:n,ur100,col=7)

legend("topright",legend=c("True precision","True posterior","Lag 50","Lag 100","Lag 150","Lag 200"),lty=rep(1,6),col=c(1,2,3,4,6,7))
