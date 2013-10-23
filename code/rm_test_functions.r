##########################################
# Functions to use within particle filters
##########################################

dllik <- function(y, x, theta) dnorm(y,x,sqrt(theta),log=TRUE)

revo <- function(x, theta) rnorm(1,x,sqrt(theta))

rprior <- function(j,a,b)
{
  mytheta = 1 / rgamma(1,a,b)
  mystate = rnorm(1,0,sqrt(mytheta))
  return(list(x=mystate,theta=mytheta))
}

rm_mcmc <- function(y, x, theta, a0, b0, n.iter)
{
  for(t in 1:n.iter)
  {
    # Sample theta from full conditional
    theta = sample.theta(y, x, theta, a0, b0)

    # Sample states by FFBS
    x = sample.states(y, theta)
  }
  return(list(state=x,theta=theta))
}

###################
# Utility functions
###################

sample.theta <- function(y, x, theta, a0, b0)
{
  K = length(y)
  a = K + 1.5 + a0 - 1
  b = .5*(sum((y - x[2:(K+1)])^2) + sum((x[2:(K+1)] - x[1:K])^2) + x[1]^2 + 2*b0)
  return(1/rgamma(1,a,b))
}

sample.states <- function(y, theta)
{
  # Initialize DLM
  nt = length(y)
  W = V = theta
  F = G = 1

  # Forward-filtering
  m = rep(NA, nt + 1)
  C = rep(NA, nt + 1)
  f = rep(NA, nt)
  Q = rep(NA, nt)
  A = rep(NA, nt)
  R = rep(NA, nt)
  m[1] = 0; C[1] = theta
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

##############################
# Functions to construct plots
##############################
pf.states_plot <- function(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab, leg = TRUE)
{
  n = length(a) - 1

  # Compute 95% credible intervals of filtered states over time
  # True posterior
  lk = qt(0.025,2*a)*sqrt(C*(b/a)) + m
  uk = qt(0.975,2*a)*sqrt(C*(b/a)) + m

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

  # Resample-move SMC - moving theta only
  require(Hmisc)
  lr1 = rep(NA,n+1)
  ur1 = rep(NA,n+1)
  for(i in 1:length(out.rm1$state))
  {
    mystates = out.rm1$state[[i]][1,,i]
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
    lr2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.025,normwt=TRUE)
    ur2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.975,normwt=TRUE)
  }

  # Plot 95% credible intervals of filtered states over time
  gmin = min(x,lk,lkd,lr1,lr2); gmax = max(x,uk,ukd,ur1,ur2)
  plot(0:n,x,ylim=c(gmin,gmax),type="l",xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=2)
  lines(0:n,lk,col=2)
  lines(0:n,uk,col=2)
  lines(0:n,lkd,col=3)
  lines(0:n,ukd,col=3)
  lines(0:n,lr1,col=6)
  lines(0:n,ur1,col=6)
  lines(0:n,lr2,col=4)
  lines(0:n,ur2,col=4)
  if(leg) legend("topright",legend=c("x","Post.","KD","RM1","RM2"),lty=rep(1,5),col=c(1,2,3,6,4),cex=1.75)
}

pf.states.zeroed_plot <- function(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab, leg = TRUE)
{
  n = length(a) - 1

  # Compute 95% credible intervals of filtered states over time
  # True posterior
  lk = qt(0.025,2*a)*sqrt(C*(b/a)) + m
  uk = qt(0.975,2*a)*sqrt(C*(b/a)) + m

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

  # Resample-move SMC - moving theta only
  require(Hmisc)
  lr1 = rep(NA,n+1)
  ur1 = rep(NA,n+1)
  for(i in 1:length(out.rm1$state))
  {
    mystates = out.rm1$state[[i]][1,,i]
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
    lr2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.025,normwt=TRUE)
    ur2[i] = wtd.quantile(mystates,out.rm2$weight[,i],probs=0.975,normwt=TRUE)
  }

  # Plot 95% credible intervals of filtered states over time
  gmin = min(x-m,lk-m,lkd-m,lr1-m,lr2-m); gmax = max(x-m,uk-m,ukd-m,ur1-m,ur2-m)
  plot(0:n,x-m,ylim=c(gmin,gmax),type="l",xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=2)
  lines(0:n,lk-m,col=2)
  lines(0:n,uk-m,col=2)
  lines(0:n,lkd-m,col=3)
  lines(0:n,ukd-m,col=3)
  lines(0:n,lr1-m,col=6)
  lines(0:n,ur1-m,col=6)
  lines(0:n,lr2-m,col=4)
  lines(0:n,ur2-m,col=4)
  if(leg) legend("bottomright",legend=c("x","Post.","KD","RM1","RM2"),lty=rep(1,5),col=c(1,2,3,6,4),cex=1.75)
}

pf.pvalues_plot <- function(m, C, a, b, out.kd, out.rm1, out.rm2, xlab, ylab, leg = TRUE)
{
  n = length(a) - 1

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
  plot(0:n,pval[,1],ylim=c(gmin,gmax),type="l",col=3,xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=2)
  lines(0:n,pval[,2],col=6)
  lines(0:n,pval[,3],col=4)
  if(leg) legend("topleft",legend=c("KD","RM1","RM2"),lty=rep(1,3),col=c(3,6,4),cex=1.75)
}

pf.precision_plot <- function(m, C, a, b, v, out.kd, out.rm1, out.rm2, xlab, ylab, leg = TRUE)
{
  n = length(a) - 1

  # Plot 95% credible intervals of unknown precision over time
  # Analytic estimates
  lk = qgamma(0.025,a,b)
  uk = qgamma(0.975,a,b)

  # Kernel density SMC
  lkd = rep(NA,n+1)
  ukd = rep(NA,n+1)
  for(i in 1:(n+1))
  {
    lkd[i] = wtd.quantile(1/exp(out.kd$theta[1,,i]),out.kd$weight[,i],probs=0.025,normwt=TRUE)
    ukd[i] = wtd.quantile(1/exp(out.kd$theta[1,,i]),out.kd$weight[,i],probs=0.975,normwt=TRUE)
  }

  # resample-move - moving theta only
  lr1 = rep(NA,n+1)
  ur1 = rep(NA,n+1)
  for(i in 1:(n+1))
  {
    lr1[i] = wtd.quantile(1/out.rm1$theta[1,,i],out.rm1$weight[,i],probs=0.025,normwt=TRUE)
    ur1[i] = wtd.quantile(1/out.rm1$theta[1,,i],out.rm1$weight[,i],probs=0.975,normwt=TRUE)
  }

  # resample-move - moving theta and states
  lr2 = rep(NA,n+1)
  ur2 = rep(NA,n+1)
  for(i in 1:(n+1))
  {
    lr2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.025,normwt=TRUE)
    ur2[i] = wtd.quantile(1/out.rm2$theta[1,,i],out.rm2$weight[,i],probs=0.975,normwt=TRUE)
  }

  burn = 1
  gmin = min(cbind(lk,v)[-(1:burn),])
  gmax = max(cbind(uk,v)[-(1:burn),])
  plot(0:n,lk,type="l",col=2,ylim=c(gmin,gmax),xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=2)
  lines(0:n,uk,col=2)
  abline(h=v)
  lines(0:n,lkd,col=3)
  lines(0:n,ukd,col=3)
  lines(0:n,lr1,col=6)
  lines(0:n,ur1,col=6)
  lines(0:n,lr2,col=4)
  lines(0:n,ur2,col=4)
  if(leg) legend("topright",legend=c("True precision","True posterior","KD","RM1","RM2"),lty=rep(1,5),col=c(1,2,3,6,4),cex=1.75)
}