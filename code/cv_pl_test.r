# Test particle filters on local level model
source("dlm_cv_functions.r")
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

cv_pl_test <- function(np, n.sim = 1, filt = c("pl","rm"), alpha = 0.05, burn = 0)
{
  # Load simulated data
  load(paste(dpath,"rw_sim.rdata",sep=""))
  y = mysims[[n.sim]]$y
  
  # Calculate true 95% CI of states and precision
  post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=1, a0=1, b0=1, m0=0, C0=1)
  lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  lp = qgamma(alpha/2,post$a,post$b)
  up = qgamma(1-alpha/2,post$a,post$b)
  
  # Run particle filters and calculate quantiles
  out = state.quant = theta.quant = list()
  if('pl' %in% filt)
  {
    source("pl_functions.r")
    source("pl.r")
    out$pl = pl(y, dlpred.ll, revo.ll, function(j) rprior.ll(), rmove.ll, smap.theta.ll, smap.state.ll, np, method="stratified",nonuniformity="none",log=F)    
    state.quant$pl = pf.quantile(out$pl$state, out$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    theta.quant$pl = pf.quantile(1/out$pl$theta, out$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  }
  if('rm' %in% filt)
  {
    source("rm_cv_functions.r")
    source("rm_pf.r")
    mydlm = list(F=1,G=1,V=1,W=1,m0=0,C0=1)
    rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0=1, b0=1, mydlm, 1)
    out$rm = rm_pf(mysims[[n.sim]]$y, dllik, revo, rprior, rmove, np, store.filter = TRUE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
    state.quant$rm = pf.quantile(out$rm$state, out$rm$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    theta.quant$rm = pf.quantile(1/out$rm$theta, out$rm$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  }

  # Plot 95% CI for states
  nt = dim(y)[2]
  if(burn > 0)
  {
    g.lk = lk[-(1:burn)]
    g.uk = uk[-(1:burn)]
    g.ls = sapply(state.quant, function(x) x[-(1:burn),1,1])
    g.us = sapply(state.quant, function(x) x[-(1:burn),1,2])
  } else {
    g.lk = lk
    g.uk = uk
    g.ls = sapply(state.quant, function(x) x[,1,1])
    g.us = sapply(state.quant, function(x) x[,1,2])
  }
  gmin = min(mysims[[n.sim]]$x[1,], g.lk, apply(g.ls,2,min))
  gmax = max(mysims[[n.sim]]$x[1,], g.uk, apply(g.us,2,max))
  pdf(paste(gpath,"pl-test-states-",n.sim,"-",np,".pdf",sep=""))
  plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position",col='gray',lwd=2)
  mtext(paste(np," particles",sep=""),side=3)
  lines(0:nt,lk)
  lines(0:nt,uk)
  for(i in 1:length(filt))
  {
    lines(0:nt,state.quant[[i]][,1,1],col=2*i)
    lines(0:nt,state.quant[[i]][,1,2],col=2*i)
  }
  legend("topleft",legend=c(filt, "True Post", "True Sim"),lty=c(rep(1,length(filt)),1,1),col=c(2*(1:length(filt)),1,'gray'),lwd=c(rep(1,length(filt)),1,2))
  title("95% CI for filtered states")
  dev.off()
  
  # Plot 95% CI for precision
  if(burn > 0)
  {
    g.lp = lp[-(1:burn)]
    g.up = up[-(1:burn)]
    g.lf = sapply(theta.quant, function(x) x[-(1:burn),1,1])
    g.uf = sapply(theta.quant, function(x) x[-(1:burn),1,2])
  } else {
    g.lp = lp
    g.up = up
    g.lf = sapply(theta.quant, function(x) x[,1,1])
    g.uf = sapply(theta.quant, function(x) x[,1,2])
  }
  gmin = min(1, g.lp, apply(g.lf,2,min))
  gmax = max(1, g.up, apply(g.uf,2,max))
  pdf(paste(gpath,"pl-test-precision-",n.sim,"-",np,".pdf",sep=""))
  plot(0:nt,lp,type="l",ylim=c(gmin,gmax),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(1/theta))
  mtext(paste(np," particles",sep=""),side=3)
  lines(0:nt,up)
  for(i in 1:length(filt))
  {
    lines(0:nt,theta.quant[[i]][,1,1],col=2*i)
    lines(0:nt,theta.quant[[i]][,1,2],col=2*i)
  }
  abline(h=1,col='gray',lwd=2)
  legend("topright",legend=c(filt, "True Post", "True Sim"),lty=c(rep(1,length(filt)),1,1),col=c(2*(1:length(filt)),1,'gray'),lwd=c(rep(1,length(filt)),1,2))
  dev.off()

  return(c(sapply(out, function(x) pf.lmarglik(x)),cv.lmarglik(y, post$f, post$Q, post$a, post$b)))
}

require(plyr)
mydata = expand.grid(np = 100, n.sim = 1, stringsAsFactors = FALSE)
maply(mydata, cv_pl_test)
