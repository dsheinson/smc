# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to construct traceplots for MCMC chains
fmri_mcmc_plot <- function(N, n, n.sim, mod, diff, dimx, n.chains, nsims, nburn, nthin, same, rm.prior)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,diff,"-",dimx,".rdata",sep=""))
  mysims = get(paste(mod,"_dat",sep=""))[[1]][[n]]
  if(strsplit(mod,"")[[1]][4] == 0) modtype = "reg" else modtype = "dlm"
  if(!same & modtype == "reg") same.lab = "" else same.lab = paste("-",same,sep="")
  if(rm.prior) prior.lab = "-prior" else prior.lab = ""
  
  # Load MCMC data
  out.all <- list()
  for(i in 1:n.chains)
  {
    file = paste(dpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,i,nsims,nburn,nthin,sep="-"),same.lab,prior.lab,".rdata",sep="")
    load(file)
    out.all[[i]] = out.est
  }
  n.sims = out.all[[1]]$out$mcmc.details$n.sims
  n.burn = out.all[[1]]$out$mcmc.details$n.burn
  n.thin = out.all[[1]]$out$mcmc.details$n.thin
  n.keep = (n.sims - n.burn) %/% n.thin

  # Traceplots for beta
  iter = (1:n.keep)*n.thin
  d = length(mysims[[n.sim]]$true.params$beta)
  file = paste(gpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chains,nsims,nburn,nthin,sep="-"),same.lab,"-traceplots-beta.pdf",sep="")
  pdf(file=file)
  par(mfrow=c(d,1), mar=c(5,6,4,2)+.1)
  mins = apply(matrix(sapply(out.all, function(x) apply(x$out$beta, 2, min)),nr=d), 1, min)
  maxs = apply(matrix(sapply(out.all, function(x) apply(x$out$beta, 2, max)),nr=d), 1, max)
  for(i in 1:d)
  {
    ylab = bquote(expression(beta[.(i-1)]))
    xlabs = c(rep("",d-1),"Iteration")
    true.beta = mysims[[n.sim]]$true.params$beta[i]
    mle.beta = out.all[[1]]$theta.mle[i]
    plot(iter,out.all[[1]]$out$beta[,i],type="l",ylim=c(min(mins[i],mle.beta,true.beta),max(maxs[i],mle.beta,true.beta)),xlab=xlabs[i],ylab=eval(ylab),cex.lab=1.5)
    ess.mcmcse = ess(rep(sapply(out.all, function(x) x$out$beta[,i]),1))
    ess.coda = 0
    for(j in 1:n.chains) ess.coda = ess.coda + effectiveSize(mcmc(out.all[[j]]$out$beta[,i]))
    mtext(paste("ESS: ",round(ess.mcmcse,2)," (mcmcse), ",round(ess.coda,2)," (coda)",sep="") ,side = 3)
    if(n.chains > 1) for(j in 2:n.chains) lines(iter,out.all[[j]]$out$beta[,i],col=2*(j-1))
    abline(h=true.beta,lty=2)
    abline(h=mle.beta,col="gray",lwd=3)
    mtext(true.beta, side = 4, at = true.beta)
  }
  dev.off()
  
  # Traceplots for sigma2m
  if(modtype == "dlm")
  {
    file = paste(gpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chains,nsims,nburn,nthin,same,sep="-"),"-traceplots-sigma2m.pdf",sep="")
    pdf(file=file)
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$out$sigma2m)))
    ymax = max(sapply(out.all, function(x) max(x$out$sigma2m)))
    true.sigma2m = mysims[[n.sim]]$true.params$V[1,1]
    mle.sigma2m = out.all[[1]]$theta.mle[length(out.all[[1]]$theta.mle)]
    plot(iter,out.all[[1]]$out$sigma2m,type="l",ylim=c(min(true.sigma2m,mle.sigma2m,ymin),max(true.sigma2m,mle.sigma2m,ymax)),xlab="Iteration",ylab=expression(sigma[m]^2),cex.lab=1.5)
    ess.mcmcse = ess(rep(sapply(out.all, function(x) x$out$sigma2m),1))
    ess.coda = 0
    for(j in 1:n.chains) ess.coda = ess.coda + effectiveSize(mcmc(out.all[[j]]$out$sigma2m))
    mtext(paste("ESS: ",round(ess.mcmcse,2)," (mcmcse), ",round(ess.coda,2)," (coda)",sep="") ,side = 3)
    if(n.chains > 1) for(j in 2:n.chains) lines(iter,out.all[[j]]$out$sigma2m,col=2*(j-1))
    abline(h=true.sigma2m,lty=2)
    abline(h=mle.sigma2m,col="gray",lwd=3)
    mtext(true.sigma2m, side = 4, at = true.sigma2m)
    dev.off()
  }
  
  # Traceplots for phi   
  P = length(out.all[[1]]$out$phi)
  p.sum = 1
  for(i in 1:P)
  {
    if(i > 1 & same) break
    p = dim(out.all[[1]]$out$phi[[i]])[2]
    file = paste(gpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chains,nsims,nburn,nthin,sep="-"),same.lab,"-traceplots-phi-",i,".pdf",sep="")
    pdf(file=file)
    par(mfrow=c(p,1), mar=c(5,6,4,2)+.1)
    mins = apply(matrix(sapply(out.all, function(x) apply(x$out$phi[[i]], 2, min)),nr=p), 1, min)
    maxs = apply(matrix(sapply(out.all, function(x) apply(x$out$phi[[i]], 2, max)),nr=p), 1, max)
    xlabs = c(rep("",p-1),"Iteration")
    for(j in 1:p)
    {
      true.phi = mysims[[n.sim]]$true.params$G[p.sum:(p.sum+j-1),p.sum]
      mle.phi = out.all[[1]]$theta.mle[d+p.sum+j-1]
      ylab = bquote(expression(phi[.(j)]))
      plot(iter,out.all[[1]]$out$phi[[i]][,j],type="l",ylim=c(min(true.phi,mle.phi,mins[j]),max(true.phi,mle.phi,maxs[j])),xlab=xlabs[j],ylab=eval(ylab),main=paste(i,"th AR component",sep=""),cex.lab=1.5)
      ar = mean(sapply(out.all, function(x) mean(x$out$accept.phi[[i]])))
      ess.mcmcse = ess(rep(sapply(out.all, function(x) x$out$phi[[i]][,j]),1))
      ess.coda = 0
      for(k in 1:n.chains) ess.coda = ess.coda + effectiveSize(mcmc(out.all[[k]]$out$phi[[i]][,j]))
      mtext(paste("ESS: ",round(ess.mcmcse,2)," (mcmcse), ",round(ess.coda,2)," (coda) Accept Rate: ",round(ar,3),sep="") ,side = 3)
      if(n.chains > 1)  for(k in 2:n.chains) lines(iter,out.all[[k]]$out$phi[[i]][,j],col=2*(k-1))
      abline(h=true.phi,lty=2)
      abline(h=mle.phi,col="gray",lwd=3)
      mtext(true.phi, side = 4, at = true.phi)
    }
    dev.off()
    p.sum = p.sum + p
  }
  
  # Traceplots for sigma2s   
  ps = 1
  for(i in 1:P)
  {
    if(i > 1 & same) break
    p = dim(out.all[[1]]$out$phi[[i]])[2]
    file = paste(gpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chains,nsims,nburn,nthin,sep="-"),same.lab,"-traceplots-sigma2s-",i,".pdf",sep="")
    pdf(file=file)
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$out$sigma2s[,i])))
    ymax = max(sapply(out.all, function(x) max(x$out$sigma2s[,i])))
    true.sigma2s = mysims[[n.sim]]$true.params$W[ps,ps]
    mle.sigma2s = out.all[[1]]$theta.mle[d+p.sum]
    plot(iter,out.all[[1]]$out$sigma2s[,i],type="l",ylim=c(min(true.sigma2s,mle.sigma2s,ymin),max(true.sigma2s,mle.sigma2s,ymax)),xlab="Iteration",ylab=expression(sigma[s]^2),main = paste(i,"th AR component",sep=""), cex.lab=1.5)
    ess.mcmcse = ess(rep(sapply(out.all, function(x) x$out$sigma2s[,i]),1))
    ess.coda = 0
    for(j in 1:n.chains) ess.coda = ess.coda + effectiveSize(mcmc(out.all[[j]]$out$sigma2s[,i]))
    mtext(paste("ESS: ",round(ess.mcmcse,2)," (mcmcse), ",round(ess.coda,2)," (coda)",sep="") ,side = 3)
    if(n.chains > 1) for(j in 2:n.chains) lines(iter,out.all[[j]]$out$sigma2s[,i],col=2*(j-1))
    abline(h=true.sigma2s,lty=2)
    abline(h=mle.sigma2s,col="gray",lwd=3)
    mtext(true.sigma2s, side = 4, at = true.sigma2s)
    dev.off()
    p.sum = p.sum + 1
    ps = ps + p
  }
  
  # Traceplots for x_t (9 points)
  if(modtype == "dlm")
  {
    tt = dim(out.all[[1]]$out$x)[3]
    npts = floor(seq(1,tt,len=9))
    p.sum = 1
    for(i in 1:P)
    {
      p = dim(out.all[[1]]$out$phi[[i]])[2]
      file = paste(gpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chains,nsims,nburn,nthin,same,sep="-"),"-traceplots-x-",i,".pdf",sep="")
      pdf(file=file)
      par(mfrow=c(3,3), mar=c(5,6,2,2)+.1)
      ymin = min(sapply(out.all, function(a) min(a$out$x[,p.sum,npts])))
      ymax = max(sapply(out.all, function(a) max(a$out$x[,p.sum,npts])))
      for(k in 1:length(npts))
      {
        ylab = bquote(expression(x[.(npts[k]-1)]))
        if(k == 1) xlab = "Iteration" else xlab = ""
        true.x = mysims[[n.sim]]$x[p.sum,npts[k]]
        mle.x = out.all[[1]]$x.mle[p.sum,npts[k]]
        plot(iter,out.all[[1]]$out$x[,p.sum,npts[k]],type="l",ylim=c(min(true.x,mle.x,ymin),max(true.x,mle.x,ymax)),xlab=xlab,ylab=eval(ylab),cex.lab=1.5)
        ess.mcmcse = ess(rep(sapply(out.all, function(x) x$out$x[,p.sum,npts[k]]),1))
        ess.coda = 0
        for(j in 1:n.chains) ess.coda = ess.coda + effectiveSize(mcmc(out.all[[j]]$out$x[,p.sum,npts[k]]))
        mtext(paste("ESS: ",round(ess.mcmcse,2)," (mcmcse), ",round(ess.coda,2)," (coda)",sep="") ,side = 3, cex = 0.5)
        if(n.chains > 1) for(j in 2:n.chains) lines(iter,out.all[[j]]$out$x[,p.sum,npts[k]],col=2*(j-1))
        abline(h=true.x,lty=2)
        abline(h=mle.x,col="gray",lwd=3)
        mtext(round(true.x,3), cex = 0.9, side = 4, at = true.x)
      }
      dev.off()
      p.sum = p.sum + p
    }
  }
}

require(coda)
require(mcmcse)
require(plyr)
data1 = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M010","M020"),diff="",dimx=2,n.chains=3,nsims=11000,nburn=1000,nthin=10,same=FALSE,rm.prior=FALSE,stringsAsFactors=FALSE)
data2 = expand.grid(N=20,n=6,n.sim=1:8,mod=c("M011","M101"),diff="",dimx=2,n.chains=3,nsims=11000,nburn=1000,nthin=10,same=FALSE,rm.prior=c(FALSE,TRUE),stringsAsFactors=FALSE)
mydata = rbind(data1,data2)
m_ply(mydata, fmri_mcmc_plot)