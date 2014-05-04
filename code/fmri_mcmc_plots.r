# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to construct traceplots for MCMC chains
fmri_mcmc_plot <- function(N, n, n.sim, mod, diff, dimx, n.chains, nsims, nburn, nthin, same)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,diff,"-",dimx,".rdata",sep=""))
  mysims = get(paste(mod,"_dat",sep=""))[[1]][[n]]
  if(strsplit(mod,"")[[1]][4] == 0) modtype = "reg" else modtype = "dlm"
  if(!same & modtype == "reg") same.lab = "" else same.lab = paste("-",same,sep="")
  
  # Load MCMC data
  out.all <- list()
  for(i in 1:n.chains)
  {
    file = paste(dpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,i,nsims,nburn,nthin,sep="-"),same.lab,".rdata",sep="")
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
data1 = expand.grid(N=20,n=6,n.sim=1:10,mod=c("M011","M101"),diff="",dimx=3,n.chains=3,nsims=11000,nburn=1000,nthin=1,same=FALSE,stringsAsFactors=FALSE)
data2 = expand.grid(N=20,n=6,n.sim=1:10,mod="M101",diff="-diff",dimx=3,n.chains=3,nsims=11000,nburn=1000,nthin=1,same=FALSE,stringsAsFactors=FALSE)
data3 = expand.grid(N=20,n=4,n.sim=1:10,mod="M101",diff="-diff",dimx=3,n.chains=3,nsims=11000,nburn=1000,nthin=1,same=TRUE,stringsAsFactors=FALSE)
data4 = expand.grid(N=20,n=6,n.sim=1:10,mod=c("M010","M020"),diff="",dimx=3,n.chains=3,nsims=11000,nburn=1000,nthin=1,same=FALSE,stringsAsFactors=FALSE)
mydata = rbind(data1,data2,data3,data4)
m_ply(mydata, fmri_mcmc_plot)

# # Function to plot histograms of posterior means (or medians) among multiple simulations
# fmri_mcmc_hist <- function(N, n, mod, n.chains, nsims, nburn, nthin, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
# {
#   # Load simulated data
#   load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
#   mysims = get(paste(mod,"_dat",sep=""))[[1]][[n]]
#   if(strsplit(mod,"")[[1]][4] == 0)
#   {
#     modtype = "reg"
#     vars = paste(beta,phi,sigma2s,sep="-")
#   } else {
#     modtype = "dlm"
#     vars = paste(x,beta,sigma2m,phi,sigma2s,sep="-")
#   }
#   
#   # Load MCMC data and calculate medians
#   for(j in 1:N)
#   {
#     out.all = list()
#     for(i in 1:n.chains)
#     {
#       load(paste(dpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,j,mod,i,nsims,nburn,nthin,vars,sep="-"),".rdata",sep=""))
#       out.all[[i]] = out
#     }
#     n.sims = out.all[[1]]$mcmc.details$n.sims
#     n.burn = out.all[[1]]$mcmc.details$n.burn
#     n.thin = out.all[[1]]$mcmc.details$n.thin
#     n.keep = (n.sims - n.burn) %/% n.thin
#     
#     # Calculate median for beta
#     if(beta)
#     {
#       d = length(mysims[[j]]$true.params$beta)
#       if(j == 1) med.beta = matrix(NA, nr=nsims, nc = d)
#       for(k in 1:d) med.beta[j,k] = median(sapply(out.all, function(x) x$beta[,k]))
#     }
#     
#     # Calculate median for sigma2m
#     if(sigma2m)
#     {
#       if(j == 1) med.sigma2m = rep(NA, nsims)
#       med.sigma2m[j] = median(sapply(out.all, function(x) x$sigma2m))
#     }
#     
#     # Calculate median for phi
#     if(phi)
#     {
#       p = dim(mysims[[j]]$true.params$G)[1]
#       if(j == 1) med.phi = matrix(NA, nr=nsims, nc = p)
#       for(k in 1:p) med.phi[j,k] = median(sapply(out.all, function(x) x$phi[,k]))
#     }
#     
#     # Calculate median for sigma2s
#     if(sigma2s)
#     {
#       if(j == 1) med.sigma2s = rep(NA, nsims)
#       med.sigma2s[j] = median(sapply(out.all, function(x) x$sigma2s))
#     }
#     
#     # Calculate medians for x (at 9 ponits)
#     if(x)
#     {
#       tt = dim(out.all[[1]]$x)[3]
#       npts = floor(seq(1,tt,len=9))
#       if(j == 1) med.x = matrix(NA, nr=nsims, nc = length(npts))
#       for(k in 1:length(npts)) med.x[j,k] = median(sapply(out.all, function(a) a$x[,1,npts[k]]))
#     }
#   }
#   
#   # Construct histograms
#   if(beta)
#   {
#     pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-beta.pdf",sep=""),width=5*d,height=5)
#     par(mfrow=c(1,d))
#     for(i in 1:d)
#     {
#       hist(med.beta[,i],xlab=eval(bquote(expression(beta[.(i-1)]))),main="")
#       abline(v=mysims[[j]]$true.params$beta[i],lwd=2,col=2)
#     }
#     dev.off()
#   }
#   if(sigma2m)
#   {
#     pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-sigma2m.pdf",sep=""))
#     hist(med.sigma2m,xlab=expression(sigma[m]^2),main="")
#     abline(v=mysims[[j]]$true.params$V[1,1],lwd=2,col=2)
#     dev.off()
#   }
#   if(phi)
#   {
#     pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-phi.pdf",sep=""),width=5*p,height=5)
#     par(mfrow=c(1,p))
#     for(i in 1:p)
#     {
#       hist(med.phi[,i],xlab=eval(bquote(expression(phi[.(i)]))),main="")
#       abline(v=mysims[[j]]$true.params$G[i,1],lwd=2,col=2)
#     }
#     dev.off()
#   }
#   if(sigma2s)
#   {
#     pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-sigma2s.pdf",sep=""))
#     hist(med.sigma2s,xlab=expression(sigma[s]^2),main="")
#     abline(v=mysims[[j]]$true.params$W[1,1],lwd=2,col=2)
#     dev.off()
#   }
#   if(x)
#   {
#     tt = dim(out.all[[1]]$x)[3]
#     npts = floor(seq(1,tt,len=9))
#     pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-x.pdf",sep=""))
#     par(mfrow=c(3,3))
#     for(k in 1:length(npts))
#     {
#       hist(med.x[,k],xlab=eval(bquote(expression(x[.(npts[k]-1)]))),main="")
#       abline(v=mysims[[j]]$x[1,npts[k]],lwd=2,col=2)
#     }
#     dev.off()
#   }
# }
# 
# require(plyr)
# data1 = expand.grid(N=20,n=6,mod=c("M011","M101"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=1,sigma2m=1,stringsAsFactors=FALSE)
# data2 = expand.grid(N=20,n=6,mod=c("M010","M020"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=0,sigma2m=0,stringsAsFactors=FALSE)
# mydata = rbind(data1,data2)
# m_ply(mydata, fmri_mcmc_hist)