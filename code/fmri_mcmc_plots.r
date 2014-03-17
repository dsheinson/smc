# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to construct traceplots for MCMC chains
fmri_mcmc_plot <- function(N, n, n.sim, mod, n.chains, nsims, nburn, nthin, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
  mysims = get(paste(mod,"_dat",sep=""))[[1]][[n]]
  if(strsplit(mod,"")[[1]][4] == 0)
  {
    modtype = "reg"
    vars = paste(beta,phi,sigma2s,sep="-")
  } else {
    modtype = "dlm"
    vars = paste(x,beta,sigma2m,phi,sigma2s,sep="-")
  }
  
  # Load MCMC data
  out.all <- list()
  for(i in 1:n.chains)
  {
    load(paste(dpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,n.sim,mod,i,nsims,nburn,nthin,vars,sep="-"),".rdata",sep=""))
    out.all[[i]] = out
  }
  n.sims = out.all[[1]]$mcmc.details$n.sims
  n.burn = out.all[[1]]$mcmc.details$n.burn
  n.thin = out.all[[1]]$mcmc.details$n.thin
  n.keep = (n.sims - n.burn) %/% n.thin

  # Traceplots for beta
  if(beta)
  {
    # MLE
    if(modtype == "dlm") beta.mle = out.all[[1]]$initial$theta$beta else beta.mle = out.all[[1]]$initial$beta

    iter = (1:n.keep)*n.thin
    d = length(mysims[[n.sim]]$true.params$beta)
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,n.sim,mod,n.chains,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-traceplots-beta.pdf",sep=""))
    par(mfrow=c(d,1), mar=c(5,6,4,2)+.1)
    mins = apply(matrix(sapply(out.all, function(x) apply(x$beta, 2, min)),nr=d), 1, min)
    maxs = apply(matrix(sapply(out.all, function(x) apply(x$beta, 2, max)),nr=d), 1, max)

    for(i in 1:d)
    {
      ylab = bquote(expression(beta[.(i-1)]))
      plot(iter,out.all[[1]]$beta[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=eval(ylab),cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$beta[,i],col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$beta[i])
      abline(h=beta.mle[i],col="gray",lwd=3)
      mtext(mysims[[n.sim]]$true.params$beta[i], side = 4, at = mysims[[n.sim]]$true.params$beta[i])
    }
    dev.off()
  }
  
  # Traceplots for sigma2m
  if(sigma2m)
  {
    # MLE
    if(modtype == "dlm") sigma2m.mle = out.all[[1]]$initial$theta$sigma2m else stop("reg model should not have sigma2m")
    
    iter = (1:n.keep)*n.thin
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,n.sim,mod,n.chains,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-traceplots-sigma2m.pdf",sep=""))
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$sigma2m)))
    ymax = max(sapply(out.all, function(x) max(x$sigma2m)))
    plot(iter,out.all[[1]]$sigma2m,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[m]^2),cex.lab=1.5)
    if(n.chains > 1)
    {
      for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2m,col=2*(j-1))
    }
    abline(h=mysims[[n.sim]]$true.params$V[1,1])
    abline(h=sigma2m.mle,col="gray",lwd=3)
    mtext(mysims[[n.sim]]$true.params$V[1,1], side = 4, at = mysims[[n.sim]]$true.params$V[1,1])
    dev.off()
  }
  
  # Traceplots for phi
  if(phi)
  {
    # MLE
    if(modtype == "dlm") phi.mle = out.all[[1]]$initial$theta$phi else phi.mle = out.all[[1]]$initial$phi
    
    iter = (1:n.keep)*n.thin
    p = length(mysims[[n.sim]]$true.params$G[,1])
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,n.sim,mod,n.chains,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-traceplots-phi.pdf",sep=""))
    par(mfrow=c(p,1), mar=c(5,6,4,2)+.1)
    mins = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, min)),nr=p), 1, min)
    maxs = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, max)),nr=p), 1, max)
    
    for(i in 1:p)
    {
      ylab = bquote(expression(phi[.(i)]))
      plot(iter,out.all[[1]]$phi[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=eval(ylab),cex.lab=1.5)
      ar = mean(sapply(out.all, function(x) mean(x$accept.phi)))
      mtext(paste("Acceptance rate: ",round(ar,3),sep="") ,side = 3)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$phi[,i],col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$G[i,1])
      abline(h=phi.mle[i],col="gray",lwd=3)
      mtext(mysims[[n.sim]]$true.params$G[i,1], side = 4, at = mysims[[n.sim]]$true.params$G[i,1])
    }
    dev.off()
  }
  
  # Traceplots for sigma2s
  if(sigma2s)
  {
    # MLE
    if(modtype == "dlm") sigma2s.mle = out.all[[1]]$initial$theta$sigma2s else sigma2s.mle = out.all[[1]]$initial$sigma2s
    
    iter = (1:n.keep)*n.thin
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,n.sim,mod,n.chains,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-traceplots-sigma2s.pdf",sep=""))
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$sigma2s)))
    ymax = max(sapply(out.all, function(x) max(x$sigma2s)))
    plot(iter,out.all[[1]]$sigma2s,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[s]^2),cex.lab=1.5)
    if(n.chains > 1)
    {
      for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2s,col=2*(j-1))
    }
    abline(h=mysims[[n.sim]]$true.params$W[1,1])
    abline(h=sigma2s.mle,col="gray",lwd=3)
    mtext(mysims[[n.sim]]$true.params$W[1,1], side = 4, at = mysims[[n.sim]]$true.params$W[1,1])
    dev.off()
  }
  
  # Traceplots for x_t (9 points)
  if(x)
  {
    tt = dim(out.all[[1]]$x)[3]
    npts = floor(seq(1,tt,len=9))
    iter = (1:n.keep)*n.thin
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,n.sim,mod,n.chains,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-traceplots-x.pdf",sep=""))
    par(mfrow=c(3,3), mar=c(5,6,2,2)+.1)
    ymin = min(sapply(out.all, function(a) min(a$x[,,npts])))
    ymax = max(sapply(out.all, function(a) max(a$x[,,npts])))
    for(k in 1:length(npts))
    {
      ylab = bquote(expression(x[.(npts[k]-1)]))
      if(k == 1) xlab = "Iteration" else xlab = ""
      plot(iter,out.all[[1]]$x[,1,npts[k]],type="l",ylim=c(ymin,ymax),xlab=xlab,ylab=eval(ylab),cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$x[,1,npts[k]],col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$x[1,npts[k]])
      abline(h=out.all[[1]]$initial$x[1,npts[k]],col="gray",lwd=3)
      mtext(round(mysims[[n.sim]]$x[1,npts[k]],3), cex = 0.9, side = 4, at = mysims[[n.sim]]$x[1,npts[k]])
    }
    dev.off()
  }
}

require(plyr)
data1 = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M011","M101"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=1,sigma2m=1,stringsAsFactors=FALSE)
data2 = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M010","M020"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=0,sigma2m=0,stringsAsFactors=FALSE)
mydata = rbind(data1,data2)
m_ply(mydata, fmri_mcmc_plot)

# Function to plot histograms of posterior means (or medians) among multiple simulations
fmri_mcmc_hist <- function(N, n, mod, n.chains, nsims, nburn, nthin, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
  mysims = get(paste(mod,"_dat",sep=""))[[1]][[n]]
  if(strsplit(mod,"")[[1]][4] == 0)
  {
    modtype = "reg"
    vars = paste(beta,phi,sigma2s,sep="-")
  } else {
    modtype = "dlm"
    vars = paste(x,beta,sigma2m,phi,sigma2s,sep="-")
  }
  
  # Load MCMC data and calculate medians
  for(j in 1:N)
  {
    out.all = list()
    for(i in 1:n.chains)
    {
      load(paste(dpath,"fmri_",modtype,"_mcmc_test-",paste(N,n,j,mod,i,nsims,nburn,nthin,vars,sep="-"),".rdata",sep=""))
      out.all[[i]] = out
    }
    n.sims = out.all[[1]]$mcmc.details$n.sims
    n.burn = out.all[[1]]$mcmc.details$n.burn
    n.thin = out.all[[1]]$mcmc.details$n.thin
    n.keep = (n.sims - n.burn) %/% n.thin
    
    # Calculate median for beta
    if(beta)
    {
      d = length(mysims[[j]]$true.params$beta)
      if(j == 1) med.beta = matrix(NA, nr=nsims, nc = d)
      for(k in 1:d) med.beta[j,k] = median(sapply(out.all, function(x) x$beta[,k]))
    }
    
    # Calculate median for sigma2m
    if(sigma2m)
    {
      if(j == 1) med.sigma2m = rep(NA, nsims)
      med.sigma2m[j] = median(sapply(out.all, function(x) x$sigma2m))
    }
    
    # Calculate median for phi
    if(phi)
    {
      p = dim(mysims[[j]]$true.params$G)[1]
      if(j == 1) med.phi = matrix(NA, nr=nsims, nc = p)
      for(k in 1:p) med.phi[j,k] = median(sapply(out.all, function(x) x$phi[,k]))
    }
    
    # Calculate median for sigma2s
    if(sigma2s)
    {
      if(j == 1) med.sigma2s = rep(NA, nsims)
      med.sigma2s[j] = median(sapply(out.all, function(x) x$sigma2s))
    }
    
    # Calculate medians for x (at 9 ponits)
    if(x)
    {
      tt = dim(out.all[[1]]$x)[3]
      npts = floor(seq(1,tt,len=9))
      if(j == 1) med.x = matrix(NA, nr=nsims, nc = length(npts))
      for(k in 1:length(npts)) med.x[j,k] = median(sapply(out.all, function(a) a$x[,1,npts[k]]))
    }
  }
  
  # Construct histograms
  if(beta)
  {
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-beta.pdf",sep=""),width=5*d,height=5)
    par(mfrow=c(1,d))
    for(i in 1:d)
    {
      hist(med.beta[,i],xlab=eval(bquote(expression(beta[.(i-1)]))),main="")
      abline(v=mysims[[j]]$true.params$beta[i],lwd=2,col=2)
    }
    dev.off()
  }
  if(sigma2m)
  {
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-sigma2m.pdf",sep=""))
    hist(med.sigma2m,xlab=expression(sigma[m]^2),main="")
    abline(v=mysims[[j]]$true.params$V[1,1],lwd=2,col=2)
    dev.off()
  }
  if(phi)
  {
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-phi.pdf",sep=""),width=5*p,height=5)
    par(mfrow=c(1,p))
    for(i in 1:p)
    {
      hist(med.phi[,i],xlab=eval(bquote(expression(phi[.(i)]))),main="")
      abline(v=mysims[[j]]$true.params$G[i,1],lwd=2,col=2)
    }
    dev.off()
  }
  if(sigma2s)
  {
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-sigma2s.pdf",sep=""))
    hist(med.sigma2s,xlab=expression(sigma[s]^2),main="")
    abline(v=mysims[[j]]$true.params$W[1,1],lwd=2,col=2)
    dev.off()
  }
  if(x)
  {
    tt = dim(out.all[[1]]$x)[3]
    npts = floor(seq(1,tt,len=9))
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(N,n,mod,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),"-hist-x.pdf",sep=""))
    par(mfrow=c(3,3))
    for(k in 1:length(npts))
    {
      hist(med.x[,k],xlab=eval(bquote(expression(x[.(npts[k]-1)]))),main="")
      abline(v=mysims[[j]]$x[1,npts[k]],lwd=2,col=2)
    }
    dev.off()
  }
}

require(plyr)
data1 = expand.grid(N=20,n=6,mod=c("M011","M101"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=1,sigma2m=1,stringsAsFactors=FALSE)
data2 = expand.grid(N=20,n=6,mod=c("M010","M020"),n.chains=3,nsims=11000,nburn=1000,nthin=10,x=0,sigma2m=0,stringsAsFactors=FALSE)
mydata = rbind(data1,data2)
m_ply(mydata, fmri_mcmc_hist)