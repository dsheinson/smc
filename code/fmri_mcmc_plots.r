# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

fmri_mcmc_plot <- function(n.sim, mod, n.chains, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-20-",mod,".rdata",sep=""))
  
  # Load MCMC data
  out.all <- list()
  for(i in 1:n.chains)
  {
    load(paste(dpath,"fmri_mcmc_test-",paste(n.sim,mod,i,beta,sigma2m,phi,sigma2s,sep="-"),".rdata",sep=""))
    out.all[[i]] = out
  }
  n.sims = out.all[[1]]$mcmc.details$n.sims
  n.burn = out.all[[1]]$mcmc.details$n.burn
  n.thin = out.all[[1]]$mcmc.details$n.thin
  n.keep = (n.sims - n.burn) %/% n.thin

  # Traceplots for beta
  if(beta)
  {
    iter = (1:n.keep)*n.thin
    d = length(mysims[[n.sim]]$true.params$beta)
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-beta.pdf",sep=""))
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
    }
    dev.off()
  }
  
  # Traceplots for sigma2m
  if(sigma2m)
  {
    iter = (1:n.keep)*n.thin
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-sigma2m.pdf",sep=""))
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$sigma2m)))
    ymax = max(sapply(out.all, function(x) max(x$sigma2m)))
    plot(iter,out.all[[1]]$sigma2m,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[m]^2),cex.lab=1.5)
    if(n.chains > 1)
    {
      for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2m,col=2*(j-1))
    }
    abline(h=mysims[[n.sim]]$true.params$V[1,1])
    dev.off()
  }
  
  # Traceplots for phi
  if(phi)
  {
    iter = (1:n.keep)*n.thin
    p = length(mysims[[n.sim]]$true.params$G[,1])
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-phi.pdf",sep=""))
    par(mfrow=c(p,1), mar=c(5,6,4,2)+.1)
    mins = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, min)),nr=p), 1, min)
    maxs = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, max)),nr=p), 1, max)
    
    for(i in 1:p)
    {
      ylab = bquote(expression(phi[.(i)]))
      plot(iter,out.all[[1]]$phi[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=eval(ylab),cex.lab=1.5)
      ar = sum(sapply(out.all, function(x) x$accept.phi)) / (n.chains*n.sims)
      mtext(paste("Acceptance rate: ",round(ar,3),sep="") ,side = 3)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$phi[,i],col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$G[i,1])
    }
    dev.off()
  }
  
  # Traceplots for sigma2s
  if(sigma2s)
  {
    iter = (1:n.keep)*n.thin
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-sigma2s.pdf",sep=""))
    par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
    ymin = min(sapply(out.all, function(x) min(x$sigma2s)))
    ymax = max(sapply(out.all, function(x) max(x$sigma2s)))
    plot(iter,out.all[[1]]$sigma2s,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[s]^2),cex.lab=1.5)
    if(n.chains > 1)
    {
      for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2s,col=2*(j-1))
    }
    abline(h=mysims[[n.sim]]$true.params$W[1,1])
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(n.sim = 1:20, mod = "M101", n.chains = 3, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1, stringsAsFactors=FALSE)
m_ply(mydata, fmri_mcmc_plot)

# Function to calculate medians
fmri_mcmc_medians <- function(nsims, mod, n.chains, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-20-",mod,".rdata",sep=""))
  
  # Load MCMC data and calculate medians
  for(j in 1:nsims)
  {
    out.all = list()
    for(i in 1:n.chains)
    {
      load(paste(dpath,"fmri_mcmc_test-",paste(j,mod,i,beta,sigma2m,phi,sigma2s,sep="-"),".rdata",sep=""))
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
      for(k in 1:d) med.beta[j,d] = median(sapply(out.all, function(x) x$beta[,d]))
    }
    # Traceplots for sigma2m
    if(sigma2m)
    {
      iter = (1:n.keep)*n.thin
      pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-sigma2m.pdf",sep=""))
      par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
      ymin = min(sapply(out.all, function(x) min(x$sigma2m)))
      ymax = max(sapply(out.all, function(x) max(x$sigma2m)))
      plot(iter,out.all[[1]]$sigma2m,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[m]^2),cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2m,col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$V[1,1])
      dev.off()
    }
    
    # Traceplots for phi
    if(phi)
    {
      iter = (1:n.keep)*n.thin
      p = length(mysims[[n.sim]]$true.params$G[,1])
      pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-phi.pdf",sep=""))
      par(mfrow=c(p,1), mar=c(5,6,4,2)+.1)
      mins = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, min)),nr=p), 1, min)
      maxs = apply(matrix(sapply(out.all, function(x) apply(x$phi, 2, max)),nr=p), 1, max)
      
      for(i in 1:p)
      {
        ylab = bquote(expression(phi[.(i)]))
        plot(iter,out.all[[1]]$phi[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=eval(ylab),cex.lab=1.5)
        ar = sum(sapply(out.all, function(x) x$accept.phi)) / (n.chains*n.sims)
        mtext(paste("Acceptance rate: ",round(ar,3),sep="") ,side = 3)
        if(n.chains > 1)
        {
          for(j in 2:n.chains) lines(iter,out.all[[j]]$phi[,i],col=2*(j-1))
        }
        abline(h=mysims[[n.sim]]$true.params$G[i,1])
      }
      dev.off()
    }
    
    # Traceplots for sigma2s
    if(sigma2s)
    {
      iter = (1:n.keep)*n.thin
      pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-sigma2s.pdf",sep=""))
      par(mfrow=c(1,1), mar=c(5,6,4,2)+.1)
      ymin = min(sapply(out.all, function(x) min(x$sigma2s)))
      ymax = max(sapply(out.all, function(x) max(x$sigma2s)))
      plot(iter,out.all[[1]]$sigma2s,type="l",ylim=c(ymin,ymax),xlab="Iteration",ylab=expression(sigma[s]^2),cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$sigma2s,col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$W[1,1])
      dev.off()
    }
  }
}
