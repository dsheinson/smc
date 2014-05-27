# Construct plots for particle filter runs
source("dlm_cv_functions.r")
source("pf_mc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

cv_pl_quantiles <- function(lambda, np, nrtot, nruns, nsims, alpha = 0.05, burn = 0)
{
  # Plot quantiles for each sim for each # particles
  for(n.sim in 1:nsims)
  {
    # Load simulated data
    load(paste(dpath,"rw_sim.rdata",sep=""))
    y = mysims[[n.sim]]$y
    
    # Calculate true 95% CI of states and precision
    post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
    lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
    uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
    lp = qgamma(alpha/2,post$a,post$b)
    up = qgamma(1-alpha/2,post$a,post$b)
    
    for(j in 1:length(np))
    {
      # Load data
      load(paste(dpath,"cv_pl-",lambda,"-",np[j],"-",nrtot,"-",n.sim,".rdata",sep=""))
      attach(pf.out,warn.conflicts=F)
      filt = rownames(lmarglik)
      
      # Plot 95% CI for states
      nt = dim(y)[2]
      if(burn > 0)
      {
        g.lk = lk[-(1:burn)]
        g.uk = uk[-(1:burn)]
        g.ls = sapply(state.quant, function(x) sapply(x, function(w) w[-(1:burn),1,1]))
        g.us = sapply(state.quant, function(x) sapply(x, function(w) w[-(1:burn),1,2]))
      } else {
        g.lk = lk
        g.uk = uk
        g.ls = sapply(state.quant, function(x) sapply(x, function(w) w[,1,1]))
        g.us = sapply(state.quant, function(x) sapply(x, function(w) w[,1,2]))
      }
      gmin = min(mysims[[n.sim]]$x[1,], g.lk, apply(g.ls,2,min))
      gmax = max(mysims[[n.sim]]$x[1,], g.uk, apply(g.us,2,max))
      pdf(paste(gpath,"cv-pl-states-",10*lambda,"-",n.sim,"-",nrtot,"-",nruns,"-",np[j],".pdf",sep=""))
      plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position",col='gray',lwd=3)
      mtext(paste(np[j]," particles",sep=""),side=3)
      for(i in 1:length(filt))
      {
        for(k in 1:nruns)
        {
          lines(0:nt,state.quant[[k]][[i]][,1,1],col=i+1)
          lines(0:nt,state.quant[[k]][[i]][,1,2],col=i+1)
        }
      }
      lines(0:nt,mysims[[n.sim]]$x[1,],lwd=3,col="gray")
      lines(0:nt,lk,lwd=3)
      lines(0:nt,uk,lwd=3)
      legend("topleft",legend=c(filt, "True Post", "True Sim"),lty=c(rep(1,length(filt)),1,1),col=c(2:(length(filt)+1),1,'gray'),lwd=c(rep(1,length(filt)),1,2))
      title("95% CI for filtered states")
      dev.off()
      
      # Plot 95% CI for precision
      if(burn > 0)
      {
        g.lp = lp[-(1:burn)]
        g.up = up[-(1:burn)]
        g.lf = sapply(theta.quant, function(x)  sapply(x, function(w) w[-(1:burn),1,1]))
        g.uf = sapply(theta.quant, function(x)  sapply(x, function(w) w[-(1:burn),1,2]))
      } else {
        g.lp = lp
        g.up = up
        g.lf = sapply(theta.quant, function(x)  sapply(x, function(w) w[,1,1]))
        g.uf = sapply(theta.quant, function(x)  sapply(x, function(w) w[,1,2]))
      }
      gmin = min(1, g.lp, apply(g.lf,2,min))
      gmax = max(1, g.up, apply(g.uf,2,max))
      pdf(paste(gpath,"cv-pl-precision-",10*lambda,"-",n.sim,"-",nrtot,"-",nruns,"-",np[j],".pdf",sep=""))
      plot(0:nt,lp,type="l",lwd=3,ylim=c(gmin,gmax),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(1/theta))
      mtext(paste(np[j]," particles",sep=""),side=3)
      for(i in 1:length(filt))
      {
        for(k in 1:nruns)
        {
          lines(0:nt,theta.quant[[k]][[i]][,1,1],col=i+1)
          lines(0:nt,theta.quant[[k]][[i]][,1,2],col=i+1)
        }
      }
      lines(0:nt,lp,lwd=3)
      lines(0:nt,up,lwd=3)
      abline(h=1,col='gray',lwd=3)
      legend("topright",legend=c(filt, "True Post", "True Sim"),lty=c(rep(1,length(filt)),1,1),col=c(2:(length(filt)+1),1,'gray'),lwd=c(rep(1,length(filt)),1,2))
      dev.off()
    }
  }
}

require(plyr)
mydata = expand.grid(lambda = c(.5,1,2), np=c(100,500,1000,5000), nrtot=20, nruns=5,nsims=1,stringsAsFactors=FALSE)
m_ply(mydata, cv_pl_quantiles)

cv_pl_loglik_wsim = function(np, nruns, nsims, lambda = c(.5, 1, 2), filt = c('pl','plrb','kd','rm'))
{
  # Plot marginal likelihoods for each sim for each # particles
  for(n.sim in 1:nsims)
  {
    # Load simulated data
    load(paste(dpath,"rw_sim.rdata",sep=""))
    y = mysims[[n.sim]]$y
    
    # Calculate true log marginal likelihoods
    true_lmarglik = function(lambda)
    {
      post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
      return(cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b))
    }
    true.lmarglik = maply(data.frame(lambda=lambda), true_lmarglik)
    
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = array(dim=c(length(filt),nruns,length(lambda)))
      for(k in 1:length(lambda))
      {
         load(paste(dpath,"cv_pl-",lambda[k],"-",np[j],"-",nruns,"-",n.sim,".rdata",sep=""))
         attach(pf.out,warn.conflicts=F)
         lmargliks[[j]][,,k] = lmarglik[which(rownames(lmarglik) %in% filt),]
      }
    }

    # Plot kernel density estimates of of log marginal likelihoods for pfs compared with truth
    dmax = max(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) max(density(a)$y,na.rm=TRUE))))
    bmax = max(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) max(density(a)$x,na.rm=TRUE))))   
    bmin = min(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) min(density(a)$x,na.rm=TRUE))))
    cols = rainbow(length(filt))
    pdf(file=paste(gpath,"cv_pl_loglik-",paste(10*lambda,sep="",collapse="-"),"-",n.sim,"-",nruns,".pdf",sep=""),width=5*length(lambda),height=5*length(np))
    par(mfrow=c(length(np),length(lambda)),mar=c(5,6,4,2)+0.1)
    for(j in 1:length(np))
    {
      for(i in 1:length(lambda))
      {
        if(j == 1 & i == 1)
        {
          ylab = paste("J = ",np[j],sep="")
          xlab = "Log marginal likelihood"
          main = substitute(paste(lambda," = ",aa),list(aa=lambda[i]))
        } else if(j == 1) {
          xlab = ""
          ylab = ""
          main = substitute(paste(lambda," = ",aa),list(aa=lambda[i]))
        } else if(i == 1) {
          xlab = ""
          ylab = paste("J = ",np[j],sep="")
          main = ""
        }
        plot(density(lmargliks[[j]][1,,i]),lwd=2,col=cols[1],main=main,xlab=xlab,ylab=ylab,xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.5,cex.lab=1.75,cex.main=2)
        if(length(filt > 1)) for(k in 2:length(filt)) lines(density(lmargliks[[j]][k,,i]),lwd=2,col=cols[k])
        abline(v=true.lmarglik[i],lwd=1,col=1)
        for(l in which(!(1:length(lambda) %in% i))) abline(v=true.lmarglik[l],lwd=1,col=1,lty=2)
        mtext(round(true.lmarglik[i],2),side=1,at=true.lmarglik[i],cex=0.75)
        if(i == 1 & j == 1) legend("topleft",c(filt,"Truth","Truth (other models)"),lty=c(rep(1,length(filt)),1,2),lwd=c(rep(2,length(filt)),1,1),col=c(cols,1,1),cex=1.5)
      }
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(nruns = 20, nsims = 1)
m_ply(mydata, function(nruns, nsims) cv_pl_loglik_wsim(c(100,500,1000,5000), nruns, nsims, filt = c('pl','plrb','rm')))

cv_pl_comp_wsim = function(np, nruns, nsims, lambda = c(.5, 1, 2), filt = c('pl','plrb','kd','rm'))
{
  require(compositions)
  
  # Plot quantiles for each sim for each # particles
  for(n.sim in 1:nsims)
  {
    # Load simulated data
    load(paste(dpath,"rw_sim.rdata",sep=""))
    y = mysims[[n.sim]]$y
    
    # Calculate true log marginal likelihoods
    true_lmarglik = function(lambda)
    {
      post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
      return(cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b))
    }
    true.lmarglik = maply(data.frame(lambda=lambda), true_lmarglik)
    
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = array(dim=c(length(filt),nruns,length(lambda)))
      for(k in 1:length(lambda))
      {
        load(paste(dpath,"cv_pl-",lambda[k],"-",np[j],"-",nruns,"-",n.sim,".rdata",sep=""))
        attach(pf.out,warn.conflicts=F)
        lmargliks[[j]][,,k] = lmarglik[which(rownames(lmarglik) %in% filt),]
      }
    }
    
    # Ternary diagrams of posterior model probabilities
    # Compute true posterior model probabilities
    true.postModProbs = postModProbs(true.lmarglik, rep(1/length(lambda), length(lambda)))
    
    # Compute rm_pf approximate posterior model probabilities
    pmargliks = list()
    for(i in 1:length(np))
    {
      pmargliks[[i]] = apply(lmargliks[[i]], 1:2, function(a) postModProbs(a, rep(1/length(lambda), length(lambda))))
      pmargliks[[i]] = aperm(pmargliks[[i]], c(3,1,2))
    }

    # Ternary diagrams of posterior model probabilities
    require(compositions)
    wlabels = rep(NA,length(lambda))
    for(i in 1:length(lambda)) wlabels[i] = as.expression(bquote(paste(lambda," = ",sep="") ~ .(lambda[i])))
    pdf(file=paste(gpath,"cv_pl_ternary-",paste(10*lambda,sep="",collapse="-"),"-",n.sim,"-",nruns,".pdf",sep=""),width=10,height=10)
    par(mfrow=c(2,2))
    for(i in 1:length(np))
    {
      plot(acomp(pmargliks[[i]][,,1]),labels=wlabels, lwd=2)
      if(length(filt) > 1) for(j in 2:length(filt)) plot(acomp(pmargliks[[i]][,,j]),add=TRUE, col = j,lwd=2)
      plot(acomp(true.postModProbs),add=TRUE,col='gray',lwd=3)
      mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
      if(i==1) legend("topleft",legend=c(filt,"True Posterior"),pch=rep(1,length(filt)+1),col=c(1:length(filt),'gray'),pt.lwd=c(rep(2,length(filt)),3))
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(nruns = 20, nsims = 1)
m_ply(mydata, function(nruns, nsims) cv_pl_comp_wsim(c(100,500,1000,5000), nruns, nsims, filt = c('pl','plrb','kd','rm')))
