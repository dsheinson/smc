# Construct plots for particle filter runs
source("dlm_cv_functions.r")
source("pf_mc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

cv_pl_quantiles <- function(lambda, np, nrtot, nruns, nsims, filt, alpha = 0.05, burn = 0, ymin, ymax)
{
  # Plot quantiles for each sim for each # particles
  for(n.sim in nsims)
  {
    # Load simulated data
    load(paste(dpath,"rw_sim.rdata",sep=""))
    y = mysims[[n.sim]]$y
    nt = dim(y)[2]
    
    # Calculate true 95% CI of states and precision
    post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
    lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
    uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
    lp = qgamma(alpha/2,post$a,post$b)
    up = qgamma(1-alpha/2,post$a,post$b)
    
    # Find ymax and ymin if missing
    if(missing(ymin) | missing(ymax))
    {
      ymin = c(Inf, Inf)
      ymax = c(-Inf, -Inf)
      for(j in 1:length(np))
      {
        load(paste(dpath,"cv_pl-",lambda,"-",np[j],"-",nrtot,"-",n.sim,".rdata",sep=""))
        attach(pf.out,warn.conflicts=F)
        if(missing(filt))
        {
          filt = rownames(lmarglik)
          filt.ind = 1:length(filt)
        } else {
          filt.ind = which(rownames(lmarglik) %in% filt)
        }
        for(k in filt.ind)
        {
          if(burn > 0)
          {
            g.lk = lk[-(1:burn)]
            g.uk = uk[-(1:burn)]
            g.ls = sapply(state.quant, function(x) x[[k]][-(1:burn),1,1])
            g.us = sapply(state.quant, function(x) x[[k]][-(1:burn),1,2])
            g.lp = lp[-(1:burn)]
            g.up = up[-(1:burn)]
            g.lf = sapply(theta.quant, function(x)  x[[k]][-(1:burn),1,1])
            g.uf = sapply(theta.quant, function(x)  x[[k]][-(1:burn),1,2])     
          } else {
            g.lk = lk
            g.uk = uk
            g.ls = sapply(state.quant, function(x) x[[k]][,1,1])
            g.us = sapply(state.quant, function(x) x[[k]][,1,2])
            g.lp = lp
            g.up = up
            g.lf = sapply(theta.quant, function(x)  x[[k]][,1,1])
            g.uf = sapply(theta.quant, function(x)  x[[k]][,1,2])
          }
          ymin[1] = min(ymin[1], mysims[[n.sim]]$x[1,], g.lk, apply(g.ls,2,min))
          ymin[2] = min(ymin[2], 1, g.lp, apply(g.lf,2,min))
          ymax[1] = max(ymax[1], mysims[[n.sim]]$x[1,], g.uk, apply(g.us,2,max))
          ymax[2] = max(ymax[2], 1, g.up, apply(g.uf,2,max))
        }
      }
    }
    
    pdf(paste(gpath,"cv-pl-quant-",10*lambda,"-",n.sim,"-",nrtot,"-",nruns,".pdf",sep=""),width=10,height=5*length(np))
    par(mfrow=c(length(np),2), mar=c(5,9,7,2)+0.1,mgp=c(6,1,0))
    for(j in 1:length(np))
    {
      # Load data
      load(paste(dpath,"cv_pl-",lambda,"-",np[j],"-",nrtot,"-",n.sim,".rdata",sep=""))
      attach(pf.out,warn.conflicts=F)
      if(missing(filt))
      {
        filt = rownames(lmarglik)
        filt.ind = 1:length(filt)
      } else {
        filt.ind = which(rownames(lmarglik) %in% filt)
      }
      
      # If only one rm in filt, delete lag from label
      filt.names = filt
      firsttwo = sapply(strsplit(filt,""), function(x) paste(x[1:2],sep="",collapse=""))
      ind = which(firsttwo == 'rm')
      if(length(ind) == 1) filt.names[ind] = 'rm'
      
      if(j == 1)
      {
        xlab = c(expression(t),"")
        ylab = c(paste("J = ",np[j],sep=""),"")
        main = paste("95% CI for Filtered ", c("States","Precision"), sep="") 
      } else {
        xlab = c("","")
        ylab = c(paste("J = ",np[j],sep=""),"")
        main = c("","")
      }
      
      # Plot 95% CI for states
      if(j == 1){
        plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(ymin[1],ymax[1]),type="l",main=main[1],xlab=xlab[1],ylab=ylab[1],col='gray',lwd=3, cex.lab = 3, cex.main = 2.65, cex.axis = 1.9)
      } else {
        plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(ymin[1],ymax[1]),type="l",axes=FALSE,main=main[1],xlab=xlab[1],ylab=ylab[1],col='gray',lwd=3, cex.lab = 3, cex.main = 2.65, cex.axis = 1.9)
        box()
      }
      for(i in 1:length(filt.ind))
      {
        for(k in 1:nruns)
        {
          lines(0:nt,state.quant[[k]][[filt.ind[i]]][,1,1],col=i+1)
          lines(0:nt,state.quant[[k]][[filt.ind[i]]][,1,2],col=i+1)
        }
      }
      lines(0:nt,mysims[[n.sim]]$x[1,],lwd=3,col="gray")
      lines(0:nt,lk,lwd=3)
      lines(0:nt,uk,lwd=3)
      if(j == 1)
      {
        mtext(expression(x[t]),side=2,cex=2,line=2.75)
        mtext(expression(t),side=1,cex=2,line=3.5)
      }
      
      # Plot 95% CI for precision
      plot(0:nt,lp,type="l",lwd=3,axes=ifelse(j==1,TRUE,FALSE),ylim=c(ymin[2],ymax[2]),main=main[2],xlab=xlab[2],ylab=ylab[2], cex.lab = 3, cex.main = 2.65, cex.axis = 1.9)
      if(j != 1) box()
      for(i in 1:length(filt.ind))
      {
        for(k in 1:nruns)
        {
          lines(0:nt,theta.quant[[k]][[filt.ind[i]]][,1,1],col=i+1)
          lines(0:nt,theta.quant[[k]][[filt.ind[i]]][,1,2],col=i+1)
        }
      }
      lines(0:nt,lp,lwd=3)
      lines(0:nt,up,lwd=3)
      if(j == 1) mtext(expression(1/theta),side=2,cex=2,line=2.75)
      abline(h=1,col='gray',lwd=3)
      if(j == 1) legend("topright",legend=c(filt.names, "True Post", "True Sim"),lty=c(rep(1,length(filt.names)),1,1),col=c(2:(length(filt.names)+1),1,'gray'),lwd=c(rep(1,length(filt.names)),3,3), cex=1.75)
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(lambda = c(.5,1,2), nrtot=20, nruns=5, stringsAsFactors=FALSE)
m_ply(mydata, function(lambda, nrtot, nruns) cv_pl_quantiles(lambda, np = c(100,500,1000,5000), nrtot, nruns, 1, c('pl', 'kd', 'rm100'), burn = 5))

cv_pl_loglik_wsim = function(np, nruns, nsims, lambda = c(.5, 1, 2), filt = c('pl','plrb','kd','rm'), blim, dmax)
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
    
    # If only one rm in filt, delete lag from label
    firsttwo = sapply(strsplit(filt,""), function(x) paste(x[1:2],sep="",collapse=""))
    ind = which(firsttwo == 'rm')
    if(length(ind) == 1) filt[ind] = 'rm'
    
    # Plot kernel density estimates of of log marginal likelihoods for pfs compared with truth
    if(missing(dmax)) dmax = max(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) max(density(a)$y,na.rm=TRUE))))
    if(missing(blim))
    {
      bmax = max(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) max(density(a)$x,na.rm=TRUE))))   
      bmin = min(sapply(lmargliks, function(x) apply(x, c(1,3), function(a) min(density(a)$x,na.rm=TRUE))))
    } else {
      bmin = blim[1]
      bmax = blim[2]
    }
    cols = rainbow(length(filt))
    pdf(file=paste(gpath,"cv_pl_loglik-",paste(10*lambda,sep="",collapse="-"),"-",n.sim,"-",nruns,".pdf",sep=""),width=5*length(lambda),height=5*length(np))
    par(mfrow=c(length(np),length(lambda)),mar=c(8,9,7,2)+0.1,mgp=c(5.5,1,0))
    for(j in 1:length(np))
    {
      for(i in 1:length(lambda))
      {
        if(j == 1 & i == 1)
        {
          ylab = paste("J = ",np[j],sep="")
          xlab = expression(log(p(y[1:T])))
          main = substitute(paste(lambda," = ",aa),list(aa=lambda[i]))
        } else if(j == 1) {
          xlab = ""
          ylab = ""
          main = substitute(paste(lambda," = ",aa),list(aa=lambda[i]))
        } else if(i == 1) {
          xlab = ""
          ylab = paste("J = ",np[j],sep="")
          main = ""
        } else {
          xlab = ylab = main = ""
        }
        plot(density(lmargliks[[j]][1,,i]),axes=ifelse(j==1&i==1,TRUE,FALSE),lwd=2,col=cols[1],main=main,xlab=xlab,ylab=ylab,xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.9,cex.lab=4,cex.main=4)
        if(!(j == 1 & i == 1)) box()
        if(j == 1 & i == 1) mtext("Density",side=2,line=2.5,cex=2)
        if(length(filt > 1)) for(k in 2:length(filt)) lines(density(lmargliks[[j]][k,,i]),lwd=2,col=cols[k])
        abline(v=true.lmarglik[i],lwd=1,col=1)
        for(l in which(!(1:length(lambda) %in% i))) abline(v=true.lmarglik[l],lwd=1,col=1,lty=2)
        mtext(round(true.lmarglik[i],2),side=3,at=true.lmarglik[i],cex=1.5)
        if(i == 1 & j == 1) legend("topleft",c(filt,"Truth","Truth (others)"),lty=c(rep(1,length(filt)),1,2),lwd=c(rep(2,length(filt)),1,1),col=c(cols,1,1),cex=1.6)
      }
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(nruns = 20, nsims = 1)
m_ply(mydata, function(nruns, nsims) cv_pl_loglik_wsim(c(100,500,1000,5000), nruns, nsims, filt = c('pl','kd','rm100'), blim = c(-200,-190)))

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
    
    # If only one rm in filt, delete lag from label
    firsttwo = sapply(strsplit(filt,""), function(x) paste(x[1:2],sep="",collapse=""))
    ind = which(firsttwo == 'rm')
    if(length(ind) == 1) filt[ind] = 'rm'
    
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
    par(mfrow=c(2,2),mar=c(5,7,4,4)+0.1,cex.lab=2)
    for(i in 1:length(np))
    {
      plot(acomp(pmargliks[[i]][,,1]),labels=wlabels, col=2, lwd=2)
      if(length(filt) > 1) for(j in 2:length(filt)) plot(acomp(pmargliks[[i]][,,j]),add=TRUE, col = j+1,lwd=2)
      plot(acomp(true.postModProbs),add=TRUE,col='gray',lwd=3)
      mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
      if(i==1) legend("topleft",legend=c(filt,"True Posterior"),pch=rep(1,length(filt)+1),col=c(2:(length(filt)+1),'gray'),pt.lwd=c(rep(2,length(filt)),3))
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(nruns = 20, nsims = 1)
m_ply(mydata, function(nruns, nsims) cv_pl_comp_wsim(c(100,500,1000,5000), nruns, nsims, filt = c('pl','kd','rm100')))
