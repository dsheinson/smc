# pf_functions - functions to help calculate quantiles of particle filtered distributions and construct plots
#
# pf.quantile - function that calculates quantiles of filtered distributions; returns a 3-D array with dimensions number of observations, number of parameters and/or states, and number of quantiles desired
# Arguments:
# out - an 3-D dimensional array of particles with dimensions number of parameters and/or states, number of particles, and number of observations; this can be the state or theta component of the output list from functions bf, apf, or kd_pf
# wts - a matrix of wts with rows corresponding to particles and columns corresponding to observations; can be the weight component of the output list from functions bf, apf, or kd_pf; number of rows must match the 2nd dimension of out and number of columns must match the third dimension of out
# ftheta - a function that maps a vector from one scale to another, e.g. exp(theta) to map theta from log(scale) to original scale; can vary with parameter
# probs - a vector of quantiles to be calculated
# normwt - TRUE or FALSE corresponding to the normwt argument in wtd.quantile
pf.quantile = function(out, wts, ftheta, probs=.5, normwt=TRUE)
{
  require(Hmisc)
  nq = length(probs)
  np = dim(out)[1]
  tt = dim(out)[3]
  if(dim(wts)[1] != dim(out)[2] | dim(wts)[2] != tt) stop("dimensions of out and wts don't match")
  quantiles = array(NA,dim=c(tt,np,nq))
  for(i in 1:tt)
  {
    for(j in 1:np)
    {
      quantiles[i,j,] = wtd.quantile(ftheta(out[j,,i], j), wts[,i], normwt=normwt, probs=probs)
    }
  }
  return(quantiles)
}

# pf_plot - function that creates a grid of plots of credible intervals over time of filtered distributions based on particle filter output; number of particles along the rows, parameter along the columns and different pf's for which to compare within plot panels
# Arguments:
# n - integer vector of number of particles (rows)
# params - character vector of expression vector of parameters (columns)
# filt - character vector of particle filter types (different lines within plot panels)
# n.sim - integer, simulation number for which to load particle filter approximations
# probs - length-2 vector giving indices corresponding to which quantiles to use as lower and upper bounds (third dimension of output from pf.quantile)
# cols - vector of colors of lines in plots (should be same length as 'filt')
# create.label - character label for output file
# load.label - function that takes elements of n, filt, and n.sim as arguments and returns the file name from which to load particle filter output
# states - boolean, if TRUE, credible intevals for states are plotted, and if FALSE credible intervals for unknown parameters are plotted
# ymins, ymaxs - vector of minimum and maximum values of plot window within columns (should be same length as params); may be left missing
# out.avg - a np by tt by 2 matrix of average quantiles for which the area outside is shaded gray, with np number of parameters, tt number of time points, and third dimension corresponding to lower and upper quantiles, respectively; if missing, no area is shaded
# cex.lab, cex.main, cex.axis, cex.leg - expansion factors for plot labels - same as those in functions plot() and legend()
# pic.fac - factor by which to multiply the length of n and params to get the height and width, respectively, of output pdf file
# burn - vector of length equal to params, how many beginning time points to ignore when finding ymins and ymaxs (only used if ymins and ymaxs are missing)
pf_plot <- function(n, params, filt, n.sim, probs, cols, create.label, load.label, states = FALSE, ymins, ymaxs, out.avg, lwd = 4, cex.lab = 6, cex.main = 7, cex.axis = 4, cex.leg = 4, pic.fac = 10, burn = 0)
{
  # If missing ymins and ymaxs, find values to make plot windows consistent across columns (parameters)
  if(missing(ymins))
  {
    if(length(burn) != length(params)) burn = rep(0,length(params))
    mins = matrix(nr=0,nc=length(params))
    for(i in 1:length(n)){
    for(j in 1:length(filt)){
      load(load.label(filt[j], n[i], n.sim))
      if(states) out = pf.quant.out$state.quant else out = pf.quant.out$theta.quant
      min.k <- rep(NA, length(params))
      for(k in 1:length(params)) if(burn[k] > 0) min.k[k] = min(out[-(1:burn[k]),k,probs[1]]) else min.k[k] = min(out[,k,probs[1]])
      mins = rbind(mins, min.k)
    }}
    ymins = apply(mins, 2, min)
  }
  if(missing(ymaxs))
  {
    maxs = matrix(nr=0,nc=length(params))
    for(i in 1:length(n)){
    for(j in 1:length(filt)){
      load(load.label(filt[j], n[i], n.sim))
      if(states) out = pf.quant.out$state.quant else out = pf.quant.out$theta.quant
      max.k <- rep(NA, length(params))
      for(k in 1:length(params)) if(burn[k] > 0) max.k[k] = max(out[-(1:burn[k]),k,probs[2]]) else max.k[k] = max(out[,k,probs[2]])
      maxs = rbind(maxs, max.k)
    }}
    ymaxs = apply(maxs, 2, max)
  }

  # Construct plots
  pdf(create.label, width = pic.fac*length(params), height = pic.fac*length(n))
  par(mfrow=c(length(n),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
  for(i in 1:length(n))
  {
    for(k in 1:length(params))
    {
      for(j in 0:length(filt))
      {
        if(j > 0)
        {
          load(load.label(filt[j], n[i], n.sim))
          if(states) out = pf.quant.out$state.quant else out = pf.quant.out$theta.quant
          tt = dim(out)[1]; nt = tt - 1
          quant = out[,k,]
          x = 0:nt
        }
        
        if(j == 0) # call plot function
        {
          if(!missing(out.avg))
          {
            x = 0:(dim(out.avg)[2]-1)
            if(k == 1 & i == 1) # label y axis and title
            {
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            } else if(k == 1 & i == length(n)) { # label x and y axes
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            } else if(k == 1) { # label y axis only
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            } else if(i == 1) { # label title only
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            } else if(i == length(n)) { # label x axis only
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            } else { # label nothing
              plot(x,out.avg[k,,1],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col="gray",xlab="",ylab="",cex.axis=cex.axis)
              lines(x,out.avg[k,,2],col="gray",lwd=lwd)
            }          
            y = out.avg[k,,1]
            polygon(c(x[length(x)],x[1],x[1],x,x[length(x)]),c(ymins[k],ymins[k],y[1],y,y[length(y)]),col="gray",border=NA)
            y = out.avg[k,,2]
            polygon(c(x[length(x)],x[1],x[1],x,x[length(x)]),c(ymaxs[k],ymaxs[k],y[1],y,y[length(y)]),col="gray",border=NA)
          } else {} # do nothing
        } else if(j == 1 & missing(out.avg)) { # call plot function
          if(k == 1 & i == 1) # label y axis and title
          {
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          } else if(k == 1 & i == length(n)) { # label x and y axes
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          } else if(k == 1) { # label y axis only
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          } else if(i == 1) { # label title only
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          } else if(i == length(n)) { # label x axis only
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          } else { # label nothing
            plot(x,quant[,probs[1]],type="l",lwd=lwd,ylim=c(ymins[k],ymaxs[k]),col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
            lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
          }
        } else { # lines only
          lines(x,quant[,probs[1]],col=cols[j],lwd=lwd)
          lines(x,quant[,probs[2]],col=cols[j],lwd=lwd)
        }
      }
      
      # Add truth
      load("../data/sim-orig.rdata")
      if(states)
      {
        if(k < 3) lines(x, mysims[[n.sim]]$sim$x[k,], col="gray47") else lines(x, 1 - mysims[[n.sim]]$sim$x[1,] - mysims[[n.sim]]$sim$x[2,], col="gray47")
      } else {
        abline(h=mysims[[n.sim]]$true.params$theta[k],col="gray47",lwd=6)
      }
      
      # Whiten borders
      if(!missing(out.avg))
      {
        edgex = 0
        edgey = 1
        polygon(c(x[length(x)],x[1],x[1],x[length(x)]),c(ymins[k],ymins[k],ymins[k]-edgey,ymins[k]-edgey),col="white",border=NA) # bottom border
        polygon(c(x[length(x)],x[1],x[1],x[length(x)]),c(ymaxs[k],ymaxs[k],ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # top border
        polygon(c(x[length(x)],x[length(x)]+edgex,x[length(x)]+edgex,x[length(x)]),c(ymins[k]-edgey,ymins[k]-edgey,ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # right border
        polygon(c(x[1],x[1]-edgex,x[1]-edgex,x[1]),c(ymins[k]-edgey,ymins[k]-edgey,ymaxs[k]+edgey,ymaxs[k]+edgey),col="white",border=NA) # left border
        box()
      }
      
      # Add legend
      if(k == 1 & i == 1)
      {
        legend("topright",legend=c("Truth",filt),col=c("gray47",cols),lty=c(1,rep(1,length(filt))),lwd=c(6,rep(1,length(filt))),bg="white",cex=cex.leg)
      }
    }
  }
  dev.off()
}

# pf_coverage - function to calculate the proportion of times credible intervals cover the truth (for each parameter, at each time point)
# Arguments:
# n.sims - integer, number of simulated data sets over which to calculate proportions
# n - integer, number of particles
# filt - character string, type of particle filter
# probs - length-2 vector giving indices corresponding to which quantiles to use as lower and upper bounds (third dimension of output from pf.quantile)
# load.label - function that takes filt, n, and simulation number as arguments and returns file name from which to load particle filter data
# states - TRUE or FALSE, indicating whether or not coverage probabilities for the states or unknown parameters should be calculated
pf_coverage <- function(n.sims, n, filt, probs, load.label, states = FALSE, mod = "orig")
{
  if(mod == "orig") load("../data/sim-orig.rdata")
  if(mod == "ext") load("../data/sim-ext.rdata")
  load(load.label(filt, n, 1))
  if(states) out = pf.quant.out$state.quant else out = pf.quant.out$theta.quant
  tt = dim(out)[1]
  n.params = dim(out)[2]
  covered = array(NA, dim = c(n.sims, n.params, tt))
  for(i in 1:n.sims)
  {
    if(mod == "orig") theta = mysims[[i]]$true.params$theta
    if(mod == "ext") theta = c(mysims[[i]]$true.params$theta, mysims[[i]]$true.params$b, mysims[[i]]$true.params$varsigma, mysims[[i]]$true.params$sigma, mysims[[i]]$true.params$eta)
    
    load(load.label(filt, n, i))
    if(states){
      out = pf.quant.out$state.quant
      for(k in 1:n.params)
      { 
        if(k < 3) covered[i,k,] = out[,k,probs[1]] < mysims[[i]]$sim$x[k,] & out[,k,probs[2]] > mysims[[i]]$sim$x[k,] else covered[i,k,] = out[,k,probs[1]] < (1 -  mysims[[i]]$sim$x[1,] - mysims[[i]]$sim$x[2,]) & out[,k,probs[2]] > (1 -  mysims[[i]]$sim$x[1,] - mysims[[i]]$sim$x[2,])
      }
    } else {
      out = pf.quant.out$theta.quant
      for(k in 1:n.params)
      {
        covered[i,k,] = out[,k,probs[1]] < theta[k] & out[,k,probs[2]] > theta[k]
      }
    }
  }
  return(apply(covered, c(2, 3), mean))
}

# pf_coverage_plot - function to plot coverage probabilities over time for particle filters for each parameter (columns) with increasing number of particles (rows)
# Arguments:
# coverage - 4-dimensional array where first two dimensions are different number of particles and types of particle filter, and last two dimensions are different parameters and total time points
# alpha - scalar, nominal proportion of time truth should be covered, to be plotted as horizontal gray line
# n.sim - scalar, how many simulations for each dimension of coverage were run?
# params - vector of parameter names, should be of length equal to third dimension of coverage
# cols - vector of colors of lines in plots (should be same length as 'filt')
# create.label - character label for output file
# ymins, ymaxs - vector of minimum and maximum values of plot window within columns (should be same length as params); may be left missing
# cex.lab, cex.main, cex.axis, cex.leg - expansion factors for plot labels - same as those in functions plot() and legend()
# pic.fac - factor by which to multiply the length of n and params to get the height and width, respectively, of output pdf file
# burn - vector of length equal to params, how many beginning time points to ignore when finding ymins and ymaxs (only used if ymins and ymaxs are missing)
pf_coverage_plot <- function(coverage, alpha, n.sim, params, cols, create.label, ymins, ymaxs, lwd = 4, cex.lab = 6, cex.main = 7, cex.axis = 4, cex.leg = 4, leg.location = "topright", pic.fac = 10, burn = 0)
{
  n <- dimnames(coverage)[[1]]
  filt <- dimnames(coverage)[[2]]

  # If missing ymins and ymaxs, find values to make plot windows consistent across columns (parameters)
  if(missing(ymins))
  {
    if(length(burn) != length(params)) burn = rep(0,length(params))
    mins = matrix(nr=0,nc=length(params))
    for(i in 1:length(n)){
    for(j in 1:length(filt)){   
      min.k <- rep(NA, length(params))
      for(k in 1:length(params)) if(burn[k] > 0) min.k[k] = min(coverage[i,j,k,-(1:burn[k])]) else min.k[k] = min(coverage[i,j,k,])
      mins = rbind(mins, min.k)
    }}
    ymins = apply(mins, 2, min)
  }
  if(missing(ymaxs))
  {
    maxs = matrix(nr=0,nc=length(params))
    for(i in 1:length(n)){
    for(j in 1:length(filt)){
      max.k <- rep(NA, length(params))
      for(k in 1:length(params)) if(burn[k] > 0) max.k[k] = max(coverage[i,j,k,-(1:burn[k])]) else max.k[k] = max(coverage[i,j,k,])
      maxs = rbind(maxs, max.k)
    }}
    ymaxs = apply(maxs, 2, max)
  }

  # Construct plots
  pdf(create.label, width = pic.fac*length(params), height = pic.fac*length(n))
  par(mfrow=c(length(n),length(params)),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
  for(i in 1:length(n))
  {
    for(k in 1:length(params))
    {
      for(j in 1:length(filt))
      {
        tt = dim(coverage)[4]; nt = tt - 1
        cov = coverage[i,j,k,]
        x = 0:nt
        if(j == 1) # call plot function
        {
          if(k == 1 & i == 1) # label y axis and title
          {
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
          } else if(k == 1 & i == length(n)) { # label x and y axes
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="Time (days)",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          } else if(k == 1) { # label y axis only
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="",ylab=paste("J = ",n[i],sep=""),cex.lab=cex.lab,cex.axis=cex.axis)
          } else if(i == 1) { # label title only
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
          } else if(i == length(n)) { # label x axis only
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="Time (days)",ylab="",cex.lab=cex.lab,cex.axis=cex.axis)
          } else { # label nothing
            plot(x,cov,type="b",cex=2.5,ylim=c(ymins[k],ymaxs[k]),lwd=lwd,col=cols[j],xlab="",ylab="",cex.axis=cex.axis)
          }
        } else { # lines only
          lines(x,cov,type="b",cex=2.5,lwd=lwd,col=cols[j])
        }
      }

      # Add nominal coverage level line
      abline(h = alpha, col = "gray", lwd = 6)
      me = 1.96*sqrt(alpha*(1-alpha)/n.sim)
      abline(h = alpha - me, col = "gray", lwd = 6, lty = 2)
      abline(h = alpha + me, col = "gray", lwd = 6, lty = 2)

      if(k == 1 & i == 1) # add legend
      {
        legend(leg.location,legend=c("Nominal level",filt),col=c("gray47",cols),lty=c(1,rep(1,length(filt))),lwd=c(6,rep(1,length(filt))),bg="white",cex=cex.leg)
      }
    }
  }
  dev.off()
}

# pf.scat - function to resample particles, compute maximum/minimum values, and compute correlations between two parameters at specified time points; returns a list with components xrw1 and xrw2 that are M x length(cutoff) matrices of resampled particles (row = particle, column = time point), xlim and ylim that are each two element vectors that are the minimum and maximum values over all time points of xrw1 and xrw2, respectively, and r that is a length(cutoff) length vector of correlations between xrw1 and xrw2 at each time point
# Arguments:
# out - 2 by n by tt array of particles (assumed on original scale)
# wts - n by tt matrix of particle weights
# cutoff - k length vector of times at which to plot histograms
# M - number of samples to use to resample particles to get equal weights for constructing histograms
# xlim - length 2 vector giving xlim argument in plot(); if any elements NA, are calculated by default
# ylim - length 2 vector giving ylim argument in plot(); if any elements NA, are calculated by default
# seed - if not null, sets the random seed prior resampling particles
# ... - additional arguments passed to resample()
pf.scat = function(out, wts, cutoff, M=500, xlim=NA, ylim=NA, seed = NULL, ...)
{
  require(smcUtils)

  # Check dimensions
  n = dim(wts)[1]
  tt = dim(wts)[2]
  if(n != dim(out)[2] | tt != dim(out)[3]) stop("dimensions of out and wts don't match")
  nd = length(cutoff)
  if(nd > tt) stop("Too many cutoff points")
  if(dim(out)[1] != 2) stop("Incorrect first dimension of out")

  # Set seed
  if(!is.null(seed)) set.seed(seed)
  
  # Resample particles to get equal weights
  tmps = matrix(NA,nr=M,nc=length(cutoff))
  for(i in 1:length(cutoff))
  {
    tmps[,i] = resample(wts[,cutoff[i]], M, nonuniformity="none", ...)$indices
  }

  # Find minimum and maximum values so scatterplots on same scale
  xrw1 = matrix(NA,nr=M,nc=length(cutoff))
  xrw2 = matrix(NA,nr=M,nc=length(cutoff))
  if(!all(!is.na(xlim)))
  {
    xlim = c(Inf,-Inf)
    for(i in 1:length(cutoff))
    {
      xrw1[,i] = out[1,tmps[,i],cutoff[i]]
      xlim[1] = min(xlim[1],xrw1[,i])
      xlim[2] = max(xlim[2],xrw1[,i])
    }
  } else {
    for(i in 1:length(cutoff))
    {
      xrw1[,i] = out[1,tmps[,i],cutoff[i]]
    }
  }
  if(!all(!is.na(ylim)))
  {
    ylim = c(Inf,-Inf)
    for(i in 1:length(cutoff))
    {
      xrw2[,i] = out[2,tmps[,i],cutoff[i]]
      ylim[1] = min(ylim[1],xrw2[,i])
      ylim[2] = max(ylim[2],xrw2[,i])
    }
  } else {
    for(i in 1:length(cutoff))
    {
      xrw2[,i] = out[2,tmps[,i],cutoff[i]]
    }
  }

  # Compute correlation between particles at every cutoff point
  r = rep(NA,length(cutoff))
  for(i in 1:length(cutoff))
  {
    r[i] = cov.wt(cbind(out[1,,cutoff[i]],out[2,,cutoff[i]]),wts[,cutoff[i]],cor=TRUE)$cor[1,2]
  }  

  # Return re-weighted particles, correlations, and min/max values of particles
  return(list(xrw1=xrw1,xrw2=xrw2,xlim=xlim,ylim=ylim,r=r))
}

# pf.hist - function to plot histograms of unknown parameters at specified time points for particle filtered samples
# Arguments:
# out - np by n by tt array of particles
# wts - n by tt matrix of particle weights
# cutoff - k length vector of times at which to plot histograms
# M - number of samples to use to resample particles to get equal weights for constructing histograms
# ftheta - function that maps an n by tt matrix of particles to original scale and takes an additional integer argument to optionally depend on parameter
# plabs - vector of length np giving labels for parameters
# truth - np length vector of true parameter values
# tsize - size of subtitle displaying true parameter values
# mr - number of rows in plot panel window
# mc - number of columns in plot panel window
# mar - corresponds to mar argument in par() for resizing plot windows
# msize - size of titles of upper-left histogram
# labsize - size of axis labels on upper-left histogram
# file - np length character string vector of names of outputted pdf files
# width, height - arguments to pdf()
# ... - additional arguments passed to resample()
pf.hist = function(out,wts,cutoff,ftheta,plabs,truth,file,M=10000,tsize=.65,mr=1,mc=1,mar=c(4,5.2,3,.5)+.1,msize=1.5,labsize=1.5,width=10,height=5,...)
{
  require(smcUtils)

  # Check dimensions
  n = dim(wts)[1]
  tt = dim(wts)[2]
  np = dim(out)[1]
  if(n != dim(out)[2] | tt != dim(out)[3]) stop("dimensions of out and wts don't match")
  nd = length(cutoff)
  if(nd > tt) stop("Too many cutoff points")
  if(length(plabs) != np) stop("labels should be length equal to number of parameters")

  # Resample particles to get equal weights
  tmps = matrix(NA,nr=M,nc=length(cutoff))
  for(i in 1:length(cutoff))
  {
    tmps[,i] = resample(wts[,cutoff[i]], M, nonuniformity="none",...)$indices
  }

  # Plot histograms of unknown parameters over time
  for(j in 1:np)
  {
    x = ftheta(out[j,,cutoff],j)
    xrw = matrix(NA,nr=M,nc=length(cutoff))
    for(i in 1:length(cutoff)) xrw[,i] = x[tmps[,i],i]
    xmin = min(apply(xrw,2,function(x) min(hist(x,plot=FALSE)$breaks)))
    xmax = max(apply(xrw,2,function(x) max(hist(x,plot=FALSE)$breaks)))
    pdf(file[j],width=width,height=height)
    par(mfrow=c(mr,mc),mar=mar) # dimensions should depend on cutoff
    for(i in 1:length(cutoff))
    {
      if(i==1){
        hist(xrw[,i],xlim=c(xmin,xmax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab=plabs[j],cex.main=msize,cex.lab=labsize)
        abline(v=truth[j],col=2,lwd=2)
        mtext(paste("Truth = ",truth[j],sep=""),side=3,cex=tsize)
      } else {
        hist(xrw[,i],xlim=c(xmin,xmax),freq=FALSE,main=paste("t = ",cutoff[i]-1,sep=""),xlab="",ylab="",cex.main=msize)
        abline(v=truth[j],col=2,lwd=2)
      }
    }
    dev.off()
  }
}
