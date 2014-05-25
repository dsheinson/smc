# pl - a function to perform particle learning, optionally with Rao-Blackwellization
# Inputs:
# y a matrix with no, dimension of each observation, rows and nt, number of observation time points, columns
# dlpred a function to evaluate the logarithm of the conditional predictive likelihood of y given the current state x (or sufficient stat suff.x) and fixed parameters theta
# revo a function to propagate the state given the current observation y, current state x (or sufficient stat suff.x), and fixed parameters theta
# rprior a function to sample from the prior and initialize the conditional sufficient statistics for the state and fixed parameters given the particle index j (needed if initial prior samples are already available); returns a list with elements x, theta, suff.x, and suff.theta
# rmove a function to sample the fixed parameters theta given previous theta and the sufficient statistics suff.theta
# smap.theta a function to update the sufficient statistics for theta given the current values of the sufficient statistics suff.theta, current observation y, new state x.new, and current state x.curr
# smap.state a function to update the sufficient statistics for the state given the current values of the sufficient statistics suff.x, current observation y, and fixed parameters theta
# n the number of particles
# progress indicating if a progress bar should be displayed
# ... additional arguments passed into the 'resample' function

pl = function(y, dlpred, revo, rprior, rmove, smap.theta, smap.state, n, progress = TRUE, ...)
{
  require(smcUtils)
  
  if (!is.matrix(y)) y = matrix(y, 1)
  no = nrow(y) # not currently used
  nt = ncol(y)
  
  # Find dimension of state, parameters, sufficient statistics
  current.seed = .Random.seed
  tmp = rprior()
  ns = length(tmp$x)
  np = length(tmp$theta)
  nfs = length(tmp$suff.x)
  nft = length(tmp$suff.theta)
  .Random.seed = current.seed
  
  # Set up initial state, parameters, sufficient statistics
  state = array(NA, dim=c(ns,n,nt+1))
  theta = array(NA, dim=c(np,n,nt+1))
  suff.state = array(NA, dim=c(nfs,n,nt+1))
  suff.theta = array(NA, dim=c(nft,n,nt+1))
  for (j in 1:n) 
  {
    tmp = rprior(j)
    state[,j,1] = tmp$x
    theta[,j,1] = tmp$theta
    suff.state[,j,1] = tmp$suff.x
    suff.theta[,j,1] = tmp$suff.theta
  }
  
  # Initialize normalized and incremental weights
  weight = matrix(NA, n, nt+1)
  increment = matrix(NA, n, nt)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)
  
  # Run particle filter
  p.weights = numeric(n)
  if(progress) pb = txtProgressBar(0,nt,style=3)
  for (i in 1:nt) 
  {
    if(progress) setTxtProgressBar(pb,i)
    
    # Calculate weights
    for(j in 1:n)
    {
      increment[j,i] = dlpred(y[,i],state[,j,i],suff.state[,j,i],theta[,j,i])
      p.weights[j] = log(weight[j,i]) + increment[j,i]
    }
    
    p.weights = renormalize(p.weights, log=T)
    
    # Resample
    tmp = resample(p.weights, ...)
    kk = tmp$indices
    weight[,i+1] = tmp$weights
    did.resample = !(all(kk == 1:n))
    
    # Propagate state and update sufficient statistics
    for(j in 1:n)
    {
      state[,j,i+1] = revo(y[,i],state[,kk[j],i],suff.state[,kk[j],i],theta[,kk[j],i])
      suff.theta[,j,i+1] = smap.theta(suff.theta[,kk[j],i], y[,i], state[,j,i+1], state[,kk[j],i])
    }

    for(j in 1:n) suff.state[,j,i+1] = smap.state(suff.state[,kk[j],i],y[,i],theta[,kk[j],i])
    
    # Move particles
    if(did.resample)
    {
      for (j in 1:n) theta[,j,i+1] = rmove(theta[,kk[j],i],suff.theta[,j,i+1])
    } else {
      theta[,,i+1] = theta[,kk,i]
    }
    
    # Track parent particles
    parent[,i+1] = kk
  }
  
  return(list(state=state, suff.state=suff.state, theta=theta, suff.theta=suff.theta, weight=weight, increment=increment, parent=parent))
}