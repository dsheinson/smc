#' A function to perform the resample-move particle filter
#'
#' @param y a matrix with no, dimension of each observation, rows and nt, number of observation time points, columns
#' @param dllik a function to evaluate the logarithm of the likelihood given y the current observation, x the current state, and theta the fixed parameters
#' @param revo a function to propagate the state given the current state x and the fixed parameters theta
#' @param rprior a function to sample from the prior for the state and fixed parameters; takes an argument j particle index in case function loads previously drawn prior samples; returns a list with elements x and theta
#' @param rmove a function to regenerate values of the state history and/or fixed parameters with arguments y, x, and theta; y is an no by k < nt matrix of observations up to the current time point k, x is an ns (dimension of state) by k < nt+1 matrix of states up to the current time point, and theta is the fixed parameters; returns a list with components state, an ns by k matrix of (possibly) regenerated states, and theta, the (possibly) regenerated values of the fixed parameters
#' @param n the number of particles
#' @param store.all should state histories for every time point be saved?
#' @param store.filter should only filtered states be saved?
#' @param progress a boolean to display progress bar if TRUE
#' @param ... arguments passed on to resample
#' @return a list containing, if store.all = TRUE, an (nt+1)-length list of state histories (each an ns by n by k matrix), an n by (nt+1) matrix of normalized particle weights, an np by n by (nt+1) array of theta draws, an n by nt matrix of unnormalized particle weights (increments), and an n by nt parent matrix. If store.all = FALSE, only state trajectories at final time point are saved, and only filtered states at each time point are saved if store.filter = TRUE.
#' @references Berzuini, C. and Gilks, W. Following a Moving Target-Monte Carlo Inference for Dynamic Bayesian Models. Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 63, No. 1 (2001), pp. 127-146
#' @seealso \code{\link{resample}}
#'
rm_pf = function(y, dllik, revo, rprior, rmove, n, store.all = FALSE, store.filter = TRUE, progress = TRUE, ...)
{
  require(smcUtils)

  if (!is.matrix(y)) y = matrix(y, 1)
  no = nrow(y) # not currently used
  nt = ncol(y)

  # Find dimension of state
  current.seed = .Random.seed
  tmp = rprior(1)
  ns = length(tmp$x)
  np = length(tmp$theta)
  .Random.seed = current.seed

  # Set up initial state
  state = array(NA, dim=c(ns,n,nt+1))
  theta = array(NA, dim=c(np,n,nt+1))
  for (j in 1:n) 
  {
    tmp = rprior(j)
    state[,j,1] = tmp$x
    theta[,j,1] = tmp$theta
  }
  if(store.all)
  {
    state.all = list()
    state.all[[1]] = array(state[,,1], dim=c(ns,n,1))
  } else if(store.filter) {
    state.filter = array(NA, dim = c(ns,n,nt+1))
    state.filter[,,1] = state[,,1]
  }
    
  # Initialize weights
  weight = matrix(NA, n, nt+1)
  increment = matrix(NA, n, nt)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)

  # Run particle filter
  if(progress) pb = txtProgressBar(0,nt,style=3)
  for (i in 1:nt) 
  {
    if(progress) setTxtProgressBar(pb,i)
    # Augmentation and update weights
    for(j in 1:n)
    {
      state[,j,i+1] = revo(state[,j,i],theta[,j,i])
      increment[j,i] = dllik(y[,i],state[,j,i+1],theta[,j,i])
      weight[j,i+1] = log(weight[j,i]) + increment[j,i]
    }
    
    # Resample particles
    tmp = resample(renormalize(weight[,i+1],log=TRUE), ...)
    weight[,i+1] = tmp$weight
    kk = tmp$indices
    did.resample = !(all(kk == 1:n))

    # Move particles
    if(did.resample)
    {
      for (j in 1:n) 
      {
        tmp2 = rmove(y[,1:i],state[,kk[j],1:(i+1)],theta[,kk[j],i])
        state[,j,1:(i+1)] = tmp2$state
        theta[,j,i+1] = tmp2$theta
      }
    } else {
      theta[,,i+1] = theta[,,i]
    }

    if(store.all)
    { 
      state.all[[i+1]] = array(state[,,1:(i+1)], dim=c(ns,n,i+1))
    } else if(store.filter) {
      state.filter[,,i+1] = state[,,i+1]
    }
    
    parent[,i+1] = kk
  }

  # Return entire state trajectories, filtered states, or smoothed states
  if(store.all)
  {
    state = state.all
  } else if(store.filter) {
    state = state.filter
  }
  
  return(list(state=state, theta=theta, weight=weight, increment=increment, parent=parent))
}
