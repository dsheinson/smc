#' A function to perform the resample-move particle filter
#'
#' @param y a matrix with no, dimension of each observation, rows and nt, number of observation time points, columns
#' @param dllik a function to evaluate the logarithm of the likelihood with arguments y, x, and theta; y is no-length vector giving the current data point, x is an ns (dimension of state) by k < nt+1 matrix giving the history of the state up to the current time point, and theta is an np (dimension of parameter) length vector giving the fixed parameters
#' @param revo a function to propagate the state with arguments x, the state history matrix, and theta, the fixed parameters
#' @param rprior a function to sample from the prior for the state and fixed parameters; returns a list with elements x and theta; takes an integer argument corresponding to the particle number to give the user the option to load already sampled prior draws
#' @param rmove a function to regenerate values of the state history and/or fixed parameters with arguments y, x, and theta; y is an no by k < nt data matrix, x is the state history matrix, and theta is the fixed parameters; returns a list with components state and theta, the regenerated values of the state history and fixed parameters
#' @param n the number of particles
#' @param progress a boolean to display progress bar if TRUE
#' @param ... arguments passed on to resample
#' @return a list containing (nt+1)-length list of state histories (each an ns x n x k matrix), an n x (nt+1) matrix of normalized particle weights, an np x n x (nt+1) array of theta draws, and an n x nt parent matrix
#' @references Berzuini, C. and Gilks, W. Following a Moving Target-Monte Carlo Inference for Dynamic Bayesian Models. Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 63, No. 1 (2001), pp. 127-146
#' @seealso \code{\link{resample}}
#'
rm_pf = function(y, dllik, revo, rprior, rmove, n, progress = TRUE, ...)
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
  state = list(); length(state) = nt + 1
  for(i in 1:(nt + 1)) state[[i]] = array(NA, dim=c(ns,n,i))
  theta = array(NA, dim=c(np,n,nt+1))
  for (j in 1:n) 
  {
    tmp = rprior(j)
    state[[1]][,j,] = tmp$x
    theta[,j,1] = tmp$theta
  }

  # Initialize weights
  weight = matrix(NA, n, nt+1)
  weight[,1] = 1/n
  
  # Initialize parent
  parent = matrix(NA, n, nt+1)

  # Run particle filter
  if(progress) pb = txtProgressBar(0,nt,style=3)
  for (i in 1:nt) 
  {
    if(progress) setTxtProgressBar(pb,i)

    # Augmentation and update weights
    state[[i+1]][,,1:i] = state[[i]]
    for(j in 1:n)
    {
      state[[i+1]][,j,i+1] = revo(state[[i]][,j,],theta[,j,i])
      weight[j,i+1] = log(weight[j,i]) + dllik(y[,i],state[[i+1]][,j,],theta[,j,i])
    }
    
    # Resample particles
    weight[,i+1] = renormalize(weight[,i+1],log=TRUE)
    tmp = resample(weight[,i+1],...)
    kk = tmp$indices
    did.resample = !(all(kk == 1:n))

    # Move particles
    if(did.resample)
    {
      weight[,i+1] = tmp$weight
      for (j in 1:n) 
      {
        tmp2 = rmove(y[,1:i],state[[i+1]][,kk[j],],theta[,kk[j],i])
        state[[i+1]][,j,] = tmp2$state
        theta[,j,i+1] = tmp2$theta
      }
    } else {
      theta[,,i+1] = theta[,,i]
    }

    parent[,i+1] = kk
  }

  return(list(state=state, theta=theta, weight=weight, parent=parent))
}
