# Function to approximate the log marginal likelihood using particle filtering given a list returned by rm_pf()
pf.lmarglik <- function(out)
{
  nt = dim(out$increment)[2]
  log.marglik <- 0
  for(i in 1:nt) log.marglik = log.marglik + log(sum(exp(out$increment[,i]+log(out$weight[,i]))))
  return(log.marglik)
}

# Function to calculate posterior model probabilities given a set of log marginal likelihoods and prior model probabilities
postModProbs <- function(lmarglik, priorModProbs)
{
  stopifnot(length(lmarglik) == length(priorModProbs))
  postModProbs.log = lmarglik + log(priorModProbs) - log(sum(exp(lmarglik)*priorModProbs))
  return(exp(postModProbs.log))
}

# Function to calculate log marginal likelihood of particle sample given a different prior
pf.lmarglik.prior <- function(out, dlprior.old, dlprior.new)
{
  tt = dim(out$weight)[2]
  nt = tt - 1
  np = dim(out$weight)[1]
  weights = matrix(NA, nr = np, nc = tt)
  for(j in 1:tt)  for(i in 1:np) weights[i,j] = out$weight[i,j] + dlprior.new(out$state[[j]][1,i,1],out$theta[,i,1]) - dlprior.old(out$state[[j]][1,i,1],out$theta[,i,1])
  weights = apply(weights, 2, function(x) renormalize(x,log=TRUE))
  return(pf.lmarglik(list(weight = weights, increment = out$increment)))
}