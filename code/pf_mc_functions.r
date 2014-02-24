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