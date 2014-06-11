# Function to approximate the log marginal likelihood using particle filtering given a list returned by rm_pf()
pf.lmarglik <- function(out, unnorm=F)
{
  nt = dim(out$weight)[2]-1
  if(unnorm) log.marglik = sum(log(apply(exp(out$n.weights),2,mean)*apply(exp(out$p.weights)*out$weight[,1:nt],2,sum))) else log.marglik = sum(log(apply(exp(out$increment)*out$weight[1:nt],2,sum)))
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
  pb = txtProgressBar(0,tt*np,style=3)
  for(j in 1:tt)
  {
    for(i in 1:np)
    {
      setTxtProgressBar(pb,(j-1)*np+i)
      weights[i,j] = out$weight[i,j] + dlprior.new(out$state[[j]][1,i,1],out$theta[,i,1]) - dlprior.old(out$state[[j]][1,i,1],out$theta[,i,1])
    }
  }
  weights = apply(weights, 2, function(x) renormalize(x,log=TRUE))
  return(pf.lmarglik(list(weight = weights, increment = out$increment)))
}