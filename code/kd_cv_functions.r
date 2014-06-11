# Functions to run kernel density particle filter for local level DLM
# theta is on the log scale

dllik.kd <- function(y, x, theta, lambda = 1)
{
  dnorm(y, x, sqrt(exp(theta)), log=T)
}

revo.kd <- function(x, theta, lambda = 1)
{
  rnorm(1, x, sqrt(exp(theta)*lambda))
}

pstate.kd <- function(x, theta, lambda = 1)
{
  x
}

rprior.kd <- function(a0=1, b0=1)
{
  theta = 1 / rgamma(1, a0, b0)
  x0 = rnorm(1, 0, sqrt(theta))
  return(list(x=x0,theta=log(theta)))
}