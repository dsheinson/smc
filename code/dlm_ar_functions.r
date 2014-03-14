# Functions to help construct parameters of a DLM formulation of regression with AR(p) errors and measurement noise

makeG <- function(phi) {
  m = length(phi)
  G <- diag(0,m)
  G[,1] <- phi
  G[-m,-1] <- diag(1,m-1)
  return(G)
}

is.stationary <- function(phi)
{
  return(all(Mod(polyroot(c(1,-phi)))>1))
}

makeC0 <- function(phi)
{
  p <- length(phi)
  G <- makeG(phi)
  f <- c(1,rep(0,p-1))
  vv <- solve(diag(p^2)-G%x%G)%*%as.numeric(f%*%t(f))
  nr <- sqrt(length(vv))
  return(matrix(vv,nr,nr))
}

bdiag <- function(mylist)
{
  if(length(mylist) > 0)
  {
    ps = sapply(mylist,dim)
    nrow = sum(ps[1,])
    ncol = sum(ps[2,])
    C0 = matrix(0,nr=nrow,nc=ncol)
    row.index = 1; col.index = 1
    for(i in 1:length(mylist))
    {
      if(all(ps[,i] != 0))
      {
        C0[row.index:(row.index + ps[1,i] - 1),col.index:(col.index + ps[2,i] - 1)] = mylist[[i]]
        row.index = row.index + ps[1,i]
        col.index = col.index + ps[2,i]
      }
    }
    return(C0)
  } else { return(matrix(nr=0,nc=0)) }
}