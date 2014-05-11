# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load MLEs
load(paste(dpath,"dlm_arwn.rdata",sep=""))

# True values?
theta.truth = c(900, 5, .8, 1, 1)

# Plot histograms of MLEs
hist_arwn <- function(nt, x)
{
  ind1 = which(dimnames(mle.all)[[1]] == as.character(nt))
  ind2 = which(dimnames(mle.all)[[2]] == as.character(x))
  pdf(file = paste(gpath,"dlm_arwn-",nt,"-",x,".pdf",sep=""), width = 15, height = 10)
  par(mfrow=c(2,3), mar = c(5,6,4,2) +0.1)
  xlab = expression(beta[0], beta[1], phi, sigma[s]^2, sigma[m]^2)
  ylab = c("Frequency",rep("",4))
  for(i in 1:5)
  {
    hist(mle.all[ind1,ind2,,i], freq = TRUE, xlab = xlab[i], ylab = ylab[i], main="", cex.lab=1.5)
    abline(v = theta.truth[i], col = 2, lwd = 3)
  }
  dev.off()
}

require(plyr)
mydata = expand.grid(nt = c(500, 1000, 5000, 10000), x = c("conv","norm","none"), stringsAsFactors = FALSE)
m_ply(mydata, hist_arwn)
