# Function to construct plots for particle filters
pf_plots <- function(np, alpha, beta, xlab, ylab)
{
  source("rm_test_functions.r")
  load("../data/rm_test-truth.rdata")

  pdf("../graphs/rm_test-states.pdf", width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:1)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(x)
    } else {
      xlab = ""
      ylab = ""
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.states_plot(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab)
  }
  dev.off()

  pdf("../graphs/rm_test-states-zeroed.pdf", width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:1)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(x)
    } else {
      xlab = ""
      ylab = ""
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.states.zeroed_plot(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab)
  }
  dev.off()

  pdf("../graphs/rm_test-pvalues.pdf", width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:1)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = "P-value"
    } else {
      xlab = ""
      ylab = ""
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.pvalues_plot(m, C, a, b, out.kd, out.rm1, out.rm2, xlab, ylab)
  }
  dev.off()

  pdf("../graphs/rm_test-precision.pdf", width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:1)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(phi)
    } else {
      xlab = ""
      ylab = ""
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.precision_plot(m, C, a, b, v, out.kd, out.rm1, out.rm2, xlab, ylab)
  }
  dev.off()
}

# Construct plots
mydata = expand.grid(np=c(100,1000), alpha=1, beta=c(1,.25), stringsAsFactors=FALSE)
require(plyr)
m_ply(.data = mydata[1,], .fun = pf_plots)