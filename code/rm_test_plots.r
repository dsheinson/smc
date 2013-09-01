# Function to construct plots for particle filters
pf_plots <- function(np, alpha, beta)
{
  source("rm_test_functions.r")
  load("../data/rm_test-truth.rdata")
  pdf("../graphs/rm_test-states.pdf")
  par(mfrow=c(3,3))
  for(label in 1:1)
  {
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.states_plot(m, C, a, b, x, out.kd, out.rm1, out.rm2)
  }
  dev.off()
}

# Construct plots
mydata = expand.grid(np=c(100,1000), alpha=1, beta=c(1,.25), stringsAsFactors=FALSE)
require(plyr)
m_ply(.data = mydata[1,], .fun = pf_plots)
