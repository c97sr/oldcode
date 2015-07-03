# THIS FILE NOW SUPERCEEDED BY abp... FUNCTIONS IN MAIN FUNCTION FILE
# THIS FILE LEFT "JUST IN CASE"

# Script to make a figure that illustrates multiannual properties of flu
rm(list=ls(all=TRUE))
options(error=NULL)
require(odesolve)

# Edit below as required. 
# source("http://idsource.googlecode.com/svn/trunk/R/stevensRfunctions.R")
source("stevensRfunctions.R")

# Next: put in a few alternate parameters from Jeff and Mark's paper

# Plot the comparative chart
abp.fig.illustrate()

# Edit below to avoid recalculating the vector, needs to be wrapped into a function
# This needs to be a genuine latin hyper cube to see if there are any real results
# I'm sure I've made a hypercube before in R, but I'm not sure what I did with it

dfParams = data.frame(params=c("gamma_R","dur_seas"),log=c(0,0),lb=c(1/20,7),ub=c(1/1,7*26))
system.time(tmp <- abp.hyper(nohypersamples=5,dfParams))
													
matmin <- tmp$matmin
hypersquare <- tmp$hypersquare
nosamples <- dim(matmin)[1]

# Think that below here is old or broken

plot(tabparam[1],)

pdf(paste("~/Dropbox/tmp/AB_persist_fig",".pdf",sep=""))
plot(rep(vecAveDurRec,mintakenof),matEachAnnMin,ylim=c(1,10000),log="y",main=paste("Amp Seasonality",vecSeasonAmp[j]))
dev.off()
write.csv(matEachAnnMin,file=paste("~/Dropbox/tmp/AB_persist_out",".csv",sep=""))
