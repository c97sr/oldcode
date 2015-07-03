# Assumes that runs_prep_data_MCMC.R has been run 
# therefore doesn't contain a rm(list=ls(all=TRUE))  - watch out!
# Next in this script:
#	- Tidy up axes
#	- Add axes labels (only one for x and one for y)
#	- Add legend

# Load up the parameters and make the 
source ("./runs_prep_data_MCMC.R")
source("./template/baseline_HK_Hosp_Params.R")
# source("../funcs.R")

# Use this for max like from an MCMC run
mlenames <- c("r","t0","p_H_base","phi_1","phi_2","phi_4","delta_sm")
optres <- read.csv("~/Dropbox/results/100511/trace.csv")
mle_loc <- which.max(optres$lnlike)
mlevals <- as.vector(as.numeric(optres[mle_loc,mlenames])) 

# mlevals <- c(3.748516e-02,1.776780e+04,1.335953e-01,5.660341e+00,6.715860e+00,2.968270e-01,9.657545e-01)
mlevals <- c(3e-02,1.776780e+04,1.335953e-01,5.660341e+00,6.715860e+00,2.968270e-01,1.4)

tmp <- SirModelFourAgeClasses(	
		pname=c(mlenames),
		pvals=c(mlevals),
		casemix=c(-5,-40,-40,-40),
		vp=baseHK4 )	

solution 	<- as.data.frame(tmp$sol)
endparams	<- as.data.frame(tmp$par)

# plot(solution$time,solution$dIv,type="l",col="black")
agmodel <- extFourAgeICU(solution,wkBins,agcats)

# Set up the multipart chart
# Total widths
xw <- 15
yw <- 10

# Set up y axes limits
ymaxA <- 800
ymaxB <- 200
ymaxC <- 400
ymaxD <- 200
ymaxE <- 50
xaxeA <- c(as.date("1May2009"),as.date("31Jan2010"))
xaxeBCDE <- c(as.date("7Jul2009"),as.date("01Nov2009"))

# Ratios of plottable areas
yr <- c(ymaxB,ymaxC,ymaxD,ymaxE)
xr 	<- c(xaxeA[2]-xaxeA[1],xaxeBCDE[2]-xaxeBCDE[1])

# Gaps in cms divded by total widths
xg1 <- 2.0 / xw
xg2 <- 1.6 / xw
xg3 <- 1.0 / xw
yg1 <- 1.5 / yw
yg2 <- 1.0 / yw
yg3 <- 0.5 / yw

# Sizes of plottable areas
ysA <- (1-yg2-yg1) 
ysB <- yr / sum(yr) * (1-yg3*3-yg1-yg2) 
xs	<- xr  / sum(xr) * (1-xg1-xg2-xg3)

# Define plottable areas
posA <- c(xg1,xg1+xs[1],yg1,yg1+ysA)
posB <- c(xg1+xg2+xs[1],1-xg3,yg1+ysB[4]+3*yg3+ysB[3]+ysB[2],yg1+ysB[4]+3*yg3+ysB[3]+ysB[2]+ysB[1])
posC <- c(xg1+xg2+xs[1],1-xg3,yg1+ysB[4]+2*yg3+ysB[3],yg1+ysB[4]+2*yg3+ysB[3]+ysB[2])
posD <- c(xg1+xg2+xs[1],1-xg3,yg1+ysB[4]+yg3,yg1+ysB[4]+yg3+ysB[3])
posE <- c(xg1+xg2+xs[1],1-xg3,yg1,yg1+ysB[4])

# Set up the pdf file
pdf("../figs/fit_HK_hosp.pdf",height=yw/cm(1),width=xw/cm(1))

# Standard margin zeroing and setting axis label formats
par(mai=(c(0,0,0,0)))
par(mgp=c(2.5,1,0))
vecColors <- c("red","green","blue","magenta")

# Part A - the complete data set
par(fig=posA)
plot(1:2,type="n",axes=TRUE,xlim=xaxeA,ylim=c(0,ymaxA),ylab="")
plotDiscInc(agdata,wkBins,vecColors)

# Part B - young children
par(fig=posB,new=TRUE)
plot(1:2,type="n",axes=TRUE,xlim=xaxeBCDE,ylim=c(0,ymaxB),ylab="")
plotHalfHalf(agmodel[agcats[1],],agdata[agcats[1],],wkBins,color=vecColors[1])

# Part C - school age children
par(fig=posC,new=TRUE)
plot(1:2,type="n",axes=TRUE,xlim=xaxeBCDE,ylim=c(0,ymaxC),ylab="")
plotHalfHalf(agmodel[agcats[2],],agdata[agcats[2],],wkBins,color=vecColors[2])

# Part D - young adults
par(fig=posD,new=TRUE)
plot(1:2,type="n",axes=TRUE,xlim=xaxeBCDE,ylim=c(0,ymaxD),ylab="")
plotHalfHalf(agmodel[agcats[3],],agdata[agcats[3],],wkBins,color=vecColors[3])

# Part E - older adults
par(fig=posE,new=TRUE)
plot(1:2,type="n",axes=TRUE,xlim=xaxeBCDE,ylim=c(0,ymaxE),ylab="")
plotHalfHalf(agmodel[agcats[4],],agdata[agcats[4],],wkBins,color=vecColors[4])

dev.off()
