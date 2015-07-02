# Next things to do:
# - figure out why the two solutions are so different 

# Other things to do
# - need to read in the ascii file to get the correct bits and bobs for the spatial plotting

# Clear all pbjects and set the error behavior
# setwd("/Users/sriley/Dropbox/svneclipse/globalflu/r-scripts")

rm(list=ls(all=TRUE))
options(error=NULL)

require("sp")

#source("../../R Source/SR_Misc.r")
source("funcs_global.R")

# Some globals
dtm <- 0.25
noreals <- 5
ylims <- c(1,1e14)

# Set up some comparative time series for the mass action version
mass <- gmm.stochSIRS(	nrs=noreals,dt=0.25,dtm=dtm,tmax=180,
						N=1257860,R0=1.6,D_I=2.6,D_R=1440,beta_s=0.00000,I0=100000)
dim_mass <- dim(mass)
notpmass <- dim_mass[2]
norealsmass <- dim_mass[1]
massincy <- vector(mode="numeric",length=(notpmass))
massincx <- (0:(notpmass-1))*dtm
plot(1:2,type="n",xlim=c(0,notpmass*dtm),ylim=c(0,50000))
for (i in 1:noreals) {
	massincy[1] <- 0
	for (j in 2:(notpmass)) massincy[j] = mass[i,j]-mass[i,j-1]
	points(massincx,(massincy+1)/dtm,type="l",col="red")
}

# Load in and compare the output from the spatial model
targdir <- "../working/"
filestem <- "testout"
x <- readSpacialInc(targdir,filestem)
# save(x,file="tmp_rtn_array.Rdata")

for (i in 1:noreals) {
	inc <- genSingInc(i,x,0.25)
	points(inc$x,(inc$y+1)/dtm,type="l",col="blue")
}
