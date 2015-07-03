rm(list=ls())
options(error=NULL) # NULL for normal behaviour, recvoer for debugging
source("SR_misc.r")
source("spatial/spatial_int_funcs.r.R")
require("mapdata")
require("rgdal")

# Set flags and key parameters to generate objects (and save them) or load them
flagRegenerate <- FALSE
flagTestLikelihood <- FALSE

# To dos and notes
# Switch to writing single sheet pdfs for all charts
# May be that the approximation for within square moves doesn't work
# so well for the larger scale

# names(fsData)
# [1] "household.id" "member.id"    "live.lat"     "live.long"    "work.status" 
# [6] "work.lat"     "work.long" 

# Files requried
strLandscanFile <- "../../data/landscan/asia03/asia03/w001001.adf"
fsData <- read.csv("../../projects/influenza/GZ_Travel_Study/data/export.fs.csv")
gzData <- read.csv("../../projects/influenza/GZ_Travel_Study/data/export_SR_cleaned_081223.gz.csv")
hkData <- read.csv("../../projects/influenza/GZ_Travel_Study/data/export.hk.csv")

# Files written
fnGZLandScanAscii <- "../spatialsim/data/GZLandScanAscii.txt"
fnFSLandScanAscii <- "../spatialsim/data/FSLandScanAscii.txt"
fnHKLandScanAscii <- "../spatialsim/data/HKLandScanAscii.txt"
fnWideData <- "../spatialsim/data/WideLandScanAscii.txt"
fnPopdata = "../spatialsim/data/prd_sim_ascii_n_44384655.txt"
fnToydata = "../spatialsim/data/toy_ascii_grid.in"
fnAuxTmp <- "spatial/aux_hk_fs_gz_wide.Rdata"

# Some housekeeping for later analyses
listData <- list(fsData,gzData,hkData)
noData <- length(listData)
nameData <- c("Foshan","Guangzhou","Hong Kong")
noData <- length(listData)
noHh <- c(max(fsData$household.id),max(gzData$household.id),max(hkData$household.id))

fsData <- cbind(city=rep("fs",dim(fsData)[1]),fsData)
gzData <- cbind(city=rep("gz",dim(gzData)[1]),gzData)
hkData <- cbind(city=rep("hk",dim(hkData)[1]),hkData)

nofsData <- dim(fsData)[1] 
nogzData <- dim(gzData)[1] 
nohkData <- dim(hkData)[1] 

hkData$dist <- decLongLatDist(hkData$live.long,hkData$live.lat,hkData$work.long,hkData$work.lat,translate=TRUE)
gzData$dist <- decLongLatDist(gzData$live.long,gzData$live.lat,gzData$work.long,gzData$work.lat,translate=TRUE)
fsData$dist <- decLongLatDist(fsData$live.long,fsData$live.lat,fsData$work.long,fsData$work.lat,translate=TRUE)

combData <- rbind(fsData,gzData,hkData)
combData <- combData[(1:length(combData$work.status))[!is.na(combData$work.status)],]
tabWorkStatus <- table(combData$work.status,combData$city)
tabNoInds <- table(combData$city)
tabAveHhSize <- tabNoInds / noHh

noInd <- c(dim(fsData)[1],dim(gzData)[1],dim(hkData)[1])
aveHhSize <- noInd*1.0 / noHh

# Get max and min values for home data
hbounds = data.frame(minx=rep(0,noData),maxx=rep(0,noData),miny=rep(0,noData),maxy=rep(0,noData),row.names=nameData)
for (i in 1:noData) {
	tmpdat <- listData[[i]]
	hbounds$minx[i] <- min(tmpdat$live.long,na.rm=TRUE)
	hbounds$maxx[i] <- max(tmpdat$live.long,na.rm=TRUE)
	hbounds$miny[i] <- min(tmpdat$live.lat,na.rm=TRUE)
	hbounds$maxy[i] <- max(tmpdat$live.lat,na.rm=TRUE)
}

# Get max and min values for reported workplaces
wbounds = data.frame(minx=rep(0,noData),maxx=rep(0,noData),miny=rep(0,noData),maxy=rep(0,noData),row.names=nameData)
for (i in 1:noData) {
	tmpdat <- listData[[i]]
	wbounds$minx[i] <- min(tmpdat$work.long,na.rm=TRUE)
	wbounds$maxx[i] <- max(tmpdat$work.long,na.rm=TRUE)
	wbounds$miny[i] <- min(tmpdat$work.lat,na.rm=TRUE)
	wbounds$maxy[i] <- max(tmpdat$work.lat,na.rm=TRUE)
}

wide_xlim=c(floor(min(wbounds$minx)),ceiling(max(wbounds$maxx)))
wide_ylim=c(floor(min(wbounds$miny)),ceiling(mq()ax(wbounds$maxy)))

# prepare objects which require some computation time 
if (flagRegenerate) {
	ExcerptLStoASCII(strLandscanFile,fnWideData,wide_xlim[1],wide_xlim[2],wide_ylim[1],wide_ylim[2])
	ExcerptLStoASCII(strLandscanFile,fnGZLandScanAscii,hbounds["Guangzhou","minx"],hbounds["Guangzhou","maxx"],hbounds["Guangzhou","miny"],hbounds["Guangzhou","maxy"])
	ExcerptLStoASCII(strLandscanFile,fnFSLandScanAscii,hbounds["Foshan","minx"],hbounds["Foshan","maxx"],hbounds["Foshan","miny"],hbounds["Foshan","maxy"])
	ExcerptLStoASCII(strLandscanFile,fnHKLandScanAscii,hbounds["Hong Kong","minx"],hbounds["Hong Kong","maxx"],hbounds["Hong Kong","miny"],hbounds["Hong Kong","maxy"])
	popwide <- read.asciigrid(fnWideData,as.image=TRUE)
	popgz <- read.asciigrid(fnGZLandScanAscii,as.image=TRUE)
	popfs <- read.asciigrid(fnFSLandScanAscii,as.image=TRUE)
	pophk <- read.asciigrid(fnHKLandScanAscii,as.image=TRUE)
	popwide <- downSampleAscii(popwide,10)
	auxWideHk <- genDistMatrix(popwide,hkData$live.long,hkData$live.lat,loud=TRUE)
	auxWideGz <- genDistMatrix(popwide,gzData$live.long,gzData$live.lat,loud=TRUE)
	auxWideFs <- genDistMatrix(popwide,fsData$live.long,fsData$live.lat,loud=TRUE)
	auxNarrowHk <- genDistMatrix(pophk,hkData$live.long,hkData$live.lat,loud=TRUE)
	auxNarrowGz <- genDistMatrix(popgz,gzData$live.long,gzData$live.lat,loud=TRUE)
	auxNarrowFs <- genDistMatrix(popfs,fsData$live.long,fsData$live.lat,loud=TRUE)
	save(popwide,popgz,popfs,pophk,auxWideHk,auxWideGz,auxWideFs,auxNarrowHk,auxNarrowGz,auxNarrowFs,file=fnAuxTmp)
} else {
	load(fnAuxTmp)
}

if (flagTestLikelihood) {
	lnlikegz <- lnlike_kern_new(c(10,-2.2),
			gzData$live.long,gzData$live.lat,
			gzData$work.long,gzData$work.lat,
			popgz,auxNarrowGz,popwide,auxWideGz,
			drclose=1,drwide=10)	

	lnlikehk <- lnlike_kern_new(c(10,-2.2),
			hkData$live.long,hkData$live.lat,
			hkData$work.long,hkData$work.lat,
			pophk,auxNarrowHk,popwide,auxWideHk,
			drclose=1,drwide=10)	

}

opt_out_gz <- 	optim(	c(10,-1),
					lnlike_kern_new,
					method="L-BFGS-B",
					lower=c(0.001,-10,1),
					upper=c(100,-2,100),
					control=list(trace=0,fnscale=-1,maxit=10000,REPORT=1),
					hx=gzData$live.long,hy=gzData$live.lat,
					wx=gzData$work.long,wy=gzData$work.lat,
					smallpopgrid=popgz,smallpopdist=auxNarrowGz,
					bigpopgrid=popwide,bigpopdist=auxWideGz
			)

# 2.951422 -3.542913			
			
opt_out_hk <- 	optim(	c(10,-1),
					lnlike_kern_new,
					method="L-BFGS-B",
					lower=c(0.001,-10,1),
					upper=c(100,-2,100),
					control=list(trace=0,fnscale=-1,maxit=10000,REPORT=1),
					hx=hkData$live.long,hy=hkData$live.lat,
					wx=hkData$work.long,wy=hkData$work.lat,
					smallpopgrid=pophk,smallpopdist=auxNarrowHk,
					bigpopgrid=popwide,bigpopdist=auxWideHk
			)

# 0.001	-2.000
		
# Plot distribution of work locations (Figure 1)
windows(width=7.5/cm(1),height=15/cm(1),rescale="R")
par(mai=c(0.3,0.3,0.01,0.01))
par(cex.axis=0.75)
par(mgp=c(3,0.5,0))
par(fig=sr_chart_pos(1,2,1,2))
plot(	c(0,1),xlim=wide_xlim,ylim=wide_ylim, 
		xlab="", ylab="",type="n")
mtext("A",side=2,at=wide_ylim[2],line=-1,las=1)
map('china', xlim=c(100, 122),ylim=c(20,42),add=TRUE)
points(fsData$work.long,fsData$work.lat,col="red")
points(gzData$work.long,gzData$work.lat,col="blue")
points(hkData$work.long,hkData$work.lat,col="green")

par(fig=sr_chart_pos(1,1,1,2),new=TRUE)
plot(1:2,xlim=c(112,116),ylim=c(21,25), xlab="Longitude", ylab="Latitude",
		type="n")
mtext("B",side=2,at=25,line=-1,las=1)
map('china', xlim=c(112, 116),ylim=c(21,25),add=TRUE)
# points(fsData$work.long,fsData$work.lat,col="red")
# points(gzData$work.long,gzData$work.lat,col="blue")
# points(hkData$work.long,hkData$work.lat,col="green")
points(gzData$live.long,gzData$live.lat,col="green")

# Plot fitted distance kernels
windows(width=7.5/cm(1),height=7.5/cm(1),rescale="R")
par(mai=c(0.5,0.5,0.01,0.01))
par(cex.axis=0.75)
par(mgp=c(1.25,0.5,0))
y1 <- lapply(x,fnKernHaz,d=5.713,a=-4.9)
y2 <- lapply(x,fnKernHaz,d=2.12,a=-2.39)
plot(x,y1,log="xy",type="l",col="red",xlab = "Distance", ylab = "Probability",lwd=2)
points(x,y2,col="blue",type="l",lwd=2)
		