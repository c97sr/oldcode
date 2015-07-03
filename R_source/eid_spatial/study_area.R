rm(list=ls(all=TRUE))
require("rgdal")
require("RArcInfo")

strLandscanFile <- "D:\\files\\data_ns\\landscan\\asia03\\asia03\\w001001.adf"

strShapeDir <- "D:\\files\\data_ns\\worldsf"
strShapeLayer <- "world_adm0"

# world <- readOGR(strShapeDir,layer=strShapeLayer)

# detail easting, detail westing..., wide eastin, wide westing...
de <- 114+5/60
dw <- 112+57/60
dn <- 23+56/60
ds <- 22+49/60

we <- 123 +  4/60
ww <- 108 +  7/60
wn <- 25  + 17/60
ws <- 18  +  3/60

spd <- 	120 	# Squares per degree
mn  <- 	85 		# map north
mw  <- 	34 		# map west

info <- GDALinfo(strLandscanFile)
wide <- readGDAL(strLandscanFile,region.dim=c((wn-ws)*spd,(we-ww)*spd),offset=c((mn-wn)*spd,(ww-mw)*spd))
detail <- readGDAL(strLandscanFile,region.dim=c((dn-ds)*spd,(de-dw)*spd),offset=c((mn-dn)*spd,(dw-mw)*spd))

f3w <- 10/cm(1)
widew <- 3/cm(1)	# Detail width
det_border <- 0.5/cm(1)
f3h <- (f3w-det_border)*(dn-ds)/(de-dw)+det_border
posdet <- c(0.01,0.99-det_border/f3w,0.01,0.99-det_border/f3h) 
poswide <- c(0.99-widew/f3w,0.99,0.99-widew*(wn-ws)/(we-ww)/f3h,0.99)

windows(width=f3w,height=f3h)
par(	mai=c(0,0,0,0), 		# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,
		mgp=c(0,0,0),
		par(cex=1))

par(fig=posdet)
image(detail,col=grey(0.9-0.8*0:100/100))
# image(detail,col="grey")
for (i in 1:4) axis(i,at=c(0,130),tck=0,col="blue",lwd=2)
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
polygon(c(0,1,1,0),c(0,0,1,1),col="white",border=NA)
par(new=TRUE,fig=poswide)
#image(wide,col="grey")
image(wide,col=grey(0.9-0.8*0:100/100))
for (i in 1:4) axis(i,at=c(0,130),tck=0)
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
lines(c(dw-ww,de-ww,de-ww,dw-ww,dw-ww)/(we-ww),c(ds-ws,ds-ws,dn-ws,dn-ws,ds-ws)/(wn-ws),col="blue",lwd=2)

writeGDAL(detail,"c:\tmp\test",drivername="AAIGrid")


