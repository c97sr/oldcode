# assumes that severity model has been run

# Set up the figure
xaxeA <- c(as.date("1May2009"),as.date("31Jan2010"))
vecColors <- c("red","green","blue","magenta")
ymaxA <- 800

# Setup a multipart figure to show the study data and the HK hosptialization data
# Prepare a different version of this for what _could_ have been done
# Keep these in landscape for slide presentations

posA <- c(0.1,0.95,0.55,0.95)
posB <- c(0.1,0.95,0.125,0.525)

head(base)
noparts <- dim(base)[1]

# Set up the pdf file
pdf("../figs/model_severity.pdf",height=9/cm(1),width=15/cm(1))

# Standard margin zeroing and setting axis label formats
par(mai=(c(0,0,0,0)))
par(mgp=c(2.5,1,0))

# Make the initial plot
par(fig=posA)
plot(1:2,type="n",axes=TRUE,xlim=xaxeA,ylim=c(0,ymaxA),ylab="")
plotDiscInc(matHOS,wkBins,vecColors)

par(fig=posB,new=TRUE)
plot(1:2,type="n",axes=FALSE,xlim=xaxeA,ylim=c(0,ymaxA),ylab="")
plotDiscInc(matHOS,wkBins,vecColors)

# Close the pdf file
dev.off()
