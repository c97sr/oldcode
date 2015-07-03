# Make sensitivity figure
# Needs script probably called "ventilation.R" to ahve been run prior to this script being run

cwidth 	<- 13

g1 <- 2
g2 <- 2
g3 <- 0.125
g4 <- 0.25

B  <- 100

# Spaces for the plots
A <- (cwidth - g1 - 2* g3) * B / (1 + B)

cheight <- g2 + A + g4
pdf("./figs/Sensetivity.pdf", width=cwidth/cm(1), height=cheight/cm(1))

# Positions of the different parts
posA <- c(g1/cwidth,(g1+A)/cwidth,(g2)/cheight,(g2+A)/cheight)
posB <- c((g1+A+g3)/cwidth,(g1+A+g3+A/B)/cwidth,g2/cheight,(g2+2 * A/3)/cheight)

# Parameter settings for all parts of the figure
par(cex=0.8)
par(mgp=c(2.5,1,0))
# par(tcl=0.25)
par(mai=(c(0,0,0,0)))
par(cex.axis=0.9)

yax_loc <- log(c(0.005,0.01,0.02,0.05),base=10)
yax_lab <- c("0.5%","1.0%","2.0%","5.0%")
xax_loc <- c(0.05,0.1,0.15,0.2,0.25,0.3)
xax_lab <- c("0.05 (13.8)","0.10 (6.9)","0.15 (4.6)","0.20 (3.5)","0.25 (2.8)","0.3 (2.3)")
# Part A 
par(fig=posA)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=x_lims,
		ylim=y_lims)

axis(2,las=1,at=yax_loc,labels=yax_lab)
axis(1,las=1,at=xax_loc,labels=xax_lab)

mtext(substitute(paste("Initial hospitalization rate per infection ","",sep="")),side=2,line=3.5)
mtext(substitute(paste("Growth rate, ",days^{-1}," (doubling time, days)",sep="")),side=1,line=2.5)
contlev <- c(1,2,5,10,20,50,100)
image(x_mids,y_mids,chtData,add=TRUE)
contour(	x_mids,y_mids,chtData,
		add=TRUE,
		levels=contlev,
		labels=as.character(contlev)
)

points(c(0.26),c(log(0.0037,base=10)))

dev.off()
