# This script assumes that "runs_ilustrate_peak_ar.R" has been run
axticks <- as.date(c(
				"1April2009","1May2009","1Jun2009","1Jul2009","1Aug2009"))

axstart <- min(axticks)
axend	<- max(axticks)
axnames	<- c("0","1","2","3","4")
# ytickvA	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175)
# yticknA	<- c("0.0%","","0.5%","","1.0%","","1.5%","")
# ymaxA 	<- max(ytickvA)
# ytickvB	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)
# yticknB	<- c("0.0%","","0.5%","","1.0%","","1.5%")
# ymaxB 	<- max(ytickvB)
# ytickvC	<- c(0,0.0025,0.005,0.0075,0.01,0.0125)
# yticknC	<- c("0.0%","","0.5%","","1.0%","")
# ymaxC 	<- max(ytickvC)

ytickvA	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175)
yticknA	<- c("0.0%","","0.5%","","1.0%","","1.5%","")
ymaxA 	<- max(ytickvA)
ytickvB	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)
yticknB	<- c("0.0%","","0.5%","","1.0%","","1.5%")
ymaxB 	<- max(ytickvB)
ytickvC	<- c(0,0.0025,0.005,0.0075,0.01,0.0125)
yticknC	<- c("0.0%","","0.5%","","1.0%","")
ymaxC 	<- max(ytickvC)

# Total widths
xw <- 15
yw <- 13

# Ratios of plottable areas
yra <- c(ymaxC,ymaxB,ymaxA)
yrb <- c(1,1)
xr 	<- c(5,5,0.5)

# Gaps in cms divded by total widths
yg1 <- 1.5 / yw
yg2 <- 0.5 / yw
yg2b <- 1.5 / yw
xg1 <- 2.0 / xw
xg2 <- 1.6 / xw
xg3 <- 1.0 / xw
xg4 <- 0.25 / xw

# Sizes of plottable areas
ysa <- yra / sum(yra) * (1-length(yra)*yg2-yg1) 
ysb <- yrb / sum(yrb) * (1-yg2b-yg1-yg2) 
xs	<- xr  / sum(xr) * (1-xg2-xg3-xg1-xg4)

# Define plottable areas
posA <- c(xg1,xg1+xs[1],yg1+ysa[1]+yg2+ysa[2]+yg2,yg1+ysa[1]+yg2+ysa[2]+yg2+ysa[3])
posB <- c(xg1,xg1+xs[1],yg1+ysa[1]+yg2,yg1+ysa[1]+yg2+ysa[2])
posC <- c(xg1,xg1+xs[1],yg1,yg1+ysa[1])
posD <- c(xg1+xs[1]+xg2,xg1+xs[1]+xg2+xs[2],yg1+ysb[1]+yg2b,yg1+ysb[1]+yg2b+ysb[2])
posE <- c(xg1+xs[1]+xg2,xg1+xs[1]+xg2+xs[2],yg1,yg1+ysb[1])
posF <- c(xg1+xs[1]+xg2+xs[2]+xg3,xg1+xs[1]+xg2+xs[2]+xg3+xs[3],yg1+ysb[1]+yg2b,yg1+ysb[1]+yg2b+ysb[2])
posG <- c(xg1+xs[1]+xg2+xs[2]+xg3,xg1+xs[1]+xg2+xs[2]+xg3+xs[3],yg1,yg1+ysb[1])

pdf("~/Dropbox/projects/influenza/age_peak_attack/Current/figs/illustrate_hetero.pdf",height=yw/cm(1),width=xw/cm(1),useDingbats=FALSE)

# Settings for all parts of the figure
# Settings for the default spaces around the plot areas
par(mai=(c(0,0,0,0)))
# Not sure
par(mgp=c(2.5,1,0))
# Scaling of axis fonts
nonaxis_cex <- 11/12
axis_cex <- 10/12
par(cex.axis=axis_cex)

# First chart setup
par(fig=posA)
plot(1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxA),ylab="")

# Put the axes on
axis(2,las=1,at=ytickvA,labels=yticknA,hadj=0.8)
text(axstart,max(ytickvA),"A",font=2,adj=c(0,1))

# Plot the lines
sol <- sol_part_A
# To plot daily incidence as a percentage of population 
# plotstack(sol,1,"red",var="dI",mult=1/popsize)
# plotstack(sol,2,"blue",var="dI",mult=1/popsize)
# plotstack(sol,3,"blue",var="dI",mult=1/popsize)
# points(sol[,"time"],sol[,"dS"]/popsize,type="l",col="black",lwd=2,lty=1)

# srg.plotstack(sol,1,"red",var="dI",mult=1/popsize)
srg.plotstack.dash(sol,"red",var="dI1",mult=1/popsize)
# srg.plotstack(sol,2,"blue",var="dI",mult=1/popsize)
srg.plotstack.dash(sol,"blue",var="dI2",varbelow=c("dI1"),mult=1/popsize)
# srg.plotstack(sol,3,"blue",var="dI",mult=1/popsize)
srg.plotstack.dash(sol,"blue",var="dI3",varbelow=c("dI1","dI2"),mult=1/popsize)
points(sol[,"time"],sol[,"dS"]/popsize,type="l",col="black",lwd=2,lty=1)

# Add a legend
legend(as.date("07Jun2009"),max(ytickvA)*3.75/4,
		c(	"Children",
			"Adults"),
		lty=c(1,1),
		lwd=c(8,8),
		col=c("red","blue"),
		cex=axis_cex*0.8, 
		bty="n",
)

# Second chart setup
par(fig=posB,new=TRUE)
plot(1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxB),ylab="")

# Put the axes on
axis(2,las=1,at=ytickvB,labels=yticknB,hadj=0.8)
text(axstart,max(ytickvB),"B",font=2,adj=c(0,1))
mtext("Daily incidence",side=2,line=2.7,cex=nonaxis_cex)

# Plot the lines
sol <- sol_part_B
srg.plotstack(sol,1,"red",var="dI",mult=1/popsize)
srg.plotstack(sol,2,"blue",var="dI",mult=1/popsize)
srg.plotstack(sol,3,"blue",var="dI",mult=1/popsize)
points(sol[,"time"],sol[,"dS"]/popsize,type="l",col="black",lwd=2,lty=1)

# Third chart setup
par(fig=posC,new=TRUE)
plot(	1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxC),ylab="")

# Put the axes on
axis(2,las=1,at=ytickvC,labels=yticknC,hadj=0.8)
axis(1,at=axticks,labels=axnames,line=0,padj=-0.8)
text(axstart,max(ytickvC),"C",font=2,adj=c(0,1))
mtext("Time (months)",side=1,line=1.5,cex=nonaxis_cex)

# Plot the lines
sol <- sol_part_C
srg.plotstack(sol,1,"red",var="dI",mult=1/popsize)
srg.plotstack(sol,2,"blue",var="dI",mult=1/popsize)
srg.plotstack(sol,3,"blue",var="dI",mult=1/popsize)
points(sol[,"time"],sol[,"dS"]/popsize,type="l",col="black",lwd=2,lty=1)

# Plot the heat chart
par(fig=posD,new=TRUE)

zlimits <- c(0.1,0.6)
zlimits2 <- c(0.0,0.5)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=x_lims,
		ylim=y_lims)

axis(2,las=1,hadj=0.3)
axis(1,las=1,padj=-0.8)

cscale1 <- rev(heat.colors(130))

image(x_mids,y_mids,chtData[,,1],add=TRUE,zlim=zlimits,col=cscale1)

cntlevs <- c(0.2,0.25,0.3,0.4,0.5)
contour(x_mids,y_mids,chtData[,,1],
		add=TRUE,
		levels=cntlevs,
		labels=as.character(cntlevs)
)


mtext("Susceptibility",side=2,line=1.5,cex=nonaxis_cex)
mtext("Mixing",side=1,line=1.5,cex=nonaxis_cex)
text(x_lims[1],y_lims[2],"D",font=2,adj=c(0,1))

# Plot the legend the heat chart
legend_resolution <- 50
start <- zlimits[1]
end <- zlimits[2]
width <- end - start
binsize <- width / legend_resolution
zmids <- start + binsize / 2 + binsize * (0:(legend_resolution-1))
legend_data <- matrix(zmids,nrow=legend_resolution,ncol=2)

par(fig=posF,new=TRUE)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,1),
		ylim=zlimits)

image(c(0.25,0.75),zmids,t(legend_data),add=TRUE,zlim=zlimits,col=cscale1)
axis(2,las=1,hadj=0.7,pos=-0.3)

# Make the scatter plot of the dependency on over-representation 
par(fig=posE,new=TRUE)
plot(1:2,type="n",axes=FALSE,xlim=rep_bnds,ylim=zlimits2)
axis(1,las=1,padj=-1.0)
axis(2,las=1,hadj=0.7)
mtext("Cumulative attack rate",side=2,line=2.0,cex=nonaxis_cex)
mtext("Over-representation",side=1,line=1.5,cex=nonaxis_cex)
text(rep_bnds[1],zlimits2[2],"E",font=2,adj=c(0,1))

# Plot the legend to the scatterplot
colres <- 100
cscale2 <- cm.colors(colres)
legend_resolution <- 50
start <- 0
end <- 2
width <- end - start
binsize <- width / legend_resolution
zmids <- start + binsize / 2 + binsize * (0:(legend_resolution-1))
legend_data <- matrix(zmids,nrow=legend_resolution,ncol=2)

vec_col <- as.vector(1:no_samples)
for (i in 1:no_samples) vec_col[i] <- cscale2[round((out_phi[i]-start)/(end-start)*colres)]
for (i in 1:noRs) {
	points(hv_rep,out_pi[,i]*100000,col=vec_col,pch=19)
	points(hv_rep,out_pi[,i]*100000,pch=1)
}

# Plot the legend for the scatter plot
par(fig=posG,new=TRUE)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,1),
		ylim=c(start,end))

image(c(0.25,0.75),zmids,t(legend_data),add=TRUE,zlim=c(start,end),col=cscale2)
axis(2,las=1,hadj=0.7,pos=-0.3)

dev.off()
