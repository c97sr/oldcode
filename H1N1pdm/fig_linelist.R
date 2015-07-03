# Make figure 1
# Assumes that linelists.r has been run

# See May 12 2009 entry in notebook for the key
# Units of length in cm

cwidth 	<- 20
cheight <- 10

pdf(file="./figs/Linelist.pdf",width=cwidth/cm(1),height=cheight/cm(1),version="1.4")

# Gaps
g1 <- 2 
g2 <- 2.5
g3 <- 2.0
g4 <- 0.5 
g5 <- 0.25

# Spaces for the plots
A <- (cwidth - g1 - g3 - g4 - 2*g5) / 14 
B <- (cheight - g2 - g5 - g4) / 2 

# Positions of the different parts
posA <- c(g1/cwidth,(g1+6*A)/cwidth,g2/cheight,(cheight-g4)/cheight)
posB1 <- c((g1+6*A+g3)/cwidth,(g1+7*A+g3)/cwidth,(g2+B+g5)/cheight,(g2+2*B+g5)/cheight)
posB2 <- c((g1+7*A+g3+g5)/cwidth,(g1+8*A+g3+g5)/cwidth,(g2+B+g5)/cheight,(g2+2*B+g5)/cheight)
posB3 <- c((g1+8*A+g3+2*g5)/cwidth,(g1+14*A+g3+2*g5)/cwidth,(g2+B+g5)/cheight,(g2+2*B+g5)/cheight)
posC1 <- c((g1+6*A+g3)/cwidth,(g1+7*A+g3)/cwidth,(g2)/cheight,(g2+B)/cheight)
posC2 <- c((g1+7*A+g3+g5)/cwidth,(g1+8*A+g3+g5)/cwidth,(g2)/cheight,(g2+B)/cheight)
posC3 <- c((g1+8*A+g3+2*g5)/cwidth,(g1+14*A+g3+2*g5)/cwidth,(g2)/cheight,(g2+B)/cheight)

axstart <- as.date("30Mar2009")
noweeks <- 6
axticks <- 0:(noweeks)*7+axstart
axlabels	<- c("30 Mar","6 Apr","13 Apr","20 Apr","27 Apr","4 May","11 May")

# Parameter settings for all parts of the figure
par(cex=0.8)
par(mgp=c(2.5,1,0))
# par(tcl=0.25)
par(mai=(c(0,0,0,0)))
par(cex.axis=0.9)

# Part A 
par(fig=posA)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axstart+noweeks*7),
		ylim=c(1,200),
		log="y")

axis(2,las=1)
mtext("Day of onset",side=1,line=4)
mtext("Incidence (+1)",side=2,line=3.0)
#axis(1,at=axticks,labels=axticks)
axis(1,at=axticks,las=3,labels=axlabels)
points((freq1$ub+freq1$lb)/2,freq1$f+1,type="p",col="Black")
points((freq2$ub+freq2$lb)/2,freq2$f+1,type="p",col="Red")
# points((freq1$ub+freq1$lb)/2,freq1$model+1,type="l",lwd=1.0,col="Black")
points((freq2$ub+freq2$lb)/2,freq2$model+1,type="l",lwd=1.0,col="Red")

legend(axstart,100,
		c("8 May linelist","13 May linelist"),
		pch=c(1,1),
		col=c("Black","Red"),
		cex=1.0, 
		bty="n",
)
text(axstart,200,"A",font=2,adj=c(0,1))

# Make the part C plots
# text(axstart,1,"B",font=2)

# Part B1
par(fig=posB1,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,7),
		ylim=c(lbrates,1),
		log="y")
axis(2,las=1)
mtext("Hosp Admission",side=2,line=3)
text(1,1,"B",font=2,adj=c(0,1))
points(3.5,hospital_all_pt,col="Black")
lines(c(3.5,3.5),c(hospital_all_ub,hospital_all_lb),col="Black")

# Part B2
par(fig=posB2,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,7),
		ylim=c(lbrates,1),
		log="y")
points(3.5,hosp_no_onset_pt,col="Black")
lines(c(3.5,3.5),c(hosp_no_onset_ub,hosp_no_onset_lb),col="Black")

# Part B3
par(fig=posB3,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axstart+noweeks*7),
		ylim=c(lbrates,1),
		log="y")

points(vecOnsetDatesWk+3.5,pltHospWeek$pt,col="Blue")
 for (i in 1:noOnsetDatesWk) 
	lines(c(vecOnsetDatesWk[i]+3.5,vecOnsetDatesWk[i]+3.5),c(pltHospWeek$ub[i],pltHospWeek$lb[i]),col="Blue")

stday <- fb-2 # for diff between freq table and vecOnset dates
endday <- lb-2
faintred <- rgb(1,0,0,alpha=0.5)
points(vecOnsetDates[stday:endday]+0.5,pltHospDay$pt[stday:endday],col=faintred)
for (i in stday:endday) 
	lines(c(vecOnsetDates[i]+0.5,vecOnsetDates[i]+0.5),c(pltHospDay$ub[i],pltHospDay$lb[i]),col=faintred)
#abline(x.glm)

# Part C1
par(fig=posC1,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,7),
		ylim=c(lbrates,1),
		log="y")
axis(2,las=1)
axis(1,at=c(0,7),las=3,labels=c("",""))
mtext("All",side=1,line=2.5)
mtext("ICU | Hosp",side=2,line=3)
text(1,1,"C",font=2,adj=c(0,1))
points(3.5,ICU_all_pt,col="Black")
lines(c(3.5,3.5),c(ICU_all_ub,ICU_all_lb),col="Black")

# Part C2
par(fig=posC2,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(0,7),
		ylim=c(lbrates,1),
		log="y")
axis(1,at=c(0,7),las=3,labels=c("",""))
mtext("No",side=1,line=2)
mtext("onset",side=1,line=3)
points(3.5,ICU_no_onset_pt,col="Black")
lines(c(3.5,3.5),c(ICU_no_onset_ub,ICU_no_onset_lb),col="Black")

# Part C3
par(fig=posC3,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axstart+noweeks*7),
		ylim=c(lbrates,1),
		log="y")
axis(1,at=axticks,las=3,labels=axlabels)
mtext("Week of onset",side=1,line=4)
points(vecOnsetDatesWk+3.5,pltICUWeek$pt,col="Blue")
for (i in 1:noOnsetDatesWk) 
	lines(c(vecOnsetDatesWk[i]+3.5,vecOnsetDatesWk[i]+3.5),c(pltICUWeek$ub[i],pltICUWeek$lb[i]),col="Blue")

dev.off()
