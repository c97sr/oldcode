# Make figure 2B
# - Assumes that ventilation.r has been run
# - who made the 

cwidth 	<- 8
cheight <- 13

pdf("~/Dropbox/projects/influenza/age_peak_attack/Current/figs/Baseline.pdf", width=cwidth/cm(1), height=cheight/cm(1))

# Gaps
g1 <- 2.5 
g2 <- 2
g3 <- 2
g3a <- 0.25
g4 <- 0.25
g5 <- 0.25

# Spaces for the plots
A <- (cwidth - g1 - g3)  
B <- (cheight - g2 - g3a - 2*g4) / 3 

# Positions of the different parts
posA <- c(g1/cwidth,(g1+A)/cwidth,(g2+2*B+2*g4)/cheight,(g2+3*B+2*g4)/cheight)
posB <- c(g1/cwidth,(g1+A)/cwidth,(g2+B+g4)/cheight,(g2+2*B+g4)/cheight)
posC <- c(g1/cwidth,(g1+A)/cwidth,g2/cheight,(g2+B)/cheight)

axticks <- as.date(c(
				"1Apr2009","1May2009","1Jun2009","1Jul2009","1Aug2009","1Sep2009",
				"1Oct2009","1Nov2009"))

axstart <- min(axticks)
axend	<- max(axticks)
axnames	<- c("1","2","3","4","5","6","7","8")
ytickv	<- c(0,0.005,0.01,0.015,0.02)
ytickn	<- c("0.0%","0.5%","1.0%","1.5%","2.0%")
ymax 	<- max(ytickv)
ytickvB	<- c(0,100,200,300)
yticknB	<- c("0","100","200","300")
ymaxB 	<- max(ytickvB)
ytickvB2	<- c(0,0.01,0.02,0.03,0.04)
yticknB2	<- c("0%","1%","2%","3%","4%")
ymaxB2 	<- max(ytickvB2)
ytickvC	<- c(0,10,20,30,40)
yticknC	<- c("0","10","20","30","40")
ymaxC 	<- max(ytickvC)
ytickvC2	<- c(0,0.001,0.002,0.003,0.004,0.005)
yticknC2	<- c("0.0%","0.1%","0.2%","0.3%","0.4%","0.5%")
ymaxC2 	<- max(ytickvC2)
lineylab <- 3.75

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
		xlim=c(axstart,axend),
		ylim=c(0,ymax))

axis(2,las=1,at=ytickv,labels=ytickn)
mtext("Incidence",side=2,line=lineylab)
axis(1,at=axticks,labels=rep("",length(axticks)))
text(axstart,max(ytickv),"A",font=2,adj=c(0,1))

pvars 	<- c("dS","dS")
pcols 	<- c("red","blue")
solvec 	<- list(sol_high_Tg_high_r_age,sol_high_Tg_low_r_age)
ltype	<- c(1,1)

for (i in 1:(length(pvars))) {
	sol <- solvec[[i]]
	#	points(sol[,"time"],sol[,pvars[i]]/popsize,type="l",col=pcols[i],lwd=2,lty=ltype[i])
}

# plotstack(sol_med_Tg_high_r_age,1,"red",var="dI",mult=1/popsize)
# plotstack(sol_med_Tg_high_r_age,2,"green",var="dI",mult=1/popsize)
# plotstack(sol_med_Tg_high_r_age,3,"blue",var="dI",mult=1/popsize)

# plotstack(sol_med_Tg_low_r_age,1,"red",var="dI",mult=1/popsize)
# plotstack(sol_med_Tg_low_r_age,2,"green",var="dI",mult=1/popsize)
# plotstack(sol_med_Tg_low_r_age,3,"blue",var="dI",mult=1/popsize)

# from http://www.nationmaster.com/graph/hea_hos_bed-health-hospital-beds

legend(as.date("15Jun2009"),max(ytickv)*3/4,
		c(	"> 54",
				"19 - 54",
				"< 19"),
		lty=c(ltype,ltype,ltype),
		lwd=c(8,8,8),
		col=c("blue","green","red"),
		cex=0.75, 
		bty="n",
)

# Part B 
par(fig=posB,new=TRUE)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axend),
		ylim=c(0,ymaxB))

axis(2,las=1,at=ytickvB,labels=yticknB)
mtext("Prevalence",side=2,line=lineylab)
mtext("(per 100,000)",side=2,line=lineylab-1)
axis(1,at=axticks,labels=rep("",length(axticks)),las=3)

text(axstart,max(ytickvB),"B",font=2,adj=c(0,1))

# plotstack(sol_med_Tg_high_r_age,1,"red",var="Ih",mult=100000/popsize)
# plotstack(sol_med_Tg_high_r_age,2,"green",var="Ih",mult=100000/popsize)
# plotstack(sol_med_Tg_high_r_age,3,"blue",var="Ih",mult=100000/popsize)

# plotstack(sol_med_Tg_low_r_age,1,"red",var="Ih",mult=100000/popsize)
# plotstack(sol_med_Tg_low_r_age,2,"green",var="Ih",mult=100000/popsize)
# plotstack(sol_med_Tg_low_r_age,3,"blue",var="Ih",mult=100000/popsize)

# Part B2 
par(fig=posB,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axend),
		ylim=c(0,ymaxB2))
axis(4,las=1,at=ytickvB2,labels=yticknB2)
mtext("Hospitalization",side=4,line=2.5)
mtext("rate per infection",side=4,line=3.5)
# points(sol_med_Tg_high_r_age[,"time"],hosprate_low,type="l")
# points(sol_med_Tg_high_r_age[,"time"],hosprate_high,type="l",col="red")

# Part C
par(fig=posC,new=TRUE)

plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axend),
		ylim=c(0,ymaxC))

axis(2,las=1,at=ytickvC,labels=yticknC)
mtext("Prevalence",side=2,line=lineylab)
mtext("(per 100,000)",side=2,line=lineylab-1)
mtext("Month since establishment",side=1,line=2)
axis(1,at=axticks,labels=axnames,las=1)

text(axstart,max(ytickvC),"C",font=2,adj=c(0,1))

# plotstack(sol_med_Tg_high_r_age,1,"red",mult=100000/popsize)
# plotstack(sol_med_Tg_high_r_age,2,"green",mult=100000/popsize)
# plotstack(sol_med_Tg_high_r_age,3,"blue",mult=100000/popsize)

# plotstack(sol_med_Tg_low_r_age,1,"red",mult=100000/popsize)
# plotstack(sol_med_Tg_low_r_age,2,"green",mult=100000/popsize)
# plotstack(sol_med_Tg_low_r_age,3,"blue",mult=100000/popsize)

# Part C2 
par(fig=posC,new=TRUE)
plot(	1:2,
		type="n",
		axes=FALSE,
		xlim=c(axstart,axend),
		ylim=c(0,ymaxC2))
axis(4,las=1,at=ytickvC2,labels=yticknC2)
mtext("ICU rate",side=4,line=3.25)
mtext("per infection",side=4,line=3.5)
# points(sol_med_Tg_high_r_age[,"time"],icurate_low,type="l")
# points(sol_med_Tg_high_r_age[,"time"],icurate_high,type="l",col="red")

dev.off()
