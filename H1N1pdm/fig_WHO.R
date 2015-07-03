rm(list=ls())
source("./funcs.r")

# Next : include the illustration and organize the text.

cwidth 	<- 10
cheight <- 10

illParams <- function(){
	
	c(		N 			= 1,
			seed 		= 1/7000000,
			r 			= log(2)/5,
			Tg 			= 2.6,
			phi_1		= 1.0,
			phi_2		= 1.0,
			p_R 		= 0.86,
			p_H_base 		= 0.01,
			p_H_base_1 		= 0.01,
			p_H_base_2 		= 0.01,
			p_H_base_3 		= 0.01,
			p_I 		= 0.15,
			p_I_1 		= 0.15,
			p_I_2 		= 0.15,
			p_I_3 		= 0.15,
			gamma_v 	= 1/5.2,						# Gowardman
			gamma_h 	= 1/7.0,						# Feagan		
			t0 			= as.date("1Jun2009"),			
			p_1			= 0.2,										
			p_2			= 0.5,
			mixmatindex	= 1,
			ngmindex	= 1,
			fitsus		= 1
	)
	
}

illParams2 <- function(){
	
	c(		N 			= 1,
			seed 		= 1/7000000,
			r 			= log(2)/5,
			Tg 			= 2.6,
			phi_1		= 1.0,
			phi_2		= 1.0,
			phi_3		= 1.0,
			phi_4		= 1.0,
			p_R 		= 0.86,
			p_H_base 	= 0.01,
			p_H_base_1 	= 0.01,
			p_H_base_2 	= 0.01,
			p_H_base_3 	= 0.01,
			p_H_base_4 	= 0.01,
			p_I 		= 0.15,
			p_I_1 		= 0.15,
			p_I_2 		= 0.15,
			p_I_3 		= 0.15,
			p_I_4 		= 0.15,
			gamma_v 	= 1/5.2,							# NEJM ANZ piece 10.1056/NEJMoa0908481
			gamma_h 	= 1/7.0,							# Feagan
			t0 			= as.numeric(as.date("1Jun2009")),
			tf			= as.numeric(as.date("31Dec2009")),
			p_1			=  0.2,	
			p_2			=  0.5,								# Weighted probability factors
			p_3			= 0.3,
			p_4 		=  0.001,
			mixmatindex	= 1,									# Select 
			fitsus		= 1
	)	
}


xdates <- c("1Aug2009","1Sep2009","1Oct2009","31Oct2009")
xtcklab <- 3:6
xdatesv <- as.date(xdates)
xdatest <- xdates 
y1vals <- c(0,0.005,0.01,0.015)
y1labs <- c("0.0%","0.5%","1.0%","1.5%")
y2vals <- (0:5)*2
y2labs <- y2vals

posA <- c(	0.2,	0.8,	0.2,	0.95)

loopfiles <- c("./figs/WHO_sero_AB.pdf","./figs/WHO_sero_CD.pdf")
loopcasekids <- c(20,80)
loopcaseelderly <- c(40,10)
loop_phc <- c(0.01,0.001)
loop_pha <- c(0.01,0.01)
loop_phe <- c(0.01,0.05)
loop_mm <- c(1,1)

for (i in 1:2) {
	
	# Make the model solutiuon
	baseline <- (SirModel3AgeClasses(
			pname=c("p_H_base_1","p_H_base_2","p_H_base_3","mixmatindex"),
			pvals=c(loop_phc[i],loop_pha[i],loop_phe[i],loop_mm[i]),
			casemix=c(loopcasekids[i],40,40),
			vp=illParams2()))$sol
	
	# t(baseline[1:3,])
	# sum(baseline[,"dS"]) / illParams()["p_R"]
	
	pdf(file=loopfiles[i],width=cwidth/cm(1), height=cheight/cm(1))
	
	# Set formats
	par(cex=1.0)
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	par(mai=(c(0,0,0,0)))
	par(cex.axis=0.9)
	
	par(fig=posA)
	
	plot(	1:2,
			type="n",xlab="",ylab="",main="",log="",
			axes=FALSE,xlim=c(min(xdatesv),max(xdatesv)),
			ylim=c(min(y1vals),max(y1vals)))
	axis(1,at=xdatesv,labels=xtcklab)
	# axis(2,las=1,at=y1vals,labels=y1labs)
	
	plotstack(baseline,1,"red",var="dI",mult=1,const=0)
	plotstack(baseline,2,"green",var="dI",mult=1,const=0)
	plotstack(baseline,3,"green",var="dI",mult=1,const=0)
	
	par(fig=posA,new=TRUE)
	plot(	1:2,
			type="n",xlab="",ylab="",main="",log="",
			axes=FALSE,xlim=c(min(xdatesv),max(xdatesv)),
			ylim=c(min(y2vals),max(y2vals)))
	#axis(4,las=1,at=y2vals,labels=y2labs)
	points(baseline[0:30*5,"time"],100000*baseline[0:30*5,"Iv"],type="l",col="black")
	
	dev.off()

}