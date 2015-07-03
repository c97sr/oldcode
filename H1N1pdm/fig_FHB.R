rm(list=ls(all=TRUE))
options(error=NULL)

source("./funcs.r")

# Load up the data
hkdata <- LoadCleanCHPData("~/files/projects/influenza/swineflu/data/chp/cleaned_sr.csv")
usdata 	<- postProcLineList(read.csv("~/files/projects/influenza/swineflu/data/cdc/linelists13may_new.csv"))
save(file="~/tmp_nh1n1_objects.R",usdata)
# load("~/tmp_nh1n1_objects.R")

# Set general parameters
today 		<- max(hkdata$conf_date,na.rm=TRUE)
lastImpData	<- today - 4

# Plot the age distribution
age_freq_us <- AddZeros1DTable(table(round(usdata$cleanAge)))
age_freq_hk <- AddZeros1DTable(table(round(hkdata$age)))
total_school_age_us <- sum(age_freq_us$vals[match(0,age_freq_us$days):match(19,age_freq_us$days)])
total_school_age_hk <- sum(age_freq_hk$vals[match(0,age_freq_hk$days):match(19,age_freq_hk$days)])
pdf("./figs/FHB_age_dist.pdf",height=10/cm(1),width=15/cm(1))
plot(	age_freq_hk$days[1:20]-0.1,age_freq_hk$vals[1:20]/total_school_age_hk,
		xlim=c(-0.3,19.3),ylim=c(0,0.2),col="red",type="h",lwd=2,
		xlab="age",ylab="Density")
mtext(paste("School-age cases (US, n=",total_school_age_us,", blue; HK, n=",total_school_age_hk,", red)"))
points(age_freq_us$days[2:21]+0.1,age_freq_us$vals[2:21]/total_school_age_us,col="blue",type="h",lwd=2)
# write.csv(cbind(age=0:19,prob=age_freq_us$vals[2:21]/total_school_age_us),file="outtmp.csv",row.names=FALSE)
dev.off()

# Calculate the best fit for the inputations - a little over the top. Mean would do this!
test <- LikePoissonOnsetReport(c(3),hkdata$onset,hkdata$conf_date,hkdata$imponset)
fit <- optimize(LikePoissonOnsetReport,c(0,10),maximum=TRUE,vecOnsets=hkdata$onset,vecReports=hkdata$conf_date,vecImputes=hkdata$imponset)

# Estimate the proportion
recent_simpimp <- hkdata[hkdata$conf_date %in% c(lastImpData,lastImpData-1,lastImpData-2),]$simpimp
p_imp <- table(recent_simpimp)[1] / sum(table(recent_simpimp))

# Impute the data
hkdata$tmponset <- ImputeOnsets(fit$maximum,hkdata$imponset,hkdata$conf_date,hkdata$onset)
hkdata$tmpimport <- ImputeImport(p_imp,hkdata$impimp,hkdata$simpimp)

# Plot the imputed distribution
pdf("./figs/FHB_fig_impute.pdf",height=10/cm(1),width=15/cm(1))
plot(	table(hkdata[hkdata$imponset==FALSE,"conf_date"]-hkdata[hkdata$imponset==FALSE,"onset"]),
		xlab="Delay",ylab="Count",main="")
points(0:15,dpois(0:15,fit$maximum)*(table(hkdata$imponset)[1]),type="l")
dev.off()

# Generate local transmission cases
local_inc <- AddZeros1DTable(table(hkdata[hkdata$tmpimport=="NotImported",]$tmponset))
plot(local_inc$days,local_inc$vals)

pvec <- as.vector(mode="numeric",data.frame(1,0.5))
start_time <- as.date("4Jun2009")
end_time <- today-7
names(pvec) <- c("A","r")
testlike 	<- HKLikeInc(pvec,local_inc$vals[1:(length(local_inc$vals)-6)])
mle1 		<- optim(	pvec,HKLikeInc,control=list(trace=0,fnscale=-1,maxit=1000),
						inc=local_inc$vals[1:(end_time-start_time+1)])
mle2 		<- optim(	pvec,HKLikeInc,control=list(trace=0,fnscale=-1,maxit=1000),
						inc=local_inc$vals[1:(end_time-start_time+1-1)])
mle3 		<- optim(	pvec,HKLikeInc,control=list(trace=0,fnscale=-1,maxit=1000),
						inc=local_inc$vals[1:(end_time-start_time+1+1)])

1+mle1$par["r"]*2.7
1+mle2$par["r"]*2.7
1+mle3$par["r"]*2.7

# Plot the growth characteristics
pdf("./figs/FHB_growth.pdf",height=10/cm(1),width=15/cm(1))
xdatesv <- as.date("1Jun2009")+(0:4)*7
xdatest <- c("1 Jun","8 Jun","15 Jun","22 Jun","29 Jun") 
plot(	table(hkdata[hkdata$tmpimport=="NotImported",]$tmponset),
		type="h",xlab="Date",ylab="Incidence",main="",
		axes=FALSE,xlim=c(min(xdatesv),max(xdatesv)),ylim=c(0,70))
axis(1,at=xdatesv,labels=xdatest)
axis(2,las=1)
points(start_time:(end_time),mle1$par["A"]*exp((start_time:(end_time)-start_time)*mle1$par["r"]),type="l",col="red")
points(start_time:(end_time-1),mle2$par["A"]*exp((start_time:(end_time-1)-start_time)*mle2$par["r"]),type="l",col="blue")
points(start_time:(end_time+1),mle3$par["A"]*exp((start_time:(end_time+1)-start_time)*mle3$par["r"]),type="l",col="green")
legend(as.date("1Jun2009"),60,
		c(	paste("End day ",end_time," Td=",round(log(2)/mle1$par["r"],2)),
			paste("End day ",end_time-1," Td=",round(log(2)/mle2$par["r"],2)),
			paste("End day ",end_time+1," Td=",round(log(2)/mle3$par["r"],2))),
		lty=c(1,1,1),
		lwd=c(1,1,1),
		col=c("red","blue","green"),
		cex=0.75, 
		bty="n",
)
dev.off()

# Calculate age specific rates for Hong Kong, for age groups
cm_tmp <- table(hkdata[hkdata$tmpimp=="NotImported","ageclass2"])
cm_tmp_per <- cm_tmp / sum(cm_tmp)
sum(cm_tmp)
casem <- c(cm_tmp["y"],cm_tmp["o"],cm_tmp["a"])
hkmm <- MixingCCA()																# Mixing habits of different age groups as per the polymod data
#hkmm <- matNull()
hkN  <- c(12/15*0.129,3/15*0.123+4/10*0.130,6/10*0.130+1-0.129-0.130)			# Proportion of population in young children, older children adn adults
ps <- c(phi_1=1,phi_2=1)

# Hong kong epidemiology
hk_ngm <- function(mix,params,vecN) {
	rtn <- 1/params["gamma"]*t(array(c(	
							params["phi_1"]*mix[1,1],					params["phi_1"]*(vecN[1]*mix[1,2]/vecN[2]),	params["phi_1"]*(vecN[1]*mix[1,3]/vecN[3]),
							params["phi_2"]*(vecN[2]*mix[2,1]/vecN[1]),	params["phi_2"]*mix[2,2],					params["phi_2"]*(vecN[2]*mix[2,3]/vecN[3]),
							vecN[3]*mix[3,1]/vecN[1],					vecN[3]*mix[3,2]/vecN[2],					mix[3,3]
					), dim=c(3,3)))
	rtn
}


likeSusEigen(c(2,1),ps,hkN,hkmm,casem,hk_ngm)
hkSusFit <- optim(	c(2,1),likeSusEigen,control=list(trace=0,fnscale=-1,maxit=10000),
					params=ps, vecN=hkN,mix=hkmm,cm=casem,ngm=hk_ngm)

mG 			<- hk_ngm(hkmm,c(phi_1=hkSusFit$par[1],phi_2=hkSusFit$par[2],gamma=1),hkN)
ev  		<- eigen(mG)$vector[,1]
tmp 		<- ev * mG
ratios		<- c(p_y=sum(tmp[,1]),p_o=sum(tmp[,2]),p_2=sum(tmp[,3]))/sum(tmp)
segmentedR0	<- c(sum(mG[,1]),sum(mG[,2]),sum(mG[,3]))/sum(mG)*1.5

# Projections

hkParams <- function(){
	c(		N 			= 7000000,
			seed 		= 1,
			r 			= 0.26,
			Tg 			= 2.7,
			phi_1		= 2.0,
			phi_2		= 1.0,
			p_R 		= 0.86,
			p_H_base 		= 0.031,
			p_H_base_1 		= 0.02489,
			p_H_base_2 		= 0.02489,
			p_H_base_3 		= 0.03575,
			p_I 		= 1.0,
			p_I_1 		= 1.0,
			p_I_2 		= 1.0,
			p_I_3 		= 1.0,
			gamma_v 	= 1/5.2,									# Gowardman
			gamma_h 	= 1/7.0,									# Feagan		
			t0 			= as.date("4Jun2009"),
			p_1			= 12/15*0.129,			
			p_2			= 3/15*0.129+4/10*0.130,
			mixmatindex	= 2,
			ngmindex	= 2
	)
}

# baseline <- (SirModel3AgeClasses(	pname=c("r","seed"),pvals=c(0.2,1),
#							casemix=casem,vp=hkParams(),NGM=hk_ngm))$sol
popmix <- c(0.1032,0.0778,0.819)
excasemix <- c(37,228,179)
baseline <- (SirModel3AgeClasses(	pname=c("r","seed","mixmatindex","p_1","p_2"),pvals=c(0.2,1,2,0.33,0.33),
							casemix=c(1,1,1),vp=hkParams(),NGM=hk_ngm))$sol

# Put in the case mix according to the population size					
baseline <- (SirModel3AgeClasses(	pname=c("r","seed","mixmatindex"),pvals=c(0.18,1,2),
							casemix=c(37,228,179),vp=hkParams(),NGM=hk_ngm))$sol

attack_rate <- (baseline[365,"dS1"]+baseline[365,"dS2"]+baseline[365,"dS3"]) / 7000000
max(baseline[,"dS"])
					
# Plot the projections
pdf("./figs/FHB_Projection.pdf",height=10/cm(1),width=15/cm(1))
xdates <- c("1Jun2009","1Jul2009","1Aug2009","31Aug2009","30Sep2009")
xdatesv <- as.date(xdates)
xdatest <- xdates 
y1vals <- c(0,5000,10000,15000,20000,25000,30000)
plot(	1:2,
		type="n",xlab="",ylab="",main="Projected daily incidence of symptomatic cases",log="",
		axes=FALSE,xlim=c(min(xdatesv),max(xdatesv)),
		ylim=c(min(y1vals),max(y1vals)))
axis(1,at=xdatesv,labels=xdatest)
axis(2,las=1,at=y1vals)
plotstack(baseline,1,"red",var="dI",mult=1,const=0)
plotstack(baseline,2,"green",var="dI",mult=1,const=0)
plotstack(baseline,3,"blue",var="dI",mult=1,const=0)
dev.off()

# Number for the text are below and are being edited by Danny
cat(paste("analyses are based on a database of: ",dim(hkdata)[1]))
cat(paste("to and including 23 June: ","look in data file for this"))

 