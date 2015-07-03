# Remove all objects 
rm(list=ls(all=TRUE))
options(error=recover)
#options(error=NULL)
require("odesolve")
setwd("C:/eclipse/R Source/gz_syphilis")

seirsModel <- function(t,y,p) {
	rtn <- array(0,c(5))
	names(rtn) <- c("S","E","I","R","dS")
	names(y) <- c("S","E","I","R","dS")
	lambda <- p["beta"]*y["I"]/(y["S"]+y["E"]+y["I"]+y["R"])
	rtn["S"] <- -lambda*y["S"] + p["omega"]*y["R"]+(1-p["p_i"])*p["gamma"]*y["I"]
    rtn["E"] <- lambda*y["S"]-p["theta"]*y["E"]
	rtn["I"] <-  p["theta"]*y["E"]- p["gamma"]*y["I"]
	rtn["R"] <- p["p_i"]*p["gamma"]*y["I"] - p["omega"]*y["R"]
	rtn["dS"] <- lambda*y["S"]
	list(rtn)
}

seirs_LH_Model <- function(t,y,p) {
	rtn <- array(0,c(10))
	names(rtn) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
	names(y) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dH","dL") 
	lambda_L <- p["beta"]*(p["m_LL"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HL"]*y["I_H"]/p["N_H"])
	lambda_H <- p["beta"]*(p["m_LH"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HH"]*y["I_H"]/p["N_H"])
	rtn["S_L"] <- -lambda_L*y["S_L"] + p["omega"]*y["R_L"]+(1-p["p_i"])*p["gamma"]*y["I_L"]
	rtn["E_L"] <- lambda_L*y["S_L"]-p["theta"]*y["E_L"]
	rtn["I_L"] <-  p["theta"]*y["E_L"]- p["gamma"]*y["I_L"]
	rtn["R_L"] <- p["p_i"]*p["gamma"]*y["I_L"] - p["omega"]*y["R_L"]
	rtn["S_H"] <- -lambda_H*y["S_H"] + p["omega"]*y["R_H"]+(1-p["p_i"])*p["gamma"]*y["I_H"]
	rtn["E_H"] <- lambda_H*y["S_H"]-p["theta"]*y["E_H"]
	rtn["I_H"] <-  p["theta"]*y["E_H"]- p["gamma"]*y["I_H"]
	rtn["R_H"] <- p["p_i"]*p["gamma"]*y["I_H"] - p["omega"]*y["R_H"]
	rtn["dH"] <- lambda_H*y["S_H"]
	rtn["dL"] <- lambda_L*y["S_L"]
	list(rtn)
}

fnPoissonLike <- function(x,y) {
	# y are model incidences
	# x are data
	# output is log likelihood of x if y mean of independent poisson variables
	noSamples <- length(y)
	lnlike <- 0
	mintol <- 1e-20
	for (i in 1:noSamples) {
		if (y[i] < mintol && x[i] > mintol) {
			lnlike <- lnlike - 1000000
		} else {
			lnlike <- lnlike + dpois(x[i],y[i],log=TRUE)
		}
	}
	lnlike
}

fnRepPs <- function(allp,pfn,pfv) {
	nofitted <- length(pfv)
	for (i in 1:nofitted) {
		allp[pfn[i]]<-pfv[i]
	}
	allp
}

fnModelLike <- function(v,pnamesfitted,allps,alldata_h,alldata_l,report=TRUE,plot=TRUE) {

	# v <- c(1000,1)
	# pnamesfitted <- c("N","seed")
	# allps <- seirsParams
	# alldata_h <- incData_H
	# alldata_l <- incData_L
			
	allps <- fnRepPs(allps,pnamesfitted,v)
	
	if (	allps["seed"] > allps["N"] 	||
			allps["seed"] < 1 			||
			FALSE) {
		like <- -10000000
	} else {
		ics <- c(allps["N"]-allps["N_H"],0,0,0,allps["N_H"]-allps["seed"],0,allps["seed"],0,0,0)
		names(ics) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
		no_years <- 14
		time_points <- 0:no_years
		timestep <- 1
		solution <- lsoda(ics,time_points,seirs_LH_Model,allps)
		yL <- array(-100,no_years)
		yH <- array(-100,no_years)
		for (i in 1:no_years) {
			yL[i] <- solution[i+1,"dL"]-solution[i,"dL"]
			yH[i] <- solution[i+1,"dH"]-solution[i,"dH"]
		}
		like <- fnPoissonLike(alldata_l,yL/(allps["N"]-allps["N_H"])*allps["births"]*allps["p_cgs"])
		like <- like + fnPoissonLike(alldata_h,yH)
	}
	
	if (report) {
		cat(like,allps,names(allps),"\n")
	}
	if (plot) { 
		# To be completed
	}
	flush.console()
	like
}

seirsParams <- c(
		N = 5680000,					# Size of the population
		N_H = 10000,					# Number in high risk group
		seed = 1,					# Number initially infectious
		omega = 1/(26/52),			# Duration immune stage
        theta = 1/(2/52),			# Duration of latent period
		gamma = 1/2,				# Duration of infectiousness
		p_i = 0.2,					# Proportion gaining immunity
		beta = 1,					# Infectiousness
		m_LH = 0,					# Relative infectivity of low risk to high risk 
		m_HH = 1,					# Relative infectivity of hight risk to high risk
		m_LL = 0,					# Relatvie infectivity of low risk to low risk
		m_HL = 1,					# Relative infectivity of high risk to low risk
		births = 70000,				# 10 year approx average from www.gzfb.gov.cn
		p_cgs = 0.5		# guess
		)

tabData <- read.csv("data\\syphilis_gz_1994_to_2007.csv",header=TRUE)
incData_H <- tabData$cases.male + tabData$cases.female
incData_L <- tabData$congenital

# Test of likelihood function
fnModelLike(c(10,1),c("seed","beta"),seirsParams,incData_H,incData_L)

psToFit <- c("beta","N_H","seed","m_HL","p_cgs")
psInitial <- c(1,10000,10,1,0.5)

fit <- optim(	psInitial,
				fn=fnModelLike,
				control=list(trace=0,fnscale=-1,maxit=1000),
				pnamesfitted=psToFit,
				allps = seirsParams, 
				alldata_h = incData_H,
				alldata_l = incData_L)

# c(1.268885e+00,1.091744e+04,1.740912e+02,2.963224e-03)
# c("beta","N_H","seed","m_HL")
		
# 1.269379e+00 1.094424e+04 1.704990e+02 5.189741e+00 1.804359e-03
# c("beta","N_H","seed","m_HL","p_cgs")

no_years <- 14
time_points <- 0:no_years
timestep <- 1

plotPs <- fnRepPs(seirsParams,psToFit,fit$par)
initial_conditions <- c(plotPs["N"]-plotPs["N_H"],0,0,0,plotPs["N_H"]-plotPs["seed"],0,plotPs["seed"],0,0,0)
names(initial_conditions) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
solution <- lsoda(initial_conditions,time_points,seirs_LH_Model,plotPs)
modelPnS <- solution[2:(no_years+1),"dH"]-solution[1:(no_years),"dH"]
modelCon <- solution[2:(no_years+1),"dL"]-solution[1:(no_years),"dL"]

# Command to make figure 3.1
windows(width=5.75,height=4)
plot(	tabData$year,
		incData_H,
		ylab="Incidence",
		xlab="Year",
		main="Primary and secondary",
		ylim=c(0,max(incData_H,modelPnS)),
		type="p",
		col="red",
		lwd=2)
points(	tabData$year,
		modelPnS,
		type="l",
		col="green",
		lwd=2)

windows(width=5.75,height=4)
plot(	tabData$year,
		incData_L,
		ylab="Incidence",
		xlab="Year",
		main="Congenital",
		ylim=c(0,max(incData_L,modelCon)),
		type="p",
		col="red",
		lwd=2)
points(	tabData$year,
		modelPnS,
		type="l",
		col="green",
		lwd=2)

