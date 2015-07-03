# Script to support paper ""
# S Riley
# sr@stevenriley.net
# last revised : 22nd March 2006, 1034
# Next: multinomial likelihood for the n and p in the PNAS paper

# This script may only be distributed in its original form
rm(list=ls(all=TRUE))
setwd("Z:\\projects\\sars\\amoy\\eclipse\\R Source")

datafile_ave = "D:\\files\\projects\\influenza\\ChinaPoultry\\pnas_data_average.in"
datafile_peak = "D:\\files\\projects\\influenza\\ChinaPoultry\\pnas_data_peak.in"
#datafile_ave = "C:\\tmp\\synch_messups\\pnas_data_average.in"
#datafile_peak = "C:\\tmp\\synch_messups\\pnas_data_peak.in"

require("odesolve")
require("grDevices")fitted

updatePVector <- function(values,names,pvec) {
	no_updates <- length(values)
	for (i in 1:no_updates) pvec[names[i]]=values[i]
	pvec
}

poultryNGM <- function(p) {

	D_A_dash = 1 / (1/p["D_A"] + 1/p["L_A"])
	D_C_dash = 1 / (1/p["D_C"] + 1/p["L_C"])

	ng_A_to_A = (1-p["v"])*D_A_dash*p["alpha_A"]
	ng_C_to_A = (1-p["v"])*D_C_dash*p["m"]
	ng_A_to_C = D_A_dash*p["m"]*p["alpha_A"]
	ng_C_to_C = D_C_dash

	M_NG = c(ng_A_to_A,  ng_C_to_A, ng_A_to_C, ng_C_to_C)

	dim(M_NG) <- c(2,2)

	M_NG

}

chinaPoultry <- function(t,y,p) {

	rtn <- c(0,0,0)
	names(rtn) <- c("S_A","I_A","I_C")
	names(y) <- c("S_A","I_A","I_C")

	R0_base = max(eigen(poultryNGM(p))$values)
	beta = p["R0"]/R0_base
	# beta = p["beta"]
	
	lambda_C = beta*(y["I_C"]/p["N_C"] + p["m"]*p["alpha_A"]*y["I_A"]/p["N_A"])
	lambda_A = beta*(p["m"]*y["I_C"]/p["N_C"] + p["alpha_A"]*y["I_A"]/p["N_A"])

	rtn["S_A"] = (1-p["v"])/p["L_A"]*p["N_A"] - (lambda_A + 1/p["L_A"])*y["S_A"]
	rtn["I_A"] = lambda_A*y["S_A"] - (1/p["L_A"]+1/p["D_A"])*y["I_A"]
	rtn["I_C"] = lambda_C*(p["N_C"]-y["I_C"]) - (1/p["L_C"]+1/p["D_C"])*y["I_C"]

	aux <- c(0,0)
	names(aux) <- c("prop_A","prop_C")
	aux["prop_A"] <- y["I_A"] / p["N_A"]
	aux["prop_C"] <- y["I_C"] / p["N_C"]

	list(rtn,aux)
}

fnSteadyState <- function(ics,params,model,tol,dt,max_it) {
	current_tol <- tol+1
	current_it <- 0
	no_vars <- length(ics)
	current_state <- ics
	names(current_state) <- c("S_A","I_A","I_C")
	while ((current_tol > tol) && (current_it < max_it)) {
		new_state <- lsoda(current_state,c(0,dt),model,params,rtol=1e-4)[2,2:(2+no_vars-1)]
		current_tol <- sqrt(sum((new_state-current_state)^2)/no_vars)
		current_state <- new_state
		current_it = current_it+1
	}
	rtn = list(current_state,current_it,current_it>max_it)
	names(rtn) <- c("state","iterations","max_it_exceeded")
	rtn
}

chinaPoultryLike <- function(fitted_ps,fitted_names,all_ps,all_initials,data) {
	fwdth=5
	goodvalue=TRUE	
	no_fitted = length(fitted_ps)
	no_data = dim(data)[1]
	for(i in 1:no_fitted) {
		all_ps[fitted_names[i]]=fitted_ps[i]
		if (fitted_ps[i] < 0) goodvalue=FALSE
	}

	# Ad hoc constraints here
	if (all_ps["m"]>1) goodvalue=FALSE

	if (goodvalue) {
		lnlike <- 0
		out <- fnSteadyState(all_initials,all_ps,chinaPoultry,0.1,5,1000)

		if (out$max_it_exceeded == TRUE) stop("maximum iterations exceeded")
		out <- out$state
		
		model_ps_A = c(all_ps["N_A"]-out["I_A"],out["I_A"])/all_ps["N_A"]
		model_ps_C = c(all_ps["N_C"]-out["I_C"],out["I_C"])/all_ps["N_C"]

		data_ns_A = c(data["A","N"]-data["A","n"],data["A","n"])
		data_ns_C = c(data["C","N"]-data["C","n"],data["C","n"])

		lnlike <- lnlike 	+ dmultinom(data_ns_A,prob=model_ps_A,log=TRUE)
		lnlike <- lnlike	+ dmultinom(data_ns_C,prob=model_ps_C,log=TRUE)

	} else {
		lnlike <- -1000000
	}

	report <- format(lnlike,digits=fwdth)
	report <- paste(report,"c(")
	for (i in 1:no_fitted) {
		report <- paste(report,format(fitted_ps[i],digits=fwdth),sep="")
		if (i != no_fitted) report <- paste(report,",",sep="")
	}
	report <- paste(report,")",sep="")
	for (i in 1:no_fitted) report <- paste(report,fitted_names[i])
	cat("test",report,"\n")
	flush.console()	
	lnlike
}

parVector = 
	c(
		R0 = 1.1536,
		m = 0.0071967,
		alpha_A = 1.1988,
		L_A = 49,
		L_C = 49,
		D_A = 4.5, 			# Sturm ramirez 2004
		D_C = 3.5, 			# Webster 2006 virology
		N_C = 374,
		N_A = 146,
		v = 0,
		seed_A = 1,
		seed_C = 1
	)

poultryNGM <- function(p) {
	D_A_dash = 1 / (1/p["D_A"] + 1/p["L_A"])
	D_C_dash = 1 / (1/p["D_C"] + 1/p["L_C"])
	ng_A_to_A = (1-p["v"])*D_A_dash*p["alpha_A"]
	ng_A_to_C = (1-p["v"])*D_C_dash*p["m"]
	ng_C_to_A = D_A_dash*p["m"]*p["alpha_A"]
	ng_C_to_C = D_C_dash
	M_NG = c(ng_A_to_A, ng_A_to_C, ng_C_to_A, ng_C_to_C)
	dim(M_NG) <- c(2,2)
	M_NG
}

initials = c(parVector["N_A"]-parVector["seed_A"],parVector["seed_A"],parVector["seed_C"])
names(initials) = c("S_A","I_A","I_C")

eigenResults <- eigen(poultryNGM(parVector))
max(eigen(poultryNGM(parVector))$values)
eigenResults$vector

ps_to_fit_names = c("R0","m","alpha_A")
ps_start_lb = c(0,0,0)
ps_start_ub = c(4,1,2)
no_fitted = length(ps_to_fit_names)
start_vals = rep(0,no_fitted)
no_opt_runs = 5

# poultryData = read.table(datafile_peak)
poultryData = read.table(datafile_ave)

A_bounds = qpois(c(0.025,0.975),poultryData["A","n"])
C_bounds = qpois(c(0.025,0.975),poultryData["C","n"])

poultryDataUB = poultryData
poultryDataUB[,"n"] <- c(A_bounds[2],C_bounds[2])

poultryDataLB = poultryData
poultryDataLB[,"n"] <- c(A_bounds[1],C_bounds[1])

fittedvals <- array(dim=c(no_opt_runs,length(ps_to_fit_names)+1),dimnames=list(NULL,c(ps_to_fit_names,"lnlike")))
for (i in 1:no_opt_runs) {
	for (j in 1:no_fitted) start_vals[j] <- runif(1,ps_start_lb[j],ps_start_ub[j])[1]
	poultryEst <- optim(	start_vals, chinaPoultryLike, fitted_names=ps_to_fit_names, all_ps = parVector,
					all_initials = initials, data = poultryDataLB, control=list(trace=0,fnscale=-1,maxit=10000) )
	fittedvals[i,1:no_fitted] <- poultryEst$par
	fittedvals[i,no_fitted+1] <- poultryEst$value
}

fittedvals_PI <- fittedvals

SetChoice = 3
PVecFitted <- updatePVector(fittedvals[SetChoice,1:no_fitted],ps_to_fit_names,parVector)
ts_output <- lsoda(initials,(0:10)*1000,chinaPoultry,PVecFitted,rtol=1e-4)
ts_output

chinaPoultryLike(c(1.137885,0.090010062,4.273606),ps_to_fit_names,parVector,initials,poultryData)

test_vals = c(1,1,1)
PVecTest <- updatePVector(test_vals,ps_to_fit_names,parVector)
ts_output <- lsoda(initials,(0:10)*1000,chinaPoultry,PVecTest,rtol=1e-4)
ts_output

# Old values
# c(0.15733,0.0073151,1.1803) c("beta","m","alpha_A")
# test -6.8899 c(1.3658,0.015535,1.2903) R0 m alpha_A 
# test -6.8899 c(1.2085,0.008858,1.2153) R0 m alpha_A
# test -6.8899 c(1.1536,0.0071967,1.1988) R0 m alpha_A

# Plotting stuff below here - everything needed
n_points = 51
v_samples = (0:(n_points-1))/(n_points-1)*1
delta_m_samples = (0:(n_points-1))/(n_points-1)*1
R0_field = rep(0,n_points*n_points)
dim(R0_field) <- c(n_points,n_points)

R0_with_CI = rep(0,n_points*3)
dim(R0_with_CI) <- c(n_points,3)

#ps_to_fit_names = c("R0","m","alpha_A")
chart_vals = c(1.241568,0.06081517,1.859271)
chart_vals_lb = c(1.182015,0.06981120,2.222470)
chart_vals_ub = c(1.341568,0.06081517,1.859271)

PVecChart <- updatePVector(chart_vals,ps_to_fit_names,parVector)
PVecChart_lb <- updatePVector(chart_vals_lb,ps_to_fit_names,parVector)
PVecChart_ub <- updatePVector(chart_vals_ub,ps_to_fit_names,parVector)

R0_base = PVecChart["R0"]
lambda_base = max(eigen(poultryNGM(PVecChart))$values)
m_base = PVecChart["m"]

R0_base_lb = PVecChart_lb["R0"]
lambda_base_lb = max(eigen(poultryNGM(PVecChart_lb))$values)

R0_base_ub = PVecChart_ub["R0"]
lambda_base_ub = max(eigen(poultryNGM(PVecChart_ub))$values)

for (i in 1:n_points) {
	for (j in 1:n_points) {
		tmpPV = updatePVector(c(v_samples[i],(1-delta_m_samples[j])*m_base),c("v","m"),PVecChart)
		lambda_dash = max(eigen(poultryNGM(tmpPV))$values)
		R0_field[i,j] = R0_base * lambda_dash / lambda_base
	}
}

for (i in 1:n_points) {

		tmpPV = updatePVector(c(v_samples[i]),c("v"),PVecChart)
		lambda_dash = max(eigen(poultryNGM(tmpPV))$values)
		R0_with_CI[i,1] = R0_base * lambda_dash / lambda_base

		tmpPV = updatePVector(c(v_samples[i]),c("v"),PVecChart_lb)
		lambda_dash = max(eigen(poultryNGM(tmpPV))$values)
		R0_with_CI[i,2] = R0_base_lb * lambda_dash / lambda_base_lb

		tmpPV = updatePVector(c(v_samples[i]),c("v"),PVecChart_ub)
		lambda_dash = max(eigen(poultryNGM(tmpPV))$values)
		R0_with_CI[i,3] = R0_base_ub * lambda_dash / lambda_base_ub

}

# R0_with_CI

# Bit for Prof Yuan's paper from here ...
windows(width=5,height=4)
par(fig=c(0,1,0,1))
par(mai=(c(0.75,0.75,0.1,0.1)))
par(mgp=c(1.5,0.75,0))
x_limits=c(0,1)

plot(	1:2,
	type="n",
	axes=FALSE, 
	xlab="Effective coverage of hatchling aquatic birds",
	ylab="Basic reproductive number",
	ylim=c(0,1.5),
	xlim=c(0,1)
)

axis(2,pos=0,at=(0:3)/3*1.5, labels=c("0.0","0.5","1.0","1.5"),las=1)
axis(1,pos=0,at=(0:5)/5, labels=c("0%","20%","40%","60%","80%","100%"),las=1)

points(v_samples,R0_with_CI[,1],type="l")
lines(v_samples,R0_with_CI[,2],lty=2)
lines(v_samples,R0_with_CI[,3],lty=2)
lines(0:100/100,rep(1,101),lty=3)

# ... to here

plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("c",side=2,las=1,line=4,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_QIA$S_365[1:max_i]/popsize,cex=ps,col="cyan",pch=20)

# Plot routines for figure 4
windows(width=5,height=4)
par(fig=c(0,0.5,0,0.8))
par(mai=(c(0.85,0.85,0.5,0.5)))
par(mgp=c(2,1,0))

filled.contour(v_samples,delta_m_samples,R0_field,nlevels=50,color = heat.colors,
	xlab="Vaccine coverage (v)",ylab = "Effectiveness of reduced mixing \n (dm, where m -> {1-dm}m)",
	key.title="R0")
contour(v_samples,delta_m_samples,R0_field,levels=c(1),add=TRUE)	
# filled.contour(v_samples,delta_m_samples,R0_field,nlevels=50,col = rgb(150:200/200,0,0))



