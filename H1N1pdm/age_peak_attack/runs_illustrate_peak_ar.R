# Script to make a figure that illustrates the imprtance of heterogeneity
# its likely you might want to run this change directory

# Needs to be switched to peak incidence of ICU admission per day of onset
# Set this working directory if needed setwd("~/Dropbox/svneclipse/H1N1pdm/age_peak_attack")

rm(list=ls(all=TRUE))
source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
options(error=NULL)
require("date")
require("odesolve")

# Define the illustrative parameters
ill_params <- 
		c(		N 			= 7000000,
				seed 		= 100,
				trickle		= 1.0,
				r 			= (1.4-1.0)/2.6,
				Tg 			= 2.6,
				phi_1		= 1.0,
				phi_2		= 1.0,
				p_R 		= 1.0,
				p_H_base 	= 0.01,
				p_H_base_1 	= 1.0,
				p_H_base_2 	= 1.0,
				p_H_base_3 	= 1.0,
				p_I 		= 1.0,
				p_I_1 		= 1.0,
				p_I_2 		= 1.0,
				p_I_3 		= 1.0,
				gamma_v 	= 1/5.2,									# Gowardman
				gamma_h 	= 1/7.0,									# Feagan		
				t0 			= as.date("1Apr2009"),
				p_1			= 0.2,			
				p_2			= 0.4,
				p_3			= 0.4,
				mixmatindex	= 4,										# CAE from polymod
				mixsense	= 0,
				fitsus		= 0,
				noyears		= 1,
				nodts		= 52*7,
				dur_seas	= 7*13,
				amp_seas	= 0.10,
				gamma_R		= 1/99999999,
				aging_on	= 0,
				mu_1		= 1/364/19, 
				mu_2		= 1/364/(65-19),
				mu_3		= 1/364/50									# Set to give very slow overall population growth
		)

illRateICU <- 10/100000

p 	 <- apa.setup.params.spq()
s2 	 <- apa.setup.state.vector.2state(p)

# Mass action and uniform susceptibility
sol_part_A <- (SirModel3AgeClasses(
					pname=c("phi_1","mixsense","fitsus"),
					pvals=c(1,0,0),
					casemix=c(20,40,40),
					vp=ill_params
			))$sol

sol_part_A_dash <- lsoda(s2,0:180,apa.spa.2state.mix,p)
sol_part_B_dash <- lsoda(s2,0:180,apa.spa.2state.mix,p)
sol_part_C_dash <- lsoda(s2,0:180,apa.spa.2state.mix,p)

# Mass action and uniform susceptibility
sol_part_B <- (SirModel3AgeClasses(
					pname=c("phi_1","mixsense","fitsus"),
					pvals=c(2,0,0),
					casemix=c(20,40,40),
					vp=ill_params
			))$sol

# Mass action and uniform susceptibility
sol_part_C <- (SirModel3AgeClasses(
					pname=c("phi_1","mixsense","fitsus"),
					pvals=c(2,1,0),
					casemix=c(20,40,40),
					vp=ill_params
			))$sol

popsize <- ill_params["N"]
sum(sol_part_A[,"dS"])/popsize
sum(sol_part_B[,"dS"])/popsize
sum(sol_part_C[,"dS"])/popsize

# Runs for the sensitivity analysis figure
x_lims 		<- c(0,1.2)		# r
y_lims 		<- c(0.5,5)		# p_R * p_I * p_H_base
no_bins_x	<- 4
no_bins_y	<- 4			# often either 2 or 10
x_size 		<- (x_lims[2]-x_lims[1])/no_bins_x
y_size		<- (y_lims[2]-y_lims[1])/no_bins_y
x_bounds	<- x_lims[1]+(0:no_bins_x) * x_size
y_bounds	<- y_lims[1]+(0:no_bins_y) * y_size
x_mids		<- x_lims[1] + (1:no_bins_x) * x_size - x_size / 2
y_mids		<- y_lims[1] + (1:no_bins_y) * y_size - y_size / 2

# Make the array to be plotted
rvals <- c((1.2-1.0)/2.6,(1.4-1.0)/2.6,(1.8-1.0)/2.6)
# rvals <- c((1.4-1.0)/2.6)
noRs <- length(rvals)
chtData <- array(dim=c(no_bins_x,no_bins_y,noRs))
chtData_pi <- array(dim=c(no_bins_x,no_bins_y,noRs))
# XXXX Up to here. Next - make the modified part E
for (i in 1:no_bins_x) {
	for (j in 1:no_bins_y) {
		# generate solution
		for (k in 1:noRs) {
			sol_tmp <- (SirModel3AgeClasses(
							pname=c("r","N","phi_1","mixsense","fitsus"),
							pvals=c(rvals[k],popsize,y_mids[j],x_mids[i],0),
							casemix=c(20,40,40),
							vp=ill_params))$sol
			# q[i,j,k] <- sum(sol_tmp[,"dS"])/popsize
			chtData[i,j,k] <- illRateICU*max(sol_tmp[,"dS"])/popsize*100000
		}
		cat(paste("Completed",j + (i-1)*no_bins_y," of ",no_bins_x*no_bins_y,"\n"))
		flush.console()
	}
}

# Generate the latin hypercube samples for the scatter plot
rep_bnds <- c(1,5)
mix_bnds <- c(0,5)
no_samples <- 10
set.seed(1234)
hv_rep <- hyper_vector(no_samples,rep_bnds[1],rep_bnds[2])
hv_mix <- hyper_vector(no_samples,mix_bnds[1],mix_bnds[2])
out_ar <- matrix(nrow=no_samples,ncol=noRs)
out_pi <- matrix(nrow=no_samples,ncol=noRs)
out_phi <- matrix(nrow=no_samples,ncol=noRs)

# Set up the output variables
for (i in 1:no_samples) {
	for (j in 1:noRs) {
		sol_tmp <- SirModel3AgeClasses(
				pname=c("r","N","phi_1","mixsense","fitsus"),
				pvals=c(rvals[j],popsize,1,hv_mix[i],1),
				casemix=c(20*hv_rep[i],40,40),
				vp=ill_params)
		out_ar[i,j] <- sum((sol_tmp$sol)[,"dS"])/popsize
		out_pi[i,j] <- illRateICU*max((sol_tmp$sol)[,"dS"])/popsize
		out_phi[i,j] <- sol_tmp$par["phi_1"]
	}
	cat(paste("Completed",i," of ",no_samples,"\n"))
	flush.console()
}

apa.gen.plot.rosens <- function(
		nsamps = 100,
		R0s = c(1.1,1.4,1.8),
		rangeSus = c(0.01,1),
		logSus = TRUE,
		rangeDelat = c(0.1,10),
		logDelta=TRUE,
		outputFile="~/Dropbox/") {
	
	# Setup data structures for the plotting
	noR0s <- length(R0s)
	params <- array()
	
	# Run the models to calculate peak attack rate and over-representation
	
	# PLot to a pdf
} 