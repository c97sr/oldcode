rm(list=ls(all=TRUE))
options(error=recover)

source("./stevensRfunctions.R")

# Runs for the baseline figure

popsize <- 8274527	# new york from http://www.census.gov/popest/cities/tables/SUB-EST2007-01.csv
# popsize <- 300000000	# new york from http://www.census.gov/popest/cities/tables/SUB-EST2007-01.csv

# Next check that the revised force of infection doesn't make much difference
st 		<- modelSIR(pname=c("r","N","seed"),pvals=c(0.26,1e9,3))

stn 	<- (SirModel3AgeClasses(pname=c("r","N","seed"),pvals=c(+0.00001,1e9,3),
						casemix=c(1492,724,71),vp=usParams(),NGM=MakeNgm))$sol

	# plot(stn[,"dS"],ylim=c(0.9,1.1),type="l")

max(st[,"Iv"]) 
max(stn[,"Iv"]) # 0.02419522

sol_low_Tg 	<- modelSIR(			pname=c("r","p_H_base","p_R","Tg","N"),
									pvals=c(0.26,0.031,0.86,1.8,popsize))
							
sol_high_Tg <- modelSIR(			pname=c("r","p_H_base","p_R","Tg","N"),
									pvals=c(0.26,0.031,0.86,3.5,popsize))

sol_low_Tg_age 			<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.1,0.031,0.86,1.8,popsize),casemix=c(1492,724,71)))$sol														
							
sol_low_Tg_low_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.125,0.031,0.86,1.8,popsize),casemix=c(1492,724,71)))$sol														

sol_low_Tg_high_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.26,0.031,0.86,1.8,popsize),casemix=c(1492,724,71)))$sol			
									
sol_high_Tg_low_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.125,0.031,0.86,3.5,popsize),casemix=c(1492,724,71)))$sol													

sol_high_Tg_high_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.26,0.031,0.86,3.5,popsize),casemix=c(1492,724,71)))$sol												
																		
sol_med_Tg_low_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.125,0.031,0.86,2.7,popsize),casemix=c(1492,724,71)))$sol														
									
sol_med_Tg_high_r_age 	<- (SirModel3AgeClasses(	pname=c("r","p_H_base","p_R","Tg","N"),
											pvals=c(0.26,0.031,0.86,2.7,popsize),casemix=c(1492,724,71)))$sol		
									
ps <- usParams()
multicu <- c(	ps["p_R"]*ps["p_H_base_1"]*ps["p_I_1"],
				ps["p_R"]*ps["p_H_base_2"]*ps["p_I_2"],
				ps["p_R"]*ps["p_H_base_3"]*ps["p_I_3"])

multhosp <- c(	ps["p_R"]*ps["p_H_base_1"],
				ps["p_R"]*ps["p_H_base_2"],
				ps["p_R"]*ps["p_H_base_3"])
		
hosprate_high <- agerate(sol_med_Tg_high_r_age,mult=multhosp)
hosprate_low <- agerate(sol_med_Tg_low_r_age,mult=multhosp)
icurate_high <- agerate(sol_med_Tg_high_r_age,mult=multicu)
icurate_low <- agerate(sol_med_Tg_low_r_age,mult=multicu)
																		
# Runs for the sensitivity analysis figure
x_lims 		<- c(0.05,0.30)									# r
y_lims 		<- c(log(0.005,base=10),log(0.05,base=10))		# p_R * p_I * p_H_base
no_bins		<- 2											# often either 2 or 10
x_size 		<- (x_lims[2]-x_lims[1])/no_bins
y_size		<- (y_lims[2]-y_lims[1])/no_bins
x_bounds	<- x_lims[1]+(0:no_bins) * x_size
y_bounds	<- y_lims[1]+(0:no_bins) * y_size
x_mids		<- x_lims[1] + (1:no_bins) * x_size - x_size / 2
y_mids		<- y_lims[1] + (1:no_bins) * y_size - y_size / 2

# Make the array to be plotted
chtData <- matrix(nrow=no_bins,ncol=no_bins)
for (i in 1:no_bins) {
	for (j in 1:no_bins) {
		# generate solution
		sol <- (SirModel3AgeClasses(
					pname=c("r","p_H_base_1","p_H_base_2","p_H_base_3","p_R","p_I","Tg","N"),
					pvals=c(x_mids[j],0.02489*10^(y_mids[i])/0.031,0.03575*10^(y_mids[i])/0.031,0.09338*10^(y_mids[i])/0.031,1,0.13,2.7,popsize),
					casemix=c(1492,724,71)))$sol
		chtData[i,j] <- max(sol[,"Iv"])/popsize * 100000
		cat(paste("Completed",j + (i-1)*no_bins," of ",no_bins*no_bins,"\n"))
		flush.console()
	}
}

icu_cap <- 59162/296410404*100000
