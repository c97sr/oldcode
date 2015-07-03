read_study_data <- function(fileTests,fileHouseChar,fileIndChar) {

	# Takes the filename of the laboratory tests, household characteristics
	# and individual characteristics, loads the files, and returns a table
	# of individual charateristics and household characteristics
	
	# Read in the three files
	hc <- read.csv(fileTests,header=TRUE)
	hh <- read.csv(fileHouseChar,header=TRUE)
	ic <- read.csv(fileIndChar,header=TRUE)
	
	# Define household IDs and household member indices, check for consistancy
	familysize <- hh$familysize
	hhIDseq <- as.character(rep(hh$hhID,familysize))
	intNumberHouses <- dim(hh)[1]
	hhIndexSeq <- as.numeric(rep(1:intNumberHouses,familysize))
	noPeople <- length(hhIDseq)
	noSamples <- dim(hc)[1]
	if (noPeople != dim(ic)[1]) stop("mismatch between hh data and individual data")
	
	# Define the main output variable
	hculture <- data.frame(	hhID = hhIDseq,
							hhIndex = hhIndexSeq,
							member = rep(0,noPeople),
							hh_size = rep(familysize, familysize), #rep(999,noPeople),
							index_first = rep(cumsum(familysize) - familysize + 1, familysize),	#rep(999,noPeople),
							int_g = rep(hh$intervention,hh$familysize),
							age = ic$age,
							t_S = ic$clinical_onset,
							T0 = rep(NA,noPeople),
							T1 = rep(NA,noPeople),
							T2 = rep(NA,noPeople),
							T3 = rep(NA,noPeople),
							T4 = rep(NA,noPeople),
							V0 = rep(NA,noPeople),
							V1 = rep(NA,noPeople),
							V2 = rep(NA,noPeople),
							V3 = rep(NA,noPeople),
							V4 = rep(NA,noPeople),
							i_g = rep(NA,noPeople),			# 0 sus never infected, 1 inf, 2 immune
							t_I = rep(-2,noPeople))
	
	# Give indiviuals the right index within the household				
	for(i in 2:noPeople){
		if(hculture$hhID[i]==hculture$hhID[i-1]) hculture$member[i] <- hculture$member[i-1]+1
	}

	# Make sure we can find the start of each household from any hh member
	hh$index_first <- cumsum(familysize) - familysize + 1
	
	# Give everyone the same household visit days in this study
	hculture$T0 <- rep(rep(0,length(familysize)), familysize) 			# so that day=0 is clinic day not symptom onset day
    hculture$T1 <- rep(hh$v1_day, familysize)	#(lab data uses *day.from.clinic*)
	hculture$T2 <- rep(hh$v2_day, familysize)	#(lab data uses *day.from.clinic*)
	hculture$T3 <- rep(hh$v3_day, familysize)	#(lab data uses *day.from.clinic*)
	hculture$T4 <- rep(hh$v4_day, familysize)	#(lab data uses *day.from.clinic*)
	
	# Put the correct test results into the hculture object, using R's logical dataframe 
	# referencing
	for (i in 1:noSamples) {
		hhIDtmp <- as.character(hc$hhID[i])
		membertmp <- hc$member[i]
		visittmp <- hc$visit[i]
		culturetmp <- hc$culture[i]
		hculture[hculture$hhID==hhIDtmp & hculture$member==membertmp, paste("V",visittmp,sep="")] <- 
				as.character(culturetmp)
	}
	
	# Return individual and household tables
	list(indDb=hculture,hhDb=hh)

}

calc_ln <-function(rps,bps,np,nt,vecT,vecI,matT,matV,ft,hhlookup,t_factor=1,v_factor=1,s_factor=1) {
	
	# Take the run parameters, biological parameters, number of tests, vector of times, matrix of test
	# times, matrix of test results, 

	nohh <- dim(hhlookup)[1]
	rtn <- as.vector(rep(0,nohh),mode="numeric")
	for (i in 1:nohh) {
		rtn[i] <- calc_ln_aux(	rps,bps,hhlookup[i,"start"],hhlookup[i,"end"],
								nt,vecT,vecI,matT,matV,ft,
								t_factor=t_factor,v_factor=v_factor,s_factor=s_factor	)
	}
	rtn
}

calc_ln_hh <-function(hh,rps,bps,np,nt,vecT,vecI,matT,matV,ft,hhlookup,t_factor=1,v_factor=1,s_factor=1) {
	rtn <- calc_ln_aux(rps,bps,hhlookup[hh,"start"],hhlookup[hh,"end"],
				nt,vecT,vecI,matT,matV,ft,
				t_factor=t_factor,v_factor=v_factor,s_factor=s_factor	)
	rtn
}

calc_ln_aux <-function(rps,bps,indstart,indfinish,nt,vecT,vecI,matT,matV,ft,t_factor=t_factor,v_factor=v_factor,s_factor=s_factor) {
	i <- 1
	lnlike <- 0
	for (i in indstart:indfinish) {
		lnlike <- 	lnlike +
					t_factor * calc_ind_l_T(i,rps,bps,vecT,vecI,ft) + 
					v_factor * calc_ind_l_V_cult(i,rps,bps,nt,vecT,vecI,matT,matV,ft)  +
					s_factor * calc_ind_l_S(i,rps,bps,vecT,vecI,ft)
	}
	lnlike
}

calc_ind_l_S <- function(i,rps,bps,vecT,vecI,ft) {

	# Calculate the log likelihood for the occurrence of symptoms 
	# given times of infection and infection status.
	# i is index of individual
	# rps are run parameters
	# bps are biological parameters
	# vecT is a vector of times of infection
	# ft is a table of indiviual level properties

	rtn <- 0 # Initialize return value
	p1_symp <- bps["p1_symp"] # scale parameter
	p2_symp <- bps["p2_symp"] # shape parameter
	t_S <- ft$t_S[i] # time of symtoms
	t_I <- vecT[i] # time of infection
	i_g <- vecI[i] # infection group
	t_max <- rps["t_max"] # maximum time considered
	t_end_S <- min(t_S,t_max)

	# Index cases 
	if (ft$member[i]==0) {
		
		# Log likelihood of symptoms occurring at time t_end_S given infection at t_I
		lnl_dense_index <- fnDIncubation(t_end_S-t_I,p1_symp,p2_symp,log=TRUE) 
		
		# Log likelihood of symptoms not appearing between t_I and t_end 
		lnl_inthaz_index <- -fnPIncubation(t_end_S-t_I,p1_symp,p2_symp) 
		
		rtn <- lnl_dense_index + lnl_inthaz_index

	# Non-index cases
	} else {
		
		# Likelihood (not log) of survival from infection due only to background symptom risk
		l_inthaz_bck_only <- exp(-bps["back_symp"]*(t_end_S-rps["t_min"]))
		
		# Likelihood (not log) of symptoms due only to background symptom risk at a particular time
		l_dense_bck_only <- bps["back_symp"] 
		
		# If the time of symptoms is prior to the time of infection, but the time of infection is 
		# less than the end of the observation window (i.e. the person was infected)
		if (t_S < t_I && t_I < t_max) {
			
			# Return the log of the product of survival and infection
			rtn <- log(l_inthaz_bck_only*l_dense_bck_only)
		
		# If the time of symptoms is greater than the time of infection or the time of infection
		# is greater than the max time?
		} else {
			
			# Likelihood (not log) of survival from symptoms due to infection only
			l_inthaz_inf_only <- exp(-fnPIncubation(t_end_S-t_I,p1_symp,p2_symp))
			
			# Likelihood (not log) of survival due to background and infection
			l_inthaz_both <-  l_inthaz_inf_only * l_inthaz_bck_only
			
			# Likelihood (not log) of symptoms at a particular time due to being infected
			l_dens_inf <- fnDIncubation(t_end_S-t_I,p1_symp,p2_symp,log=FALSE)
			
			# If symtoms come after infection but before the max time
			if (t_S >= t_I && t_S < t_max) {
				
				# If the individual is infected
				if (i_g == 1) {
					
					# The exponent of the return value is assumed to be due to symptoms if the 
					# infection caused symptoms, but could still have happened if the infection
					# did not cause symptoms
					exprtn <- (bps["p_symp"])*(l_dens_inf+l_dense_bck_only)*l_inthaz_both +
										(1-bps["p_symp"])*l_dense_bck_only*l_inthaz_bck_only			
					
				} else {
					
					# If the person was not infected, then symptoms can only have occurred due to 
					# the background presence of symptoms
					exprtn <- l_dense_bck_only*l_inthaz_bck_only
					
				}
			
			# If symtoms came after the max time
			} else {
				
				# If the individual was infected 
				if (i_g == 1) {
					exprtn <- (bps["p_symp"])*l_inthaz_both +
							(1-bps["p_symp"])*l_inthaz_bck_only
				
				# Otherwise, the individual avoided only background symptoms
				} else {
					exprtn <- l_inthaz_bck_only
				} 
				
			}
			
			rtn <- log(exprtn)
			
		}
		
	}
	
	rtn
	
}

gen_symp <- function(ps,rps,t_i,i_g) {
	
	# Generates a time of symptoms for an individual
	# Takes biological and run parameters, time of infection and infection group as arguments
	# Returns the simulated time of symptoms
	
	rtn <- 999
	
	# For those individuals infected during the study
	if (i_g==1) {
		
		# Test to see if they have a syptomatic infection
		if (runif(1) > ps["p_symp"]) rtn <- 999
		
		# If they do, draw an infection time from the appropriate distribution
		else rtn <- t_i + fnRIncubation(1,ps["p1_symp"],ps["p2_symp"])
		
	}
	
	# Set t_end to either the last observed time or the time of symptoms
	t_end <- min(rtn,rps["t_max"])
	
	# Test to see if symptoms occur spontaneously prior to infection generated symptoms
	if (runif(1) < (1-exp(-ps["back_symp"]*(t_end-rps["t_min"])))) {
		
		# If they do, choose the time of infection
		rtn <- runif(1,min=rps["t_min"],max=t_end)
		
	}
	
	rtn
	
}

calc_ind_l_T <- function(i,rps,bps,vecT,vecI,ft) {														
	
	# ith ind, tab culture table, biological parameters
	rtn <- 0
	if (vecI[i]!=2) {
		if (vecI[i]==0) t_e <- rps["t_max"]
		else t_e <- vecT[i]
		set_others <- setdiff(ft$index_first[i]:(ft$index_first[i]+ft$hh_size[i]-1),i)
		rtn <- rtn - (t_e - rps["t_min"]) * bps["h_B"]
		for (j in set_others) rtn <- rtn - int_ptp_lambda(t_e-vecT[j],0,bps,i,ft)
		if (vecI[i]==1) {
			h_inf <- bps["h_B"]
			for (j in set_others) h_inf <- h_inf + dens_ptp_lambda(t_e-vecT[j],bps,i,ft)
			rtn <- rtn + log(h_inf)
		}
	}
	rtn
}

calc_ind_l_V_cult <- function(i,rps,bps,nt,vecT,vecI,matT,matV,ft) {
	
	# Calculate the likelihood of an individual viral culture result
	# Takes index, run parameters, biological parameters, number of tests, vector of times of infection,
	# vector of infection states (0 - never, 1 - infected during study, 2 - immune at outset), a matrix
	# of test times and a matrix of test values and a table of individual characteristics
	
	# Initialize the return variable
	rtn <- 0
	
	# Loop through the number of tests
	for (j in 1:nt) {
		
		# Current result is the jth test of individual i
		res <- matV[i,j]
		
		# If the result is a non-negative number
		if (res >= 0) {
			
			# If the individual was never infected or immune to start with ...
			if (vecI[i] == 0 || vecI[i] == 2) {
				
				# Likelihood of a positive is one minus the specificity
				if (res == 1) rtn <- rtn + log(1-bps["spec"])
				
				# Likelihood of a negative is the specificity
				else rtn <- rtn + log(bps["spec"]) 
			
			# ... otherwise, the individual was infected at some point during the study
			} else {
				
				# The time of the jth test of the ith individual
				t <- matT[i,j]
				
				# The time since infection at the jth test
				tau <- t - vecT[i]
				
				# The hazard if infection at the time of the test
				haz <- dens_ptp_lambda(tau,bps,i,ft)
				
				# If the hazard is greater than the detectable threshold...
				if (haz > bps["h_detect"]) {
					
					# Likelihood of positive is the sensitivity of the test
					if (res==1) rtn <- rtn + log(bps["sens"]) 
					
					# Likelihood of a negative is one minus the sensitivity of the test
					else rtn <- rtn + log(1-bps["sens"]) 
				
				# ... otherwise, the hazard is less than the detectable threshold
				} else {
					
					# Likelihood of a positive is one minus the specificity
					if (res==1) rtn <- rtn + log(1-bps["spec"])
					
					# Likelihood of a negative is the specificity
					else rtn <- rtn + log(bps["spec"]) 			
					
				}
			}
		}
	}
	
	rtn
	
}

dens_ptp_lambda_aux <- function(t,ps,int=FALSE) {
	rtn <- t*0
	if (int) ptp <- ps["q_int_1"]
	else ptp <- ps["q_base"] 
	if (t >= ps["t_off_ptp"] ) {
		rtn <- rtn + log(1/(1-ptp)) * 	
				(fnDInfection(t-ps["t_off_ptp"],ps["mu_ptp"]/ps["alpha_ptp"],ps["alpha_ptp"],log=FALSE))				
	}
	rtn
}

dens_ptp_lambda <- function(tau,ps,i,ft) {
	
	# Give the hazard of person to person tranmsission at time tau for the ith person in table
	# ft with parameters ps
	
	# If the person is in a no intervention group thene the time of intervention is set to a very large number
	if (ft$int_g[i]==0) tau_int <- 9999	# No intervention must be a zero code
	
	# Otherwise it is set to minus the time of infection. Remember that chronological time zero is the 
	# time that the interventions are implemented
	else tau_int <- -ft$t_I[i]
	
	# If the intervention has already started, call the aux function with intervention turned on
	# otherwise call the aux function without the intervention turned on.
	if (tau_int < tau) rtn <- dens_ptp_lambda_aux(tau,ps,int=TRUE)
	else rtn <- dens_ptp_lambda_aux(tau,ps,int=FALSE)
	
	rtn
}

int_ptp_lambda_aux <- function(t,ts,ps,int=FALSE,intgroup=0) {
	
	# A auxilliary function for int_ptp_lamnda
	
	t <- t
	ts <- ts
	rtn <- 0
	
	if (int) ptp <- ps["q_int_1"]
	else ptp <- ps["q_base"] 

	# if (intgroup == 0) ptp <- ps["q_base"]
	# else if (intgroup == 1) ptp <- ps["q_int_1"]
	# else if (intgroup == 2) ptp <- ps["q_int_2"]

	if (t >= ps["t_off_ptp"] && ts >= ps["t_off_ptp"]) {
		rtn <- rtn + 1-exp(-log(1/(1-ptp))*(
							fnPInfection(t-ps["t_off_ptp"],ps["mu_ptp"]/ps["alpha_ptp"],ps["alpha_ptp"]) -
							fnPInfection(ts-ps["t_off_ptp"],ps["mu_ptp"]/ps["alpha_ptp"],ps["alpha_ptp"]) ))}

	else if (t >= ps["t_off_ptp"]) {
		rtn <- rtn + 1-exp(-log(1/(1-ptp))*(	
							fnPInfection(t-ps["t_off_ptp"],ps["mu_ptp"]/ps["alpha_ptp"],ps["alpha_ptp"]) ))	}
	
	rtn

}

int_ptp_lambda <- function(tau,ts,ps,i,ft) {
	
	# Integrated person to person hazard
	# Takes tau time since infection, ts start time tau since infection, ps biological parameters
	# i index of the individual and ft table of features
	# Returns the integrated hazard of person to person infection
	
	# If the person is not in an intervention group, the time since infection
	ig <- ft$int_g[i]
	if (ig==0) t_int <- 9999
	
	# If they are in an intervention group, the time of the intervention is equal to minus 1 times t_I
	else t_int <- -ft$t_I[i]
	
	# If the tau time of intervention is greater than the end tau time, then
	# return the full infectiousness regardless of infection group
	if (t_int > tau) rtn <- int_ptp_lambda_aux(tau,ts,ps,int=FALSE, intgroup = 0)
	
	# Else if the time of intervention is less than the start tau time time, then
	# apply the intervention reduction to the whole period
	else if (t_int < ts) rtn <- int_ptp_lambda_aux(tau,ts,ps,int=TRUE, intgroup = ig)
	
	# Else apply the intervention at some point during the time
	else {
		rtn <- 	int_ptp_lambda_aux(t_int,ts,ps,int=FALSE, intgroup = 0) +
				int_ptp_lambda_aux(tau,t_int,ps,int=TRUE, intgroup = ig)
	}
	
	rtn
	
}

opt_int_ptp_lambda <- function(t,p,ts,ps,i,ft) {
	
	# An auxilliary function using the integrated person to person hazard function (intptpt).
	# Returns the difference between the intptp argument p the target probability.
	# Used for during lazy simulation of times to next event.
	
	rtn <- p-int_ptp_lambda(t,ts,ps,i,ft)
	
	rtn
	
}

fnPropTSUpdates <- function(rps,np,vecT,vecI,vecTM) {

	# Takes run parameters, number of people, vector of current infection
	# times, vector of current infection states, and a vector of max infection
	# times and returns a vector or proposed infection times
	# Three different proposal steps are implemented below
	# Cyclical random walk is probably the best solution overall

	rtn <- vecT
	
	for (i in 1:np) {
	
		if (vecI[i] == 1) {
		
			# Naive uniform timestep proposal from entire range
			# rtn[i] <- rps["t_min"] + (vecTM[i]-rps["t_min"]) * runif(1)
			
			# Random walk update step
			rtn[i] <- rtn[i] + (runif(1)-0.5)*rps["t_U_max"]*2
			
			# Bouncing boundary conditions for random walk => non-uniform prior
			# if (rtn[i] < rps["t_min"]) rtn[i] <- rps["t_min"] + (rps["t_min"] - rtn[i])  
			# if (rtn[i] > vecTM[i]) rtn[i] <- vecTM[i] - (rtn[i] - vecTM[i])  
			
			# Cyclical boundary conditions for random walk => uniform prior
			if (rtn[i] < rps["t_min"]) rtn[i] <- vecTM[i] - (rps["t_min"] - rtn[i])
			if (rtn[i] > vecTM[i]) rtn[i] <- rps["t_min"] + (rtn[i] - vecTM[i])

		}

	}

	rtn
	
}

fnPropTSUpdatesSingle <- function(rps,hhlu,nohh,vecT,vecI,vecTM) {
	
	# Takes run parameters, hhlookup value, number of households, vector of current infection
	# times, vector of current infection states, and a vector of max infection
	# times and returns the index of the household, the index of the indiviuals and the new time
	
	# Changes may be needed here XXXX - currently assumes that at least one individual in
	# the household is infectious. Should probably stick with that for now.
	# For expediency, should probably just delete the non-infectious index cases
	
	# Choose household and setup aux variables
	intChosenHouse <- as.numeric(floor(runif(1)*nohh + 1))
	szHouse <- hhlu[intChosenHouse,"size"]
	stHouse <- hhlu[intChosenHouse,"start"]
	# debug_check <- 0 # Used to escape livelock in while loop below
	
	# Make a list of the infectious indiviuals in the household
	intNumberInfectious <- 0
	for (i in 1:szHouse) {
		if (vecI[stHouse+i-1]==1) {
			intNumberInfectious <- intNumberInfectious + 1 
		}
	}	
	
	# Choose an infectious indiviual at random using only one random number
	intChosenFromInfectious <- floor(runif(1)*intNumberInfectious)+1
	boolNotYetChosen <- TRUE
	intChosenPerson <- stHouse
	intLoopIndex <- 0
	while (boolNotYetChosen) {
		if (vecI[stHouse+intLoopIndex] == 1) {
			if (intChosenFromInfectious == 1) {
				intChosenPerson <- stHouse + intLoopIndex
				boolNotYetChosen <- FALSE
			} else {
				intChosenFromInfectious <- intChosenFromInfectious - 1
			}
		}
		intLoopIndex <- intLoopIndex + 1
	}
		
	# Random walk update step
	newT <- vecT[intChosenPerson]
	newT <- newT + (runif(1)-0.5)*rps["t_U_max"]*2
		
	# Cyclical boundary conditions for random walk => uniform prior
	if (newT < rps["t_min"]) newT <- vecTM[intChosenPerson] - (rps["t_min"] - newT)
	if (newT > vecTM[intChosenPerson]) newT <- rps["t_min"] + (newT - vecTM[intChosenPerson])
	
	# Return a list, maybe change if slow
	list(person=intChosenPerson,house=intChosenHouse,time=newT)
	
}

fnPropImmuneUpdatesSingle <- function(nononind,mapnonind,vecI) {
	
	# Propose an update to a single non-index individual
	# Takes the number of non indexes, a map of non indexes 
	# a list of the current immune states and a list of attributes as input
	# Returns the new immune state, the index of the person
	# and the household of the person
	
	toUpdate <- floor(runif(1)*nononind)+1
	index <- mapnonind[toUpdate]
	curState <- vecI[index]
	if (curState==0) rtnstate <- 1
	else rtnstate <- 0
	
	list(person=index,state=rtnstate)
	
}

fnResT <- function(str) {
	rtn <- FALSE
	if (!is.na(str)) if (str=="A" || str=="B") rtn <- TRUE
	rtn
}

MH_sym_ap <- function(delta) {
	if (delta > 0) rtn <- TRUE
	else if (delta < -150) rtn <- FALSE
	else {
		test <- runif(1)
		p <- exp(delta)
		if (p > test) rtn <- TRUE
		else rtn <- FALSE
	}
	rtn
}

fnSetupOutputTab <- function(bp, fmask, tis, igs, nss, fno) {
	
	# parameters modified
	# time steps recorded
	# immune states recorded
	# number of samples stored
	# file name for output
	
	no_ps <- length(fmask)
	no_ts <- length(tis)
	no_gs <- length(igs)
	
	ot <- data.frame(
			index = rep(NA,nss),
			lnl = rep(NA,nss),
			noInf = rep(NA,nss),
			noSus = rep(NA,nss),
			noImm = rep(NA,nss),
			pAcT = rep(NA,nss),
			pAcP = rep(NA,nss),
			pAcI = rep(NA,nss)
			)
		
	offset <- dim(ot)[2]+1
	
	for (i in 1:no_ps) ot <- cbind(ot,rep(NA,nss))
	names(ot)[offset:(offset+no_ps-1)] <- names(bp)[fmask]
	
	for (i in 1:no_ts) ot <- cbind(ot,rep(NA,nss))
	names(ot)[(offset+no_ps):(offset+no_ps+no_ts-1)] <- paste("t_",1:no_ts,sep="")
	
	for (i in 1:no_gs) ot <- cbind(ot,rep(NA,nss))
	names(ot)[(offset+no_ps+no_ts):(offset+no_ps+no_ts+no_gs-1)] <- paste("i_",1:no_gs,sep="")
	
	write(names(ot),file=fno,sep="\t",ncolumns=length(ot))
	
	ot
	
}

fnCountStates <- function(vecI,size) {
	i <- 1
	rtn = vector("numeric",3)
	while (i <= size) {
		x <- vecI[i]
		rtn[x+1] <- rtn[x+1] + 1
		i <- i+1
	}
	rtn
}

fnExtract <- function(otr,ind,lnl,bp,nf,tss,igs,vecT,vecI,vecAccPs) {

	# Writes selected output to the global table output_tab
	# Returns nothing.
	# Should be rewritted to return a vector of numeric to go into the 
	# table at a higher level of scope.
	
	no_ts <- length(tss)
	no_is <- length(igs)
	no_psm <- length(nf)
	no_people <- length(vecI)
	veccounts <- fnCountStates(vecI,no_people)
	output_tab$index[otr] <<- ind
	output_tab$lnl[otr] <<- lnl	
	output_tab$noInf[otr] <<- veccounts[2]
	output_tab$noSus[otr] <<- veccounts[1]
	output_tab$noImm[otr] <<- veccounts[3]
	output_tab$pAcT[otr] <<- sum(vecAccPs[,1])
	output_tab$pAcP[otr] <<- sum(vecAccPs[,2])
	output_tab$pAcI[otr] <<- sum(vecAccPs[,3])
	
	for (p in nf) {
		output_tab[otr,p] <<- bp[p]
	}
	
	for (i in 1:no_ts) {
		output_tab[otr,8+no_psm+i] <<- vecT[tss[i]]
	}
	
	for (i in 1:no_is) {
		output_tab[otr,8+no_psm+no_ts+i] <<- vecI[igs[i]]
	}

}

fnProposeParamUpdates <- function(bps,fmask,ptab,nofit) {

	# Takes biological parameters, a mask which is a list of the parameters
	# being fitted, a table of the fitted parameters and their ranges and
	# the number of parameters to be fitted. 
	# Returns a proposed vector of parameters.
	# Comment in and out bouncing and cyclical boundary conditions
	# for the random walk
	
	rtn <- bps
	
	for (i in 1:nofit) {
	
		# Set up and transform to unit scale
		rv <- runif(1)
		rv <- (rv-0.5)* ptab[i,3]
		x <- bps[fmask[i]]
		x <- SR_to_unit(x,min=ptab[i,1],max=ptab[i,2],logflag=ptab[i,4])
		x <- x + rv
		
		# Bouncing boundary conditons
		# if (x < 0) x <- -x	
		# if (x > 1) x <- 2 - x
		
		# Cyclical boundary conditions
		if (x < 0) x <- 1 + x	
		if (x > 1) x <- x - 1
		
		# Test for errors and return to originl scales
		if (x < 0 || x > 1) stop("problem here")		
		rtn[fmask[i]] <- SR_from_unit(x,min=ptab[i,1],max=ptab[i,2],logflag=ptab[i,4])

	}

	rtn
}

fnStdOutputAnalysis <- function(od,ps_to_plot = c("lnl"),plot=TRUE,burnin=0,burninprop=TRUE) {
	
	ns <- dim(od)[1]
	nops <- length(ps_to_plot)
	if (burninprop) burnin <- round(ns*burnin)
	
	if (plot) {
		wlist <- dev.list()
		for (w in wlist) dev.off(w)
		for (i in 1:nops) {
			
			# Traces
			windows(width=20/2.5,height=10/2.5)
			plot(	od[(burnin+1):ns,1],
					od[(burnin+1):ns,ps_to_plot[i]],
					ylab=ps_to_plot[i],type="l")
			
			# Histogram
			windows(width=10/2.5,height=10/2.5)
			hist(od[(burnin+1):ns,ps_to_plot[i]],main=ps_to_plot[i])
			
		}
	}

	tmpvec <- rep(-999,nops)
	crint <- data.frame(mean=tmpvec,median=tmpvec,lb=tmpvec,ub=tmpvec)

	for (i in 1:nops) {
		curpar <- ps_to_plot[i]
		crint$mean[i] <- mean(od[(burnin+1):ns,curpar])
		crint$median[i] <- quantile(od[(burnin+1):ns,curpar],c(0.5),type=7)
		crint$lb[i] <- quantile(od[(burnin+1):ns,curpar],c(0.025),type=1)
		crint$ub[i] <- quantile(od[(burnin+1):ns,curpar],c(0.975),type=1)
		row.names(crint)[i] <- curpar
	}
	
	list(vals=od,sumtab=crint)

}

fnLotsOfTiPlots <- function(file=fnOutput, ts_to_plot=1:20, burnin=0.1,	burninprop=TRUE) {
	# function to plot a grid of histograms of samples from each T_i
	#
	od <- read.table(file,header=TRUE)
	ns <- dim(od)[1]
	if (burninprop) burnin <- round(ns*burnin)
	
	windows(width=8, height=7)
	layout(matrix(ts_to_plot, nrow=round(length(ts_to_plot)/4,0), byrow=TRUE))
	par(mar=c(4,3,1,1), las=1)
	for(i in ts_to_plot) hist(od[(burnin+1):ns,paste("t_", i, sep="")], main=paste("t_", i, sep=""), xlab="")
	return()
}

fnConvCulture <- function(x) {
	if (is.na(x) || x==-999) rtn <- -999
	else if (x=="A" || x=="B" || x==1) rtn <- 1
	else if (x=="0" || x==0) rtn <- 0
	else stop("Problem with the conversion in tmpFn.")
	rtn
}

fnMakeHHLookup <- function(hhdb) {
	noHH <- dim(hhdb)[1]
	rtn <- as.matrix(data.frame(start=rep(-1,noHH),end=rep(-1,noHH),size=rep(-1,noHH)),mode="numeric")
	currentrow <- 1
	for (i in 1:noHH) {
		rtn[i,"start"] <- currentrow
		sizetmp <- hhdb[i,"familysize"]
		rtn[i,"end"] <- currentrow + sizetmp - 1
		rtn[i,"size"] <- sizetmp
		currentrow <- currentrow + sizetmp
	}
	rtn
}

fnRunMcmc <- function(	bp,rp,p_table,indDb,hhlookup,
						fnOutputStem = "mcmc",
						fnSnapshot = "NULL",
						t_I_s_stored = 1:2,
						i_g_s_stored = 1:2,
						db_prop = FALSE,
						db_like = FALSE,
						db_freq=1) {			
				
	# Function to take samples from the markov chain and output them to a file
	# Takes bio params, run params, a table of fitted bio params, a table of individual
	# level data, a household lookup table, a filename stem for the output, a filename for the 
	# saved state of the chain
	# if it is to be started from a previous run, a list of the indexes of individuals 
	# who are to hhave their times of infections stored, similar for infection statuses and a 
	# flag for whether to do an update debug run. In update debug runs, all proposed updates are
	# accepted, regardless of likelihood.
	# Returns no values, but does output its state as fnOutputStem+"*snapshot*" which can be
	# used as an argument to start the chain again if required.
					
	# Assigne all the filenames
	fnOutputTrace 		<- paste(fnOutputStem,"_trace.out",sep="")
	fnOutputTable 		<- paste(fnOutputStem,"_table.out",sep="")
	fnOutputBiops 		<- paste(fnOutputStem,"_biops.out",sep="")
	fnOutputPsrange 	<- paste(fnOutputStem,"_psrange.out",sep="")
	fnOutputRunps 		<- paste(fnOutputStem,"_runps.out",sep="")
	fnOutputSnapshot	<- paste(fnOutputStem,"_snapshot.out",sep="")
	fnOutputSnapshotDB	<- paste(fnOutputStem,"_db_snapshot.out",sep="")
	
	# Setup vectors from the input dataframe if a new run
	# Load a previous napshot if not a new run.
	if (fnSnapshot == "NULL") {
		noTests <- 5
		matT <- data.matrix(indDb[,c("T0","T1","T2","T3","T4")])
		tmpMatV <- as.matrix(indDb[,c("V0","V1","V2","V3","V4")])
		matV <- apply(tmpMatV,c(1,2),fnConvCulture)
		vecT <- as.vector(indDb$t_I,mode="numeric")
		vecI <- as.vector(indDb$i_g,mode="numeric")
		noPeople <- dim(indDb)[1]
		# Set baseline parameter values
		loglike <- calc_ln(rp,bp,noPeople,noTests,vecT,vecI,matT,matV,indDb,hhlookup)
	} else {
		load(fnSnapshot)
	}
	
	# Make a list of the indices of fitted parameters and aux structures
	fmask <- match(row.names(p_table),names(bp))
	nofitted <- dim(p_table)[1]
	namesfitted <- names(bp)[fmask]
	
	# Set up acceptance rate calculators. p for probability matrix and c for the count vector
	# Also some sux variables
	p_1000 <- array(0,c(1000,3))
	c_1000 <- array(1,c(3))
	intNumberHouses <- dim(hhlookup)[1]
	ps_mod <- row.names(p_table)
	
	# Set up vector of max update times. Non-index cases will be different to index cases
	vec_t_I_max <- as.vector(indDb$t_S,mode="numeric")
	for (i in 1:noPeople) if (indDb$member[i] != 0) vec_t_I_max[i] <- rp["t_max"]
	
	# Set up the weight vectors for the updates steps
	# One is redundent, but its much easier to have three parameters
	tmpvec <- data.frame(p=rp["weightParams"],t=rp["weightTime"],i=rp["weightImmune"])	# Min length For timestep, parameters, immune status
	weights <- as.vector(tmpvec,mode="numeric")	# Min length For timestep, parameters, immune status
	names(weights) <- names(tmpvec)
	weights <- cumsum(weights)/sum(weights)
	
	# Set up to be able to update only the infection status of non-index cases
	nononind <- noPeople-intNumberHouses
	mapnonind <- array(-1,c(nononind))
	nonindcount <- 1	
	for (i in 1:noPeople) {
		if (indDb$member[i] != 0) {
			mapnonind[nonindcount] <- i
			nonindcount <- nonindcount + 1
		}
	}
	
	# Preconditions for the main loop
	index <- 0
	sampCurStored <- 0
	output_tab <<- fnSetupOutputTab(bp,fmask,t_I_s_stored,i_g_s_stored,rp["noSampStored"],fnOutputTrace)
	fnExtract(sampCurStored+1,index,sum(loglike),bp,namesfitted,t_I_s_stored,i_g_s_stored,vecT,vecI,p_1000)
	sampCurStored <- sampCurStored + 1 
	next_extract <- rp["sampFreq"]
	
	# Set a minimum index threshold for testing the likelihood calculation}
	next_db_thresh <- db_freq

	# Main loop for MCMC sampler
	while (index < rp["totSamples"]) {
		
		# Proposal type random number
		p_type <- runif(1)

		# Propose to update all parameters	
		if (p_type < weights["p"]) {
						
			propPs <- fnProposeParamUpdates(bp,fmask,p_table,nofitted)
			propLoglike <- calc_ln(rp,propPs,hhlookup,noTests,vecT,vecI,matT,matV,indDb,hhlookup)					
			x <- sum(propLoglike - loglike)
			
			if (db_prop || MH_sym_ap(x)) {
				bp <- propPs
				loglike <- propLoglike
				p_1000[c_1000[2],2] <- 1
			} else {
				p_1000[c_1000[2],2] <- 0
			}
			if (c_1000[2] == 1000) c_1000[2] <- 1
			else c_1000[2] <- c_1000[2] + 1
		
		# Propose to update all timesteps
		} else if (p_type < weights["t"]) {
			
			# Update all infection times at one go...
			# propVecT <- fnPropTSUpdates(rp,noPeople,vecT,vecI,vec_t_I_max)
			# propLoglike <- calc_ln(rp,bp,noPeople,noTests,propVecT,vecI,matT,matV,indDb,hhlookup)		
			# x <- sum(propLoglike - loglike)
		
			# or update only a single infection time
			propVecT <- vecT
			updates <- fnPropTSUpdatesSingle(rp,hhlookup,intNumberHouses,vecT,vecI,vec_t_I_max)
			propVecT[updates$person] <- updates$time
			propLoglike <- calc_ln_hh(updates$house,rp,bp,noPeople,noTests,propVecT,vecI,matT,matV,indDb,hhlookup)
			x <- propLoglike - loglike[updates$house]
			
			# If debug or accept prob, accept the update step
			if (db_prop || MH_sym_ap(x)) {
				
				# Update all infection times
				# vecT <- propVecT
				# loglike <- propLoglike
		
				# or update only a single infection time
				vecT[updates$person] <- updates$time
				loglike[updates$house] <- propLoglike
				
				# If debug like is true, and the debug frequency threshold has been passed
				# check for self consistency and dump state if problem
				if (db_like && index > next_db_thresh) {		
					
					# Increment the minimum debug threshold
					next_db_thresh <- index + db_freq
					
					# Test to see if the piecemeal likelihood is the same as the overall recalculated likelihood
					if (abs(sum(calc_ln(rp,bp,noPeople,noTests,vecT,vecI,matT,matV,indDb,hhlookup))-sum(loglike)) > 1e-6) {
						browser()
						save(	file=fnOutputSnapshotDB,
								noTests,
								matT,
								matV,
								vecT,
								vecI,
								bp,
								indDb,
								hhlookup,
								loglike,
								noPeople
						)
						stop("problem with the likelihood",index-1,p_type,sum(loglike),sum(calc_ln(rp,bp,noPeople,noTests,vecT,vecI,matT,matV,indDb,hhlookup)))
					}
				}
				
				# Increment the accept probabilities
				p_1000[c_1000[1],1] <- 1
				
			# Decline the update step
			} else {
				p_1000[c_1000[1],1] <- 0
			}
			
			# Manage the accept probabilities
			if (c_1000[1] == 1000) c_1000[1] <- 1
			else c_1000[1] <- c_1000[1] + 1
			
		# Propose to update immune states	
		} else {
			
			# Propose only a single non-index case
			propupdate <- fnPropImmuneUpdatesSingle(nononind,mapnonind,vecI)
			propvecI <- vecI
			propvecI[propupdate$person] <- propupdate$state
			propHH <- indDb$hhIndex[propupdate$person]
			propLoglike <- calc_ln_hh(propHH,rp,bp,noPeople,noTests,vecT,propvecI,matT,matV,indDb,hhlookup)
			x <- propLoglike - loglike[propHH]	
			
			# If debig or accept prob, accept the update
			if (db_prop || MH_sym_ap(x)) {

				# Update all immune states
				# vecI <- propvecI
				# loglike <- propLoglike
		
				# Or update only a single immune state
				vecI[propupdate$person] <- propupdate$state
				loglike[propHH] <- propLoglike
				
				# Manage the accept probabilities
				p_1000[c_1000[3],3] <- 1
				
			# Decline the update	
			} else {
				p_1000[c_1000[3],3] <- 0
			}
			
			# Manage the accept probabilities
			if (c_1000[3] == 1000) c_1000[3] <- 1
			else c_1000[3] <- c_1000[3] + 1 
			
		}
		
		# Increment the main loop variable
		index <- index + 1
				
		# If the sample is to be sampled, then sample
		if (index == next_extract) {
			fnExtract(	sampCurStored+1,index,sum(loglike),bp,
						namesfitted,t_I_s_stored,i_g_s_stored,vecT,vecI,p_1000)
			sampCurStored <- sampCurStored + 1
			next_extract <- next_extract + rp["sampFreq"]
			if (sampCurStored == rp["noSampStored"]) {
				write.table(output_tab,file=fnOutputTrace,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
				sampCurStored <- 0
			}
		}
		
	}
	
	# Check to see if there are samples left over from the main loop
	if (sampCurStored > 0) {
		write.table(output_tab[1:sampCurStored,],file=fnOutputTrace,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
	}
	
	# Save the state of the chain as a binary
	save(	file=fnOutputSnapshot,
			noTests,
			matT,
			matV,
			vecT,
			vecI,
			bp,
			indDb,
			hhlookup,
			loglike,
			noPeople
	)
	
}

npi_sim <- function(mnTable,params,runps,loud=FALSE) {

# Function to simulate an outbreak in a single household
# Takes table of household characteristics
	
	back_dur = 1/params["h_B"] # Averge duration between background infections
	max_ID <- dim(mnTable)[1] # Number of individuals in the household
	mnTable$i_g[1] <- 1 # Set the index to have been infected at some point
			
	# Generate time of infection for index
	mnTable$t_I[1] <- mnTable$t_S[1] - fnRIncubation(1,params["p1_symp"],params["p2_symp"]) # draw from gamma distribution
	if (mnTable$t_I[1] < runps["t_min"]) mnTable$t_I[1] = runps["t_min"] 					# if infection time is less that start time, make it sart time
	
	# Loop through from current time to max time
	# need to change the line below to reflect the start of the simulation
	t 		<- runps["t_min"] 
	max_t 	<- runps["t_max"]
	
	while (t < max_t) {
		
		# Set up variables for the loop
		if (t < mnTable$t_I[1]) shortest_delay <- mnTable$t_I[1] - t
		else shortest_delay <- max_t - t
		inf_ID <- 1 
		index_infectee <- -999			
		index_infector <- -999
		event <- 0
		
		# Person to person transmission
		while (inf_ID <= max_ID) {
			if (mnTable$t_I[inf_ID] <= t) {
				sus_ID <- 1 
				while (sus_ID <= max_ID) {
					if (mnTable$i_g[sus_ID] == 0 && sus_ID != inf_ID) {
						
						# Lazy evaluation of more complex inverse cumulative probabilities
						udev <- runif(1)
						test_cum_prob <- int_ptp_lambda(t-mnTable$t_I[inf_ID]+shortest_delay,
														t-mnTable$t_I[inf_ID],
														params,
														i=inf_ID,
														ft=mnTable)
						if (udev < test_cum_prob) {
							opt_res <- uniroot(	opt_int_ptp_lambda,
												interval=c(t-mnTable$t_I[inf_ID],
												t-mnTable$t_I[inf_ID]+shortest_delay),
												p=udev,
												ps=params,
												ts=t-mnTable$t_I[inf_ID],
												i=inf_ID,
												ft=mnTable)
							shortest_delay <- opt_res$root+mnTable$t_I[inf_ID]-t
							event <- 1
							index_infectee <- sus_ID
							index_infector <- inf_ID
						}
						
					}
					sus_ID <- sus_ID + 1
				}
			}
			inf_ID <- inf_ID + 1
		} 
		
		# Background transmission
		sus_ID <- 1
		while (sus_ID <= max_ID) {
			if (mnTable$i_g[sus_ID] == 0) {
				udev = runif(1)
				delay <- -back_dur*log(udev)
				if (delay < shortest_delay) {
					event <- 3
					shortest_delay <- delay
					index_infectee <- sus_ID
					index_infector <- -888
				}
			}
			sus_ID <- sus_ID + 1
		}
		
		t <- t + shortest_delay
		if (t < max_t && index_infectee > 0)  {
			mnTable$t_I[index_infectee] <- t
			mnTable$i_g[index_infectee] <- 1
			if (loud) {
				cat("r\t",-999,"\t",index_infector,"\tinfected\t",index_infectee,"\tat\t",t,"\tdays.\n")
				flush.console()
			}
		}

	}
	
	# Generate times of symptoms for non-index cases
	# (time of symptoms for index cases is assumed known)
	for (i in 2:max_ID) {
		mnTable$t_S[i] <- gen_symp(params,runps,mnTable$t_I[i],mnTable$i_g[i])
	}
	
	mnTable

}

simStudy <- function(	bps,rps,tabChar,tabRes,
						inf_times_known=TRUE,
						immunity_known=TRUE,
						sym_times_known=TRUE) {

	hh <- tabChar
	reald <- tabRes
	simd <- tabRes
	intNumberHouses <- dim(hh)[1]
	max_hh_size <- max(hh$familysize)
	
	maxhouse <- as.data.frame(array(max_hh_size,c(0,7)))
	names(maxhouse) <- c("ID","Number","t_I","i_g","int_g","age","t_S")
	for (i in 1:max_hh_size) maxhouse[i,] <- c(i,1,999,0,-1,-1,999)
	cur_ind <- 1
	
	# Cycle through each household and simulate the outbreak
	for (i in 1:intNumberHouses) {
		tmpsize <- hh$familysize[i]
		tmpindst <- hh$index_first[i]
		tmphouse <- maxhouse[1:tmpsize,]
		tmphouse$int_g <- tabRes$int_g[tmpindst:(tmpindst+tmpsize-1)]
		tmphouse$age <- tabRes$age[tmpindst:(tmpindst+tmpsize-1)]
		tmphouse$t_S <- tabRes$t_S[tmpindst:(tmpindst+tmpsize-1)]
		op <- npi_sim(tmphouse,bps,rps)
		for (j in 1:tmpsize) {
			simd$t_S[tmpindst+j-1] <- op$t_S[j]
			if (immunity_known) simd$i_g[tmpindst+j-1] <- op$i_g[j]
			if (inf_times_known) simd$t_I[tmpindst+j-1] <- op$t_I[j]
			for (k in 0:4) {	
				if (!is.na(reald[tmpindst+j-1,paste("V",k,sep="")])) {
					time <- reald[tmpindst+j-1,paste("T",k,sep="")] 
					tau <-  time - op$t_I[j]
					simd[tmpindst+j-1,paste("V",k,sep="")] <- 0
					if (dens_ptp_lambda(tau,bps,tmpindst+j-1,tabRes) > bps["h_detect"]) {
						if (runif(1)< bps["spec"]) simd[tmpindst+j-1,paste("V",k,sep="")] <- "A"
					} else {
						if (runif(1)> bps["sens"]) simd[tmpindst+j-1,paste("V",k,sep="")] <- "A"
					}
				}
			}
		}	
	}
	simd
}

fnDebugLike <- function(bp,rp,indDb,hhlookup,t_factor=1,v_factor=1,s_factor=1) {			
	
	# New bits below
	matT <- data.matrix(indDb[,c("T0","T1","T2","T3","T4")])
	tmpMatV <- as.matrix(indDb[,c("V0","V1","V2","V3","V4")])
	matV <- apply(tmpMatV,c(1,2),fnConvCulture)
	vecT <- as.vector(indDb$t_I,mode="numeric")
	vecI <- as.vector(indDb$i_g,mode="numeric")
	
	noTests <- 5
	noPeople <- dim(indDb)[1]
	loglike <- calc_ln(rp,bp,noPeople,noTests,vecT,vecI,matT,matV,indDb,hhlookup,
						t_factor=t_factor,v_factor=v_factor,s_factor=s_factor)
	
	rtn <- sum(loglike)
	
	rtn
	
}

fnMargLike <- function(pname,vals,bp,rp,data,hhlookup,t_factor=1,v_factor=1,s_factor=1) {
	rtn <- vals
	tmpps <- bp
	nvals <- length(vals)
	for (i in 1:nvals) {
		tmpps[pname] <- vals[i]
		rtn[i] <- fnDebugLike(tmpps,rp,data,hhlookup,t_factor=t_factor,v_factor=v_factor,s_factor=s_factor)			
	}
	rtn
}

fnMargLikeIndexInfs <- function(bp,rp,data,vals = 1:20/20*2, t_factor=1,v_factor=1,s_factor=1) {
	rtn <- vals
	tmptis <- data$t_I
	nvals <- length(vals)
	for (i in 1:nvals) {
		data$t_I <- tmptis * vals[i]
		rtn[i] <- fnDebugLike(bp,rp,data,t_factor=t_factor,v_factor=v_factor,s_factor=s_factor)			
	}
	rtn
}

fnAssignInitialTimeValues <- function(idb,rps) {
	rtn <- idb
	norows <- dim(rtn)[1]
	for (i in 1:norows) {
		if (rtn$member[i]==0) t_I_guess <- rps["t_min"] + (min(idb$t_S) - rps["t_min"]) / 2
		else t_I_guess <- 0
		rtn$t_I[i] <- t_I_guess
	}
	rtn
}

fnDataIsPositive <- function(s) {
	rtn <- FALSE
	if (is.numeric(s)) {
		if (s > 1e-100) rtn <- TRUE
	} else {
		if (!is.na(s)) {
			if (s == "A" || s == "B") rtn <- TRUE
		}
	}
	rtn
}

fnAssignInitialInfectionStates <- function(idb) {
	rtn <- idb
	norows <- dim(idb)[1]
	for (i in 1:norows) {
		vals <- idb[i,c("V0","V1","V2","V3","V4")]
		if (any(sapply(vals,fnDataIsPositive))) {
			rtn$i_g[i] <- 1
		}
		else rtn$i_g[i] <- 0
		}
	rtn
}

fnCountApparentInfections <- function(idb) {
	rtn <- 0
	norows <- dim(idb)[1]
	for (i in 1:norows) {
		vals <- idb[i,c("V0","V1","V2","V3","V4")]
		if (any(sapply(vals,fnDataIsPositive))) {rtn <- rtn + 1}
	}
	rtn
}

fnRIncubation <- function(n,p1,p2) {
	# rgamma(1,scale=p1,shape=p2)
	rlnorm(1,meanlog=p1,sdlog=p2)
}

fnDIncubation <- function(t,p1,p2,log) {
	# dgamma(t,scale=p1,shape=p2,log=log)
	dlnorm(t,meanlog=p1,sdlog=p2,log=log)
}

fnPIncubation <- function(t,p1,p2) {
	# pgamma(t,scale=p1,shape=p2)
	plnorm(t,meanlog=p1,sdlog=p2)
}

fnDInfection <- function(t,p1,p2,log) {
	dgamma(t,scale=p1,shape=p2,log=log)
}

fnPInfection <- function(t,p1,p2) {
	pgamma(t,scale=p1,shape=p2)
}
