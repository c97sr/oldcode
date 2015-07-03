# Next - debug the code to 

rm(list=ls(all=TRUE))
options(error=recover)
#options(error=NULL)
require(date)

deliver <- function(S,E,I,R,p_di,acs,nc) {
	
	# Function to take current state of the system S susceptibles
	# E exposed, I infectious and R recovered, p_di probability
	# that delivered birds are infected, acs average cage size and 
	# nc number of cages. Returns a list of vectors of infectious delivered
	# and susceptible delivered
	
	N <- 		S+E+I+R
	N_del <- 	rep(acs,nc) - N
	I_del <- 	rbinom(nc,N_del,p_di)
	S_del <- 	N_del - I_del
	list(I=I_del,S=S_del)
	
}

lambda <- function(S,E,I,R,beta,alpha,ncages) {
	N <- S+E+I+R
	lambda <- rep(0,ncages)
	sigma_I <- sum(I)
	sigma_N <- sum(N)
	if (sigma_N > 0) {
		for (i in 1:ncages) {
			lambda[i] <- lambda[i]+alpha*(sigma_I-I[i])/(sigma_N-N[i])
			if (N[i] > 0) lambda[i] <- lambda[i]+I[i]/N[i]
		}
		lambda <- lambda*beta
	}
	lambda
}

noCagesInfected <- function(I,epsilon) {
	n <- length(I)
	rtn <- 0
	for (i in 1:n) {if (I[i] > epsilon) rtn <- rtn + 1}
	rtn
}

estParamsDpoiszero <- function(v,n) {
	browser()
	f <- sum(v)
	mean <- sum((0:n)*v)
	poi_zero_cont <- dpois(0,mean,FALSE)
	if (poi_zero_cont*f > v[0]) pzero <- v[0]/f - poi_zero_cont
	else pzero <- 0
	list(mean=mean,pzero=pzero)
}

dpoiszero <- function(x,mean,pzero) {
	if (x==0) rtn <- (1-pzero)*dpois(0,mean,FALSE) + pzero
	else rtn <- rtn <- (1-pzero)*dpois(x,mean,FALSE)
	rtn
}

simMarket <- function(bps,rps) {
	
	no_measures <- length(rps$measure_times)
	cage_prev <- array (0,c(no_measures,bps$n_c+1))
	indinc <- array(0,c(no_measures))
	max_time <- max(rps$measure_times)
	rps$measure_times <- c(rps$measure_times,99999)
	rps$wash_days <- c(rps$wash_days,99999)
	
	for (i in 1:rps$no_realizations) {
		
		# Initial conditions
		S <- rep(0,bps$n_c)
		E <- rep(0,bps$n_c)
		I <- rep(0,bps$n_c)
		R <- rep(0,bps$n_c)
		
		time <- 0
		next_day <- 0
		measure_index <- 1
		wash_index <- 1
		next_measure <- rps$measure_times[measure_index]
		next_wash <- rps$wash_days[wash_index]
		indincdt <- 0
		
		while (time < (max_time+1+rps$epsilon)) {
			
			if (time > (next_measure - rps$epsilon)) {
				# Measure cages
				noInfCages <- noCagesInfected(I,rps$epsilon)
				cage_prev[measure_index,noInfCages+1] <- cage_prev[measure_index,noInfCages+1]+1
				indinc[measure_index] <- indinc[measure_index] + indincdt
				measure_index <- measure_index + 1
				next_measure <- rps$measure_times[measure_index]
				indincdt <- 0
			}

			if (time > (next_wash - rps$epsilon)) {
				# Measure cages
				S <- rep(0,bps$n_c)
				E <- rep(0,bps$n_c)
				I <- rep(0,bps$n_c)
				R <- rep(0,bps$n_c)
				wash_index <- wash_index + 1
				next_wash <- rps$wash_days[wash_index]
			}
			
			# Morning delivery
			if (time > (next_day - rps$epsilon)) {
				delivery <- deliver(S,E,I,R,bps$p_di,bps$s_c,bps$n_c)
				S <- S + delivery$S
				I <- I + delivery$I
				next_day <- next_day + 1
			}
			
			N <- S+E+I+R
			
			# cat(i,time,N,I,"\n",sep="\t")
			# flush.console()
			
			# Define probabilities
			p_Ifn <- 1-exp(-lambda(S,E,I,R,bps$beta,bps$alpha,bps$n_c)*rps$dt)
			p_Ifs <- 1-exp(-bps$D_E*rps$dt)
			p_Rec <- 1-exp(-bps$D_I*rps$dt)
			p_Buy <- 1-exp(-bps$h_db*rps$dt)
			
			# Draw random events
			x_Ifn <- rbinom(bps$n_c,S,p_Ifn)
			x_Ifs <- rbinom(bps$n_c,E,p_Ifs)
			x_Rec <- rbinom(bps$n_c,I,p_Rec)
			
			if (is.na(x_Ifn) || is.na(x_Ifs) || is.na(x_Rec)) browser()
			
			# Update state variables for infection process
			S <- S - x_Ifn
			E <- E + x_Ifn - x_Ifs
			I <- I + x_Ifs - x_Rec
			R <- R + x_Rec
			
			# Update state variables for buying process
			S <- S - rbinom(bps$n_c,S,p_Buy)
			E <- E - rbinom(bps$n_c,E,p_Buy)		
			I <- I - rbinom(bps$n_c,I,p_Buy)		
			R <- R - rbinom(bps$n_c,R,p_Buy)		
			
			# if (sum(x_Ifn > 0)) browser()
			indincdt <- indincdt + sum(x_Ifn)
			time <- time + rps$dt
			
			# cat("                      ",time-rps$dt,time,x_Ifn,"\n",sep="\t")
			# flush.console()
		}
	}

	list(cp=cage_prev,ii=indinc)

}

# Define biological parameters
biops <- data.frame(
	R0 				= 1,
	alpha 			= 0.1,
	D_E 			= 2,
	D_I 			= 4.3,
	h_db 			= 0.5,
	p_di 			= 0.01,
	n_c 			= 10,
	s_c 			= 10
	)
biops$beta = biops$R0 / biops$D_I / (biops$alpha + (biops$s_c - 1) / biops$s_c )

# Define run parameters
runps <- list(
	dt 				= 1/4,
	epsilon 		= 1e-20,
	no_realizations = 100,
	measure_times 	= 1:30,
	wash_days 		= c(15.5)
	)

biops$R0 <- 10	
biops$h_db <- 0.00001
biops$n_c <- 10
biops$p_di <- 0.5
simMarket(biops,runps)

# Data stuff
fileData <- "hk_poultry/data/8market_raw1.csv"
matData <- read.csv(fileData)
matData <- matData[which(matData$Market==1),]
matData$Sample.Date <- as.date(as.character(matData$Sample.Date),order="dmy")
tabtmpDat <- table(matData$Sample.Date,matData$H9)
tabData <- 	data.frame(	date=as.numeric(row.names(tabtmpDat)),
		N=tabtmpDat[,"0"] + tabtmpDat[,"1"],
		n=tabtmpDat[,"1"]
		)
plot(as.date(tabData$date),tabData$n/tabData$N,ylim=c(0,1))
