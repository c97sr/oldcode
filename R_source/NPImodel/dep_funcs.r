lnlike_aux <- function(k,tab,ps) {
	back_dur <- -ps$T_B/log(1-ps$q_B)
	t_e <- min(tab$t_I[tab$index_first[k]:(tab$index_first[k]+tab$hh_size[k]-1)]) 	# Earliest household infection
	t_l <- tab$t_R[k]																# Last household observation
	set_others <- setdiff(tab$index_first[k]:(tab$index_first[k]+tab$hh_size[k]-1),k)
	rtn <- 0
	if (tab$t_I[k]>t_l) t_end <- t_l else t_end <- tab$t_I[k] 
	rtn <- rtn - (t_end - t_e) * 1 / back_dur
	for (i in set_others)
		rtn <- rtn - int_ptp_lambda(t_end-tab$t_I[i],0,ps)
	if (t_end < t_l) {
		rtn <- rtn + log(1/back_dur)
		p_inf <- 0
		for (i in set_others) p_inf <- p_inf + dens_ptp_lambda(t_end-tab$t_I[i],amPs)
		if (p_inf > 0) rtn <- rtn + log(p_inf)
	}
	rtn
}

lnlike_flu <- function(tab,ps,reval=TRUE,incindex=TRUE) {
	if (reval) for (i in 1:(dim(tab)[1])) tab$lnlike[i] <- lnlike_aux(i,tab,ps)
	rtn <- 0
	for (i in 1:(dim(tab)[1])) {
		if (tab$member[i]!=0) rtn <- rtn + tab$lnlike[i]
		else if (incindex) rtn <- rtn + tab$lnlike[i]
	}
	rtn
}

amoy_sim_old <- function (mnTable,params,outfile="NONE",seed=-1,realno=1,loud=TRUE,max_t=9) {
	
	if (seed > 0) set.seed(seed)
	back_dur = -params$T_B/log(1-params$q_B)
	t <- 0
	max_ID <- dim(mnTable)[1]
	
	# Generate onset times for individuals known to be infectious at t=0
	for (i in 1:max_ID) {
		if (mnTable$t_I[i] < max_t && mnTable$t_O[i] > max_t) {
			if (params$p_A < runif(1)) {
				gdev <- rgamma(1,params$alpha_O,scale=params$mu_O/params$alpha_O)
				mnTable$t_O[i] <- mnTable$t_I[i] + gdev
			}
		}
	}
	
	# Loop through from current time to max time
	while (t < max_t) {
		index_inf <- -999
		shortest_delay <- max_t - t
		inf_ID <- 1 
		index_infectee <- -999			
		index_infector <- -999
		event <- 0
		# Person to person transmission
		while (inf_ID <= max_ID) {
			if (mnTable$t_I[inf_ID] <= t && mnTable$t_R[inf_ID] > t) {
				sus_ID <- 1 
				while (sus_ID <= max_ID) {
					if (mnTable$t_I[sus_ID] > max_t && mnTable$t_R[sus_ID] > t && sus_ID != inf_ID) {
						# This is how to efficiently evaluate more complex inverse cumulative probabilities
						udev <- runif(1)
						# XXX have to correct time here
						test_cum_prob <- int_ptp_lambda(t-mnTable$t_I[inf_ID]+shortest_delay,t-mnTable$t_I[inf_ID],params,1,1,1,1,1,1)
						if (udev < test_cum_prob) {
							opt_res <- uniroot(opt_int_ptp_lambda,interval=c(t-mnTable$t_I[inf_ID],t-mnTable$t_I[inf_ID]+shortest_delay),p=udev,ps=params,ts=t-mnTable$t_I[inf_ID])
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

		# Seed transmission
		sus_ID <- 1
		while (sus_ID <= max_ID) {
			if (mnTable$t_I[sus_ID] > max_t && mnTable$t_R[sus_ID] > t) {
				udev <- runif(1)
				test_cum_prob <- int_seed_lambda(t+shortest_delay,t,params,1,1,1,1,1,1)
				if (udev < test_cum_prob) {
					opt_res <- uniroot(opt_int_seed_lambda,interval=c(t,t+shortest_delay),p=udev,ps=params,ts=t)
					shortest_delay <- opt_res$root-t
					event <- 2
					index_infectee <- sus_ID
					index_infector <- -999
				}
			}
			if (mnTable$t_I[sus_ID] > max_t && mnTable$t_R[sus_ID] > t) {
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
		if (t < max_t)  {
			mnTable$t_I[index_infectee] <- t
			if (loud) {
				cat("r\t",realno,"\t",index_infector,"\tinfected\t",index_infectee,"\tat\t",t,"\tdays.\n")
				flush.console()
			}
			if (params$p_A < runif(1)) {
				gdev <- rgamma(1,params$alpha_O,scale=params$mu_O/params$alpha_O)
				mnTable$t_O[index_infectee] <- mnTable$t_I[index_infectee] + gdev
			}
		} 
	}
	if (outfile != "NONE") {
		if (loud) {	
			cat("Simulation ends at ",t,"\tdays.\n")
			flush.console()
		}
		mnTable <- cbind(mnTable,array(999,max_ID))
		noField <- dim(mnTable)[2]
		names(mnTable)[noField] <- "t_I_guess"
		mnTable <- mnTable[order(mnTable$t_O,mnTable$ID),]
		mnTable$t_I_guess[1] <- 0
		for (i in 2:max_ID) if (mnTable$t_I[i] < max_t) mnTable$t_I_guess[i] <- mnTable$t_I_guess[i-1]+ (min(mnTable$t_O[i],mnTable$t_R[i])-mnTable$t_I_guess[i-1])/2
		write.table(mnTable,file=outfile,sep="\t",col.names=TRUE,row.names=FALSE)
	}
	mnTable
}


plotMCMCDensity <- function(df,variable,log=TRUE,vmin=1,vmax=1000,burnin=0,plot_points=100,nk=10) {

	e2 = dim(df)[1]
	if (burnin > e2-10) stop("Burnin range too large")
	s1 = burnin+1
	s2 = floor(burnin + (e2-burnin)/2)
	e1 = s2-1

	xvals = (0:plot_points)/plot_points
	xvals = sapply(xvals,fnUnitToLog10Scale,max=vmax,min=vmin)

	spRes1 <- logspline(df[s1:e1,variable],vmin,vmax,nknots=nk)
	spRes2 <- logspline(df[s2:e2,variable],vmin,vmax,nknots=nk)

	yvals_s1 <- sapply(xvals,dlogspline,fit=spRes1)
	yvals_s2 <- sapply(xvals,dlogspline,fit=spRes2)
	uniprior <- array(1/(vmax-vmin),c(plot_points+1))

	plot(1:2,xlim=c(vmin,vmax),ylim=c(0,max(yvals_s1,yvals_s2,uniprior)),log="x",type="n")
	points(xvals,yvals_s1,type="l",col="green")
	points(xvals,yvals_s2,type="l",col="red")
	points(xvals,uniprior,type="l",col="blue")
	
}

# Way to test pairwise function
test_N <- 1000
test_n <- 0
for (i in 1:test_N) {
	tmpsize <- 2
	tmpindst <- hh$index_first[i]
	tmphouse <- maxhouse[1:tmpsize,]
	op <- amoy_sim(tmphouse,amPs,"NONE",realno=j,loud=FALSE)
	if (op$t_I[2] < 990) test_n <- test_n + 1
}

test_n

# Useful lines for profiling
Rprof(filename="c:\\tmp\\Rprof.out")
source("C:\\eclipse\\R Source\\NPImodel\\flu_time_since_sim.r")
Rprof(NULL)
summaryRprof("c:\\tmp\\Rprof.out")


