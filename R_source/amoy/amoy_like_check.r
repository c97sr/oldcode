wrapper <- function() {

	strFName = "D:\\share\\files\\projects\\inference2005\\data\\4floors.in"
	data = read.table(strFName,header=TRUE)
	no_records = dim(data)[1]
	
	inf_log_mean = 5
	inf_log_sd = 2 
	inc_mean = 4.64
	inc_alpha = 1.353
	inc_beta = inc_mean / inc_alpha

	max_inf_time = 100

	lnlike <- 0
	for (i in 1:(no_records-1)) {
		inf_time_i <- data$Guess_Infection[i]
		rem_time_i <- data$Time_Removal[i]
		ons_time_i <- data$Time_Onset[i]
		tau_ons = ons_time_i - inf_time_i
		if (inf_time_i < max_inf_time) lnlike <- lnlike + log(dgamma(tau_ons,shape=inc_alpha,scale=inc_beta))
		for (j in ((i+1):no_records)) {
			if (inf_time_i < max_inf_time) {
				inf_time_j = data$Guess_Infection[j]
				rem_time_j = rem_time_i = data$Time_Removal[i]
				tau_inf = min(inf_time_j,rem_time_j) - min(inf_time_i,rem_time_i)
				lnlike <- lnlike - plnorm(tau_inf,inf_log_mean,inf_log_sd)
				if (inf_time_j < max_inf_time) lnlike <- lnlike + log(dlnorm(tau_inf,inf_log_mean,inf_log_sd))
				# cat(i,j,tau_inf,tau_ons,dlnorm(tau_inf,inf_log_mean,inf_log_sd),"\n")
			}
		}
	}
	lnlike
}

# undebug(wrapper)
out <- wrapper()
cat("lnlike",out,"\n")
 
xdata = 1:3000/100
plot(xdata,plnorm(xdata,inf_log_mean,inf_log_sd))