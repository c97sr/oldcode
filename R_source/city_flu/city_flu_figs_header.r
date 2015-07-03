flu_processor <- function(root,stem,leaf,policy,s,f) {
	## s for start and f for finish. Inclusive of both.
	if (s >= f) stop("f must be greater than s\n")
	colnames <- c("S_365","Max_Q","Max_I","Sum_A","Max_T","Max_CT","Max_Q_Inc","Max_I_Inc","Ave_I")
	nocols <- length(colnames)
	norows <- f-s+1
	rtnval = array(0,dim=c(norows,nocols))
	rtnval <- as.data.frame(rtnval)
	names(rtnval) <- colnames 
	for (i in s:f) {
		tab <- read.table(file=paste(root,stem,i,"_",policy,leaf,sep=""),header=TRUE,sep="\t",row.names=NULL)
		## names(tab) <- names(tab)[2:length(names(tab))]
		rtnval$S_365[i-s+1] <- mean(tab$S365)
		rtnval$Max_Q[i-s+1] <- mean(tab$maxIndividualQuarantined)
		rtnval$Max_I[i-s+1] <- mean(tab$maxIndividualIsolated)
		rtnval$Sum_A[i-s+1] <- mean(tab$DrugUsed)
		rtnval$Max_T[i-s+1] <- mean(tab$maxActiveTestServers)
		rtnval$Max_CT[i-s+1] <- mean(tab$maxActiveCTServers)
		rtnval$Max_Q_Inc[i-s+1] <- mean(tab$maxDailyQuarantined)
		rtnval$Max_I_Inc[i-s+1] <- mean(tab$maxDailyIsolated)
		rtnval$Ave_I[i-s+1] <- mean(tab$avgIndividualIsolated)
	}
	rtnval
	## tab
}

flu_processor_dash <- function(root,stem,leaf,policy,s,f) {
	## s for start and f for finish. Inclusive of both.
	if (s >= f) stop("f must be greater than s\n")
	colnames <- c("S_365","Max_Q","Max_I","Sum_A","Max_T","Max_CT","Max_Q_Inc","Max_I_Inc","Ave_I")
	nocols <- length(colnames)
	norows <- f-s+1
	rtnval = array(0,dim=c(norows,nocols))
	rtnval <- as.data.frame(rtnval)
	names(rtnval) <- colnames 
	for (i in s:f) {
		tab <- read.table(file=paste(root,stem,i,"_",policy,leaf,sep=""),header=TRUE,sep="\t",row.names=NULL)
		names(tab) <- names(tab)[2:length(names(tab))]
		rtnval$S_365[i-s+1] <- mean(tab$S365)
		rtnval$Max_Q[i-s+1] <- mean(tab$maxIndividualQuarantined)
		rtnval$Max_I[i-s+1] <- mean(tab$maxIndividualIsolated)
		rtnval$Sum_A[i-s+1] <- mean(tab$DrugUsed)
		rtnval$Max_T[i-s+1] <- mean(tab$maxActiveTestServers)
		rtnval$Max_CT[i-s+1] <- mean(tab$maxActiveCTServers)
		rtnval$Max_Q_Inc[i-s+1] <- mean(tab$maxDailyQuarantined)
		rtnval$Max_I_Inc[i-s+1] <- mean(tab$maxDailyIsolated)
		rtnval$Ave_I[i-s+1] <- mean(tab$avgIndividualIsolated)
	}
	rtnval
	## tab
}


## Colour assignment
c_ass <- function(int) {
	if (int == "Q") col="green"
	else if (int == "None") col="black"
	else if (int == "QI") col="blue"
	else if (int == "QIA") col="cyan"
	else if (int == "QIAT") col="red"
	else if (int == "QA") col="magenta"
	else if (int == "QIAC") col="orange"
	else col="grey"
	col
}

plot_int <- function(stem,type,int,colour,PI,popsize,thickness=1) {
	data_filename = paste(stem,type,"_",int,".txt",sep="")
	data_raw <- read.table(file=paste(data_filename,sep=""))
	data <- ts_mean_95_pi(data_raw)
	if (PI==TRUE) {for (i in 1:2) points((data[i+1,]/popsize),col=colour,lwd=1,type="l")}
	points((data[1,]/popsize),type="l",col=colour,lwd=thickness)
}

table_int <- function(stem,list_int) {
	no_ints <- length(list_int)
	colnames <- c("total_infs","peak_inc_inf","total_treat","peak_Q","peak_Isol","peak_prev_I")
	nocols <- length(colnames)
	rtnval <- array(0,c(no_ints,nocols))
	rtnval <- as.data.frame(rtnval)
	names(rtnval) <- colnames 
	rownames(rtnval) <- list_int
	for (i in 1:no_ints) {
		scen = list_int[i]
		data_raw_inf <- read.table(file=paste(stem,"DailyIncidence_",scen,".txt",sep=""))
		no_days = dim(data_raw_inf)[1]
		no_reals = dim(data_raw_inf)[2]
		aux = array(0,c(no_reals))
		for (j in 1:no_reals) aux[j]=max(data_raw_inf[,j])
		rtnval$peak_inc_inf[i]=quantile(aux,probs=c(0.95))[1]
		rtnval$total_infs[i]=sum(data_raw_inf)/no_reals
		data_raw_q <- read.table(file=paste(stem,"numIndividualQuarantined_",scen,".txt",sep=""))
		for (j in 1:no_reals) aux[j]=max(data_raw_q[,j])
		rtnval$peak_Q[i]=quantile(aux,probs=c(0.95))[1]
		data_raw_isol <- read.table(file=paste(stem,"numIndividualIsolated_",scen,".txt",sep=""))
		for (j in 1:no_reals) aux[j]=max(data_raw_isol[,j])
		rtnval$peak_Isol[i]=quantile(aux,probs=c(0.95))[1]
		for (j in 1:no_reals) aux[j]=max(data_raw_inf[,j])
		rtnval$peak_prev_I[i]=quantile(aux,probs=c(0.95))[1]
		data_raw_drugs <- read.table(file=paste(stem,"DrugsUsed_",scen,".txt",sep=""))
		rtnval$total_treat[i]=mean(data_raw_drugs[no_days,])[1]

	}
	rtnval
}


ts_mean_95_pi <- function(simarr) {
	no_days = dim(simarr)[1]
	no_reals = dim(simarr)[2]
	rtnmat = array(0,c(3,no_days))
	for (i in 1:no_days) {
		rtnmat[1,i]=mean(as.numeric(simarr[i,]))
		quant = quantile(as.numeric(simarr[i,]),probs=c(0.025,0.975))
		for (j in 1:2) rtnmat[j+1,i] <- quant[j]
	}
	rtnmat
}
