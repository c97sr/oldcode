rm(list=ls(all=TRUE))
datafile <- "D:\\files\\projects\\schisto\\individual_based\\dummy_humans.in"
data <- as.data.frame(read.table(datafile,header=TRUE))

nohumans <- dim(data)[1]
data <- cbind(data,array(-999,nohumans))
names(data)[dim(data)[2]]<-"t_Infed"
data <- cbind(data,array(-999,nohumans))
names(data)[dim(data)[2]]<-"t_Infous"
data <- cbind(data,array(-999,nohumans))
names(data)[dim(data)[2]]<-"t_Rec"
data <- cbind(data,array(-999,nohumans))
names(data)[dim(data)[2]]<-"NxtSt"

dt <- 0.1
max_t <- 50
beta <- 1
gamma <- 1/3
sigma <- 1/4

for (i in 1:nohumans) if (data$ISt[i]==3) {
	data$t_Infous[i]=0
}

for (t in (0:(max_t/dt))*dt) {
	noInfectious <- 0
	for (i in 1:nohumans) { 
		if (data$ISt[i]==3) noInfectious <- noInfectious + 1
		data$NxtSt[i] <- data$ISt[i]
	}
	for (i in 1:nohumans) {
		if (data$ISt[i]==1) {
			FOI <- noInfectious / nohumans * beta
			if (runif(1) < 1-exp(-FOI*dt)) {
				data$NxtSt[i] <- 2
				data$t_Infed[i] <- t+dt
			} 
		}  
		if (data$ISt[i]==2) {
			if (runif(1) < 1-exp(-gamma*dt)) {
				data$NxtSt[i] <- 3
				data$t_Infous[i] <- t+dt
			} 
		}  
		if (data$ISt[i]==3) {
			if (runif(1) < 1-exp(-sigma*dt)) {
				data$NxtSt[i] <- 4
				data$t_Rec[i] <- t+dt
			} 
		}  
	}
	for (i in 1:nohumans) data$ISt[i] <- data$NxtSt[i] 
}

outfile <- "D:\\files\\projects\\schisto\\individual_based\\dummy_out.out"
write.table(data,file=outfile,sep="\t",row.names=FALSE,col.names=TRUE)




