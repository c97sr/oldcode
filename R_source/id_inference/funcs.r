strFLU <- function(
			R0_1=1.8,
			R0_2=1.8,
			p_1=0,
			p_2=0,
			alpha_12=0.5,
			t_1=0,
			t_2=50,
			D_I=2.6,
			N=10000,
			D_S=0.25,
			tps=(0:400)/4,
			reals=1,
			dt=0.25,
			stochastic=TRUE,
			plot=TRUE) {
	
	noTPts <- length(tps) 
	epsilon <- 1e-10
	
	if (noTPts < 3) 
		stop("simSEIR: tps must be of length 2 or greater")
	if (reals > 1 && stochastic==FALSE) 
		stop("simSEIR: why make multiple realisations of a deterministic model?")
	
	rtn <- array(0,c(noTPts-1,2))
	
	for (i in 1:reals) {
		
		t_cur <- tps[1]
		ind_t <- 2
		t_next <- tps[ind_t]
		
		S_0 <- round(N*(1-p_1-p_2))
		S_1 <- round(N*p_1)
		S_2 <- N-S_0-S_1
		I_1 <- 0
		I_2 <- 0
		
		while (ind_t <= noTPts) {

			if (t_cur > t_1 - epsilon) hSed_1 <- 1/D_S/N 
			else hSed_1 <- 0
			if (t_cur > t_2 - epsilon) hSed_2 <- 1/D_S/N 
			else hSed_2 <- 0
			if (hSed_1 > 0 && hSed_2 > 0) {
				hSed_1 <- hSed_1 / 2
				hSed_2 <- hSed_2 / 2
			}
			
			lambda_1 <- R0_1*I_1/D_I/N + hSed_1
			lambda_2 <- R0_2*I_2/D_I/N + hSed_2
			
			pInf_0 	<- 1 - exp(-(lambda_1+lambda_2)*dt)
			pInf_1 	<- 1 - exp(-lambda_1*dt)
			pInf_2 	<- 1 - exp(-lambda_2*dt)
			pRec 	<- 1 - exp(-dt/D_I)
			
			if (stochastic) {
				if (lambda_1 + lambda_2 > 0) {					
					nInf_0 	<- rbinom(1,S_0,pInf_0)
					nInf_01 <- rbinom(1,nInf_0,lambda_1/(lambda_1+lambda_2))
					nInf_02 <- rbinom(1,nInf_0,lambda_2/(lambda_1+lambda_2))
				} else {
					nInf_01 <- 0
					nInf_02 <- 0
				}
				nInf_21 <- rbinom(1,S_2,pInf_1*(1-alpha_12))
				nInf_12 <- rbinom(1,S_1,pInf_2*(1-alpha_12))
				nRec_1 	<- rbinom(1,I_1,pRec)
				nRec_2 	<- rbinom(1,I_2,pRec)
				# browser()
			} else {
				if (lambda_1 + lambda_2 > 0) {					
					nInf_0 	<- S_0*pInf_0
					nInf_01 <- nInf_0*lambda_1/(lambda_1+lambda_2)
					nInf_02 <- nInf_0*lambda_2/(lambda_1+lambda_2)
				} else {
					nInf_01 <- 0
					nInf_02 <- 0
				}
				nInf_21 <- S_2*pInf_1*(1-alpha_12)
				nInf_12 <- S_1*pInf_2*(1-alpha_12)
				nRec_1 	<- I_1*pRec
				nRec_2 	<- I_2*pRec
			}
			
			S_0 <- S_0 - nInf_01 - nInf_02 
			S_1 <- S_1 + nRec_1  - nInf_12
			S_2 <- S_2 + nRec_2  - nInf_21
			I_1 <- I_1 + nInf_01 + nInf_21 - nRec_1
			I_2 <- I_2 + nInf_02 + nInf_12 - nRec_2
			
			rtn[ind_t-1,1] <- rtn[ind_t-1,1] + nInf_01 + nInf_21 
			rtn[ind_t-1,2] <- rtn[ind_t-1,2] + nInf_02 + nInf_12 
			
			t_cur <- t_cur + dt
			
			if (t_cur > (t_next - epsilon)) {
				t_next <- tps[ind_t]
				ind_t <- ind_t+1
			}
			
		}
	}
	
	rtn <- rtn / reals
	
	if (plot) {
		plot(rtn[,1],type="l",col="red",ylim=c(0,max(rtn[,1],rtn[,2])))
		points(rtn[,2],type="l",col="green")
	}

	list(ave=rtn)
	
}

approxPoissLike <- function(strData,strModel) {
	nstr <- dim(strData)[2]
	lnlike <- 0
	for (i in 1:nstr) {
		lnlike <- lnlike + sum(dpois(strData[,i],strModel[,i],log=TRUE))
	}
	lnlike
}

x <- strFLU(stochastic=TRUE,R0_2=1.8,plot=FALSE)
param_range <- 1:150/100+1
like_range <- rep(0,length(param_range))
for (i in 1:length(param_range)) {
	y <- strFLU(stochastic=FALSE,R0_2=param_range[i],plot=FALSE)
	like_range[i] <- approxPoissLike(x$ave,y$ave)
}
plot(param_range,like_range,type="l")