
# construct the data set

lab <- data.frame(ID = c(1:21))

lab$trunct <- c(rep(0,7),rep(1,13),2)
lab$timeL <- c(0,1,3,1,2,4,1,1,1,2,2,2,5,5,1,4,6,4,2,5,2)
lab$timeR <- c(3,4,7,3,6,7,3,4,4,5,5,5,9,9,4,8,10,6,4,8,6)

lab_1 <- data.frame(ID = c(1:21))
lab_1$trunct <- lab$trunct+1
lab_1$timeL <- lab$timeL+1
lab_1$timeR <- lab$timeR+1

lab_2 <- data.frame(ID = c(1:21))
lab_2$trunct <- lab$trunct-1
lab_2$trunct[lab_2$trunct<0] <- 0
lab_2$timeL <- lab$timeL-1
lab_2$timeL[lab_2$timeL<0] <- 0
lab_2$timeR <- lab$timeR-1

weibull.loglik <- function(parameters, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    alpha <- exp(parameters[1])
    beta <- exp(parameters[2])
   
    n <- length(data$trunct)  
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( (pweibull(data$timeR[i], shape=alpha, scale=beta) - pweibull(data$timeL[i], shape=alpha, scale=beta)) 
	                      / (1-pweibull(data$trunct[i],shape=alpha, scale=beta)) )
    }
    -sum(log.lik)
    #log.lik
}

serial.weibull.day0 <- optim(c(1,1), weibull.loglik, data=lab,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

serial.weibull.day1 <- optim(c(1,1), weibull.loglik, data=lab_1,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

 serial.weibull.day2 <- optim(c(1,1), weibull.loglik, data=lab_2,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

# the mean
exp(serial.weibull.day0$par[2])*gamma(1+exp(-serial.weibull.day0$par[1]))
exp(serial.weibull.day1$par[2])*gamma(1+exp(-serial.weibull.day1$par[1]))
exp(serial.weibull.day2$par[2])*gamma(1+exp(-serial.weibull.day2$par[1]))

# the variance
exp(2*serial.weibull.day0$par[2])*(gamma(1+2*exp(-serial.weibull.day0$par[1]))-(gamma(1+exp(-serial.weibull.day0$par[1])))^2)
exp(2*serial.weibull.day1$par[2])*(gamma(1+2*exp(-serial.weibull.day1$par[1]))-(gamma(1+exp(-serial.weibull.day1$par[1])))^2)
exp(2*serial.weibull.day2$par[2])*(gamma(1+2*exp(-serial.weibull.day2$par[1]))-(gamma(1+exp(-serial.weibull.day2$par[1])))^2)

weibull.log <- function(parameters, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    alpha <- exp(parameters[1])
    beta <- exp(parameters[2])
    
    n <- nrow(data) 
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dweibull(data$time[i], shape=alpha, scale=beta) / (1-pweibull(trunct[i],shape=alpha, scale=beta)) )
    }
    -sum(log.lik)
    #log.lik
}

boot.weibull <- function(reps, data){
  # function to resample the data and fit weibull models to each resample
  n <- length(data[,1])
  temp <- rep(NA, reps)
  output <- data.frame(par1=temp, par2=temp, mean=temp)
  for(i in 1:reps){
    a <- sample(1:n, replace=TRUE)
    temp.data <- data[a,]
    temp.weibull <- optim(c(1.5,1.5), weibull.loglik, data=temp.data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
    output$mean[i] <- exp(temp.weibull$par[2])*gamma(1+exp(-temp.weibull$par[1]))
    output$par1[i] <- temp.weibull$par[1]
    output$par2[i] <- temp.weibull$par[2]
  }
  output
}

bt <- 50

day0.wei.boot <- boot.weibull(bt, lab)     # reps=500,c(0.5,1.5)
quantile(day0.wei.boot$mean, c(0.5, 0.025, 0.975))



pdf("D:\\Work\\Influenza\\output NPI year 1\\serial interval\\serial_lab.pdf", width=5, height=4, version="1.4")
#windows(width=5,height=4)


par(mar=c(5,5,2,3))

plot(0, xlim=c(0,10), ylim=c(0,0.3), type="n",
  xlab="Serial interval (days)", ylab="Density", bty="n", axes=FALSE)

b <- data.frame(times=c(0:200)/10, des.prob=c(dweibull(0:200/10, rep(exp(serial.weibull.day0$par[1]),201),
   rep(exp(serial.weibull.day0$par[2]),201) )))
lines(b$times, b$des.prob, lty=1, col="black", lwd=1)

b <- data.frame(times=c(0:200)/10, des.prob=c(dweibull(0:200/10, rep(exp(serial.weibull.day1$par[1]),201),
   rep(exp(serial.weibull.day1$par[2]),201) )))
lines(b$times, b$des.prob, lty=2, col="black", lwd=1)

b <- data.frame(times=c(0:200)/10, des.prob=c(dweibull(0:200/10, rep(exp(serial.weibull.day2$par[1]),201),
   rep(exp(serial.weibull.day2$par[2]),201) )))
lines(b$times, b$des.prob, lty=3, col="black", lwd=1)

axis(1, pos=0, lwd=1, cex.axis=1, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1, las=1)


dev.off()

