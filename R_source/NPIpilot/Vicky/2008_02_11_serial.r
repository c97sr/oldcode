source("D:\\Work\\Influenza\\programs\\2008_01_03\\2008_02_turnbull.r")

visitday <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)
cdat <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_28_clinicdat_h.csv", header=TRUE)
lab <- read.csv("P:\\GROUP\\NPIpilot\\data\\lab2nd2007.csv", header=TRUE)
symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)

visitday <- merge(visitday,cdat[,c(1,26)],by="hhID",all.x=TRUE)   # hhID, onsettime

lab2nd <- data.frame(hhID=unique(lab$hhID[lab$labsedcase==1]))
sedcase <- lab[lab$labsedcase==1,c(1,2)]

onset <- merge(visitday[,c(1,21:26)],lab2nd,by="hhID",all.y=TRUE)
#onset$v2_day <- as.numeric(onset$v2_day)

#### -------------------------------------------------------------------------------------------------------------------------------------

for (i in 1:nrow(onset)){
   if(onset$onsettime[i]==5) {
          onset$clinic_day[i] <- onset$clinic_day[i]-1
          onset$v1_day[i] <- onset$v1_day[i]-1
	  onset$v2_day[i] <- onset$v2_day[i]-1
	  onset$v3_day[i] <- onset$v3_day[i]-1
	  onset$v4_day[i] <- onset$v4_day[i]-1
   }
}

symptom_2nd <- symptom[symptom$hhID==sedcase$hhID[1]&symptom$member==sedcase$member[1],]
for (i in 2:nrow(sedcase)){
     tmp <- symptom[symptom$hhID==sedcase$hhID[i]&symptom$member==sedcase$member[i],]
     symptom_2nd <- rbind(symptom_2nd,tmp)
}
symptom_2nd <- merge(symptom_2nd,onset[,c(1,2,3)],by="hhID",all.x=TRUE)         # hhID, clinic_day -> v1_day
v1_day <- symptom_2nd[symptom_2nd$day==0,c(1,24,25)]
symptom_2nd$day <- symptom_2nd$day + symptom_2nd$v1_day

# set day 0 as index symtom onset, to obtain day of contact symtom onset
symptom_2nd$fever <- 1*(symptom_2nd$bodytemp >=37.8)
for (i in 1:nrow(symptom_2nd)){
      symptom_2nd$resp_symp[i] <- sum(symptom_2nd$sthroat[i],symptom_2nd$cough[i],
                                      symptom_2nd$pmuscle[i],symptom_2nd$headache[i],symptom_2nd$fever[i],na.rm=TRUE)
}
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07075"&symptom_2nd$day==1] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07117"&symptom_2nd$day==0] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07132"&symptom_2nd$day==1] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07132"&symptom_2nd$day==2] <- 0

symptom_2nd$mark_1 <- 0
for (j in 1:nrow(sedcase)){
     i <- (j-1)*10+1
     endc <- j*10
     while(i <= endc){
        if (!is.na(symptom_2nd$fever[i])& symptom_2nd$fever[i]==1) {
	symptom_2nd$mark_1[i] <- 1 
	break
	}
        else  i <- i+1
     }
}

symptom_2nd$mark_2 <- 0
for (j in 1:nrow(sedcase)){
     i <- (j-1)*10+1
     endc <- j*10
     while(i <= endc){
        if (symptom_2nd$resp_symp[i]>=2) {
	symptom_2nd$mark_2[i] <- 1 
	break
	}
        else  i <- i+1
     }
}

symptom_2nd$mark_3 <- 0
for (j in 1:nrow(sedcase)){
     i <- (j-1)*10+1
     endc <- j*10
     while(i <= endc){
        if (symptom_2nd$resp_symp[i]>=1) {
	symptom_2nd$mark_3[i] <- 1 
	break
	}
        else  i <- i+1
     }
}


# construct the table
#sedcase$v1_day <- v1_day$v1_day
sedcase$clinic_day <- v1_day$clinic_day
sedcase_1 <- merge(sedcase,symptom_2nd[symptom_2nd$mark_1==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_1)[4] <- "day_fever"
sedcase_2 <- merge(sedcase_1,symptom_2nd[symptom_2nd$mark_2==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_2)[5] <- "day_resp_symp"
sedcase_3 <- merge(sedcase_2,symptom_2nd[symptom_2nd$mark_3==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_3)[6] <- "day_any_symp"


data.fever <- sedcase_3[,1:3]
names(data.fever)[3] <- "trunct"
data.fever$timeL <- sedcase_3$day_fever
data.fever$timeR <- sedcase_3$day_fever
data.fever$exact <- 1
data.fever <- data.fever[!is.na(data.fever$timeL),]
fever.turnbull <- make.turnbull(data.fever)

data.resp <- sedcase_3[,1:3]
names(data.resp)[3] <- "trunct"
data.resp$timeL <- sedcase_3$day_resp_symp
data.resp$timeR <- sedcase_3$day_resp_symp
data.resp$exact <- 1
data.resp <- data.resp[!is.na(data.resp$timeL),]
resp.turnbull <- make.turnbull(data.resp)

data.any <- sedcase_3[,1:3]
names(data.any)[3] <- "trunct"
data.any$timeL <- sedcase_3$day_any_symp
data.any$timeR <- sedcase_3$day_any_symp
data.any$exact <- 1
data.any <- data.any[!is.na(data.any$timeL),]
any.turnbull <- make.turnbull(data.any)

# Fit the lognormal/gamma/weibull model (allow for truncted serial interval)

gamma.loglik <- function(parameters, time,data){               # time <- sedcase_3$day_fever/day_resp_symp/day_any_symp
    if(length(parameters)!=2) stop("gamma distribution should have two parameters")
    k <- exp(parameters[1])
    lambda <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dgamma(time.g[i], shape=k,scale=1/lambda) / (1-pgamma(trunct[i], shape=k,scale=1/lambda)) )
    }
    -sum(log.lik)
    #log.lik
}

serial.gamma.fever <- optim(c(1,1), gamma.loglik, time=sedcase_3$day_fever,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.gamma.fever

serial.gamma.resp <- optim(c(1,1), gamma.loglik, time=sedcase_3$day_resp_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.gamma.resp

serial.gamma.any <- optim(c(0.5,0.5), gamma.loglik, time=sedcase_3$day_any_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.gamma.any


# the mean
exp(serial.gamma.fever$par[1])/exp(serial.gamma.fever$par[2])
exp(serial.gamma.resp$par[1])/exp(serial.gamma.resp$par[2])
exp(serial.gamma.any$par[1])/exp(serial.gamma.any$par[2])

# the variance
exp(serial.gamma.fever$par[1])/exp(2*serial.gamma.fever$par[2])
exp(serial.gamma.resp$par[1])/exp(2*serial.gamma.resp$par[2])
exp(serial.gamma.any$par[1])/exp(2*serial.gamma.any$par[2])


lognorm.loglik <- function(parameters, time, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    mean <- exp(parameters[1])
    sd <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dlnorm(time.g[i], meanlog=mean, sdlog=sd) / (1-plnorm(trunct[i],meanlog=mean, sdlog=sd)) )
    }
    -sum(log.lik)
    #log.lik
}

serial.lognorm.fever <- optim(c(1,1), lognorm.loglik, time=sedcase_3$day_fever,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.lognorm.fever

serial.lognorm.resp <- optim(c(1,0.5), lognorm.loglik, time=sedcase_3$day_resp_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.lognorm.resp

serial.lognorm.any <- optim(c(0.5,-1.5), lognorm.loglik, time=sedcase_3$day_any_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.lognorm.any


# the mean
exp(exp(serial.lognorm.fever$par[1])+exp(2*serial.lognorm.fever$par[2])/2)
exp(exp(serial.lognorm.resp$par[1])+exp(2*serial.lognorm.resp$par[2])/2)
exp(exp(serial.lognorm.any$par[1])+exp(2*serial.lognorm.any$par[2])/2)

# the variance
exp(exp(2*serial.lognorm.fever$par[2])+2*exp(serial.lognorm.fever$par[1]))*(exp(exp(2*serial.lognorm.fever$par[2]))-1)
exp(exp(2*serial.lognorm.resp$par[2])+2*exp(serial.lognorm.resp$par[1]))*(exp(exp(2*serial.lognorm.resp$par[2]))-1)
exp(exp(2*serial.lognorm.any$par[2])+2*exp(serial.lognorm.any$par[1]))*(exp(exp(2*serial.lognorm.any$par[2]))-1)




weibull.loglik <- function(parameters, time, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    alpha <- exp(parameters[1])
    beta <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    #trunct <- rep(0,length(trunct))
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dweibull(time.g[i], shape=alpha, scale=beta) / (1-pweibull(trunct[i],shape=alpha, scale=beta)) )
    }
    -sum(log.lik)
    #log.lik
}

serial.weibull.fever <- optim(c(1,1), weibull.loglik, time=sedcase_3$day_fever,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.weibull.fever

serial.weibull.resp <- optim(c(1,1), weibull.loglik, time=sedcase_3$day_resp_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.weibull.resp

serial.weibull.any <- optim(c(0.5,1.5), weibull.loglik, time=sedcase_3$day_any_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#serial.weibull.any

# the mean
exp(serial.weibull.fever$par[2])*gamma(1+exp(-serial.weibull.fever$par[1]))
exp(serial.weibull.resp$par[2])*gamma(1+exp(-serial.weibull.resp$par[1]))
exp(serial.weibull.any$par[2])*gamma(1+exp(-serial.weibull.any$par[1]))

# the variance
exp(2*serial.weibull.fever$par[2])*(gamma(1+2*exp(-serial.weibull.fever$par[1]))-(gamma(1+exp(-serial.weibull.fever$par[1])))^2)
exp(2*serial.weibull.resp$par[2])*(gamma(1+2*exp(-serial.weibull.resp$par[1]))-(gamma(1+exp(-serial.weibull.resp$par[1])))^2)
exp(2*serial.weibull.any$par[2])*(gamma(1+2*exp(-serial.weibull.any$par[1]))-(gamma(1+exp(-serial.weibull.any$par[1])))^2)

############### Plot the graph ######################

pdf("D:\\Work\\Influenza\\output NPI year 1\\serial interval\\fig.tiff\\serial.pdf", width=6, height=8, version="1.4")
#windows(width=6,height=8)
layout(matrix(1:2,ncol=1))

## any resp symp
#
#par(mar=c(5,10,2,3))
#
#plot(0, xlim=c(0,8), ylim=c(0,1), type="n",
#  xlab=" ", ylab=" ", bty="n", axes=FALSE)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dlnorm(0:200/20,
#  meanlog=rep(exp(serial.lognorm.any$par[1]),201),
#  sdlog=rep(exp(serial.lognorm.any$par[2]),201))))
#lines(b$times, b$des.prob, lty=2, col="blue", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dgamma(0:200/20,rep(exp(serial.gamma.any$par[1]),201),
#  rep(exp(serial.gamma.any$par[2]),201) )) ) 
#lines(b$times, b$des.prob, lty=1, col="green", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.any$par[1]),201),
#   rep(exp(serial.weibull.any$par[2]),201) )))
#lines(b$times, b$des.prob, lty=4, col="cyan", lwd=3)
#
#axis(1, pos=0, lwd=2, cex.axis=1.5, at=2*0:4)
#axis(2, pos=0, lwd=2, cex.axis=1.5, las=1)
#
#legend(4,0.9, c("Log normal", "Gamma", "Weibull"), col=c("blue", "green","cyan"),
#  lty=c(2,1,4), bty="n",lwd=2,cex=1.3)      
#
#mtext("Density",
# side=2, line=9, cex=1.3, adj=0, las=1, at=0.5) 
#
#mtext("Serial interval (days)",
# side=1, line=2, cex=1.3)
#mtext("(a)",side=3,cex=1.3,line=-1)
#
## 2 resp symp
#
#par(mar=c(5,10,2,3))
#
#plot(0, xlim=c(0,8), ylim=c(0,1), type="n",
#  xlab=" ", ylab=" ", bty="n", axes=FALSE)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dlnorm(0:200/20,
#  meanlog=rep(exp(serial.lognorm.resp$par[1]),201),
#  sdlog=rep(exp(serial.lognorm.resp$par[2]),201))))
#lines(b$times, b$des.prob, lty=2, col="blue", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dgamma(0:200/20,rep(exp(serial.gamma.resp$par[1]),201),
#  rep(exp(serial.gamma.resp$par[2]),201) )) ) 
#lines(b$times, b$des.prob, lty=1, col="green", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.resp$par[1]),201),
#   rep(exp(serial.weibull.resp$par[2]),201) )))
#lines(b$times, b$des.prob, lty=4, col="cyan", lwd=3)
#
#axis(1, pos=0, lwd=2, cex.axis=1.5, at=2*0:4)
#axis(2, pos=0, lwd=2, cex.axis=1.5, las=1)
#
#legend(4,0.9, c("Log normal", "Gamma", "Weibull"), col=c("blue", "green","cyan"),
#  lty=c(2,1,4), bty="n",lwd=2,cex=1.3)      
#
#mtext("Density",
# side=2, line=9, cex=1.3, adj=0, las=1, at=0.5) 
#
#mtext("Serial interval (days)",
# side=1, line=2, cex=1.3)
#mtext("(c)",side=3,cex=1.3,line=-1)
#
## fever
#
#par(mar=c(5,10,2,3))
#
#plot(0, xlim=c(0,8), ylim=c(0,1), type="n",
#  xlab=" ", ylab=" ", bty="n", axes=FALSE)
#
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dlnorm(0:200/20,
#  meanlog=rep(exp(serial.lognorm.fever$par[1]),201),
#  sdlog=rep(exp(serial.lognorm.fever$par[2]),201))))
#lines(b$times, b$des.prob, lty=2, col="blue", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dgamma(0:200/20,rep(exp(serial.gamma.fever$par[1]),201),
#  rep(exp(serial.gamma.fever$par[2]),201) )) ) 
#lines(b$times, b$des.prob, lty=1, col="green", lwd=3)
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.fever$par[1]),201),
#   rep(exp(serial.weibull.fever$par[2]),201) )))
#lines(b$times, b$des.prob, lty=4, col="cyan", lwd=3)
#
#axis(1, pos=0, lwd=2, cex.axis=1.5, at=2*0:4)
#axis(2, pos=0, lwd=2, cex.axis=1.5, las=1)
#
#legend(4,0.9, c("Log normal", "Gamma", "Weibull"), col=c("blue", "green","cyan"),
#  lty=c(2,1,4), bty="n",lwd=2,cex=1.3)     
#
#mtext("Density",
# side=2, line=9, cex=1.3, adj=0, las=1, at=0.5) 
#
#mtext("Serial interval (days)",
# side=1, line=2, cex=1.3)
#mtext("(e)",side=3,cex=1.3,line=-1)


#-------------------------------------------#

# any resp symp

par(mar=c(5.5,5,3,1))

plot(0, xlim=c(0,10), ylim=c(0,1), type="n",
  xlab="Serial interval (days)", ylab="Cumulative proportion", bty="n", axes=FALSE)

#turnbull.plot(any.turnbull, plot.S=T, plot.CI=T, line.type=1,
#line.col="red", line.width=3, max.time=25, CI.col=rgb(1,0,0,0.5))

turnbull.plot(any.turnbull, plot.S=T, plot.CI=F, line.type=3,
line.col="black", line.width=1, max.time=10)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0, plnorm(1:100/10,
  meanlog=rep(exp(serial.lognorm.any$par[1]),100),
  sdlog=rep(exp(serial.lognorm.any$par[2]),100))))
lines(b$times, b$cum.prob, lty=4, col="black", lwd=1)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0,
  pgamma(1:100/10,rep(exp(serial.gamma.any$par[1]),100),
  rep(exp(serial.gamma.any$par[2]),100) )) ) 
lines(b$times, b$cum.prob, lty=2, col="black", lwd=1)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0, 
   pweibull(1:100/10, rep(exp(serial.weibull.any$par[1]),100),
   rep(exp(serial.weibull.any$par[2]),100) )))
lines(b$times, b$cum.prob, lty=1, col="black", lwd=1)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)

#legend(12,0.9, c("Log normal", "Gamma", "Weibull"), col=c("blue", "green","cyan"),
#  lty=c(2,1,4), bty="n",lwd=2,cex=1.3)      

mtext("(a)",side=3,cex=0.9,line=0.5)

# 2 resp symp

par(mar=c(5.5,5,3,1))

plot(0, xlim=c(0,10), ylim=c(0,1), type="n",
  xlab="Serial interval (days)", ylab="Cumulative proportion", bty="n", axes=FALSE)

turnbull.plot(resp.turnbull, plot.S=T, plot.CI=F, line.type=3,
line.col="black", line.width=1, max.time=10)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0, plnorm(1:100/10,
  meanlog=rep(exp(serial.lognorm.resp$par[1]),100),
  sdlog=rep(exp(serial.lognorm.resp$par[2]),100))))
lines(b$times, b$cum.prob, lty=4, col="black", lwd=1)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0,
  pgamma(1:100/10,rep(exp(serial.gamma.resp$par[1]),100),
  rep(exp(serial.gamma.resp$par[2]),100) )) ) 
lines(b$times, b$cum.prob, lty=2, col="black", lwd=1)

b <- data.frame(times=c(0:100)/10, cum.prob=c(0, 
   pweibull(1:100/10, rep(exp(serial.weibull.resp$par[1]),100),
   rep(exp(serial.weibull.resp$par[2]),100) )))
lines(b$times, b$cum.prob, lty=1, col="black", lwd=1)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)

mtext("(b)",side=3,cex=0.9,line=0.5)

## fever
#
#par(mar=c(5.5,5,3,1))
#
#plot(0, xlim=c(0,10), ylim=c(0,1), type="n",
#  xlab="Serial interval (days)", ylab="Cumulative proportion", bty="n", axes=FALSE)
#
#turnbull.plot(fever.turnbull, plot.S=T, plot.CI=F, line.type=3,
#line.col="black", line.width=1, max.time=10)
#
#b <- data.frame(times=c(0:100)/10, cum.prob=c(0, plnorm(1:100/10,
#  meanlog=rep(exp(serial.lognorm.fever$par[1]),100),
#  sdlog=rep(exp(serial.lognorm.fever$par[2]),100))))
#lines(b$times, b$cum.prob, lty=4, col="black", lwd=1)
#
#b <- data.frame(times=c(0:100)/10, cum.prob=c(0,
#  pgamma(1:100/10,rep(exp(serial.gamma.fever$par[1]),100),
#  rep(exp(serial.gamma.fever$par[2]),100) )) ) 
#lines(b$times, b$cum.prob, lty=2, col="black", lwd=1)
#
#b <- data.frame(times=c(0:100)/10, cum.prob=c(0, 
#   pweibull(1:100/10, rep(exp(serial.weibull.fever$par[1]),100),
#   rep(exp(serial.weibull.fever$par[2]),100) )))
#lines(b$times, b$cum.prob, lty=1, col="black", lwd=1)
#
#axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
#axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)
#
#
#mtext("(c)",side=3,cex=0.9,line=0.5)


dev.off()


