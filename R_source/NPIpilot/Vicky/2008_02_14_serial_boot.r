
visitday <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)
cdat <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_28_clinicdat_h.csv", header=TRUE)
lab <- read.csv("P:\\GROUP\\NPIpilot\\data\\lab2nd.csv", header=TRUE)
symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)

visitday <- merge(visitday,cdat[,c(2,26)],by="hhID",all.x=TRUE)   # hhID, onsettime

lab2nd <- data.frame(hhID=unique(lab$hhID[lab$labsedcase==1]))
sedcase <- lab[lab$labsedcase==1,c(1,2)]

onset <- merge(visitday[,c(1,21:26)],lab2nd,by="hhID",all.y=TRUE)
onset$v2_day <- as.numeric(onset$v2_day)

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
symptom_2nd <- merge(symptom_2nd,onset[,c(1,2)],by="hhID",all.x=TRUE)         # hhID, clinic_day
clinic_day <- symptom_2nd[symptom_2nd$day==0,c(1,24)]
symptom_2nd$day <- symptom_2nd$day + symptom_2nd$clinic_day

# set day 0 as index symtom onset, to obtain day of contact symtom onset
symptom_2nd$fever <- 1*(symptom_2nd$bodytemp >=37.8)
for (i in 1:nrow(symptom_2nd)){
      symptom_2nd$resp_symp[i] <- sum(symptom_2nd$sthroat[i],symptom_2nd$cough[i],
                                      symptom_2nd$pmuscle[i],symptom_2nd$headache[i],symptom_2nd$fever[i],na.rm=TRUE)
}
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07075"&symptom_2nd$day==0] <- 0
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
sedcase$clinic_day <- clinic_day$clinic_day
sedcase_1 <- merge(sedcase,symptom_2nd[symptom_2nd$mark_1==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_1)[4] <- "day_fever"
sedcase_2 <- merge(sedcase_1,symptom_2nd[symptom_2nd$mark_2==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_2)[5] <- "day_resp_symp"
sedcase_3 <- merge(sedcase_2,symptom_2nd[symptom_2nd$mark_3==1,c(1,3,4)],by=c("hhID","member"),all.x=T)
names(sedcase_3)[6] <- "day_any_symp"



# Fit the weibull model (allow for truncted serial interval)

weibull.loglik <- function(parameters, time, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    alpha <- exp(parameters[1])
    beta <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dweibull(time.g[i], shape=alpha, scale=beta) / (1-pweibull(trunct[i],shape=alpha, scale=beta)) )
    }
    -sum(log.lik)
    #log.lik
}

serial.weibull.fever <- optim(c(1,1), weibull.loglik, time=sedcase_3$day_fever,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

serial.weibull.resp <- optim(c(1,1), weibull.loglik, time=sedcase_3$day_resp_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

serial.weibull.any <- optim(c(1,1), weibull.loglik, time=sedcase_3$day_any_symp,data=sedcase_3,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

### naive bootstrap
#boot.weibull <- function(reps, time.b, data){
#  # function to resample the data and fit weibull models to each resample
#  n <- length(data[,1])
#  temp <- rep(NA, reps)
#  output <- data.frame(par1=temp, par2=temp, mean=temp)
#  for(i in 1:reps){
#    a <- sample(1:n, replace=TRUE)
#    temp.data <- data[a,]
#    temp.weibull <- optim(c(1,1), weibull.loglik, time=time.b,data=temp.data,
#                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
#    output$mean[i] <- exp(temp.weibull$par[2])*gamma(1+exp(-temp.weibull$par[1]))
#    output$par1[i] <- temp.weibull$par[1]
#    output$par2[i] <- temp.weibull$par[2]
#  }
#  output
#}
#
#bt <- 1000
#
#resp.wei.boot <- boot.weibull(bt, sedcase_3$day_resp_symp, sedcase_3)
#quantile(resp.wei.boot$mean, c(0.5, 0.025, 0.975))
#
#fever.wei.boot <- boot.weibull(bt, sedcase_3$day_fever, sedcase_3)
#quantile(fever.wei.boot$mean, c(0.5, 0.025, 0.975))

### parametric bootstrap
pboot.weibull <- function(reps, time.b, data){
   # function to draw iid samples from fitted model and fit weibull models to each resample
  n <- length(data[,1])
  temp <- rep(NA, reps)
  output <- data.frame(mean=temp, par1=temp, par2=temp)
  fitted.weibull <- optim(c(0.5,1.5), weibull.loglik, time=time.b,data=data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
  for(i in 1:reps){
    new.data <- data.frame(rep(NA,2*n))
    new.data[,1] <- rweibull(2*n, exp(fitted.weibull$par[1]), exp(fitted.weibull$par[2]))
    names(new.data)[1] <- "day_symp"
    new.data$clinic_day <- sample(c(0,1,2),2*n,replace=TRUE,prob=c(8,24,7))             # proportion based on sample!!!!!!
    new.data <- new.data[new.data$day_symp>=new.data$clinic_day,]
    new.data <- new.data[1:n,]
    
    new.weibull <- optim(c(0.5,1.5), weibull.loglik, time=new.data$day_symp,data=new.data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
    output$mean[i] <- exp(new.weibull$par[2])*gamma(1+exp(-new.weibull$par[1]))
    output$par1[i] <- new.weibull$par[1]
    output$par2[i] <- new.weibull$par[2]
  }
  output
}


bt <- 1000

any.wei.pboot <- pboot.weibull(bt, sedcase_3$day_any_symp, sedcase_3)     # reps=1000,c(0.5,1.5)
quantile(any.wei.pboot$mean, c(0.5, 0.025, 0.975))

resp.wei.pboot <- pboot.weibull(bt, sedcase_3$day_resp_symp, sedcase_3)  # reps=1000,c(1.0,1.5)
quantile(resp.wei.pboot$mean, c(0.5, 0.025, 0.975))

fever.wei.pboot <- pboot.weibull(bt, sedcase_3$day_fever, sedcase_3)     # reps=1000,c(0.5,1.5)
quantile(fever.wei.pboot$mean, c(0.5, 0.025, 0.975))


############### Plot the graph ######################

pdf("D:\\Work\\Influenza\\output NPI year 1\\serial interval\\fig.tiff\\serialboot.pdf", width=6, height=8, version="1.4")
#windows(width=6,height=4)
layout(matrix(1:2,ncol=1))


## resp symp
#
#par(mar=c(5,10,2,3))
#
#plot(0, xlim=c(0,8), ylim=c(0,0.4), type="n",
#  xlab=" ", ylab=" ", bty="n", axes=FALSE)
#
#
#for (i in 1:100){
#  b <- data.frame(times=c(0:200)/20, des.prob=c( 
#   dweibull(0:200/20, rep(exp(resp.wei.boot$par1[i]),201),
#   rep(exp(resp.wei.boot$par2[i]),201) )))
#   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
#}
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.resp$par[1]),201),
#   rep(exp(serial.weibull.resp$par[2]),201) )))
#lines(b$times, b$des.prob, lty=1, col="cyan", lwd=1.5)
#
#axis(1, pos=0, lwd=2, cex.axis=1.5, at=2*0:4)
#axis(2, pos=0, lwd=2, cex.axis=1.5, las=1)
#
#legend(4.5,0.4, c("Weibull","Bootstrap"), col=c("cyan","gray"),lty=c(1,1), bty="n",lwd=1.5,cex=1.3) 
#
#mtext("Density",side=2, line=9, cex=1.3, adj=0, las=1, at=0.2) 
#
#mtext("Serial interval (days) for resp symptoms",side=1, line=2, cex=1.3)
#mtext("(b1)",side=3,cex=1.3,line=-1)
#
## fever
#
#par(mar=c(5,10,2,3))
#
#plot(0, xlim=c(0,8), ylim=c(0,0.4), type="n",
#  xlab=" ", ylab=" ", bty="n", axes=FALSE)
#
#for (i in 1:100){
#  b <- data.frame(times=c(0:200)/20, des.prob=c( 
#   dweibull(0:200/20, rep(exp(fever.wei.boot$par1[i]),201),
#   rep(exp(fever.wei.boot$par2[i]),201) )))
#   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
#}
#
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.fever$par[1]),201),
#   rep(exp(serial.weibull.fever$par[2]),201) )))
#lines(b$times, b$des.prob, lty=1, col="cyan", lwd=1.5)
#
#axis(1, pos=0, lwd=2, cex.axis=1.5, at=2*0:4)
#axis(2, pos=0, lwd=2, cex.axis=1.5, las=1)
#
#legend(4.5,0.4, c("Weibull","Bootstrap"), col=c("cyan","gray"),lty=c(1,1), bty="n",lwd=1.5,cex=1.3) 
#
#mtext("Density",side=2, line=9, cex=1.3, adj=0, las=1, at=0.2) 
#
#mtext("Serial interval (days) for fever",side=1, line=2, cex=1.3)
#mtext("(c1)",side=3,cex=1.3,line=-1)

##-------------parametric bootstrap------------------##

# any resp symp

par(mar=c(5.5,5,3,1))

plot(0, xlim=c(0,10), ylim=c(0,0.4), type="n",
  xlab="Serial interval (days)", ylab="Density", bty="n", axes=FALSE)


for (i in 1:100){
  b <- data.frame(times=c(0:200)/20, des.prob=c( 
   dweibull(0:200/20, rep(exp(any.wei.pboot$par1[i]),201),
   rep(exp(any.wei.pboot$par2[i]),201) )))
   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
}

b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.any$par[1]),201),
   rep(exp(serial.weibull.any$par[2]),201) )))
lines(b$times, b$des.prob, lty=1, col="black", lwd=1.5)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)

#legend(6,0.4, c("Weibull","Para-boot"), col=c("black","gray"),lty=c(1,1), bty="n",lwd=1.5,cex=1.0) 

mtext("(a)",side=3,cex=0.9,line=1)

# 2 resp symp

par(mar=c(5.5,5,3,1))

plot(0, xlim=c(0,10), ylim=c(0,0.4), type="n",
  xlab="Serial interval (days)", ylab="Density", bty="n", axes=FALSE)


for (i in 1:100){
  b <- data.frame(times=c(0:200)/20, des.prob=c( 
   dweibull(0:200/20, rep(exp(resp.wei.pboot$par1[i]),201),
   rep(exp(resp.wei.pboot$par2[i]),201) )))
   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
}

b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.resp$par[1]),201),
   rep(exp(serial.weibull.resp$par[2]),201) )))
lines(b$times, b$des.prob, lty=1, col="black", lwd=1.5)

axis(1, pos=0, lwd=1, cex.axis=1, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1, las=1)

mtext("(b)",side=3,cex=0.9,line=1)

# fever
#
#par(mar=c(5.5,5,3,1))
#
#plot(0, xlim=c(0,10), ylim=c(0,0.4), type="n",
#  xlab="Serial interval (days)", ylab="Density", bty="n", axes=FALSE)
#
#for (i in 1:100){
#  b <- data.frame(times=c(0:200)/20, des.prob=c( 
#   dweibull(0:200/20, rep(exp(fever.wei.pboot$par1[i]),201),
#   rep(exp(fever.wei.pboot$par2[i]),201) )))
#   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
#}
#
#
#b <- data.frame(times=c(0:200)/20, des.prob=c(dweibull(0:200/20, rep(exp(serial.weibull.fever$par[1]),201),
#   rep(exp(serial.weibull.fever$par[2]),201) )))
#lines(b$times, b$des.prob, lty=1, col="black", lwd=1.5)
#
#axis(1, pos=0, lwd=1, cex.axis=1, at=2*0:5)
#axis(2, pos=0, lwd=1, cex.axis=1, las=1)
#
#mtext("(c)",side=3,cex=0.9,line=1)


dev.off()


