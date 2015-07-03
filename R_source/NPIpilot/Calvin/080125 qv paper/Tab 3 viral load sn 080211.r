cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
#cdc <- read.csv("d:/r/080125 qv paper/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")
#source ("D:/r/080125 qv paper/080211_add_groups.r")

# Sn for flu A and B having qPCR result

cdc$logqPCR <- log(cdc$qPCR,10)
cdc$logqPCR[cdc$logqPCR == -Inf] <- 0

cdc$qPCRgp <-cut(cdc$logqPCR, breaks=c(-1 ,0, 4, 5 , 6, 10))
cdc$qPCRgp <- as.numeric(cdc$qPCRgp)
pcrava <- cdc[!is.na(cdc$qPCR),]

#function for outputing the result
sn <- function(gp, gpQVpos){
table(gpQVpos)
# handle only 1 case in table
 if (dim(table(gpQVpos)) <2) {
   new1 <- vector(length =2)
   new1[1] <- length(gpQVpos) - sum(gpQVpos) # total -ve 
   new1[2] <- sum(gpQVpos)    #total +ve

    sn <- binom.test(x= new1[2], dim(gp)[1]) 
 }
 
 else{
 sn <- binom.test(x=table(gpQVpos)[[2]], dim(gp)[1])
 }

 output <- data.frame(Total = dim(gp)[1], Value = c(sn[[1]] / sn[[2]]),
  lowerCI = c(sn[[4]][1]),
  upperCI = c(sn[[4]][2]),
  row.names=c("Sensitivity"))
 round(output, 2)
}

# for influenza A
pcrcultposA <- pcrava[pcrava$cultposA ==1,]
pcrcultposA1 <- pcrcultposA[pcrcultposA$qPCRgp==1,]
pcrcultposA2 <- pcrcultposA[pcrcultposA$qPCRgp==2,]
pcrcultposA3 <- pcrcultposA[pcrcultposA$qPCRgp==3,]
pcrcultposA4 <- pcrcultposA[pcrcultposA$qPCRgp==4,]
pcrcultposA5 <- pcrcultposA[pcrcultposA$qPCRgp==5,]

 a0<- sn(pcrcultposA, pcrcultposA$QVpos) # overall influenza A
 a1<- sn(pcrcultposA1, pcrcultposA1$QVpos)
 a2<- sn(pcrcultposA2, pcrcultposA2$QVpos)
 a3<- sn(pcrcultposA3, pcrcultposA3$QVpos)
 a4<- sn(pcrcultposA4, pcrcultposA4$QVpos)
 a5<- sn(pcrcultposA5, pcrcultposA5$QVpos)
 out <- rbind(a0,a1,a2,a3,a4,a5)
 row.names(out)<- c("Overall flu A PCR","Undetectable","<10^4 copies/ml",
 "10^4 to 10^5 copies/ml","10^5 to 10^6 copies/ml","> 10^6 copies/ml")
 out

#for influenza B
pcrcultposB <- pcrava[pcrava$cultposB ==1,] # overall influenza B
pcrcultposB1 <- pcrcultposB[pcrcultposB$qPCRgp==1,]
pcrcultposB2 <- pcrcultposB[pcrcultposB$qPCRgp==2,]
pcrcultposB3 <- pcrcultposB[pcrcultposB$qPCRgp==3,]
pcrcultposB4 <- pcrcultposB[pcrcultposB$qPCRgp==4,]
pcrcultposB5 <- pcrcultposB[pcrcultposB$qPCRgp==5,]

a0<- sn(pcrcultposB, pcrcultposB$QVpos)
a1<- sn(pcrcultposB1, pcrcultposB1$QVpos)
a2<- sn(pcrcultposB2, pcrcultposB2$QVpos)
a3<- sn(pcrcultposB3, pcrcultposB3$QVpos)
a4<- sn(pcrcultposB4, pcrcultposB4$QVpos)
a5<- sn(pcrcultposB5, pcrcultposB5$QVpos)
out <- rbind(a0,a1,a2,a3,a4,a5)
 row.names(out)<- c("Overall flu B PCR","Undetectable","<10^4 copies/ml",
 "10^4 to 10^5 copies/ml","10^5 to 10^6 copies/ml","> 10^6 copies/ml")
 out

#####################################
# for table 3a using Odd ratio
pcrcultposA$qPCRgp <- factor(pcrcultposA$qPCRgp, levels=c(3,1,2,4,5))

glm.fit <- glm( QVposA ~ qPCRgp, family = binomial, data=pcrcultposA) 

results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 2)

#### influenza B
pcrcultposB$qPCRgp <- factor(pcrcultposB$qPCRgp, levels=c(3,1,2,4,5))

glm.fit <- glm( QVposB ~ qPCRgp, family = binomial, data=pcrcultposB) 

results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 2)
#####################################

##########################################################
#for t test between influenza A and influenza B in the paper
#cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")

dat <-cdc[!is.na(cdc$qPCR),] 

dat$logqPCR <- log(dat$qPCR,10)
dat$logqPCR[dat$logqPCR == -Inf] <- log(900/2,10)  

#seperate into 3 groups, quickvue culture --, A-+, B-+, A++, B++
group1 <- dat[dat$QVres == 3 & dat$culture==0,] #both -ve
group2 <- dat[dat$QVres == 3 & dat$culture!=0,] # QV-ve , culture +ve
group3 <- dat[dat$QVres != 3 & dat$culture!=0,] # QV +ve , culture +ve

group2a <- group2[group2$culture =="A",] 
group2b <- group2[group2$culture =="B",]
group3a <- group3[group3$culture =="A",]
group3b <- group3[group3$culture =="B",]

groupa <-rbind (group2a,group3a)
groupb <-rbind (group2b,group3b)

## for data in result paragraph 
table(group1$qPCR)
mean (group2a$logqPCR)
mean (group3a$logqPCR)

mean (group2b$logqPCR)
mean (group3b$logqPCR)

# t test using 0 for undetectable  qPCR
groupa$QVres[groupa$QVres==3] <-0
t.test(groupa$logqPCR~groupa$QVres)

groupb$QVres[groupb$QVres==3] <-0
groupb$QVres[groupb$QVres==2] <-1
t.test(groupb$logqPCR~groupb$QVres)

#  t test using 900/2 for undetectable  qPCR
groupa$logqPCR[groupa$logqPCR == log(900/2,10) ] <- 0
groupb$logqPCR[groupb$logqPCR == log(900/2,10) ] <- 0
t.test(groupa$logqPCR~groupa$QVres)
t.test(groupb$logqPCR~groupb$QVres)

groupa$logqPCR[groupa$logqPCR == 0] <- log(900,10) 
groupb$logqPCR[groupb$logqPCR == 0] <- log(900,10) 
t.test(groupa$logqPCR~groupa$QVres)
t.test(groupb$logqPCR~groupb$QVres)