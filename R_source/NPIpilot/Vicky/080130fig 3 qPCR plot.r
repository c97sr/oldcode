# boxplot of pcr results for 3 groups, QV culture --, -+ , ++
#get data

#dat <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvPCRdat.csv")
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")

dat <-cdc[!is.na(cdc$qPCR),]

dat$logqPCR <- log(dat$qPCR,10)
dat$logqPCR[dat$logqPCR == -Inf] <- log(70,10)
dat$logqPCR[dat$logqPCR < log(900,10) ] <- log(70,10)


#seperate into 3 groups, quickvue culture --, A-+, B-+, A++, B++
group0a <- dat[dat$QVres == 1 & dat$culture==0,]
group0b <- dat[dat$QVres == 2 & dat$culture==0,]
group1 <- dat[dat$QVres == 3 & dat$culture==0,]
group2 <- dat[dat$QVres == 3 & dat$culture!=0,]
group3 <- dat[dat$QVres != 3 & dat$culture!=0,]

group2a <- group2[group2$culture =="A",]
group2b <- group2[group2$culture =="B",]
group3a <- group3[group3$culture =="A",]
group3b <- group3[group3$culture =="B",]

groupa <-rbind (group2a,group3a)
groupb <-rbind (group2b,group3b)


lower.lim <- log(900,10)


################



windows(width=10, height=7)
par(mar=c(5.5,4,1,2))
plot(0, type="n", axes=FALSE, xlim=c(-1,14), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1,c(0.5,2,3.5, 6,7.5, 10,11.5),labels=rep("", 7), cex.axis=1.5)         
axis(2, pos=0, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)
mtext("Viral load (copies/ml)",side=2,line=1)

set.seed(100)

points(jitter(rep(0.5, length(group1$logqPCR)), factor=1), pch=16,group1$logqPCR)
points(jitter(rep(2, length(group0a$logqPCR)), factor=1), pch=16,group0a$logqPCR)
points(jitter(rep(3.5, length(group0b$logqPCR)), factor=1), pch=16,group0b$logqPCR)
points(jitter(rep(6, length(group2a$logqPCR)), factor=1), pch=16,group2a$logqPCR)
points(jitter(rep(7.5, length(group3a$logqPCR)), factor=1), pch=16,group3a$logqPCR)
points(jitter(rep(10, length(group2b$logqPCR)), factor=1), pch=16,group2b$logqPCR)
points(jitter(rep(11.5, length(group3b$logqPCR)), factor=1), pch=16,group3b$logqPCR)

mtext("QuickVue result:",side=1,at=-1.5, line=1.0)
mtext("Culture result:",side=1,at=-1.7, line=2.5)
mtext("Neg" ,side=1,at=0.5, line=1)
mtext("Pos A" ,side=1,at=2, line=1)
mtext("Pos B" ,side=1,at=3.5, line=1)
mtext("Neg" ,side=1,at=6, line=1)
mtext("Pos A" ,side=1,at=7.5, line=1)
mtext("Neg" ,side=1,at=10, line=1)
mtext("Pos B",side=1,at=11.5, line=1)

mtext("Negative",  side=1, at=2, line=2.5)
mtext("Influenza A", side=1, at=6.75, line=2.5)
mtext("Influenza B", side=1, at=10.75, line=2.5)


lines(c(0,15),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=13.5, line=-6)


########################################################
#for median and IQR of each group
mediangp <-vector(length=4)
lowiqrgp<-vector(length=4)
upiqrgp <-vector(length=4)
mediangp[1] <- median(group2a$logqPCR)
mediangp[2] <- median(group2b$logqPCR)
mediangp[3] <- median(group3a$logqPCR)
mediangp[4] <- median(group3b$logqPCR)

lowiqrgp[1] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[1] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[2] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[2] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[3] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[3] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[4] <-quantile (group3b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[4] <-quantile (group3b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]



lines(c(6.35,6.45), c(mediangp[1],mediangp[1]),col="red")
lines(c(6.4,6.4),c(lowiqrgp[1],upiqrgp[1]),col="red")

lines(c(7.85,7.95), c(mediangp[3],mediangp[3]),col="red")
lines(c(7.9,7.9),c(lowiqrgp[3],upiqrgp[3]),col="red")

lines(c(10.35,10.45), c(mediangp[2],mediangp[2]),col="red")
lines(c(10.4,10.4),c(lowiqrgp[2],upiqrgp[2]),col="red")

lines(c(11.85,11.95), c(mediangp[4],mediangp[4]),col="red")
lines(c(11.9,11.9),c(lowiqrgp[4],upiqrgp[4]),col="red")



