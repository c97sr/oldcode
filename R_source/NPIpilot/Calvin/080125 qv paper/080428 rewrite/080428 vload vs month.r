#080428 rewrite analysis

cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

qPCR <- cdc[!is.na(cdc$qPCR),]
qPCR$logqPCR <- log(qPCR$qPCR,10)
qPCR$logqPCR[qPCR$logqPCR ==-Inf] <-0

ci <- function( x ){
  n <- length(x)
  stderror <- function(x){ sqrt(var(x)/length(x))}
  mean_v <- mean(x)
  stderror_v <- stderror(x)
  upperCI <- mean_v + (1.96*stderror_v)
  lowerCI <- mean_v - (1.96*stderror_v)
  y <- cbind(n, mean_v,lowerCI,upperCI)
  
  return(round(y,2))
}

pcr_v <- ci(qPCR$qPCR[qPCR$month ==2])
for (i in 3:9){
temp <- ci(qPCR$qPCR[qPCR$month ==i])
pcr_v <- rbind(pcr_v,temp)
}

#logpcr_v <- log(pcr_v ,10)

#
# 1st window - QV sn and sp
par(mar=c(4,5,3,1))
plot(0, type="n", axes=FALSE, xlim=c(0,8), ylim=c(0,10), xlab=" ", ylab="log PCR",las=1,cex.lab=1.5)
axis(1, pos=0, at=c(1:8), rep(NA, 8))
axis(2, pos=0, at=0:5*2, labels=c("0","2","4","6","8","10"), las=1,cex.axis=1.5)

lines(type= "o",  as.matrix(log(pcr_v[1:8,2],10)), pch=16,lty=1,cex=1.6)
#for (i in 1:12){
#  segments(i-0.55, sn[i,2], i-0.55, sn[i,3])
#}

#lines(type= "o", 0.55:11.55 , as.matrix(sp[1:12,1]), pch=17,lty=2,cex=1.6)
#for (i in 1:12){
#  segments(i-0.45, sp[i,2], i-0.45, sp[i,3])
#}

#legend(7.8, 1.0, cex=1.2,
#  legend=c("Specificity", "Sensitivity"), pch = c(17,16) , lty=c(2,1))

mtext("(a)", side=3,line=1,at=4.5)
#