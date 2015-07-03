# figure 2
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/070801 qv function.r")
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

sn <- data.frame(v=rep(NA,8), lower=rep(NA,8), upper=rep(NA,8))
sp <- data.frame(v=rep(NA,8), lower=rep(NA,8), upper=rep(NA,8))
ppv <- data.frame(v=rep(NA,8), lower=rep(NA,8), upper=rep(NA,8))
npv <- data.frame(v=rep(NA,8), lower=rep(NA,8), upper=rep(NA,8))

for (i in 2:8){
sn[i,1] <- snspresult (cdc, cdc$month,i) [[1]][1]
sp[i,1] <- snspresult (cdc, cdc$month,i) [[1]][2]
ppv[i,1] <- snspresult (cdc, cdc$month,i) [[1]][3]
npv[i,1] <- snspresult (cdc, cdc$month,i) [[1]][4]

sn[i,2] <- snspresult (cdc, cdc$month,i) [[2]][1]
sp[i,2] <- snspresult (cdc, cdc$month,i) [[2]][2]
ppv[i,2] <- snspresult (cdc, cdc$month,i) [[2]][3]
npv[i,2] <- snspresult (cdc, cdc$month,i) [[2]][4]

sn[i,3] <- snspresult (cdc, cdc$month,i) [[3]][1]
sp[i,3] <- snspresult (cdc, cdc$month,i) [[3]][2]
ppv[i,3] <- snspresult (cdc, cdc$month,i) [[3]][3]
npv[i,3] <- snspresult (cdc, cdc$month,i) [[3]][4]
}


# for plots side by side
windows(width=13, height=15)
layout(matrix(1:3, ncol=1))


#
# 1st window - QV sn and sp
par(mar=c(4,5,3,1))
plot(0, type="n", axes=FALSE, xlim=c(0,10), ylim=c(0,1), xlab=" ", ylab="Sensitivity/specificity",las=1,cex.lab=1.5)
axis(1, pos=0, at=c(0:9), rep(NA, 10))
axis(2, pos=0, at=0:5*0.2, labels=c("0","0.2","0.4","0.6","0.8","1.0"), las=1,cex.axis=1.5)

lines(type= "o", 0.45:11.45 , as.matrix(sn[1:12,1]), pch=16,lty=1,cex=1.6)
for (i in 1:12){
  segments(i-0.55, sn[i,2], i-0.55, sn[i,3])
}

lines(type= "o", 0.55:11.55 , as.matrix(sp[1:12,1]), pch=17,lty=2,cex=1.6)
for (i in 1:12){
  segments(i-0.45, sp[i,2], i-0.45, sp[i,3])
}

legend(7.8, 1.0, cex=1.2,
  legend=c("Specificity", "Sensitivity"), pch = c(17,16) , lty=c(2,1))

mtext("(a)", side=3,line=1,at=4.5)
#
# 2nd window - qv ppv

par(mar=c(4,5,3,2))
plot(0, type="n", axes=FALSE, xlim=c(0,10), ylim=c(0,1), xlab=" ", ylab="Predictive value",las=1,cex.lab=1.5)
axis(1, pos=0, at=c(0:9), rep(NA, 10))
axis(2, pos=0, at=0:5*0.2, labels=c("0","0.2","0.4","0.6","0.8","1.0"), las=1, cex.axis=1.5)

lines(type= "o", 0.45:11.45 , as.matrix(ppv[1:12,1]), pch=16,lty=1,cex=1.6)
for (i in 1:12){
  segments(i-0.55, ppv[i,2], i-0.55, ppv[i,3])
}

lines(type= "o", 0.55:11.55 , as.matrix(npv[1:12,1]), pch=17,lty=2,cex=1.6)
for (i in 1:12){
  segments(i-0.45, npv[i,2], i-0.45, npv[i,3])
}

legend(7.8, 1.0, cex=1.2,
  legend=c("PPV", "NPV"), pch = c(16,17) , lty=c(1,2))

mtext("(b)", side=3,line=1,at=4.5)


dat <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/rawdata/monthly qm rate.csv")
#
# 3rd window - queen mary Lab surveillance data

par(mar=c(4,5,3,1))

barplot(height=t(as.matrix(dat[1:12,4:5])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,10), ylim=c(0,0.3), xlab=" ", ylab="Isolation rate", cex.lab=1.5,
 col=c(gray(0.2),0)

#legend.text=c("Flu A +ve", "Flu B +ve")
) 

axis(1, pos=0, at=c(0:9), rep(NA, 10), , cex.axis=1.2)
axis(2, pos=0, at=0:3/10, labels=c("0","0.1","0.2","0.3"), las=1, cex.axis=1.2)

legend(7.8, 0.3, cex=1.2,
  legend=c("Influenza B", "Influenza A"), fill=c(0, gray(0.2)), lty=NULL)

mtext("(c)",  side=3,line=1,at=4.5)


mtext("J", side=1, at=0.5, line=1, cex=1)
mtext("F", side=1, at=1.5, line=1, cex=1)
mtext("M", side=1, at=2.5, line=1, cex=1)
mtext("A", side=1, at=3.5, line=1, cex=1)
mtext("M", side=1, at=4.5, line=1, cex=1)
mtext("J", side=1, at=5.5, line=1, cex=1)
mtext("J", side=1, at=6.5, line=1, cex=1)
mtext("A", side=1, at=7.5, line=1, cex=1)
mtext("S", side=1, at=8.5, line=1, cex=1)

