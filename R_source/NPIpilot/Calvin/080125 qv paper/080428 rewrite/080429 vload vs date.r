# appendix fig 2, viral load against time
require(chron)
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")
dat <-cdc
dat <-dat[!is.na(dat$qPCR),]

dat$logqPCR <- log(dat$qPCR,10)
dat$logqPCR[dat$logqPCR == -Inf] <- log(70,10)
dat$logqPCR[dat$logqPCR < log(900,10) ] <- log(70,10)

dat$startd <-dates(as.character(dat$date), format= "m/d/y") - dates("1/1/2007",format= "m/d/y")

dat.a <- dat[dat$culture =="A",]
dat.b <- dat[dat$culture =="B",]

#dat$date <- dates(as.character(dat$date), format= "m/d/y")

month <- vector(length=10)
month[1] <- dates("1/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[2] <-dates("2/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[3] <-dates("3/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[4] <-dates("4/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[5] <-dates("5/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[6] <-dates("6/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[7] <-dates("7/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[8] <-dates("8/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[9] <-dates("9/1/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")
month[10] <-dates("9/30/2007", format= "m/d/y") - dates("1/1/2007",format= "m/d/y")

#influenza A
childa <- dat.a[dat.a$isadult ==0,]
adulta <- dat.a[dat.a$isadult ==1,]

windows(width=12, height=15)
layout(matrix(1:2, nrow=2))

par(mar=c(4,7,1,4))
plot(0, type="n", axes=FALSE, xlim=c(0,272), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1, month, labels=rep(NA,10))
#
axis(2, pos=1, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)

points(adulta$startd, adulta$logqPCR, pch=19)
points(childa$startd, childa$logqPCR, pch=1)

#lines(  ksmooth(adulta$startd, adulta$logqPCR, "normal",bandwidth=10), type="l")
#lines(  ksmooth(adulta$startd, adulta$logqPCR, "normal",bandwidth=30), type="l")
#lines(  ksmooth(childa$startd, childa$logqPCR, "normal",bandwidth=30), type="l",lty=2)
#lines(  ksmooth(adulta$startd, adulta$logqPCR, "normal",bandwidth=30), type="l")

#lines(  smooth.spline(adulta$startd, adulta$logqPCR, spar=.1), type="l")
#lines(  smooth.spline(adulta$startd, adulta$logqPCR, spar=.5), type="l")
#lines(  smooth.spline(adulta$startd, adulta$logqPCR, spar=1), type="l")

mtext("Viral load (copies/ml)",side=2,at=6,line=4)

#mtext("J", side=1, at=15, line=1, cex=1)
#mtext("F", side=1, at=45, line=1, cex=1)
#mtext("M", side=1, at=75, line=1, cex=1)
#mtext("A", side=1, at=105, line=1, cex=1)
#mtext("M", side=1, at=135, line=1, cex=1)
#mtext("J", side=1, at=165, line=1, cex=1)
#mtext("J", side=1, at=195, line=1, cex=1)
#mtext("A", side=1, at=225, line=1, cex=1)
#mtext("S", side=1, at=255, line=1, cex=1)

mtext("Jan 1", side=1, at=month[1], line=1, cex=0.8)
#mtext("Feb", side=1, at=45, line=1, cex=0.8)
mtext("Mar 1", side=1, at=month[3], line=1, cex=0.8)
#mtext("Apr", side=1, at=105, line=1, cex=0.8)
mtext("May 1", side=1, at=month[5], line=1, cex=0.8)
#mtext("Jun", side=1, at=165, line=1, cex=0.8)
mtext("Jul 1", side=1, at=month[7], line=1, cex=0.8)
#mtext("Aug", side=1, at=225, line=1, cex=0.8)
mtext("Sep 1", side=1, at=month[9], line=1, cex=0.8)


lower.lim <- log(900,10)
lines(c(0,280),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=290, line=-4)
mtext("(a)",side=3)

#influenza b
childb <- dat.b[dat.b$isadult ==0,]
adultb <- dat.b[dat.b$isadult ==1,]

par(mar=c(4,7,1,4))
plot(0, type="n", axes=FALSE, xlim=c(0,272), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1, month, labels=rep("", 10))

#
axis(2, pos=1, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)

points(adultb$startd, adultb$logqPCR, pch=19)
points(childb$startd, childb$logqPCR, pch=1)

#lines(  ksmooth(adultb$startd, adultb$logqPCR, "normal",bandwidth=30), type="l")
#lines(  ksmooth(childb$startd, childb$logqPCR, "normal",bandwidth=30), type="l",lty=2)

mtext("Viral load (copies/ml)",side=2,at=6,line=4)


mtext("Jan 1", side=1, at=month[1], line=1, cex=0.8)
#mtext("Feb", side=1, at=45, line=1, cex=0.8)
mtext("Mar 1", side=1, at=month[3], line=1, cex=0.8)
#mtext("Apr", side=1, at=105, line=1, cex=0.8)
mtext("May 1", side=1, at=month[5], line=1, cex=0.8)
#mtext("Jun", side=1, at=165, line=1, cex=0.8)
mtext("Jul 1", side=1, at=month[7], line=1, cex=0.8)
#mtext("Aug", side=1, at=225, line=1, cex=0.8)
mtext("Sep 1", side=1, at=month[9], line=1, cex=0.8)


lower.lim <- log(900,10)
lines(c(0,280),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=295, line=-4)
mtext("(b)",side=3)

