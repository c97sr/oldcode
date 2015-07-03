plots <- read.csv("P:/GROUP/NPIpilot/data/surveillance/07_10_24_qvchpqmrate.csv", header=TRUE)
plots <- plots[1:40,]
#plots <-na.omit(plots)

months <- rep(0,40)
months[1] <- 1
for(i in 2:40){
    if(plots$month[i]!=plots$month[i-1]) months[i] <- 1
}
labs <- c(sort(unique(months*1:40))[-1], 40)

# for plots side by side
windows(width=8, height=7)
layout(matrix(1:3, ncol=1))
par(mar=c(3,5,3,1))

#
# 1st window - QV recruitment data

barplot(height=t(as.matrix(plots[1:40,9:10])), space=0,
 axes=FALSE, axisnames=FALSE, col=c(gray(0.2),0),
# main="No. of QV tests done per working day",
 xlim=c(0,40), ylim=c(0,12), xlab=" ", ylab="Recruitment rate", 
# legend.text=c("QV +ve", "QV -ve")
 cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,3,6,9,12) ,  labels=c("0","3","6","9","12"), las=1, cex.axis=1.2)
legend(33, 11, cex=1,
  legend=c("QV -ve", "QV +ve"), fill=c(0, gray(0.2)), lty=NULL)

mtext("(a)", cex=1)

#
# 2nd window - chp surveillance data, GOPC and GP rate

par(mar=c(3,5,3,1))

plot(as.matrix(plots[1:40,18]), axes=FALSE, type="l",
# main = "CHP surveillance data, GP rate per 1000 patients", 
 xlim=c(0,40), ylim=c(0,80), xlab=" ", ylab="ILI rate", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:4*20, labels=c("0","20","40","60","80"), las=1, cex.axis=1.2)

mtext("(b)", cex=1)
#
# 3rd window - queen mary Lab surveillance data

par(mar=c(5,5,3,1))

barplot(height=t(as.matrix(plots[1:40,15:16])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,40), ylim=c(0,0.3), xlab=" ", ylab="Isolation rate", cex.lab=1.5,
 col=c(gray(0.2),0)

#legend.text=c("Flu A +ve", "Flu B +ve")
) 
legend(32.5, 0.3, cex=1,
  legend=c("Flu B +ve", "Flu A +ve"), fill=c(0, gray(0.2)), lty=NULL)


axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:3/10, labels=c("0%","10%","20%","30%"), las=1, cex.axis=1.2)

mtext("Jan", side=1, at=2.5, line=1, cex=1)
mtext("Feb", side=1, at=7, line=1, cex=1)
mtext("Mar", side=1, at=11, line=1, cex=1)
mtext("Apr", side=1, at=15.5, line=1, cex=1)
mtext("May", side=1, at=20, line=1, cex=1)
mtext("Jun", side=1, at=24, line=1, cex=1)
mtext("Jul", side=1, at=28.5, line=1, cex=1)
mtext("Aug", side=1, at=33, line=1, cex=1)
mtext("Sep", side=1, at=37, line=1, cex=1)

mtext("(c)", cex=1)
