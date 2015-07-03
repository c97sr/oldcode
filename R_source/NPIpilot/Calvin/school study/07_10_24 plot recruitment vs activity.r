plots <- read.csv("P:\\GROUP\\NPIpilot\\data\\surveillance\\07_10_24_qvchpqmrate.csv", header=TRUE)
plots <- plots[1:week,]
#plots <-na.omit(plots)

#update how many weeks
week <- 38

months <- rep(0,week)
months[1] <- 1
for(i in 2:week){
   
    if(plots$month[i]!=plots$month[i-1]) months[i] <- 1
   
}
labs <- c(sort(unique(months*1:week))[-1], week+1)

# for plots side by side
windows(width=8, height=6)
layout(matrix(1:3, ncol=1))
par(mar=c(3,5,2,1))

#
# 1st window - QV recruitment data

barplot(height=t(as.matrix(plots[8:week,9:10])), space=0,
 axes=FALSE, axisnames=FALSE, col=c(gray(0.2),0),
# main="No. of QV tests done per working day",
 xlim=c(-7,32), ylim=c(0,12), xlab=" ", ylab="Recruitment rate", 
# legend.text=c("QV +ve", "QV -ve")
 cex.lab=1.5)
axis(1, pos=0, at=labs-8, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=-7, at=c(0,3,6,9,12) ,  labels=c("0","3","6","9","12"), las=1, cex.axis=1.2)
legend(-5, 9, cex=1,
  legend=c("QV -ve", "QV +ve"), fill=c(0, gray(0.2)), lty=NULL)


#
# 2nd window - chp surveillance data, GP rate

par(mar=c(3,5,1,1))

plot(as.matrix(plots[1:week,18]), axes=FALSE, type="l",
# main = "CHP surveillance data, GP rate per 1000 patients", 
 xlim=c(0,40), ylim=c(30,70), xlab=" ", ylab="ILI rate", cex.lab=1.5)
axis(1, pos=30, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=3:7*10, labels=c("30","40","50","60","70"), las=1, cex.axis=1.2)

#
# 3rd window - queen mary Lab surveillance data

par(mar=c(5,5,1,1))

barplot(height=t(as.matrix(plots[1:week,15:16])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,40), ylim=c(0,0.3), xlab=" ", ylab="Isolation rate", cex.lab=1.5,
legend.text=c("Flu A +ve", "Flu B +ve"))
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:3/10, labels=c("0%","10%","20%","30%"), las=1, cex.axis=1.2)

mtext("Jan", side=1, at=2, line=1, cex=1)
mtext("Feb", side=1, at=6, line=1, cex=1)
mtext("Mar", side=1, at=10.5, line=1, cex=1)
mtext("Apr", side=1, at=15, line=1, cex=1)
mtext("May", side=1, at=19, line=1, cex=1)
mtext("Jun", side=1, at=23, line=1, cex=1)
mtext("Jul", side=1, at=27.5, line=1, cex=1)
mtext("Aug", side=1, at=32, line=1, cex=1)
mtext("Sep", side=1, at=36.5, line=1, cex=1)


# for scatter
plots$qmtotalrate <- (plots[12]+plots[13])/plots[14]
plots$qvposrate <- plots[6]/plots[8]
data <- plots[c(18,17,19,6,8,20 )] 


#data <- read.csv("d:\\r\\r raw data\\07_10_24_simrate.csv")
library(lattice)
splom(data[,1:6]) 