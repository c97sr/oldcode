dat <- read.csv("P:\\GROUP\\NPIpilot\\data\\surveillance\\07_11_26_chplab_edited.csv", header=TRUE)


barplot(height=t(as.matrix(dat[1:12,8:9])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,12), ylim=c(0,0.2), xlab=" ", ylab="CHP isolation rate", cex.lab=1.5,
 col=c(gray(0.2),0)
 ) 
legend(8, 0.2, cex=1,
  legend=c("Flu B +ve", "Flu A +ve"), fill=c(0, gray(0.2)), lty=NULL)

axis(1, pos=0, at=0:12, labels=rep(NA,12), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.05,0.1,0.15,0.2), labels=c("0%","5%","10%","15%","20%"), las=1, cex.axis=1.2)

mtext("Jan", side=1, at=0.5, line=1, cex=1)
mtext("Feb", side=1, at=1.5, line=1, cex=1)
mtext("Mar", side=1, at=2.5, line=1, cex=1)
mtext("Apr", side=1, at=3.5, line=1, cex=1)
mtext("May", side=1, at=4.5, line=1, cex=1)
mtext("Jun", side=1, at=5.5, line=1, cex=1)
mtext("Jul", side=1, at=6.5, line=1, cex=1)
mtext("Aug", side=1, at=7.5, line=1, cex=1)
mtext("Sep", side=1, at=8.5, line=1, cex=1)
mtext("Oct", side=1, at=9.5, line=1, cex=1)
mtext("Nov", side=1, at=10.5, line=1, cex=1)
mtext("Dec", side=1, at=11.5, line=1, cex=1)
