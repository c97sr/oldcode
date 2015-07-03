#
#
#

flu <- read.csv("d:\\work\\influenza\\data\\qmh0507.csv", header=TRUE)

months <- rep(0,104)
months[1] <- 1
for(i in 2:104){
    if(flu$month[i]!=flu$month[i-1]) months[i] <- 1
}
labs <- c(sort(unique(months*1:104))[-1], 104)


windows(width=6, height=3.5)
par(mar=c(0,4,1,0.5))
barplot(height=t(as.matrix(flu[1:104,8:9])), space=0, axes=FALSE, axisnames=FALSE,
 xlim=c(0,103), ylim=c(-0.05,0.3), xlab=" ", ylab=" ", legend.text=c("Influenza A", "Influenza B"))
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)))
axis(2, pos=0, at=0:3/10, las=1, labels=c("0%", "10%", "20%", "30%"))
for(i in 1:12){
    text(0.5*(labs[i+ 1]+labs[i   ])-0.75, 0.005, substr(month.abb[i],1,1), cex=0.5, pos=1)
    text(0.5*(labs[i+13]+labs[i+12])-0.75, 0.005, substr(month.abb[i],1,1), cex=0.5, pos=1)
}
lines(c(  0,  0), c(0,-0.04))
lines(c( 52, 52), c(0,-0.04))
lines(c(103,103), c(0,-0.04))
mtext("2005", side=1, at=25, line=-1.5)
mtext("2006", side=1, at=77, line=-1.5)
mtext("Proportion of positive isolations", side=2, at=0.15, line=2.5)


#
# the end
#
#
