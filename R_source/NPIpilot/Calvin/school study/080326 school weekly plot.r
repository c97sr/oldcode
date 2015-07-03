#combine school data with chp qv data

sch_dat <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/abs_total for graphs/080616sch_abs_total.csv")
chpqv <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/chpqvdata/080616 chpqvdata.csv")

#Weekly plots
plots <- merge(chpqv, sch_dat, by="week", all=TRUE)

plots <- plots[1:27,]
#plots <-na.omit(plots)

months <- rep(0,27)
months[1] <- 1
for(i in 2:27){
    if(plots$month[i]!=plots$month[i-1]) months[i] <- 1
}
labs <- c(sort(unique(months*1:27))[-1], 27)


#function for plotting graphs for each school
plotgraph <- function(plots, sch_num, sch_name, level){

# for plots side by side
windows(width=10, height=7)
layout(matrix(1:3, ncol=1))
par(mar=c(3,5,3,1))

#
# 1st window - individual school absenteeism rate 
# for low , medium high absenteeism rate's school
par(mar=c(5,5,2,1))
if (level == "high"){

barplot(height=t(as.matrix(plots[1:27,sch_num+10])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,27), ylim=c(0,0.08), xlab=" ", ylab="% ILI of absentees", cex.lab=1.5,) 

axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.02,0.04,0.06,0.08), labels=c("0%","2%","4%","6%","8%"), las=1, cex.axis=1.2)
mtext(paste(sch_name, "ILI absenteeism rate"), cex=1.2)
}

else if (level == "med"){

barplot(height=t(as.matrix(plots[1:27,sch_num+10])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,27), ylim=c(0,0.04), xlab=" ", ylab="% ILI of absentees", cex.lab=1.5,) 

axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.01,0.02,0.03,0.04), labels=c("0%","1%","2%","3%","4%"), las=1, cex.axis=1.2)
mtext(paste(sch_name, "ILI absenteeism rate"), cex=1.2)
}

else {

barplot(height=t(as.matrix(plots[1:27,sch_num+10])), space=0, axes=FALSE, axisnames=FALSE,
#main = "influenza rate of qm data ", 
xlim=c(0,27), ylim=c(0,0.02), xlab=" ", ylab="% ILI of absentees", cex.lab=1.5,) 

axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.005,0.01,0.015,0.02), labels=c("0%","0.5%","1%","1.5%","2%"), las=1, cex.axis=1.2)
mtext(paste(sch_name, "ILI absenteeism rate"), cex=1.2)
}

#
# 2nd window - overall school absenteeism rate

par(mar=c(5,5,1,1))

#plot(ksmooth(plots$week[8:15] - 0.5, as.matrix(plots[8:15,10]), "normal",2) , axes=FALSE, type="l",lwd=2,
plot(plots$week[8:27] - 0.5, as.matrix(plots[8:27,10]) , axes=FALSE, type="l",lwd=2,

xlim=c(0,27), ylim=c(0,0.04), xlab=" ", ylab="% of ILI absentees", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.01,0.02,0.03,0.04), labels=c("0%","1%","2%","3%","4%"), las=1, cex.axis=1.2)


mtext("Overall average school ILI absenteeism rate", cex=1.2)


#
# 3rd window - QV recruitment data
par(mar=c(3,5,1,1))

barplot(height=t(as.matrix(plots[1:27,9])), space=0,
 axes=FALSE, axisnames=FALSE, col=c(gray(0.2),0),

xlim=c(0,27), ylim=c(0,0.3), xlab=" ", ylab="Influenza detections",  cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.1,0.2,0.3) ,  labels=c("0","10%","20%","30%"), las=1, cex.axis=1.2)


mtext("Community influenza activity", cex=1.2)


mtext("Jan", side=1, at=2, line=1, cex=1)
mtext("Feb", side=1, at=6, line=1, cex=1)
mtext("Mar", side=1, at=10.5, line=1, cex=1)
mtext("Apr", side=1, at=15, line=1, cex=1)
mtext("May", side=1, at=19.5, line=1, cex=1)
mtext("Jun", side=1, at=24, line=1, cex=1)

}

# function for printing overall data for school 
plot_overall <- function(plots){

# for plots side by side
windows(width=12, height=15)
layout(matrix(1:3, ncol=1))
par(mar=c(3,5,3,1))

#
# 1st window - overall school absenteeism rate


par(mar=c(5,5,2,1))
#plot(ksmooth(plots$week[8:20] - 0.5, as.matrix(plots[8:20,10]), "normal",2) , axes=FALSE, type="l",lwd=1.5,
plot(plots$week[8:27] - 0.5, as.matrix(plots[8:27,10]) , axes=FALSE, type="l",lwd=1.5,
xlim=c(0,27), ylim=c(0,0.04), xlab=" ", ylab="% of absentees", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.01,0.02,0.03, 0.04), labels=c("0%","1%","2%","3%","4%"), las=1, cex.axis=1.2)


mtext("Overall average school absenteeism rate", cex=1)


#
# 2nd window - chp surveillnace GP rate per 1000 patients

par(mar=c(3,5,1,1))
#plot( ksmooth(plots$week[1:20] - 0.5, as.matrix(plots[1:27,5])/1000, "normal",2) , axes=FALSE, type="l", lwd=1.5,
plot( plots$week[1:27] - 0.5, as.matrix(plots[1:27,5])/1000 , axes=FALSE, type="l", lwd=1.5,
xlim=c(0,27), ylim=c(0,0.080), xlab=" ", ylab="Consultation rate", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.02,0.04,0.06,0.08), labels=c("0","2%","4%","6%","8%"), las=1, cex.axis=1.2)

mtext("CHP influenza-like-illness consultation rate (GP)", cex=1)

#
# 3rd window - QV recruitment data
par(mar=c(3,5,1,1))

barplot(height=t(as.matrix(plots[1:27,9])), space=0,
 axes=FALSE, axisnames=FALSE, col=c(gray(0.2),0),

xlim=c(0,27), ylim=c(0,0.3), xlab=" ", ylab="Influenza detections",  cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.1,0.2,0.3) ,  labels=c("0","10%","20%","30%"), las=1, cex.axis=1)


mtext("Community influenza activity", cex=1)


mtext("Jan", side=1, at=2, line=1, cex=1)
mtext("Feb", side=1, at=6, line=1, cex=1)
mtext("Mar", side=1, at=10.5, line=1, cex=1)
mtext("Apr", side=1, at=15, line=1, cex=1)
mtext("May", side=1, at=19.5, line=1, cex=1)
mtext("Jun", side=1, at=24, line=1, cex=1)

}
setwd("D:/Calvin flu studies/School absenteeism study/Documents/Doc to eclass")
plot_overall(plots)


setwd("D:/Calvin flu studies/School absenteeism study/Documents/documents to school/080616 to participating schools")


plotgraph(plots,1, "Baptist (STW) Lui Ming Choi Primary School","high")
dev.copy(png,filename="blmc1.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,2, "Carmel Alison Lam Foundation Secondary School", "high")
dev.copy(png,filename="calf2.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,3, "CUHKFAA Chan Chun Ha Secondary School", "med")
dev.copy(png,filename="ccch3.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,4, "The Chinese Foundation Secondary School", "low")
dev.copy(png,filename="cfss4.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,5, "Po Leung Kuk Chee Jing Yin Primary School", "med")
dev.copy(png,filename="cjyp5.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,6, "Delia Memorial School (Broadway)", "high")
dev.copy(png,filename="dmsb6.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,7, "Ho Dao College (Sponsored By Sik Sik Yuen)", "high")
dev.copy(png,filename="hdcl7.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,8, "True Light Middle School Of Hong Kong", "low")
dev.copy(png,filename="hktl8.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,9, "Ma On Shan Ling Liang Primary School", "med")
dev.copy(png,filename="llps9.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,10, "Maryknoll Convent School (Secondary Section)", "med")
dev.copy(png,filename="mcss10.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,11, "Precious Blood Primary School (South Horizons)", "med")
dev.copy(png,filename="pbps11.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,12, "Pentecostal Lam Hon Kwong School", "low")
dev.copy(png,filename="plhk12.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,13, "St. Francis Of Assisi's English Primary School", "low")
dev.copy(png,filename="sfps13.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,14, "FDBWA Szeto Ho Secondary School", "high")
dev.copy(png,filename="sthsl4.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,15, "Shan Tsui Public School", "med")
dev.copy(png,filename="stpsl5.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,16, "Salesian Yip Hon Millennium Primary School", "low")
dev.copy(png,filename="syhm16.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,17, "Yan Oi Tong Tin Ka Ping Secondary School", "low")
dev.copy(png,filename="tkps17.png",height=600, width=800,bg="white")
dev.off()
plotgraph(plots,18, "Wah Yan College, Hong Kong", "med")
dev.copy(png,filename="wyhkl8.png",height=600, width=800,bg="white")
dev.off()


# plot all school's data
#for (i in (length(plots)-6):1){ 
#  plotgraph(plots,i)
#}


