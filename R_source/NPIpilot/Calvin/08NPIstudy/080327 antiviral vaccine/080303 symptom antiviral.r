
#NPI study analysis for vaccination and antiviral
###############

alldat <- read.csv("P:/GROUP/NPIstudy/data/combined/080417_av_vac.csv", header=TRUE)
dat <- alldat[alldat$member ==0,] #get index patient only

isadult <- dat[dat$isadult ==1,] #compare those child and adult with different antiviral taken
ischild <- dat[dat$isadult ==0,]
table(isadult$av)
table(ischild$av)

#antiviral precribed or not for index
in_noav <-dat[is.na(dat$av) &  dat$member==0,]
in_av <-dat[!is.na(dat$av) & dat$member==0,]

#antiviral perscribed 24hr before and after 24hrs
#dat24noav <- dat[dat$onsettime<24 & is.na(dat$av),] 
dat24av <- dat[dat$onsettime<24 & !is.na(dat$av),] 
#dat48noav <- dat[dat$onsettime>=24 & is.na(dat$av),] 
dat48av <- dat[dat$onsettime>=24 & !is.na(dat$av),] 




# function for getting the entry number by inputting the entry name
 get_entry_num<-function(data ,entry_name){
	return (match(entry_name, names(data)))
  }
#end

# function for calculating the std error
ci <- function( x ){
  stderror <- function(x){ sqrt(var(x)/length(x))}
  mean_v <- mean(x)
  stderror_v <- stderror(x)
  upperCI <- mean_v + (1.96*stderror_v)
  lowerCI <- mean_v - (1.96*stderror_v)
  y <- cbind(mean_v,lowerCI,upperCI)
  
  return(round(y,2))
}
#


#080407 for calculating daily (day0-day6) symptom score or duration with different conditions

#function for overall mean score
mean_score <-function(dat,condi,value,output){
  mean_v <- mean(na.exclude(output[ condi==value & output != 0]))  #mean duration/score of symptoms of overall ppl
  frequen <- length(na.exclude(output[ condi==value & output != 0]))
  return(round(cbind(mean_v, frequen),2))
}
#

#function for overall mean and 95% CI of score
ci_score <-function(dat,condi,value,output){
  mean_v <- ci(na.exclude(output[ condi==value & output != 0]))  #mean duration/score of symptoms of overall ppl
  frequen <- length(na.exclude(output[ condi==value & output != 0]))
  return(round(cbind(mean_v, frequen),2))
}
#

# function for average score in each day
cal_score <-function(dat, condi, value, output, start_entname){
  start_ent <- get_entry_num(dat, start_entname)
  mean_score <- data.frame()
  for (i in 1:8){
    temp <- ci(na.exclude(dat[condi == value & output != 0, i+start_ent]))
    mean_score <- rbind(mean_score, temp) 
    
  }
  mean_score <- as.data.frame(t(mean_score))
  #names(mean_score) <- c("day0","day1","day2","day3","day4","day5","day6")
  row.names(mean_score)[1] <- paste(value, "mean", sep="_")
return(mean_score)
}



#analysis

# overall average duration of symptoms 
mean_score(dat,dat$av,"none", dat$all_duration)
mean_score(dat,dat$av,"amantadine", dat$all_duration)
mean_score(dat,dat$av,"tamiflu", dat$all_duration)

#overall symptom score
mean_score(dat,dat$av,"none", dat$all_score)
mean_score(dat,dat$av,"amantadine", dat$all_score)
mean_score(dat,dat$av,"tamiflu", dat$all_score)

cal_score(dat, dat$dummy , 1,dat$all_score,"all_score") #total averge symptom score each day (index only)
table(dat$av)#analysis different antivirals taken
cal_score(dat, dat$av, "none", dat$all_score, "all_score")
cal_score(dat, dat$av, "amantadine", dat$all_score, "all_score")
cal_score(dat, dat$av, "tamiflu", dat$all_score, "all_score")

table(dat$isadult)
cal_score(dat, dat$isadult, 0, dat$all_score) #children <16
cal_score(dat, dat$isadult, 1, dat$all_score) #adult

cal_score(isadult, isadult$av, "none", isadult$all_score, "all_score")
cal_score(isadult, isadult$av, "amantadine", isadult$all_score,  "all_score")
cal_score(isadult, isadult$av, "tamiflu", isadult$all_score , "all_score")

cal_score(ischild, ischild$av, "none", ischild$all_score , "all_score")
cal_score(ischild, ischild$av, "amantadine", ischild$all_score , "all_score")
cal_score(ischild, ischild$av, "tamiflu", ischild$all_score , "all_score")
###############

#for plot average daily symptom score



#function for plotting graphs
plot_g <-  function(plot1,plot2,plot3){

windows(width=5, height=7)
layout(matrix(1:3, ncol=1))
par(mar=c(3,5,3,1))

#first window

par(mar=c(5,5,1,1))
plot( c(1:7-0.5), plot1, axes=FALSE, type="b",
xlim=c(0,7), ylim=c(0,5), xlab=" ", ylab="Symptom score ", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4,5), labels=c("0","1","2","3","4","5"), las=1, cex.axis=1.2)

#2nd window

par(mar=c(5,5,1,1))
plot( c(1:7-0.5), plot2, axes=FALSE, type="b",
xlim=c(0,7), ylim=c(0,5), xlab=" ", ylab="Symptom score", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4,5), labels=c("0","1","2","3","4","5"), las=1, cex.axis=1.2)

#3rd window 

par(mar=c(5,5,1,1))
plot( c(1:7-0.5), plot3, axes=FALSE, type="b",
xlim=c(0,7), ylim=c(0,5), xlab=" ", ylab="Symptom score", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4,5), labels=c("0","1","2","3","4","5"), las=1, cex.axis=1.2)

#x axis
mtext("Day 0", side=1, at=0.5, line=1, cex=1)
mtext("1", side=1, at=1.5, line=1, cex=1)
mtext("2", side=1, at=2.5, line=1, cex=1)
mtext("3", side=1, at=3.5, line=1, cex=1)
mtext("4", side=1, at=4.5, line=1, cex=1)
mtext("5", side=1, at=5.5, line=1, cex=1)
mtext("6", side=1, at=6.5, line=1, cex=1)

}


#function to plot 3 lines in 1 graph
plot_1g <-  function(plot1,plot2,plot3,name1,name2,name3,title){

#first window

par(mar=c(3,5,3,1))
plot( c(1:8-0.5), plot1, axes=FALSE, type="l", main = title,
xlim=c(0,8), ylim=c(0,5), xlab=" ", ylab="Symptom score ", cex.lab=1.2)
axis(1, pos=0, at=0:8, labels=rep(NA, 8), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4,5), labels=c("0","1","2","3","4","5"), las=1, cex.axis=1.2)

lines(c(1:8-0.5), plot2, col = "red")
lines( c(1:8-0.5), plot3, col = "blue")

legend(4, 4, cex=1,
  legend=c(name1,name2,name3), col=c("black","red","blue"),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
)
#x axis
mtext("Day      C", side=1, at=0, line=1, cex=1.2)
mtext("0", side=1, at=1.5, line=1, cex=1.2)
mtext("1", side=1, at=2.5, line=1, cex=1.2)
mtext("2", side=1, at=3.5, line=1, cex=1.2)
mtext("3", side=1, at=4.5, line=1, cex=1.2)
mtext("4", side=1, at=5.5, line=1, cex=1.2)
mtext("5", side=1, at=6.5, line=1, cex=1.2)
mtext("6", side=1, at=7.5, line=1, cex=1.2)

}


#plot daily average symptom score , different antivirals
#overall no antiviral, amantadine and tamiflu 's symptom score per day
plot1 <- t(cal_score(dat, dat$av, "none", dat$all_score, "all_score")[1,])
plot2 <- t(cal_score(dat, dat$av, "amantadine", dat$all_score, "all_score")[1,])
plot3 <- t(cal_score(dat, dat$av, "tamiflu", dat$all_score, "all_score")[1,])

plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Overall average symptom score 
among different antiviral percribed")
plot_all <- cbind(plot1,plot2,plot3)

#adult >=16 no antiviral, amantadine and tamiflu 's symptom score per day
plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$all_score, "all_score")[1,])
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$all_score, "all_score")[1,])
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$all_score, "all_score")[1,])
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Average symptom score among 
different antiviral percribed in adults")
plot_adult <- cbind(plot1,plot2,plot3)

#children <16
plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$all_score, "all_score")[1,])
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$all_score, "all_score")[1,])
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$all_score, "all_score")[1,])

plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Average symptom score among 
different antiviral percribed in children")
plot_child <- cbind(plot1,plot2,plot3)



#function to plot 3 lines in 1 graph
plot_2g <-  function(plot1,plot2,plot3,name1,name2,name3,title){

#first window

par(mar=c(3,5,3,1))
plot( c(1:8-0.5), plot1, axes=FALSE, type="l", main = title,
xlim=c(0,8), ylim=c(0,2.5), xlab=" ", ylab="Symptom score ", cex.lab=1.2)
axis(1, pos=0, at=0:8, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.5,1,1.5,2,2.5), labels=c("0","0.5","1","1.5","2","2.5"), las=1, cex.axis=1.2)

lines(c(1:8-0.5), plot2, col = "red")
lines( c(1:8-0.5), plot3, col = "blue")

legend(5, 1.5, cex=1,
  legend=c(name1,name2,name3), col=c("black","red","blue"),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
)
#x axis
mtext("Day      C", side=1, at=0, line=1, cex=1.2)
mtext("0", side=1, at=1.5, line=1, cex=1.2)
mtext("1", side=1, at=2.5, line=1, cex=1.2)
mtext("2", side=1, at=3.5, line=1, cex=1.2)
mtext("3", side=1, at=4.5, line=1, cex=1.2)
mtext("4", side=1, at=5.5, line=1, cex=1.2)
mtext("5", side=1, at=6.5, line=1, cex=1.2)
mtext("6", side=1, at=7.5, line=1, cex=1.2)

}


#function to plot 3 lines in 1 graph
plot_fg <-  function(plot1,plot2,plot3,name1,name2,name3,title){

#first window

par(mar=c(3,5,3,1))
plot( c(1:8-0.5), plot1, axes=FALSE, type="l", main = title,
xlim=c(0,8), ylim=c(0,4), xlab=" ", ylab="Body temperature (C)", cex.lab=1.2)
axis(1, pos=0, at=0:8, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4), labels=c("36","37","38","39","40"), las=1, cex.axis=1.2)

lines(c(1:8-0.5), plot2, col = "red")
lines( c(1:8-0.5), plot3, col = "blue")

legend(5, 2.5, cex=1,
  legend=c(name1,name2,name3), col=c("black","red","blue"),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
)
#x axis
mtext("Day      C", side=1, at=0, line=1, cex=1.2)
mtext("0", side=1, at=1.5, line=1, cex=1.2)
mtext("1", side=1, at=2.5, line=1, cex=1.2)
mtext("2", side=1, at=3.5, line=1, cex=1.2)
mtext("3", side=1, at=4.5, line=1, cex=1.2)
mtext("4", side=1, at=5.5, line=1, cex=1.2)
mtext("5", side=1, at=6.5, line=1, cex=1.2)
mtext("6", side=1, at=7.5, line=1, cex=1.2)


}



#check upper respiratory scores between different antivirals
plot1 <- t(cal_score(dat, dat$av, "none", dat$up_score,"up_score")[1,])
plot2 <- t(cal_score(dat, dat$av, "amantadine", dat$up_score,"up_score")[1,])
plot3 <- t(cal_score(dat, dat$av, "tamiflu", dat$up_score,"up_score")[1,])

plot_2g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Overall average upper repiratory symptom score 
among different antiviral percribed")
plot_upper <- cbind(plot1,plot2,plot3)

#check lower respiratory scores between different antivirals
plot1 <- t(cal_score(dat, dat$av, "none", dat$low_score,"low_score")[1,])
plot2 <- t(cal_score(dat, dat$av, "amantadine", dat$low_score,"low_score")[1,])
plot3 <- t(cal_score(dat, dat$av, "tamiflu", dat$low_score,"low_score")[1,])

plot_2g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Overall average lower repiratory symptom score 
among different antiviral percribed")
plot_lower <- cbind(plot1,plot2,plot3)

#check systmatic scores between different antivirals
plot1 <- t(cal_score(dat, dat$av, "none", dat$sys_score,"sys_score")[1,])
plot2 <- t(cal_score(dat, dat$av, "amantadine", dat$sys_score,"sys_score")[1,])
plot3 <- t(cal_score(dat, dat$av, "tamiflu", dat$sys_score,"sys_score")[1,])

plot_2g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Overall average systemic symptom score 
among different antiviral percribed")
plot_systemic <- cbind(plot1,plot2,plot3)

#check body temperature between different antivirals
plot1 <- t(cal_score(dat, dat$av, "none", dat$f_score,"f_score")[1,]) -36
plot2 <- t(cal_score(dat, dat$av, "amantadine", dat$f_score,"f_score")[1,]) -36
plot3 <- t(cal_score(dat, dat$av, "tamiflu", dat$f_score,"f_score")[1,]) -36

plot_fg(plot1,plot2,plot3,"None","Amantadine","Tamiflu", 
"Overall average Bodytemperature 
among different antiviral percribed")
plot_fever <- cbind(plot1,plot2,plot3)



