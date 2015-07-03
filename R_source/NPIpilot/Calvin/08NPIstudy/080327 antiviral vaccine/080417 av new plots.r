# new plots 080416


#for antiviral  study
symptom <-read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_symptomday0708_d.csv" , header=TRUE)


#check if it is the same day go home visit
housechar <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_housechar0708_h.csv", header=TRUE)
housechar$clinic_day[is.na(housechar$clinic_day)] <-0 # treat NA clinic day to be going home visit at same night
housechar$clin1vis <- housechar$v1_day - housechar$clinic_day
housechar <-housechar[c(1,15)]
symptom <-merge(symptom, housechar, by="hhID", all.x=TRUE)

# adjust the day according to clinic visit day 
symptom$day <- symptom$day + symptom$clin1vis
#delete those day 0
symptom <- symptom[symptom$day !=0,]
symptom <- symptom[c(1:11)]
#names(symptom)[9] <- "aches" 

#get clinic data
clinic <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_clinicdat0708_h.csv", header=TRUE)
clinic <- clinic[clinic$hhID != "h080NA" & clinic$hhID != "      ",]

#chang "clinic" data to "symptom" format
clinic$member <-0
clinic$day <-0
clinic <-clinic[c(1:2,19,20,16,8:13)]
names(clinic)[9] <- "pmuscle" 

symptom <- rbind(symptom, clinic) 

symptom <- symptom[!is.na(symptom$hhID),]
symptom <- symptom[order(symptom$hhID, symptom$member, symptom$day),]
symptom <- symptom[symptom$member ==0,]
write.csv(symptom, "P:/GROUP/NPIstudy/data/combined/2008_4_17_symwfclin0708_d.csv", row.names=FALSE)

# then run 080407av_vac-dat.r


############# 080416 plot new graphs, in 1 graph, seperate adult and children
# for plots side by side

dat <- read.csv("P:/GROUP/NPIstudy/data/combined/080417_av_vac.csv" , header=TRUE)

isadult <- dat[dat$isadult==1,]
ischild <- dat[dat$isadult==0,]


# function for getting the entry number by inputting the entry name
 get_entry_num<-function(data ,entry_name){
	return (match(entry_name, names(data)))
  }
#

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

# function for average score in each day
cal_score <-function(dat, condi, value, output, start_entname){
  start_ent <- get_entry_num(dat, start_entname)
  mean_score <- data.frame()
  for (i in 1:7){
    temp <- ci(na.exclude(dat[condi == value & output != 0, i+start_ent]))
    mean_score <- rbind(mean_score, temp) 
    
  }
  mean_score <- as.data.frame(t(mean_score))
  #names(mean_score) <- c("day0","day1","day2","day3","day4","day5","day6")
  row.names(mean_score)[1] <- paste(value, "mean", sep="_")
return(mean_score)
}
#

#function to plot 3 lines in 1 graph
plot_1g <-  function(plot1,plot2,plot3,name1,name2,name3,title){

#first window

par(mar=c(3,5,3,1))
plot( c(1:7-0.5), plot1, axes=FALSE, type="l", main = title,
xlim=c(0,7), ylim=c(0,1), xlab=" ", ylab="Symptom score ", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,0.2,0.4,0.6,0.8,1.0), labels=c("0","0.2","0.4","0.6","0.8","1.0"), las=1, cex.axis=1.2)

lines(c(1:7-0.5), plot2, lty=2)
lines( c(1:7-0.5), plot3, lty=3)

#legend(4, 4, cex=1,
#  legend=c(name1,name2,name3), col=c("black","red","blue"),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
#)
#x axis
mtext("Day  0", side=1, at=0, line=1, cex=1)
mtext("1", side=1, at=1.5, line=1, cex=1)
mtext("2", side=1, at=2.5, line=1, cex=1)
mtext("3", side=1, at=3.5, line=1, cex=1)
mtext("4", side=1, at=4.5, line=1, cex=1)
mtext("5", side=1, at=5.5, line=1, cex=1)
mtext("6", side=1, at=6.5, line=1, cex=1)
}
#

#function to plot 3 lines in 1 graph
plot_fg <-  function(plot1,plot2,plot3,name1,name2,name3,title){

#first window
par(mar=c(3,5,3,1))
plot( c(1:7-0.5), plot1, axes=FALSE, type="l", main = title,
xlim=c(0,7), ylim=c(0,4), xlab=" ", ylab="Body temperature (C)", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,1,2,3,4), labels=c("36","37","38","39","40"), las=1, cex.axis=1.2)

lines(c(1:7-0.5), plot2, lty=2)
lines( c(1:7-0.5), plot3, lty=3)

#legend(5, 2.5, cex=1,
#  legend=c(name1,name2,name3), lty=c(1,2,3),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
#)
#x axis
mtext("Day  0", side=1, at=0, line=1, cex=1)
mtext("1", side=1, at=1.5, line=1, cex=1)
mtext("2", side=1, at=2.5, line=1, cex=1)
mtext("3", side=1, at=3.5, line=1, cex=1)
mtext("4", side=1, at=4.5, line=1, cex=1)
mtext("5", side=1, at=5.5, line=1, cex=1)
mtext("6", side=1, at=6.5, line=1, cex=1)
}
#

cal_score(dat, dat$dummy , 1,dat$all_score,"all_score") #total averge symptom score each day (index only)
table(dat$av)#analysis different antivirals taken
cal_score(dat, dat$av, "none", dat$all_score, "all_score")
cal_score(dat, dat$av, "amantadine", dat$all_score, "all_score")
cal_score(dat, dat$av, "tamiflu", dat$all_score, "all_score")


windows(width=24, height=12)
layout(matrix(1:6, nrow=2))

#plot daily average symptom score , different antivirals , adult and children
plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$all_score, "all_score")[1,]) /7
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$all_score, "all_score")[1,])/7
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$all_score, "all_score")[1,])/7
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "Adult symptom score")

plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$up_score, "up_score")[1,]) /2
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$up_score, "up_score")[1,])/2
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$up_score, "up_score")[1,])/2
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "UR symptom score")

plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$sys_score, "sys_score")[1,]) /3
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$sys_score, "sys_score")[1,])/3
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$sys_score, "sys_score")[1,])/3
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "SYS symptom score")

plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$low_score, "low_score")[1,]) /2
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$low_score, "low_score")[1,])/2
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$low_score, "low_score")[1,])/2
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "LR symptom score")

plot1 <- t(cal_score(isadult, isadult$av, "none", isadult$f_score, "f_score")[1,]) -36
plot2 <- t(cal_score(isadult, isadult$av, "amantadine", isadult$f_score, "f_score")[1,]) -36
plot3 <- t(cal_score(isadult, isadult$av, "tamiflu", isadult$f_score, "f_score")[1,]) -36
plot_fg(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "Body Temperature") 


#########################################

windows(width=24, height=12)
layout(matrix(1:6, nrow=2))

#plot daily average symptom score , different antivirals , adult and children
plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$all_score, "all_score")[1,]) /7
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$all_score, "all_score")[1,])/7
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$all_score, "all_score")[1,])/7
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "Children symptom score")

plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$up_score, "up_score")[1,]) /2
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$up_score, "up_score")[1,])/2
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$up_score, "up_score")[1,])/2
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "UR symptom score")

plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$sys_score, "sys_score")[1,]) /3
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$sys_score, "sys_score")[1,])/3
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$sys_score, "sys_score")[1,])/3
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "SYS symptom score")

plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$low_score, "low_score")[1,]) /2
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$low_score, "low_score")[1,])/2
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$low_score, "low_score")[1,])/2
plot_1g(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "LR symptom score")

plot1 <- t(cal_score(ischild, ischild$av, "none", ischild$f_score, "f_score")[1,]) -36
plot2 <- t(cal_score(ischild, ischild$av, "amantadine", ischild$f_score, "f_score")[1,]) -36
plot3 <- t(cal_score(ischild, ischild$av, "tamiflu", ischild$f_score, "f_score")[1,]) -36
plot_fg(plot1,plot2,plot3,"None","Amantadine","Tamiflu", "Body Temperature") 
