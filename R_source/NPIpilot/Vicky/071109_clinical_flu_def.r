#
# Factor analysis for clinical influenza definition (related home culture to symptoms)
#
library(gee)
library(ROCR)

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_27_home_culture.csv", header=TRUE)
symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_28_symptomday_d.csv", header=TRUE)
lab2nd <- read.csv("D:\\Work\\Influenza\\data\\071029_hculture.csv",header=TRUE)

lab2nd <- lab2nd[,c(2,3,12)]

# Add day_diary where the culture days are consistent with those in symptom diaries
hc <- hc[hc$visit!=0,]
for ( i in 1:nrow(hc)){          
         hc$day[i] <- hc$day.from.clinic[i]-hc$day.from.clinic[hc$hhID==hc$hhID[i]&hc$member==0&hc$visit==1]
}
  
# Define fever if bodytemp>38
for (i in 1:nrow(symptom)){
    if (!is.na(symptom$bodytemp[i]) & as.numeric(symptom$bodytemp[i])>38) {symptom$fever[i]<-1}
    else  {symptom$fever[i]<-0}
}

## Define fever if bodytemp>38 for 2-continuous days
#i<-1
#while (i < nrow(symptom)+1){
#    if (!is.na(symptom$bodytemp[i]) & as.numeric(symptom$bodytemp[i])>38 &
#         !is.na(symptom$bodytemp[i+1]) & as.numeric(symptom$bodytemp[i+1])>38 ) {symptom$fever[i]<-1}
#      else  {symptom$fever[i]<-0}
#    i<- i+1
#    for (j in 1:8){
#      if (!is.na(symptom$bodytemp[i]) & as.numeric(symptom$bodytemp[i])>38 &
#        ( (!is.na(symptom$bodytemp[i-1]) & as.numeric(symptom$bodytemp[i])>38) 
#	  | (!is.na(symptom$bodytemp[i+1]) & as.numeric(symptom$bodytemp[i])>38) )) {symptom$fever[i]<-1}
#      else  {symptom$fever[i]<-0}
#      i<- i+1
#    }
#    if (!is.na(symptom$bodytemp[i]) & as.numeric(symptom$bodytemp[i])>38 &
#         !is.na(symptom$bodytemp[i-1]) & as.numeric(symptom$bodytemp[i])>38 ) {symptom$fever[i]<-1}
#      else  {symptom$fever[i]<-0}
#    i<- i+1
#}

cflu <- hc[,c(1,2,6,8:12)]

new1 <- merge(cflu,symptom)
new <- merge(new1,lab2nd,all.x=T)

for (i in 1:nrow(new)){
    if (as.character(new$culture[i]) == "A" | as.character(new$culture[i]) == "B" |
         (!is.na(new$AqPCR[i])&new$AqPCR[i]>0)|(!is.na(new$AqPCR2[i])&new$AqPCR2[i]>0)|
	 (!is.na(new$AqPCR3[i])&new$AqPCR3[i]>0)) {new$labflu[i] <- 1}
        else{new$labflu[i] <- 0}
}

new <- new[,-c(4:10)]
#write.csv(new,"D:\\Work\\Influenza\\data\\071030cflu.csv")

#cflu.gee <- gee(labflu~fever+headache+lnode+rnose+hvoice+sthroat+cough+phlegm+sbreath
#               +cpain+apain+vomit+diarrhoea+sjoint+pmuscle+dizzy+tired+chill+nsweat,
#               id=factor(hhID), data=new, family="binomial")
#
#summary(cflu.gee)

#for (i in 1:nrow(new)){
#    new$ind[i]<-sum(new$fever[i],new$headache[i],new$rnose[i],new$hvoice[i],new$sthroat[i],
#                    new$cough[i],new$phlegm[i],new$sbreath[i],new$pmuscle[i],new$tired[i],new$nsweat[i])
#    if(new$ind[i]>=3) new$cflu[i]<-1
#    else new$cflu[i]<-0
#}

#####################################################################
#                                                                   #
# May run the following scripts which are independent of the above  #
#                                                                   #
#####################################################################

symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_08_symptomday_d.csv", header=TRUE)

# Lab secondary cases
lab2nd <- read.csv("D:\\Work\\Influenza\\data\\071029_hculture.csv",header=TRUE)
lab2nd <- lab2nd[,c(2,3,12)]

#########
# Clinic secondary cases (need to run 071018SAR_lab_clinic.r first)

clinic2nd <- hclinic[,c(1,2,14)]

symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>38) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}

new <- merge(symptom,clinic2nd,all.x=T)
new <- merge(new,lab2nd,all.x=T)
new.contact <- new[new$member>0&new$day>0,-c(3:5)]
new.contact$hhID <- as.numeric(substr(new.contact$hhID,3,6))
i<-1
j<-1
fq <- rep(NA,23)
freq <- data.frame()
while (i < nrow(new.contact)+1){
      while(new.contact$hhID[i+1]==new.contact$hhID[i]&new.contact$member[i+1]==new.contact$member[i]){i<-i+1}
       fq[1] <- new.contact$hhID[j]
       fq[2] <- new.contact$member[j]
       for (m in 3:23){
         if (sum(na.exclude(new.contact[j:i,m]))>0) fq[m] <- 1
	 else if (!is.na(sum(new.contact[j:i,m]))) fq[m] <- 0
	 else fq[m] <- NA
       }
      freq<-rbind(freq,fq)
      i<-i+1
      j<-i
}
names(freq)<-c("hhID","member","headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath","cpain",
             "apain","vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chill","nsweat","fever","clinic2nd","lab2nd") 

# plot the data

# For clinical 2nd cases but not lab-confirmed 
freq_2nd <- freq[freq$clinic2nd==1&freq$lab2nd==0,]
for (i in 1:nrow(freq_2nd)){
freq_2nd$asy[i] <- sum(freq_2nd[i,c(3:21)])           #Asymptomatic
}
symp_freq_2nd <- rep(NA,20)
for (i in 2:20-1){
     symp_freq_2nd[i] <- sum(na.exclude(freq_2nd[,c(i+2)]))/nrow(freq_2nd) 
    # symp_freq_2nd[i] <- sum(na.exclude(freq_2nd[,i+2]))/length(na.exclude(freq_2nd[,i+2]))  
}
symp_freq_2nd[20] <- length(freq_2nd$asy[!is.na(freq_2nd$asy)&freq_2nd$asy==0])/nrow(freq_2nd)

names(symp_freq_2nd) <- c("headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath","cpain","apain",
                      "vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chill","nsweat","fever","asymptomatic" )

# For clinical 2nd cases and lab-confirmed 
freq_2nd_lab <- freq[freq$clinic2nd==1&freq$lab2nd==1,]
for (i in 1:nrow(freq_2nd_lab)){
freq_2nd_lab$asy[i] <- sum(freq_2nd_lab[i,c(3:21)])           #Asymptomatic
}
symp_freq_2nd_lab <- rep(NA,20)
for (i in 2:20-1){
     symp_freq_2nd_lab[i] <- sum(na.exclude(freq_2nd_lab[,c(i+2)]))/nrow(freq_2nd_lab)  
    # symp_freq_2nd_lab[i] <- sum(na.exclude(freq_2nd_lab[,i+2]))/length(na.exclude(freq_2nd_lab[,i+2]))  
}
symp_freq_2nd_lab[20] <- length(freq_2nd_lab$asy[!is.na(freq_2nd_lab$asy)&freq_2nd_lab$asy==0])/nrow(freq_2nd_lab)

names(symp_freq_2nd_lab) <- c("headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath","cpain","apain",
                      "vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chill","nsweat","fever","asymptomatic" )

od <- order(symp_freq_2nd_lab)
plotdata<-data.frame()
for (i in 1:20){
plotdata[1,i] <- symp_freq_2nd_lab[od[i]]
plotdata[2,i] <- symp_freq_2nd[od[i]]
names(plotdata)[i] <- names(symp_freq_2nd_lab[od[i]])
}

plotdata <- plotdata[,-2]  # delete asymptomatic column

windows(width=5, height=4)
par(mar=c(5,7,2,2), xaxs="i", yaxs="i")
barplot(as.matrix(rev(plotdata)), beside=T,names.arg=FALSE, horiz=T, xlab="% of symptoms",
        axes=FALSE,xlim=(c(0,1.0)),col=rep(c("blue","red")))
axis(1, cex.axis=1.5)
axis(2, at=19:1*3.0-1.0,labels=names(plotdata),las=1.5)
legend(0.6,45,"lab-confirmed secondary case (n=11)",fill="blue",bty="n")
legend(0.6,50,"lab-confirmed non-secondary case (n=26)",fill="red",bty="n")
title("Symptoms ever record during 9 days of follow-up for 37 clinic confirmed secondary cases")



##########################################################################################################
# barchart: Symptoms ever recorded during ~9 days of follow-up for 22 laboratory-confirmed secondary cases in NPI pilot study.
# The following scripts could be run independent of the above parts.
##########################################################################################################

symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_08_symptomday_d.csv", header=TRUE)

# Lab secondary cases
lab2nd <- read.csv("D:\\Work\\Influenza\\data\\071029_hculture.csv",header=TRUE)
lab2nd <- lab2nd[,c(2,3,12)]

symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>38) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}

new <- merge(symptom,lab2nd,all.x=T)
new.contact <- new[new$member>0&new$day>0,-c(3:5)]
new.contact$hhID <- as.numeric(substr(new.contact$hhID,3,6))

i<-1
j<-1
fq <- rep(NA,22)
freq <- data.frame()
while (i < nrow(new.contact)+1){
      while(new.contact$hhID[i+1]==new.contact$hhID[i]&new.contact$member[i+1]==new.contact$member[i]){i<-i+1}
       fq[1] <- new.contact$hhID[j]
       fq[2] <- new.contact$member[j]
       for (m in 3:22){
         if (sum(na.exclude(new.contact[j:i,m]))>0) fq[m] <- 1
	 else if (!is.na(sum(new.contact[j:i,m]))) fq[m] <- 0
	 else fq[m] <- NA
       }
      freq<-rbind(freq,fq)
      i<-i+1
      j<-i
}
names(freq)<-c("hhID","member","headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath","cpain",
             "apain","vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chill","nsweat","fever","lab2nd") 

# plot the data
freq_2nd_lab <- freq[freq$lab2nd==1,]
for (i in 1:nrow(freq_2nd_lab)){
freq_2nd_lab$asy[i] <- sum(freq_2nd_lab[i,c(3:21)])           #Asymptomatic
}
symp_freq_2nd_lab <- rep(NA,20)
for (i in 2:20-1){
     #symp_freq_2nd_lab[i] <- sum(na.exclude(freq_2nd_lab[,c(i+2)]))/nrow(freq_2nd_lab)  
     symp_freq_2nd_lab[i] <- sum(na.exclude(freq_2nd_lab[,i+2]))/length(na.exclude(freq_2nd_lab[,i+2]))  
}
symp_freq_2nd_lab[20] <- length(freq_2nd_lab$asy[!is.na(freq_2nd_lab$asy)&freq_2nd_lab$asy==0])/nrow(freq_2nd_lab)

names(symp_freq_2nd_lab) <- c("headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath","cpain","apain",
                      "vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chill","nsweat","fever","asymptomatic" )

plotdata <- sort(symp_freq_2nd_lab[1:19])
plotdata[20] <- symp_freq_2nd_lab[20]
names(plotdata)[20] <- "asymptomatic"

windows(width=5, height=4)
par(mar=c(5,7,2,2), xaxs="i", yaxs="i")
barplot(rev(plotdata), names.arg=FALSE, horiz=T, xlab="% of symptoms",axes=FALSE,xlim=(c(0,0.8)),col=rep(c("gray","blue"),c(1,19)))
axis(1, cex.axis=1.5)
axis(2, at=20:1*1.2-0.5,labels=names(plotdata),las=1.5)
title("Proportion of various symptoms under lab-confirmed secondary cases")
legend(0.3,20,"lab 2nd case",fill="blue",bty="n")
