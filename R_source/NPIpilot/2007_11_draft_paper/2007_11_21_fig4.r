
################################################
#                                              #
#  Proportion of subjects reporting symptom    #
#                                              #
################################################

symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)
lab2nd <- read.csv("D:\\Work\\Influenza\\output NPI year 1\\2007_11_draft\\lab2nd.csv",header=TRUE) # wrote in table4.r
lab2nd <- lab2nd[,c(1,2,12)]  # hhID,member,labsedcase

symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=38) 
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
names(freq)<-c("hhID","member","Headache","Swollen lymph nodes","Runny nose","Hoarseness of voice","Sore throat",
               "Cough","Phlegm","Shortness of breath","Chest pain","Abdominal pain","Vomiting","Diarrhoea",
	       "Stiffness in muscles or joints","Aches or pains in muscles or joints","Dizziness","Fatigue/tiredness",
	       "Chills","Night sweats","Fever","lab2nd") 

# plot the data
freq_2nd_lab <- freq[freq$lab2nd==1,]
for (i in 1:nrow(freq_2nd_lab)){
freq_2nd_lab$asy[i] <- sum(freq_2nd_lab[i,c(3:21)])           #Asymptomatic
}
symp_freq_2nd_lab <- rep(NA,20)
for (i in 2:20-1){
     symp_freq_2nd_lab[i] <- sum(na.exclude(freq_2nd_lab[,c(i+2)]))/length(na.exclude(freq_2nd_lab[,c(i+2)]))  
}
symp_freq_2nd_lab[20] <- length(freq_2nd_lab$asy[!is.na(freq_2nd_lab$asy)&freq_2nd_lab$asy==0])/nrow(freq_2nd_lab)

names(symp_freq_2nd_lab) <- c("Headache","Swollen lymph nodes","Runny nose","Hoarseness of voice","Sore throat",
               "Cough","Phlegm","Shortness of breath","Chest pain","Abdominal pain","Vomiting","Diarrhoea",
	       "Stiffness in muscles or joints","Aches or pains in muscles or joints","Dizziness","Fatigue/tiredness",
	       "Chills","Night sweats","Fever","Asymptomatic" )

plotdata <- sort(symp_freq_2nd_lab[1:19])
plotdata[20] <- symp_freq_2nd_lab[20]
names(plotdata)[20] <- "Asymptomatic"

windows(width=5, height=4)
par(mar=c(5,15,2,2), xaxs="i", yaxs="i")
barplot(rev(plotdata), names.arg=FALSE, horiz=T, xlab="Proportion of subjects reporting symptom",
        axes=FALSE,xlim=(c(0,0.8)),col=rep(c("gray","blue"),c(1,19)))
axis(1, cex.axis=1.5)
axis(2, at=20:1*1.2-0.5,labels=names(plotdata),las=1.5)


