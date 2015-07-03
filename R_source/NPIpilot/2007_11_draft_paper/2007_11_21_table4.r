
########################################################################################################
#                                                                                                      #
# ROC analysis for different clinical definitions on influenza vs gold standard (viral culture & qPCR) #
#                                                                                                      #
########################################################################################################

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv", header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_20_baseflu_m.csv", header=TRUE)

####################################################### Lab-confirmed secondary cases ####################################################

hculture <- data.frame(hhID = baseflu$hhID)
hculture$member <- 0
for(i in 2:nrow(hculture)){
  if(hculture$hhID[i]==hculture$hhID[i-1]) hculture$member[i] <- hculture$member[i-1]+1
}

j <- 1
i <- 1
while(i<nrow(hculture)+1){  
  while(as.character(hculture$hhID[i])==as.character(hc$hhID[j])){
       if(hculture$member[i]==hc$member[j]){
         if( (!is.na(hc$culture[j])&as.character(hc$culture[j])>0) | (!is.na(hc$AqPCR[j])&as.character(hc$AqPCR[j])>0) |
             (!is.na(hc$AqPCR2[j])&as.character(hc$AqPCR2[j])>0) | (!is.na(hc$AqPCR3[j])&as.character(hc$AqPCR3[j])>0) |
	     (!is.na(hc$BqPCR[j])&as.character(hc$BqPCR[j])>0) )
	     { hculture[i,hc$visit[j]+3] <- 1}
	     else { hculture[i,hc$visit[j]+3] <- 0}          
       j <- j+1
       } 
       else{i<-i+1}            
  }
  while(as.character(hculture$hhID[i])!=as.character(hc$hhID[j])){
       while(as.character(hculture$hhID[i])<as.character(hc$hhID[j])){i<-i+1}
       while(as.character(hculture$hhID[i])>as.character(hc$hhID[j])){j<-j+1}
  }
}
names(hculture) <- c("hhID","member","V0","V1","V2","V3","V4")

## exd_index: none of V0/V1 culture is A/B; Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
				  | (is.na(hculture$V1[i]) & hculture$V0[i]==0) | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
              {hculture$exd_index[i]<-1}    else {hculture$exd_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]==1) ))
              {hculture$exd_contact[i]=1}   else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]   # for calculating secondary cases

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
dim(hculture)

for (i in 1:nrow(hculture)){
    if ( hculture$exd_index[i]==1 | hculture$exd_contact[i] ==1)
       {hculture$exclude[i] <-1}
       else  {hculture$exclude[i] <-0}
}

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exclude[i] == 0 & !is.na(hculture$V1[i]) & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) | (hculture$V4[i] !=0 & !is.na(hculture$V4[i])) ) )
              {hculture$labsedcase[i] <- 1}
	 else {hculture$labsedcase[i] <- 0}
}

####################################################### Clinic secondary cases ###########################################################


symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)

symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=38) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}

#######################################################################################################
# Case 1 (pilot protocal)
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$cough[i],symptom$sthroat[i],symptom$rnose[i],
                        symptom$tired[i],symptom$headache[i],symptom$sjoint[i],symptom$pmuscle[i])
    if( (!is.na(symptom$fever[i]) & symptom$fever[i]==1)
          | (!is.na(symptom$indicate[i]) & symptom$indicate[i]>=2) ) symptom$flu[i] <- 1
       else if ( !is.na(symptom$fever[i]) & symptom$fever[i]==0
                 & !is.na(symptom$indicate[i]) & symptom$indicate[i]<2 ) symptom$flu[i] <- 0
}

# Case 2 (Monto 2002)
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$fever[i],symptom$cough[i],symptom$headache[i],symptom$sthroat[i],symptom$pmuscle[i])
    if( !is.na(symptom$indicate[i]) & symptom$indicate[i]>=2 ) symptom$flu[i] <- 1
       else if ( !is.na(symptom$indicate[i]) & symptom$indicate[i]<2 ) symptom$flu[i] <- 0
}

# Case 3 (WHO)
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$cough[i],symptom$sthroat[i])
    if( (!is.na(symptom$fever[i]) & symptom$fever[i]==1)
          & (!is.na(symptom$indicate[i]) & symptom$indicate[i]>=1) ) symptom$flu[i] <- 1
       else if ( (!is.na(symptom$fever[i]) & symptom$fever[i]==0)
                 | (!is.na(symptom$indicate[i]) & symptom$indicate[i]==0) ) symptom$flu[i] <- 0
}

# Case 4 (Monto's with fever >=37.8)
symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=37.8) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$fever[i],symptom$cough[i],symptom$headache[i],symptom$sthroat[i],symptom$pmuscle[i])
    if( !is.na(symptom$indicate[i]) & symptom$indicate[i]>=2 ) symptom$flu[i] <- 1
       else if ( !is.na(symptom$indicate[i]) & symptom$indicate[i]<2 ) symptom$flu[i] <- 0
}

# Case 5 (Hayden 2004)
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$cough[i],symptom$rnose[i])
    if( (!is.na(symptom$fever[i]) & symptom$fever[i]==1)
          & (!is.na(symptom$indicate[i]) & symptom$indicate[i]>=1) ) symptom$flu[i] <- 1
       else if ( (!is.na(symptom$fever[i]) & symptom$fever[i]==0)
                 | (!is.na(symptom$indicate[i]) & symptom$indicate[i]==0) ) symptom$flu[i] <- 0
}

# Case 6 (Welliver 2001)
symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=37.2) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$ind_1[i] <- sum(symptom$cough[i],symptom$rnose[i],symptom$sthroat[i])
    symptom$ind_2[i] <- sum(symptom$headache[i],symptom$tired[i],symptom$pmuscle[i],symptom$chill[i],symptom$nsweat[i])
    if( (!is.na(symptom$fever[i]) & symptom$fever[i]==1) & (!is.na(symptom$ind_1[i]) & symptom$ind_1[i]>=1)
          & (!is.na(symptom$ind_2[i]) & symptom$ind_2[i]>=1) ) symptom$flu[i] <- 1
       else if ( (!is.na(symptom$fever[i]) & symptom$fever[i]==0) | (!is.na(symptom$ind_1[i]) & symptom$ind_1[i]==0)
                 | (!is.na(symptom$ind_2[i]) & symptom$ind_2[i]==0) ) symptom$flu[i] <- 0
}

# Case 7 (fever)
symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=37.8) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    if( !is.na(symptom$fever[i]) & symptom$fever[i]==1 ) symptom$flu[i] <- 1
       else if ( !is.na(symptom$fever[i]) & symptom$fever[i]==0 ) symptom$flu[i] <- 0
}

# Case 8 (Monto's another definition)
symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=37.8)
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$fever[i],symptom$chill[i],symptom$headache[i],symptom$sthroat[i],
                            symptom$cpain[i],symptom$apain[i],symptom$pmuscle[i])
    if( (!is.na(symptom$cough[i]) & symptom$cough[i]==1)
          & (!is.na(symptom$indicate[i]) & symptom$indicate[i]>=1) ) symptom$flu[i] <- 1
       else if ( (!is.na(symptom$cough[i]) & symptom$cough[i]==0)
                 | (!is.na(symptom$indicate[i]) & symptom$indicate[i]==0) ) symptom$flu[i] <- 0
}

#######################################################################################################

hclinic <- data.frame(hhID = baseflu$hhID)
hclinic$member <- 0
for(i in 2:nrow(hclinic)){
  if(hclinic$hhID[i]==hclinic$hhID[i-1]) hclinic$member[i] <- hclinic$member[i-1]+1
}

j <- 1
i <- 1

## Read from dataset, one member per row
while(i<nrow(hclinic)+1){  
  while(as.character(hclinic$hhID[i])==as.character(symptom$hhID[j])){
       if(hclinic$member[i]==symptom$member[j]){
       hclinic[i,symptom$day[j]+3] <- symptom$flu[j]
       j <- j+1
       } 
       else{i<-i+1}            
  }
  while(as.character(hclinic$hhID[i])!=as.character(symptom$hhID[j])){
       while(as.character(hclinic$hhID[i])<as.character(symptom$hhID[j])){i<-i+1}
       while(as.character(hclinic$hhID[i])>as.character(symptom$hhID[j])){j<-j+1}
  }
}


names(hclinic) <- c("hhID","member","day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")

hclinic$exd_index <- hculture$exd_index
hclinic$exclude <- hculture$exclude

## Define secondary cases
for (i in 1:nrow(hclinic)){
    if ( hclinic$member[i] != 0 & hclinic$exclude[i] == 0 & !is.na(hclinic$day0[i]) &
         ( (!is.na(hclinic$day1[i]) & hclinic$day1[i]==1) | (!is.na(hclinic$day2[i]) & hclinic$day2[i]==1) | 
	   (!is.na(hclinic$day3[i]) & hclinic$day3[i]==1) | (!is.na(hclinic$day4[i]) & hclinic$day4[i]==1) | 
	   (!is.na(hclinic$day5[i]) & hclinic$day5[i]==1) | (!is.na(hclinic$day6[i]) & hclinic$day6[i]==1) | 
	   (!is.na(hclinic$day7[i]) & hclinic$day7[i]==1) | (!is.na(hclinic$day8[i]) & hclinic$day8[i]==1) | 
	   (!is.na(hclinic$day9[i]) & hclinic$day9[i]==1) ))
              {hclinic$clinicsedcase[i] <- 1}
	 else {hclinic$clinicsedcase[i] <- 0}
}
sum(hclinic$clinicsedcase)

############################################ Calculate sensitivity / specificity / Area under ROC ########################################

library(ROCR)

x <- hculture[,c(1,2,10,11)] # hhID,member, exclude, labsedcase 
y <- hclinic[,c(1,2,15)]
z <- merge(x,y)
zz <- z[z$member!=0&z$exclude==0,]

clinic <- zz$clinicsedcase
lab <- zz$labsedcase          # gold standard
p1 <- prediction(clinic, lab)
ss <- performance(p1, "sens", "spec")
round(ss@y.values[[1]][2],2)  # sensitivity
round(ss@x.values[[1]][2],2)  # specificity
p2 <- performance(p1, "auc")
round(p2@y.values[[1]],2)     # area under roc

################# Bootstrap for aur ########################

boot.auc <- function(reps,fludata){
    n <- length(fludata[,1])
    temp <- rep(NA,reps)
    output <- data.frame(auc=temp)
    for(i in 1:reps){
        a <- sample(1:n,replace=TRUE)
	temp.fludata <- fludata[a,]
        temp.clinic <- temp.fludata$clinicsedcase
	temp.lab <- temp.fludata$labsedcase
	temp.p1 <- prediction(temp.clinic,temp.lab)
	temp.p2 <- performance(temp.p1,"auc")
	output$auc[i] <- temp.p2@y.values[[1]]
    }
    output
}

flu.boot <- boot.auc(1000,zz)
round(quantile(flu.boot$auc,c(0.5,0.025,0.975)),2)

