#########################################################
#                                                       #
# SAR functions for lab-confirmed / clinical influenza   #
# Data stratified according to onset - intervention     #
#                                                       #
#########################################################

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv", header=TRUE)
housechar <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_20_baseflu_m.csv", header=TRUE)

###---------------------------------------------------- Lab-confirmed secondary cases ------------------------------------------------###

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

## exd_index: none of V0/V1 culture is A/B; d_index: both V0 and V1 culture is 0; Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
				  | (is.na(hculture$V1[i]) & hculture$V0[i]==0) | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
              {hculture$exd_index[i]<-1}    else {hculture$exd_index[i]<-0}
     if(hculture$member[i]==0 &  !is.na(hculture$V0[i]) & !is.na(hculture$V1[i]) & hculture$V0[i]==0 & hculture$V1[i]==0)
              {hculture$d_index[i]<-1}      else {hculture$d_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]==1) ))
              {hculture$exd_contact[i]=1}   else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]   # for calculating secondary cases
d_index <- hculture[hculture$member==0,c(1,9)]     # for calculating SAR

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
hculture <- merge(hculture[,-8], d_index)
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

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

hh <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_demographic_m.csv", header=TRUE)

### Add variable 'vaccine' <- 1 if the subject received vaccination in any of the past 2 years
for ( i in 1:nrow(hculture)){
    if ( (!is.na(baseflu$vaccine07[i]) & baseflu$vaccine07[i] == 1) )#| (!is.na(baseflu$vaccine06[i]) & baseflu$vaccine06[i] == 1 ))
          {hculture$vaccine[i] <-1}    else {hculture$vaccine[i] <- 0}
}

### Add variable 'indexsex' & 'indexage'
j <- 1
i <- 1
while (i < 1+ nrow(hculture)){  
  if(as.character(hculture$hhID[i])==as.character(hh$hhID[j])){
       hculture$indexsex[i] <- hh$male[j]
       hculture$indexage[i] <- hh$age[j]
       i <- i+1
       while(as.character(hculture$hhID[i])==as.character(hculture$hhID[i-1])){
       hculture$indexsex[i] <- hculture$indexsex[i-1]
       hculture$indexage[i] <- hculture$indexage[i-1]
       i <- i+1
       } 
  }
  j <- j+1
  while(as.character(hculture$hhID[i])<as.character(hh$hhID[j])){i <- i+1}
}
hculture$indexagegp <- cut(hculture$indexage, c(-0.1,15,100))

### Add variable 'intervention' and 'indexdelay' (time from symptom onset to intervention)
j<-1
i<-1
while (i<1+nrow(hculture)){
     while(as.character(hculture$hhID[i])==as.character(housechar$hhID[j])){
        hculture$arm[i] <- housechar$intervention[j]
	hculture$delay[i] <- housechar$v1_day[j]
	i <- i+1
     }
     j <- j+1
}
hculture$arm <- factor(hculture$arm, levels=1:3, labels=c("control", "mask", "hand"))
hculture$indexdelay <- 1*(hculture$delay<=1)

### Add variable 'age' and 'sex' (for household contact)
hculture$age <- hh$age
hculture$sex <- hh$male
hculture$agegp <- cut(hculture$age,c(-0.1,15,100))

#####################################  Compute SAR based on lab-confirmed results  ########################################################

library(gee)

hculture.sar.ds <- hculture[hculture$member!=0 & hculture$d_index == 0 & !is.na(hculture$indexdelay) & hculture$indexdelay==1,]
inf.gee <- gee(labsedcase~arm+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hculture.sar.ds, corstr = "exchangeable", family="binomial")
summary(inf.gee)
results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2)                 # Odds ratio and 95% CI

hculture.sar.db <- hculture[hculture$member!=0 & hculture$d_index == 0 & !is.na(hculture$indexdelay) & hculture$indexdelay==0,]
inf.gee <- gee(labsedcase~arm+agegp+sex+indexagegp+indexsex, id=factor(hhID),   # no lab2nd case received vaccination before
  data=hculture.sar.db, corstr = "exchangeable", family="binomial")
summary(inf.gee)
results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2) 

###################################### n ##################################################################################################

new <- hculture.sar.ds[!is.na(hculture.sar.ds$agegp),]  # For appendix table 1
new <- hculture.sar.db[!is.na(hculture.sar.db$agegp),]  # For appendix table 2

dim(new[new$arm=="control",])[1]
dim(new[new$arm=="mask",])[1]
dim(new[new$arm=="hand",])[1]

dim(new[new$age<=15,])[1]
dim(new[new$age>15,])[1]

dim(new[new$sex==1,])[1]
dim(new[new$sex==0,])[1]

dim(new[new$vaccine==0,])[1]
dim(new[new$vaccine==1,])[1]

length(unique(new$hhID[new$indexage<=15]))
length(unique(new$hhID[new$indexage>15]))

length(unique(new$hhID[new$indexsex==1]))
length(unique(new$hhID[new$indexsex==0]))


###----------------------------------------------------- Clinic secondary cases --------------------------------------------------------###

symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)

###########################################################################################################################################
# Case 1 (pilot protocol)
symptom$fever <- rep(NA,)
for (i in 1:dim(symptom)[1]){
 if(!is.na(symptom$bodytemp[i])){
    if(as.numeric(symptom$bodytemp[i])>=38) 
      {symptom$fever[i] <- 1} else {symptom$fever[i] <- 0}
 }
}
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

# Case 3 (WHO)
symptom$flu <- rep(NA,)
for (i in 1:dim(symptom)[1]){
    symptom$indicate[i]=sum(symptom$cough[i],symptom$sthroat[i])
    if( (!is.na(symptom$fever[i]) & symptom$fever[i]==1)
          & (!is.na(symptom$indicate[i]) & symptom$indicate[i]>=1) ) symptom$flu[i] <- 1
       else if ( (!is.na(symptom$fever[i]) & symptom$fever[i]==0)
                 | (!is.na(symptom$indicate[i]) & symptom$indicate[i]==0) ) symptom$flu[i] <- 0
}

###########################################################################################################################################

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
hclinic$d_index <- hculture$d_index

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

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################


hclinic$arm <- hculture$arm
hclinic$indexagegp <- hculture$indexagegp
hclinic$indexsex <- hculture$indexsex
hclinic$agegp <- hculture$agegp
hclinic$sex <- hculture$sex
hclinic$vaccine <- hculture$vaccine
hclinic$indexdelay <- hculture$indexdelay

#####################################  Compute SAR based on lab-confirmed results  ######################################################

library(gee)

hclinic.sar.ds <- hclinic[hclinic$member!=0 & hclinic$d_index == 0 & !is.na(hclinic$indexdelay) & hclinic$indexdelay==1,]
hclinic.sar.ds[hclinic.sar.ds$clinicsedcase==1,]
inf.gee <- gee(clinicsedcase~arm+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hclinic.sar.ds, corstr = "exchangeable", family="binomial")
summary(inf.gee)
results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2) 

hclinic.sar.db <- hclinic[hclinic$member!=0 & hclinic$d_index == 0 & !is.na(hclinic$indexdelay) & hclinic$indexdelay==0,]
hclinic.sar.db[hclinic.sar.db$clinicsedcase==1,]

inf.gee <- gee(clinicsedcase~arm+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hclinic.sar.db, corstr = "exchangeable", family="binomial")
summary(inf.gee)
results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2) 

