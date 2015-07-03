#########################################################
#                                                       #
# SAR functions for  lab-confirmed / clinical influenza #
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
# write.csv(hculture,"D:\\Work\\Influenza\\output NPI year 1\\2007_11_draft\\lab2nd.csv",row.names=FALSE)

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
	hculture$flatsize[i] <- housechar$house_size[j]
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
hculture$occup <- hh$occupation
hculture$edu <- hh$education

#####################################  Compute SAR based on lab-confirmed results  ########################################################

library(gee)
hculture.sar <- hculture[hculture$member!=0 & hculture$d_index == 0,]

inf.gee <- gee(labsedcase~arm+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hculture.sar, corstr = "exchangeable", family="binomial")

summary(inf.gee)

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2)                 # Odds ratio and 95% CI

###################################### n ##################################################################################################

new <- hculture.sar[!is.na(hculture.sar$agegp),]

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

##################################### Calculate intracluster correlation coefficient , for lab secondary case ############################

hculture.cc <- merge(hculture,housechar[,c(1,5)])
names(hculture.cc)[dim(hculture.cc)[2]]<- "hhsize"
hculture.cc <- hculture.cc[hculture.cc$member!=0&hculture.cc$d_index==0,]
hculture.cc$hhsize <- hculture.cc$hhsize -1
hculture.cc$hhsize[hculture.cc$member>1] <- 0

## Begin function 'icc' for intracluster correlation coefficient
icc <- function(data_id,data_arm,data_clustersize,data_success){
 M <- length(data_id)
 Y <- sum(data_success)
 P <- Y/M
 K <- length(unique(data_id))

 arm <- c("control","mask","hand")   # interventions
 Yij <- Mij <- Pij <- hhID <- matrix(nrow=3,ncol=length(unique(data_id[data_arm==arm[1]])))
 for (i in 1:length(arm)){
      hhID[i,1:length(unique(data_id[data_arm==arm[i]]))] <- as.character(unique(data_id[data_arm==arm[i]]))
 }

 for (i in 1:length(arm)){
      for (j in 1:length(unique(data_id[data_arm==arm[i]]))){
        Yij[i,j] <- sum(data_success[data_id==hhID[i,j]])
        Mij[i,j] <- sum(data_clustersize[data_id==hhID[i,j]])
        Pij[i,j] <- Yij[i,j]/Mij[i,j]
      }
 }

 Y_i <- M_i <- P_i <- m_Ai <- rho_i <- rep(NA,length(arm))
 output <- data.frame(Y_i,M_i,m_Ai,P_i,rho_i)
 for (i in 1:length(arm)){
     output$Y_i[i] <- sum(na.exclude(Yij[i,]))
     output$M_i[i] <- sum(na.exclude(Mij[i,]))
     output$P_i[i] <- output$Y_i[i]/output$M_i[i]
     output$m_Ai[i] <- sum(na.exclude(Mij[i,]^2))/output$M_i[i]
 }

 m_0 <- (M-sum(output$m_Ai))/(K-2)
 MSW <- sum(rowSums(Mij*Pij*(1-Pij),na.rm = TRUE)) / (M-K)
 MSC <- sum(rowSums(Mij*(Pij-output$P_i)*(Pij-output$P_i),na.rm = TRUE)) / (K-2)
 rho <- (MSC-MSW)/(MSC+(m_0-1)*MSW)
 output$rho_i <- rho
 output
}
## End of function 'icc'

## Begin bootstrap loop
idlist <- unique(hculture.cc$hhID)
h <- length(idlist)
B <- 1000
icc.boot <- rep(NA,B)
for (b in 1:B){ 
	  newidlist <- sample(idlist,replace=TRUE)
	  newdata <- hculture.cc[hculture.cc$hhID==newidlist[1],]
	  for (i in 2:h){
	       newdata <- rbind(newdata,hculture.cc[hculture.cc$hhID==newidlist[i],])            
	  }
          icc.boot[b] <- icc(newdata$hhID,newdata$arm,newdata$hhsize,newdata$labsedcase)$rho_i[1]
}
## End of bootstrap loop

icc(hculture.cc$hhID,hculture.cc$arm,hculture.cc$hhsize,hculture.cc$labsedcase)$rho_i[1]
quantile(icc.boot,c(0.5,0.025,0.975))

## Chi-square test

chi_stat <- icc(hculture.cc$hhID,hculture.cc$arm,hculture.cc$hhsize,hculture.cc$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
round(achi_lab_m <- sum(chi[c(1,2)]),2)
round(pchisq(achi_lab_m,df=1,lower.tail=FALSE),2)
round(achi_lab_h <- sum(chi[c(1,3)]),2)
round(pchisq(achi_lab_h,df=1,lower.tail=FALSE),2)

## End of chi-square test

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
length(unique(hclinic$hhID[hclinic$clinicsedcase==1]))
###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

hclinic$arm <- hculture$arm
hclinic$indexagegp <- hculture$indexagegp
hclinic$indexsex <- hculture$indexsex
hclinic$agegp <- hculture$agegp
hclinic$sex <- hculture$sex
hclinic$vaccine <- hculture$vaccine
hclinic$occup <- hculture$occup
hclinic$edu <- hculture$edu
hclinic$flatsize <- hculture$flatsize

#####################################  Compute SAR based on lab-confirmed results  ######################################################

library(gee)
hclinic.sar <- hclinic[hclinic$member!=0 & hclinic$d_index == 0,]

inf.gee <- gee(clinicsedcase~arm+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hclinic.sar, corstr = "exchangeable", family="binomial")

summary(inf.gee)

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2) 

##################################### Calculate intracluster correlation coefficient , for lab secondary case ############################

hclinic.cc <- merge(hclinic,housechar[,c(1,5)])
names(hclinic.cc)[dim(hclinic.cc)[2]]<- "hhsize"
hclinic.cc <- hclinic.cc[hclinic.cc$member!=0&hclinic.cc$d_index==0,]
hclinic.cc$hhsize[hclinic.cc$member>1] <- 0

## Begin bootstrap loop
idlist <- unique(hclinic.cc$hhID)
h <- length(idlist)
B <- 1000
icc.boot <- rep(NA,B)
for (b in 1:B){ 
	  newidlist <- sample(idlist,replace=TRUE)
	  newdata <- hclinic.cc[hclinic.cc$hhID==newidlist[1],]
	  for (i in 2:h){
	       newdata <- rbind(newdata,hclinic.cc[hclinic.cc$hhID==newidlist[i],])            
	  }
          icc.boot[b] <- icc(newdata$hhID,newdata$arm,newdata$hhsize,newdata$labsedcase)$rho_i[1]
}
## End of bootstrap loop

icc(hclinic.cc$hhID,hclinic.cc$arm,hclinic.cc$hhsize,hclinic.cc$clinicsedcase)$rho_i[1]
quantile(icc.boot,c(0.5,0.025,0.975))

## Chi-square test

chi_stat <- icc(hclinic.cc$hhID,hclinic.cc$arm,hclinic.cc$hhsize,hclinic.cc$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
round(achi_c_m <- sum(chi[c(1,2)]),2)
round(pchisq(achi_c_m,df=1,lower.tail=FALSE),2)
round(achi_c_h <- sum(chi[c(1,3)]),2)
round(pchisq(achi_c_h,df=1,lower.tail=FALSE),2)

## End of chi-square test
