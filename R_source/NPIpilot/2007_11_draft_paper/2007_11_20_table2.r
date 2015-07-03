#
# Table 2: description of household analysed
#

source ("D:\\Work\\SVNrepository\\NPIpilot\\2007_11_draft_paper\\2007_11_20_fig1.r")

demog <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_demographic_m.csv",header=TRUE)
demog <- demog[,c(1,3:5)]         # hhID, member, age, sex
analysis <- analysis[,c(1,3,26)]  # hhID, intervention, baseline
tab2 <- merge(demog,analysis)
tab2_als <- tab2[tab2$baseline==0,]

tab2_c <- tab2_als[tab2_als$intervention==1,]
tab2_m <- tab2_als[tab2_als$intervention==2,]
tab2_h <- tab2_als[tab2_als$intervention==3,]


## Index subject
#
### age group ##
#length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=15]))
#round(length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=15]))/length(na.exclude(tab2_c$age[tab2_c$member==0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=30&tab2_c$age>15]))
#round(length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=30&tab2_c$age>15]))/length(na.exclude(tab2_c$age[tab2_c$member==0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=50&tab2_c$age>30]))
#round(length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age<=50&tab2_c$age>30]))/length(na.exclude(tab2_c$age[tab2_c$member==0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age>50]))
#round(length(na.exclude(tab2_c$age[tab2_c$member==0&tab2_c$age>50]))/length(na.exclude(tab2_c$age[tab2_c$member==0])),2)
#
#length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=15]))
#round(length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=15]))/length(na.exclude(tab2_m$age[tab2_m$member==0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=30&tab2_m$age>15]))
#round(length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=30&tab2_m$age>15]))/length(na.exclude(tab2_m$age[tab2_m$member==0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=50&tab2_m$age>30]))
#round(length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age<=50&tab2_m$age>30]))/length(na.exclude(tab2_m$age[tab2_m$member==0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age>50]))
#round(length(na.exclude(tab2_m$age[tab2_m$member==0&tab2_m$age>50]))/length(na.exclude(tab2_m$age[tab2_m$member==0])),2)
#
#length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=15]))
#round(length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=15]))/length(na.exclude(tab2_h$age[tab2_h$member==0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=30&tab2_h$age>15]))
#round(length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=30&tab2_h$age>15]))/length(na.exclude(tab2_h$age[tab2_h$member==0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=50&tab2_h$age>30]))
#round(length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age<=50&tab2_h$age>30]))/length(na.exclude(tab2_h$age[tab2_h$member==0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age>50]))
#round(length(na.exclude(tab2_h$age[tab2_h$member==0&tab2_h$age>50]))/length(na.exclude(tab2_h$age[tab2_h$member==0])),2)
#
### sex ##
#sum(tab2_c$male[tab2_c$member==0], na.rm = TRUE)
#round(sum(tab2_c$male[tab2_c$member==0], na.rm = TRUE)/length(na.exclude(tab2_c$male[tab2_c$member==0])),2)
#sum(tab2_m$male[tab2_m$member==0], na.rm = TRUE)
#round(sum(tab2_m$male[tab2_m$member==0], na.rm = TRUE)/length(na.exclude(tab2_m$male[tab2_m$member==0])),)
#sum(tab2_h$male[tab2_h$member==0], na.rm = TRUE)
#round(sum(tab2_h$male[tab2_h$member==0], na.rm = TRUE)/length(na.exclude(tab2_h$male[tab2_h$member==0])),)
#
## Household contact
#
### age group ##
#length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=15]))
#round(length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=15]))/length(na.exclude(tab2_c$age[tab2_c$member>0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=30&tab2_c$age>15]))
#round(length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=30&tab2_c$age>15]))/length(na.exclude(tab2_c$age[tab2_c$member>0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=50&tab2_c$age>30]))
#round(length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age<=50&tab2_c$age>30]))/length(na.exclude(tab2_c$age[tab2_c$member>0])),2)
#length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age>50]))
#round(length(na.exclude(tab2_c$age[tab2_c$member>0&tab2_c$age>50]))/length(na.exclude(tab2_c$age[tab2_c$member>0])),2)
#
#length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=15]))
#round(length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=15]))/length(na.exclude(tab2_m$age[tab2_m$member>0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=30&tab2_m$age>15]))
#round(length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=30&tab2_m$age>15]))/length(na.exclude(tab2_m$age[tab2_m$member>0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=50&tab2_m$age>30]))
#round(length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age<=50&tab2_m$age>30]))/length(na.exclude(tab2_m$age[tab2_m$member>0])),2)
#length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age>50]))
#round(length(na.exclude(tab2_m$age[tab2_m$member>0&tab2_m$age>50]))/length(na.exclude(tab2_m$age[tab2_m$member>0])),2)
#
#length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=15]))
#round(length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=15]))/length(na.exclude(tab2_h$age[tab2_h$member>0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=30&tab2_h$age>15]))
#round(length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=30&tab2_h$age>15]))/length(na.exclude(tab2_h$age[tab2_h$member>0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=50&tab2_h$age>30]))
#round(length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age<=50&tab2_h$age>30]))/length(na.exclude(tab2_h$age[tab2_h$member>0])),2)
#length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age>50]))
#round(length(na.exclude(tab2_h$age[tab2_h$member>0&tab2_h$age>50]))/length(na.exclude(tab2_h$age[tab2_h$member>0])),2)
#
### sex ##
#sum(tab2_c$male[tab2_c$member>0], na.rm = TRUE)
#round(sum(tab2_c$male[tab2_c$member>0], na.rm = TRUE)/length(na.exclude(tab2_c$male[tab2_c$member>0])),2)
#sum(tab2_m$male[tab2_m$member>0], na.rm = TRUE)
#round(sum(tab2_m$male[tab2_m$member>0], na.rm = TRUE)/length(na.exclude(tab2_m$male[tab2_m$member>0])),2)
#sum(tab2_h$male[tab2_h$member>0], na.rm = TRUE)
#round(sum(tab2_h$male[tab2_h$member>0], na.rm = TRUE)/length(na.exclude(tab2_h$male[tab2_h$member>0])),2)
#
## vaccination history
#baseflu$vaccine <- 1*((!is.na(baseflu$vaccine07) & baseflu$vaccine07==1) | (!is.na(baseflu$vaccine06) & baseflu$vaccine06==1))
#tab2$vaccine <- baseflu$vaccine
#tab2_als_vac <- tab2[tab2$baseline==0,]
#
#tab2_c <- tab2_als_vac[tab2_als_vac$intervention==1,]
#tab2_m <- tab2_als_vac[tab2_als_vac$intervention==2,]
#tab2_h <- tab2_als_vac[tab2_als_vac$intervention==3,]
#
#sum(tab2_c$vaccine[tab2_c$member>0], na.rm = TRUE)
#round(sum(tab2_c$vaccine[tab2_c$member>0], na.rm = TRUE)/length(na.exclude(tab2_c$vaccine[tab2_c$member>0])),2)
#sum(tab2_m$vaccine[tab2_m$member>0], na.rm = TRUE)
#round(sum(tab2_m$vaccine[tab2_m$member>0], na.rm = TRUE)/length(na.exclude(tab2_m$vaccine[tab2_m$member>0])),2)
#sum(tab2_h$vaccine[tab2_h$member>0], na.rm = TRUE)
#round(sum(tab2_h$vaccine[tab2_h$member>0], na.rm = TRUE)/length(na.exclude(tab2_h$vaccine[tab2_h$member>0])),2)


############################ For secondary cases calculation ##############################################################################

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

###---------------------------------------------- Lab confirmed ------------------------------------------------------------------------###

hculture <- read.csv("P:\\GROUP\\NPIpilot\\data\\lab2nd.csv",header=TRUE)

j<-1
i<-1
while (i<1+nrow(hculture)){
     while(as.character(hculture$hhID[i])==as.character(hchar$hhID[j])){
        hculture$arm[i]<-hchar$intervention[j]
	hculture$delay[i] <- hchar$v1_day[j]
	i <- i+1
     }
     j <- j+1
}
hculture$arm <- factor(hculture$arm, levels=1:3, labels=c("control", "mask", "hand"))
hculture$indexdelay <- 1*(hculture$delay<=1)
hculture <- merge(hculture,hchar[,c(1,5)])
names(hculture)[dim(hculture)[2]]<- "hhsize"
hculture$hhsize <- hculture$hhsize -1
hculture$hhsize[hculture$member!=1] <- 0

hc_contact <- hculture[hculture$member>0&hculture$d_index==0,]

########################################################### Full data #####################################################################

c2nd_all <- dim(hc_contact[hc_contact$arm=="control",])[1]
m2nd_all <- dim(hc_contact[hc_contact$arm=="mask",])[1]
h2nd_all <- dim(hc_contact[hc_contact$arm=="hand",])[1]
c2nd_all 
m2nd_all 
h2nd_all 

c2nd <- dim(hc_contact[hc_contact$arm=="control"&hc_contact$labsedcase==1,])[1]
round(c2nd/c2nd_all,2)
round(binom.test(c2nd,c2nd_all)$conf[1:2],2)

m2nd <- dim(hc_contact[hc_contact$arm=="mask"&hc_contact$labsedcase==1,])[1]
round(m2nd/m2nd_all,2)
round(binom.test(m2nd,m2nd_all)$conf[1:2],2)

h2nd <- dim(hc_contact[hc_contact$arm=="hand"&hc_contact$labsedcase==1,])[1]
round(h2nd/h2nd_all,2)
round(binom.test(h2nd,h2nd_all)$conf[1:2],2)

## Adjusted chi-square test 
chi_stat <- icc(hc_contact$hhID,hc_contact$arm,hc_contact$hhsize,hc_contact$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

############################## Repeat above stats stratified by symptom onset-intervention <= 36 hrs #####################################

c2nd_all_s <- dim(hc_contact[hc_contact$arm=="control" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
m2nd_all_s <- dim(hc_contact[hc_contact$arm=="mask" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
h2nd_all_s <- dim(hc_contact[hc_contact$arm=="hand" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
c2nd_all_s 
m2nd_all_s 
h2nd_all_s 

c2nd <- dim(hc_contact[hc_contact$arm=="control"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
round(c2nd/c2nd_all_s,2)
round(binom.test(c2nd,c2nd_all_s)$conf[1:2],2)

m2nd <- dim(hc_contact[hc_contact$arm=="mask"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
round(m2nd/m2nd_all_s,2)
round(binom.test(m2nd,m2nd_all_s)$conf[1:2],2)

h2nd <- dim(hc_contact[hc_contact$arm=="hand"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
round(h2nd/h2nd_all_s,2)
round(binom.test(h2nd,h2nd_all_s)$conf[1:2],2)

## Adjusted chi-square test 
hc.lt <- hc_contact[!is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,] 
chi_stat <- icc(hc.lt$hhID,hc.lt$arm,hc.lt$hhsize,hc.lt$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

############################### Repeat above stats stratified by symptom onset-intervention > 36 hrs ######################################

c2nd_all_g <- dim(hc_contact[hc_contact$arm=="control" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
m2nd_all_g <- dim(hc_contact[hc_contact$arm=="mask" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
h2nd_all_g <- dim(hc_contact[hc_contact$arm=="hand" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
c2nd_all_g 
m2nd_all_g 
h2nd_all_g 

c2nd <- dim(hc_contact[hc_contact$arm=="control"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(c2nd/c2nd_all_g,2)
round(binom.test(c2nd,c2nd_all_g)$conf[1:2],2)

m2nd <- dim(hc_contact[hc_contact$arm=="mask"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(m2nd/m2nd_all_g,2)
round(binom.test(m2nd,m2nd_all_g)$conf[1:2],2)

h2nd <- dim(hc_contact[hc_contact$arm=="hand"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(h2nd/h2nd_all,2)
round(binom.test(h2nd,h2nd_all_g)$conf[1:2],2)

## Adjusted chi-square test 
hc.gt <- hc_contact[!is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,] 
chi_stat <- icc(hc.gt$hhID,hc.gt$arm,hc.gt$hhsize,hc.gt$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

###---------------------------------------------- 3 clinic definitions -----------------------------------------------------------------###

###################################### Repeat part of scprits for table 3 ################################################################

symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)

######################################################################################
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

#####################################################################################

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

hclinic$arm <- hculture$arm
hclinic$d_index <- hculture$d_index
hclinic$indexdelay <- hculture$indexdelay
hclinic$hhsize <- hculture$hhsize
clicon <- hclinic[hclinic$member>0&hclinic$d_index==0,]

########################################## End of copy of scripts ########################################################################

########################################################### Full data ####################################################################

c2nd <- dim(clicon[clicon$arm=="control"&clicon$clinicsedcase==1,])[1]
round(c2nd/c2nd_all,2)
round(binom.test(c2nd,c2nd_all)$conf[1:2],2)

m2nd <- dim(clicon[clicon$arm=="mask"&clicon$clinicsedcase==1,])[1]
round(m2nd/m2nd_all,2)
round(binom.test(m2nd,m2nd_all)$conf[1:2],2)

h2nd <- dim(clicon[clicon$arm=="hand"&clicon$clinicsedcase==1,])[1]
round(h2nd/h2nd_all,2)
round(binom.test(h2nd,h2nd_all)$conf[1:2],2)

## Adjusted chi-square test 
chi_stat <- icc(clicon$hhID,clicon$arm,clicon$hhsize,clicon$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test


#################################### Repeat above stats stratified by symptom onset-intervention <= 36 hrs ################################

c2nd <- dim(clicon[clicon$arm=="control"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==1,])[1]
round(c2nd/c2nd_all_s,2)
round(binom.test(c2nd,c2nd_all_s)$conf[1:2],2)

m2nd <- dim(clicon[clicon$arm=="mask"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==1,])[1]
round(m2nd/m2nd_all_s,2)
round(binom.test(m2nd,m2nd_all_s)$conf[1:2],2)

h2nd <- dim(clicon[clicon$arm=="hand"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==1,])[1]
round(h2nd/h2nd_all_s,2)
round(binom.test(h2nd,h2nd_all_s)$conf[1:2],2)

## Adjusted chi-square test 
clicon.lt <- clicon[!is.na(clicon$indexdelay) & clicon$indexdelay==1,] 
chi_stat <- icc(clicon.lt$hhID,clicon.lt$arm,clicon.lt$hhsize,clicon.lt$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test

###################################### Repeat above stats stratified by symptom onset-intervention > 36 hrs ###############################

c2nd <- dim(clicon[clicon$arm=="control"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==0,])[1]
round(c2nd/c2nd_all_g,2)
round(binom.test(c2nd,c2nd_all_g)$conf[1:2],2)

m2nd <- dim(clicon[clicon$arm=="mask"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==0,])[1]
round(m2nd/m2nd_all_g,2)
round(binom.test(m2nd,m2nd_all_g)$conf[1:2],2)

h2nd <- dim(clicon[clicon$arm=="hand"&clicon$clinicsedcase==1 & !is.na(clicon$indexdelay) & clicon$indexdelay==0,])[1]
round(h2nd/h2nd_all_g,2)
round(binom.test(h2nd,h2nd_all_g)$conf[1:2],2)

## Adjusted chi-square test 
clicon.gt <- clicon[!is.na(clicon$indexdelay) & clicon$indexdelay==0,] 
chi_stat <- icc(clicon.gt$hhID,clicon.gt$arm,clicon.gt$hhsize,clicon.gt$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test
