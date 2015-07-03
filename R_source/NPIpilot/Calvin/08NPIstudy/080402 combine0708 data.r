require(chron)
source ("D:/Work/SVNrepository/NPIpilot/Calvin/08NPIstudy/071008 reformat function.r")

#input and output directory and today's date
############################
today <- "2008_4_11"
input_dir07 <-"P:/GROUP/NPIpilot/data/"
input_dir08 <-"P:/GROUP/NPIstudy/data/"
output_dir <- "P:/GROUP/NPIstudy/data/combined/"
############################


#function to combine 07 and 08 data, add the year to screening number
combine0708 <-function(yr07,yr08){
  yr07$scrID <- paste("07", yr07$scrID, sep="") #add screening number
  yr08$scrID <- paste("08", yr08$scrID, sep="")
  yr08$hhID <- paste("h080", yr08$hhID, sep="")
  return (rbind(yr07,yr08))
}
# end

#For 1 houshold as 1 entry
#combine 07 and 08 clinicdat_m 
clinic07 <- read.csv(paste(input_dir07, "2007_11_28_clinicdat_h.csv", sep=""), header=TRUE)
clinic08 <- read.csv(paste(input_dir08, "2008_4_11_clinicdat_h.csv",sep=""), header=TRUE)

clinic07 <- clinic07[c(1:8,12,13,21,10,14,26:28,33,34)]
names(clinic07)[11] <-"aches"
clinic <- combine0708(clinic07,clinic08)
write.csv(clinic, paste(output_dir, today,"_clinicdat0708_h.csv", sep=""), row.names=FALSE)
#

#combine household characteristics
house07 <- read.csv(paste(input_dir07, "2007_11_28_housechar_h.csv", sep=""), header=TRUE)
house07 <-house07[c(1:3,5:11,22:25)]
names(house07)[4] <-"familysize"
house08 <- read.csv(paste(input_dir08,"2008_4_11_housechar_h.csv", sep=""), header=TRUE)
#onset08<-onset08[-c(6)]

house <- combine0708(house07,house08)
write.csv(house, paste(output_dir, today,"_housechar0708_h.csv", sep=""), row.names=FALSE)
#

#For 1 houshold member as 1 entry
#combine 07 and 08 demographics
demo07 <- read.csv(paste(input_dir07,"2007_11_27_demographic_m.csv",sep=""), header=TRUE)
demo08 <- read.csv(paste(input_dir08,"2008_4_11_demographic_m.csv",sep=""), header=TRUE)

demographic <- combine0708(demo07,demo08)

write.csv(demographic, paste(output_dir, today,"_demo0708_m.csv", sep=""), row.names=FALSE)

#combine 07 and 08 baseflu_m 
baseflu07 <- read.csv(paste(input_dir07,"2007_11_20_baseflu_m.csv", sep=""), header=TRUE)
#baseflu07 <- baseflu07[,c(1,3,7:9)]
names(baseflu07)[c(7,8,9)] <-c("vac1","vac2","vac3")

baseflu08 <- read.csv(paste(input_dir08,"2008_4_11_baseflu_m.csv", sep=""), header=TRUE)
#baseflu08 <- baseflu08[,c(1,3,7:9)]
names(baseflu08)[c(7,8,9)] <-c("vac1","vac2","vac3")

baseflu <- combine0708(baseflu07,baseflu08)
write.csv(baseflu, paste(output_dir, today,"_baseflu0708_m.csv", sep=""), row.names=FALSE)
#



#combine antiviral usage 07 and 08
av07 <- read.csv(paste(input_dir07,"080411av07.csv", sep=""), header=TRUE)
av08 <- read.csv(paste(input_dir08, "080411av08.csv", sep=""), header=TRUE)
av <- combine0708(av07,av08)
#av <- av[c(1:4)]

write.csv(av, paste(output_dir, today,"_av0708_m.csv", sep=""), row.names=FALSE)
#

# For clinically defined 2nd cases 
clin2nd07 <- read.csv(paste(input_dir07,"clinic2nd2007.csv", sep=""), header=TRUE) #last year
clin2nd07 <- clin2nd07[c(1,2)]
clin2nd08 <- read.csv(paste(input_dir08,"clinic2nd2008.csv", sep=""), header=TRUE) #this year
clin2nd08 <- clin2nd08[c(1,2)]
clin2nd <- combine0708(clin2nd07,clin2nd08)

write.csv(clin2nd, paste(output_dir, today,"_clin2nd0708_m.csv", sep=""), row.names=FALSE)
#

# For laboratory defined 2nd cases 
lab2nd07 <- read.csv(paste(input_dir07,"lab2nd2007.csv", sep=""), header=TRUE) #last year
lab2nd07 <- lab2nd07[lab2nd07$labsedcase ==1,]
lab2nd07 <- lab2nd07[c(1,2)]

lab2nd08 <- read.csv(paste(input_dir08,"lab2nd2008.csv", sep=""), header=TRUE) #this year
lab2nd08 <- lab2nd08[lab2nd08$labsedcase ==1,]
lab2nd08 <- lab2nd08[c(1,2)]
lab2nd <- combine0708(lab2nd07,lab2nd08)

write.csv(lab2nd, paste(output_dir, today,"_lab2nd0708_m.csv", sep=""), row.names=FALSE)



#For 1 houshold member in 1 day as 1 entry
#combine 07 and 08 symptomday_d 
symptom07 <- read.csv(paste(input_dir07,"2007_11_27_symptomday_d.csv", sep=""), header=TRUE) #last year
symptom08 <- read.csv(paste(input_dir08,"2008_4_11_symptomday_d.csv", sep=""), header=TRUE) #this year

#modify symptom07 to symptom08 format and combine
symptom07 <- symptom07[,c(1:6,10, 11,19, 8, 12)]

symptom <- combine0708(symptom07,symptom08)
#del column with no data
symptom <- symptom[!is.na(symptom[11]) & !is.na(symptom[10]) & !is.na(symptom[9])
& !is.na(symptom[8]) & !is.na(symptom[7]) & !is.na(symptom[6]),]



write.csv(symptom, paste(output_dir, today,"_symptomday0708_d.csv", sep=""), row.names=FALSE)
#

#combine pcr data
pcr07 <- read.csv(paste(input_dir07,"2008_01_home_PCR.csv", sep=""), header=TRUE) #last year
pcr07$flu.type <- pcr07$culture
pcr07$flu.type[(!is.na(pcr07$AqPCR) & pcr07$AqPCR !=0) | (!is.na(pcr07$AqPCR2) & pcr07$AqPCR2 !=0) ] <-"A"
pcr07$flu.type[!is.na(pcr07$BqPCR) & pcr07$BqPCR!=0] <-"B"
pcr07$qPCR <- pcr07$AqPCR
#pcr07$qPCR[is.na(pcr07$qPCR)] <-0


# For other format to combine
#combine pcr results
for (i in 1:dim(pcr07)[1]){
  if (!is.na(pcr07$AqPCR2[i]) & pcr07$AqPCR2[i] > pcr07$qPCR[i]){
    pcr07$qPCR[i] <- pcr07$AqPCR2[i]
  }
  if (!is.na(pcr07$BqPCR[i])){
    pcr07$qPCR[i] <- pcr07$BqPCR[i]    
  }
 }

#add flu type to other household members (some of them become 0, as the clinic culture is 0, need to correct them in other means
flutype <- pcr07$flu.type[1] #initial flu type
for (i in 2:dim(pcr07)[1]){
  if (pcr07$hhID[i] == pcr07$hhID[i-1]) { 
    pcr07$flu.type[i] <- flutype
  }
  else {flutype <- pcr07$flu.type[i]}
} 

#write.csv(pcr07,"d:\\r\\temp02.csv")


pcr07 <- pcr07[pcr07$visit !=0 & pcr07$visit !=4,] # delete clinic and 4th visit
pcr07 <- pcr07[c(1:3,5,11,12,7)]
names(pcr07)[4] <- "visitdate"
pcr08 <- read.csv(paste(input_dir08,"2008_04_11_home_culture.csv", sep=""), header=TRUE) #this year

flutype <- pcr08$flu.type[1] #initial flu type
for (i in 2:dim(pcr08)[1]){
  if (pcr08$hhID[i] == pcr08$hhID[i-1]) { 
    pcr08$flu.type[i] <- flutype
  }
  else {flutype <- pcr08$flu.type[i]}
} 

pcr <- combine0708(pcr07,pcr08)
pcr <- pcr[c(1:6)]
write.csv(pcr, paste(output_dir, today,"_pcr0708_o.csv", sep=""), row.names=FALSE)