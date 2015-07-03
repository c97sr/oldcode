
##########################################################################
#                                                                        #
#  Characteristics of 198 randomized index subjects by intervention arm  #
#                                                                        #
##########################################################################


q1data <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_28_clinicdat_h.csv",header=TRUE)
random <- read.csv("P:\\GROUP\\NPIpilot\\data\\randomarm_198.csv",header=TRUE)

####################################### Part 1: For randomized index subjects (n=198) #####################################################

q1 <- q1data[,c(1,4:6,8:26)]
arm <- random[,c(1:3)]
arm$hhID <- tolower(arm$hhID)
arm$scrID <- tolower(arm$scrID)

tab1 <- merge(arm,q1,by="scrID",all.x=TRUE)

tab1_c <- tab1[tab1$arm==1,]
tab1_m <- tab1[tab1$arm==2,]
tab1_h <- tab1[tab1$arm==3,]

#---------------------------------------------start: copy scripts for Part 2--------------------------------------------------------------#
####### Age group #######

dim(tab1_c[tab1_c$age<=15,])[1]
dim(tab1_c[tab1_c$age<=30&tab1_c$age>15,])[1]
dim(tab1_c[tab1_c$age<=50&tab1_c$age>30,])[1]
dim(tab1_c[tab1_c$age>50,])[1]
round(dim(tab1_c[tab1_c$age<=15,])[1]/dim(tab1_c)[1],2)
round(dim(tab1_c[tab1_c$age<=30&tab1_c$age>15,])[1]/dim(tab1_c)[1],2)
round(dim(tab1_c[tab1_c$age<=50&tab1_c$age>30,])[1]/dim(tab1_c)[1],2)
round(dim(tab1_c[tab1_c$age>50,])[1]/dim(tab1_c)[1],2)

dim(tab1_m[tab1_m$age<=15,])[1]
dim(tab1_m[tab1_m$age<=30&tab1_m$age>15,])[1]
dim(tab1_m[tab1_m$age<=50&tab1_m$age>30,])[1]
dim(tab1_m[tab1_m$age>50,])[1]
round(dim(tab1_m[tab1_m$age<=15,])[1]/dim(tab1_m)[1],2)
round(dim(tab1_m[tab1_m$age<=30&tab1_m$age>15,])[1]/dim(tab1_m)[1],2)
round(dim(tab1_m[tab1_m$age<=50&tab1_m$age>30,])[1]/dim(tab1_m)[1],2)
round(dim(tab1_m[tab1_m$age>50,])[1]/dim(tab1_m)[1],2)

dim(tab1_h[tab1_h$age<=15,])[1]
dim(tab1_h[tab1_h$age<=30&tab1_h$age>15,])[1]
dim(tab1_h[tab1_h$age<=50&tab1_h$age>30,])[1]
dim(tab1_h[tab1_h$age>50,])[1]
round(dim(tab1_h[tab1_h$age<=15,])[1]/dim(tab1_h)[1],2)
round(dim(tab1_h[tab1_h$age<=30&tab1_h$age>15,])[1]/dim(tab1_h)[1],2)
round(dim(tab1_h[tab1_h$age<=50&tab1_h$age>30,])[1]/dim(tab1_h)[1],2)
round(dim(tab1_h[tab1_h$age>50,])[1]/dim(tab1_h)[1],2)

####### Sex #######
dim(tab1_c[tab1_c$male==1,])[1]
round(dim(tab1_c[tab1_c$male==1,])[1]/dim(tab1_c)[1],2)
dim(tab1_m[tab1_m$male==1,])[1]
round(dim(tab1_m[tab1_m$male==1,])[1]/dim(tab1_m)[1],2)
dim(tab1_h[tab1_h$male==1,])[1]
round(dim(tab1_h[tab1_h$male==1,])[1]/dim(tab1_h)[1],2)

####### Symptoms #######
colSums(tab1_c[,c(6:24)],na.rm=TRUE)
round(colSums(tab1_c[,c(6:24)],na.rm=TRUE)/dim(tab1_c)[1],2)
colSums(tab1_m[,c(6:24)],na.rm=TRUE)
round(colSums(tab1_m[,c(6:24)],na.rm=TRUE)/dim(tab1_m)[1],2)
colSums(tab1_h[,c(6:24)],na.rm=TRUE)
round(colSums(tab1_h[,c(6:24)],na.rm=TRUE)/dim(tab1_h)[1],2)

####### Delay from symptom onset to randomization #######

dim(tab1_c[tab1_c$onsettime<=2,])[1]
dim(tab1_c[tab1_c$onsettime>2&tab1_c$onsettime<=4,])[1]
dim(tab1_c[tab1_c$onsettime>4&tab1_c$onsettime!=9,])[1]
round(dim(tab1_c[tab1_c$onsettime<=2,])[1]/dim(tab1_c)[1],2)
round(dim(tab1_c[tab1_c$onsettime>2&tab1_c$onsettime<=4,])[1]/dim(tab1_c)[1],2)
round(dim(tab1_c[tab1_c$onsettime>4&tab1_c$onsettime!=9,])[1]/dim(tab1_c)[1],2)

dim(tab1_m[tab1_m$onsettime<=2,])[1]
dim(tab1_m[tab1_m$onsettime>2&tab1_m$onsettime<=4,])[1]
dim(tab1_m[tab1_m$onsettime>4&tab1_m$onsettime!=9,])[1]
round(dim(tab1_m[tab1_m$onsettime<=2,])[1]/dim(tab1_m)[1],2)
round(dim(tab1_m[tab1_m$onsettime>2&tab1_m$onsettime<=4,])[1]/dim(tab1_m)[1],2)
round(dim(tab1_m[tab1_m$onsettime>4&tab1_m$onsettime!=9,])[1]/dim(tab1_m)[1],2)

dim(tab1_h[tab1_h$onsettime<=2,])[1]
dim(tab1_h[tab1_h$onsettime>2&tab1_h$onsettime<=4,])[1]
dim(tab1_h[tab1_h$onsettime>4&tab1_h$onsettime!=9,])[1]
round(dim(tab1_h[tab1_h$onsettime<=2,])[1]/dim(tab1_h)[1],2)
round(dim(tab1_h[tab1_h$onsettime>2&tab1_h$onsettime<=4,])[1]/dim(tab1_h)[1],2)
round(dim(tab1_h[tab1_h$onsettime>4&tab1_h$onsettime!=9,])[1]/dim(tab1_h)[1],2)
#---------------------------------------------End: copy scripts for Part 2--------------------------------------------------------------#

####################################### Part 2: For followed up index subjects (n=128) ###################################################

hchar <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)

hchar <- hchar[,c(1,3)] # hhID, intervention
tab1.2 <- merge(tab1,hchar,by="hhID") # Only include questionnaire records for 128 randomized index subjects

tab1_c <- tab1.2[tab1.2$intervention==1,]
tab1_m <- tab1.2[tab1.2$intervention==2,]
tab1_h <- tab1.2[tab1.2$intervention==3,]
dim(tab1_c)[1]     # n 
dim(tab1_m)[1]     # n 
dim(tab1_h)[1]     # n 

###### Following, copy of scripts above for randomized index (age group, sex, symptoms, delay) #########

#----#----#----#----#

####################################### Part 3: For contacts in 128 followed-up households ################################################

demog <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_demographic_m.csv",header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_20_baseflu_m.csv", header=TRUE)

demog <- demog[,c(1,3:5)] # hhID, member, age, sex
demog$vaccine <- baseflu$vaccine07
tab1.2 <- merge(demog,hchar,by="hhID")
contact <- tab1.2[tab1.2$member>0,]

ct_c <- contact[contact$intervention==1,]
ct_m <- contact[contact$intervention==2,]
ct_h <- contact[contact$intervention==3,]

dim(ct_c)[1]   # n
dim(ct_m)[1]   # n
dim(ct_h)[1]   # n

#### Age group ####

length(na.exclude(ct_c$age[ct_c$age<=15]))
length(na.exclude(ct_c$age[ct_c$age<=30&ct_c$age>15]))
length(na.exclude(ct_c$age[ct_c$age<=50&ct_c$age>30]))
length(na.exclude(ct_c$age[ct_c$age>50]))
round(length(na.exclude(ct_c$age[ct_c$age<=15]))/length(ct_c$age),2)
round(length(na.exclude(ct_c$age[ct_c$age<=30&ct_c$age>15]))/length(ct_c$age),2)
round(length(na.exclude(ct_c$age[ct_c$age<=50&ct_c$age>30]))/length(ct_c$age),2)
round(length(na.exclude(ct_c$age[ct_c$age>50]))/length(ct_c$age),2)

length(na.exclude(ct_m$age[ct_m$age<=15]))
length(na.exclude(ct_m$age[ct_m$age<=30&ct_m$age>15]))
length(na.exclude(ct_m$age[ct_m$age<=50&ct_m$age>30]))
length(na.exclude(ct_m$age[ct_m$age>50]))
round(length(na.exclude(ct_m$age[ct_m$age<=15]))/length(ct_m$age),2)
round(length(na.exclude(ct_m$age[ct_m$age<=30&ct_m$age>15]))/length(ct_m$age),2)
round(length(na.exclude(ct_m$age[ct_m$age<=50&ct_m$age>30]))/length(ct_m$age),2)
round(length(na.exclude(ct_m$age[ct_m$age>50]))/length(ct_m$age),2)

length(na.exclude(ct_h$age[ct_h$age<=15]))
length(na.exclude(ct_h$age[ct_h$age<=30&ct_h$age>15]))
length(na.exclude(ct_h$age[ct_h$age<=50&ct_h$age>30]))
length(na.exclude(ct_h$age[ct_h$age>50]))
round(length(na.exclude(ct_h$age[ct_h$age<=15]))/length(ct_h$age),2)
round(length(na.exclude(ct_h$age[ct_h$age<=30&ct_h$age>15]))/length(ct_h$age),2)
round(length(na.exclude(ct_h$age[ct_h$age<=50&ct_h$age>30]))/length(ct_h$age),2)
round(length(na.exclude(ct_h$age[ct_h$age>50]))/length(ct_h$age),2)

#### Sex ####

sum(ct_c$male, na.rm = TRUE)
round(sum(ct_c$male, na.rm = TRUE)/length(ct_c$male),2)
sum(ct_m$male, na.rm = TRUE)
round(sum(ct_m$male, na.rm = TRUE)/length(ct_m$male),2)
sum(ct_h$male, na.rm = TRUE)
round(sum(ct_h$male, na.rm = TRUE)/length(ct_h$male),2)

#### Vaccine history ####

sum(ct_c$vaccine, na.rm = TRUE)
round(sum(ct_c$vaccine, na.rm = TRUE)/length(ct_c$vaccine),2)
sum(ct_m$vaccine, na.rm = TRUE)
round(sum(ct_m$vaccine, na.rm = TRUE)/length(ct_c$vaccine),2)
sum(ct_h$vaccine, na.rm = TRUE)
round(sum(ct_h$vaccine, na.rm = TRUE)/length(ct_c$vaccine),2)
