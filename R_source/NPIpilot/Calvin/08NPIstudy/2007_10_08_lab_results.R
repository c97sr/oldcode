#
# Script to read raw data on lab tests, clean and restructure format
#

#
# first sort out the home cultures ...
#

home.culture <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_10_27_home_culture.csv", header=TRUE,
  colClasses=rep(c("character", "integer", "character", "integer", "character"), c(1,1,1,1,11)))

for(i in c(3,6:14)){
 home.culture[home.culture[,i]=="" & !is.na(home.culture[,i]),i] <- NA
 home.culture[home.culture[,i]=="NA " & !is.na(home.culture[,i]),i] <- NA
}

names(home.culture)[1] <- "hhID"
names(home.culture)[5] <- "scrID"

m0 <- cbind(home.culture[,c(1:5, 6)], member=0)
m1 <- cbind(home.culture[,c(1:5, 7)], member=1)
m2 <- cbind(home.culture[,c(1:5, 8)], member=2)
m3 <- cbind(home.culture[,c(1:5, 9)], member=3)
m4 <- cbind(home.culture[,c(1:5,10)], member=4)
m5 <- cbind(home.culture[,c(1:5,11)], member=5)
m6 <- cbind(home.culture[,c(1:5,12)], member=6)
m7 <- cbind(home.culture[,c(1:5,13)], member=7)
m8 <- cbind(home.culture[,c(1:5,14)], member=8)

names(m0)[6] <- names(m1)[6] <- names(m2)[6] <- names(m3)[6] <- names(m4)[6] <-
  names(m5)[6] <- names(m6)[6] <- names(m7)[6] <- names(m8)[6] <- "culture"

home.culture.res <- rbind(m0, m1, m2, m3, m4, m5, m6, m7, m8)
home.culture.res <- home.culture.res[!is.na(home.culture.res$culture),]

home.culture.res <- home.culture.res[order(home.culture.res$hhID, home.culture.res$member, home.culture.res$visit),c(1,7,2:4,6,5)]

home.culture.res[1:10,]
table(home.culture.res$culture)
sum(table(home.culture.res$culture))

#
# then sort out the clinic samples ...
#

clinic.culture <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_11_09_clinic_culture.csv", header=TRUE,
  colClasses=rep(c("character"), c(29)))

clinic.culture[1:10,]

max.n <- dim(clinic.culture)[1]
n.clinic <- dim(clinic.culture)[2] - 1

tmp <- rep(NA, max.n*n.clinic)
clinic.res <- data.frame(scrID=tmp, culture=tmp)

for(i in 1:n.clinic){
  for(j in 1:max.n){
    clinic.res$scrID[(i-1)*max.n + j] <- tolower(paste(names(clinic.culture)[i+1], clinic.culture[j,1], sep=""))
    clinic.res$culture[(i-1)*max.n + j] <- clinic.culture[j,i+1]
  }
}

clinic.res[1:20,]

dim(clinic.res)
clinic.res <- clinic.res[!(clinic.res$culture=="" | clinic.res$culture=="P"),]
dim(clinic.res)

clinic.res$culture[clinic.res$culture=="Pa2"] <- "0"     # what is Pa2 -- parainfluenza?

table(clinic.res$culture)
sum(table(clinic.res$culture))

#
# then check for double counting and remove duplicates from the _clinic_ file ...
#

home.clinic <- home.culture.res[home.culture.res$visit==0,c(7,6)]
names(home.clinic)[2] <- "culture.home"

table(merge(clinic.res, home.clinic)[,-1])

new.clinic <- merge(clinic.res, home.clinic, all.x=TRUE)
new.clinic <- new.clinic[is.na(new.clinic$culture.home),1:2]

#
# read in the first qPCR results and merge with the clinic file
#

home.qPCR <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_06_11_qPCR.csv", header=TRUE,
  colClasses=rep(c("character", "integer", "character", "integer", "character", "numeric"), c(1,1,1,1,1,10)))

names(home.qPCR)[1] <- "hhID"
names(home.qPCR)[5] <- "qscrID"
home.qPCR$qscrID <- tolower(home.qPCR$qscrID)

v0 <- cbind(home.qPCR[,c(1:5, 6,11)], visit=0)
v1 <- cbind(home.qPCR[,c(1:5, 7,12)], visit=1)
v2 <- cbind(home.qPCR[,c(1:5, 8,13)], visit=2)
v3 <- cbind(home.qPCR[,c(1:5, 9,14)], visit=3)
v4 <- cbind(home.qPCR[,c(1:5,10,15)], visit=4)

names(v0)[6] <- names(v1)[6] <- names(v2)[6] <- names(v3)[6] <- names(v4)[6] <- "AqPCR"
names(v0)[7] <- names(v1)[7] <- names(v2)[7] <- names(v3)[7] <- names(v4)[7] <- "BqPCR"

home.qPCR.res <- rbind(v0, v1, v2, v3, v4)
home.qPCR.res <- home.qPCR.res[!is.na(home.qPCR.res$AqPCR) | !is.na(home.qPCR.res$BqPCR),c(1:2,8,3:4,6:7,5)]

#names(home.qPCR.res)[4:5] <- c("qdate", "qday")
home.res1 <- merge(home.culture.res, home.qPCR.res[,-c(4:5,8)], all=TRUE)

#
# read in the tamiflu qPCR results and merge with the clinic file
#

tamiflu.qPCR <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_10_18_PCR_treated_tamiflu.csv", header=TRUE,
  colClasses=rep(c("character", "integer", "character", "numeric"), c(1,2,3,1)))

tamiflu.qPCR$AqPCR2 <- 0
for(i in 1:nrow(tamiflu.qPCR)) if(tamiflu.qPCR$PCR[i]=="A Pos" &
  !is.na(tamiflu.qPCR$Viral.load[i])) tamiflu.qPCR$AqPCR2[i] <- tamiflu.qPCR$Viral.load[i]

dim(home.res1)
home.res2 <- merge(home.res1, tamiflu.qPCR[,c(1:3,8)], all=TRUE)
dim(home.res2)
table(home.res2$AqPCR, home.res2$AqPCR2) # note there are 3 overlapping qPCR results

#
# read in the symptomatic qPCR results and merge with the clinic file
#

symptom.qPCR <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_10_18_PCR_with_symptoms.csv", header=TRUE,
  colClasses=rep(c("character", "integer", "character", "numeric"), c(2,2,1,1)))

symptom.qPCR$AqPCR3 <- 0
for(i in 1:nrow(symptom.qPCR)) if(symptom.qPCR$PCR[i]=="A Pos" &
  !is.na(symptom.qPCR$Viral.load[i])) symptom.qPCR$AqPCR3[i] <- symptom.qPCR$Viral.load[i]

dim(symptom.qPCR)
symptom.qPCR[symptom.qPCR$PCR!="no sample collected",]
dim(symptom.qPCR)

dim(home.res2)
home.res <- merge(home.res2, symptom.qPCR[,c(1,3,4,7)], all=TRUE)
dim(home.res)
sum(!is.na(home.res$AqPCR) & home.res$AqPCR!=0 & !is.na(home.res$AqPCR2) & home.res$AqPCR2!=0)
sum(!is.na(home.res$AqPCR) & home.res$AqPCR!=0 & !is.na(home.res$AqPCR3) & home.res$AqPCR3!=0)
sum(!is.na(home.res$AqPCR2) & home.res$AqPCR2!=0 & !is.na(home.res$AqPCR3) & home.res$AqPCR3!=0)

#
# then write the finalised datasets ...
#

dim(new.clinic)
dim(home.res)


write.csv(new.clinic, "P:\\GROUP\\NPIpilot\\data\\2007_11_16_clinic_culture.csv", row.names=FALSE)
write.csv(home.res,   "P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv", row.names=FALSE)

#
# the end
#
#
