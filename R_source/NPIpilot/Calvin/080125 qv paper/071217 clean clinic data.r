require("chron")
#########################################################
#Program begins

qv_cdc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_28_clinicdat_h.csv")
cliniclab_cdc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_clinic_culture.csv")
homelab_cdc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv")
homelab_cdc$scrID <-tolower(homelab_cdc$scrID)

###########################################
#071212 combine pre-pilot data

qv_pre <- read.csv("P:\\GROUP\\HKUpilot\\data\\2007_12_18_clinicdat_h.csv")
cliniclab_pre <- read.csv("P:\\GROUP\\HKUpilot\\data\\2007_12_17_clinic_culture.csv")
homelab_pre <- read.csv("P:\\GROUP\\HKUpilot\\data\\2007_12_17_home_culture.csv")

qv <- rbind(qv_cdc, qv_pre)
cliniclab <- rbind(cliniclab_cdc, cliniclab_pre)
homelab <- rbind(homelab_cdc, homelab_pre)

###########################################

lab <- homelab[!is.na(homelab$scrID) & !homelab$scrID=="",]
lab <- lab[c(7,6)]

homelab <- merge(lab, cliniclab, by="scrID",  all=TRUE)
homelab$culture <- ifelse(is.na(homelab$culture.x), as.character(homelab$culture.y),  as.character(homelab$culture.x))
lab <-homelab[c(1,4)]

#Add clinic code column
cdc <- merge(lab, qv, by="scrID", all.y=TRUE)
cdc$clinic <-substr(cdc$scrID,0,3)

# remove NA of culture, inconconclusive result of QV to -ve, 1008 cases -> 998 cases
cdc <- cdc[!is.na(cdc$culture),]
#cdc <-cdc[!cdc$QVres == 4,]
cdc$QVres[is.na(cdc$QVres)] <-0
cdc$QVres[cdc$QVres == 4] <-3
cdc$QVres[cdc$QVres ==9] <-3

#
# handle data input error for 3 wrongly inputed result, flu B is written in QVres but actually culture result is flu A
#
cdc$QVres[cdc$scrID =="hks139" | cdc$scrID =="ywm111" | cdc$scrID =="ywm118"] <-1
#
#

#080121 add qPCR result
dat_lab <- read.csv("P:/GROUP/NPIpilot/data/2008_01_clinic_lab.csv", header=TRUE)
dat_lab <- dat_lab[c(1,3)]
cdc <- merge(cdc,dat_lab,by="scrID", all.x=TRUE) 

cdc$qPCR[cdc$scrID == "dvc010"] <- 0 # pre-pilot pcr result

#080131 add log qPCR result
#cdc$logqPCR <- log(cdc$qPCR,10)
#cdc$logqPCR[cdc$logqPCR == -Inf] <- 0


write.csv(cdc, "P:\\GROUP\\NPIpilot\\data\\qv paper data\\cleandata\\qvdat.csv", row.names=FALSE)
#write.csv(cdc, "d:/r/080125 qv paper/qvdat.csv", row.names=FALSE)
