#
# Extract information for household characteristics
#
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

#location of input and output file and today's date
###########################################
input_filename <- "080519 dat.csv"
input_dir <- "P:/GROUP/NPIstudy/data/raw/"
today <- "2008_05_05"
output_dir <-"P:/GROUP/NPIstudy/data/"
###########################################


  #for 2008 data
  cdcdat<-read.csv(paste(input_dir, input_filename, sep=""))
  hh<-cdcdat[!is.na(cdcdat$q2_h_no),]
 
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==3 | hh$q2_arm==4) & !is.na(hh$q2_arm),]

  control <- hh[hh$q2_arm==1,]  # for control group only
  hand <- hh[hh$q2_arm==3,] # for mask group only
  maskhand <- hh[hh$q2_arm==4,] # for hand hygiene group only


# 1. get the clinic data from index patient
  clinicdat_entry <- c("q1_date","q1_dob","q1_sex",
   "q1_fever","q1_hdache","q1_sthroat","q1_cough","q1_aches","q1_nose",
   "q1_phlegm","q1_symonset","q1_QV","q1_tempc","q1_antivir","q1_antibi")

  clinicdat_name <- c("date","dob","male","self_fever","headache",
   "sthroat","cough","aches","rnose",
   "phlegm","onsettime", "QVres", "bodytemp","antiviral","antibiotics")

clinicdat_output <- getallhousehold_thing(cdcdat, clinicdat_name, clinicdat_entry) 
#!!!clinicdat_output <-del_miss(clinicdat_output, "household") #delete the NA entries

# get the age
clinicdat_output$dob <- dates(as.character(clinicdat_output$dob), format= "d/m/y")
#clinicdat_output$age <-round((dates("07-07-10", format= "y-m-d") - clinicdat_output$dob)/365,0)
clinicdat_output$age <-round((dates(as.character(clinicdat_output$date), format= "d/m/y") - clinicdat_output$dob)/365,0)
clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] <- 
clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] + 100

# delete dob column
num <- get_entry_num(clinicdat_output,"dob" )
clinicdat_output <- clinicdat_output[-c(num)]

# add the measured fever column (fever>37.8)
clinicdat_output$measure_fever[clinicdat_output$bodytemp < 37.8] <-0
clinicdat_output$measure_fever[clinicdat_output$bodytemp >= 37.8] <-1


#rearrange order
clinicdat_output <- as.data.frame(clinicdat_output[c(1,2,3,4,17,18,5:16)])
clinicdat_output$hhID <- as.numeric(as.character(clinicdat_output$hhID))
clinicdat_output<-clinicdat_output[order(clinicdat_output$hhID,clinicdat_output$scrID),]

write.csv(clinicdat_output, paste(output_dir, today,"_clinicdat_h.csv", sep=""), row.names=FALSE)


# 2. get generaal household characteristics
housechar_name <- c("intervention" ,"familysize", "house_size",
"onsettime", "clinic_date","v1_date", "v2_date", "v3_date")

housechar_entry <- c("q2_arm", "q2_hhsize", "q2_hhsqf", "q1_symonset",
"q1_date","q2_date","q3a_date","q4_date")

# store all demographic data into data frame "demographic_output
housechar_output <- getallhousehold_thing(hh, housechar_name, housechar_entry) 
#!!!housechar_output <-del_miss(housechar_output, "household") #delete the NA entries


# for substracting the date into days for homevisit
date_substract <- function(date1, date2){
  date1 <- dates(as.character(date1), format= "d/m/y")
  date2 <- dates(as.character(date2), format= "d/m/y")
  result<- date2 - date1
 
  return (result)
}

#add visit days since symptom onset

housechar_output$v1_day <- date_substract(housechar_output$clinic_date, housechar_output$v1_date)
housechar_output$v2_day <- date_substract(housechar_output$clinic_date, housechar_output$v2_date)
housechar_output$v3_day <- date_substract(housechar_output$clinic_date, housechar_output$v3_date)


housechar_output$onsetday[housechar_output$onsettime ==1 & !is.na(housechar_output$onsettime)]<- 0 
housechar_output$onsetday[housechar_output$onsettime ==2 & !is.na(housechar_output$onsettime)]<- 1
housechar_output$onsetday[housechar_output$onsettime ==3 & !is.na(housechar_output$onsettime)]<- 1
housechar_output$onsetday[housechar_output$onsettime ==4 & !is.na(housechar_output$onsettime)]<- 2
housechar_output$onsetday[housechar_output$onsettime ==5 & !is.na(housechar_output$onsettime)]<- 3
housechar_output$onsetday[housechar_output$onsettime ==9| 
is.na(housechar_output$onsettime)]<- NA

housechar_output$clinic_day <- housechar_output$onsetday
housechar_output$v1_day <- housechar_output$v1_day + housechar_output$onsetday
housechar_output$v2_day <- housechar_output$v2_day + housechar_output$onsetday
housechar_output$v3_day <- housechar_output$v3_day + housechar_output$onsetday


#change onsettime 1 to 12, 2 to 24, 3 to 36, 4 to 48 and 5 to 60, 9 and NA to NA
housechar_output$onsettime[housechar_output$onsettime ==1 & !is.na(housechar_output$onsettime)]<- 12 
housechar_output$onsettime[housechar_output$onsettime ==2 & !is.na(housechar_output$onsettime)]<- 24
housechar_output$onsettime[housechar_output$onsettime ==3 & !is.na(housechar_output$onsettime)]<- 36
housechar_output$onsettime[housechar_output$onsettime ==4 & !is.na(housechar_output$onsettime)]<- 48
housechar_output$onsettime[housechar_output$onsettime ==5 & !is.na(housechar_output$onsettime)]<- 60
#housechar_output$onsettime[housechar_output$onsettime ==9| 
#is.na(housechar_output$onsettime)]<- NA

# delete onsettime and onsetday column
#num <- get_entry_num(housechar_output,"onsettime" )
#housechar_output <- housechar_output[-c(num)]
num <- get_entry_num(housechar_output,"onsetday" )
housechar_output <- housechar_output[-c(num)]

#rearrange order
housechar_output <- as.data.frame(housechar_output[c(1:10,14,11:13)])

housechar_output$hhID <- as.numeric(as.character(housechar_output$hhID))
housechar_output<-housechar_output[order(housechar_output$hhID),]

write.csv(housechar_output, paste(output_dir, today,"_housechar_h.csv", sep=""), row.names=FALSE)

###########
#old data structure in 2007, modified in 2008

# 3. get the objective hand hygiene compliance
#handcom_name <- c("soap_given", "soap_remain", "alcohol_given","alcohol_remain")
#handcom_entry <- c("q4_14_11", "q4_14_12","q4_14_21","q4_14_22")

#handcom_output <- getallhousehold_thing(hand, handcom_name, handcom_entry) 
#!!!handcom_output <-del_miss(handcom_output, "household") #delete the NA entries

#handcom_output<-handcom_output[order(handcom_output$hhID),]
#write.csv(handcom_output, paste(output_dir", today,"_handcom_h.csv", sep=""), row.names=FALSE)
