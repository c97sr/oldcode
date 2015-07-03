#
# Extract information for household characteristics
#
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
today <- "2007_11_28"
#

  cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""))
  hh <- cdcdat[cdcdat$q2_h_no!="      ",]
  
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==2 | hh$q2_arm==3) & !is.na(hh$q2_arm),]
  control <- hh[hh$q2_arm==1,]  # for control group only
  mask <- hh[hh$q2_arm==2,] # for mask group only
  hand <- hh[hh$q2_arm==3,] # for hand hygiene group only


# 1. get the clinic data from index patient
  clinicdat_entry <- c("q1_date","q1_p2","q1_p4",
   "q1_4_1p","q1_4_2p","q1_4_3p","q1_4_4p","q1_4_5p","q1_4_6p",
   "q1_4_7p","q1_4_8p","q1_4_9p","q1_4_10p","q1_4_11p","q1_4_12p",
   "q1_4_13p","q1_4_14p","q1_4_15p","q1_4_16p","q1_4_17p",
   "q1_4_18p","q1_4_19p","q1_5","q1_s1","q1_v1","q1_v2","q1_v3",
   "q1_v4","q1_v5","q1_av","q1_ab")

  clinicdat_name <- c("date","dob","male","self_fever","headache",
   "lnode","rnose","hvoice","sthroat","cough",
   "phlegm","sbreath","cpain","apain","vomit",
   "diarrhoea","sjoint","pmuscle","dizzy","tired",
   "chill","nsweat","onsettime", "QVres", "bodytemp","heartrate",
   "resprate","sysBP","diaBP","antiviral","antibiotics")

clinicdat_output <- getallhousehold_thing(cdcdat, clinicdat_name, clinicdat_entry) 
#!!!clinicdat_output <-del_miss(clinicdat_output, "household") #delete the NA entries

# get the age
clinicdat_output$dob <- dates(as.character(clinicdat_output$dob), format= "y-m-d")
#clinicdat_output$age <-round((dates("07-07-10", format= "y-m-d") - clinicdat_output$dob)/365,0)
clinicdat_output$age <-round((dates(as.character(clinicdat_output$date), format= "y-m-d") - clinicdat_output$dob)/365,0)
clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] <- 
clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] + 100

# delete dob column
num <- get_entry_num(clinicdat_output,"dob" )
clinicdat_output <- clinicdat_output[-c(num)]

# add the measured fever column
clinicdat_output$measure_fever[clinicdat_output$bodytemp < 38] <-0
clinicdat_output$measure_fever[clinicdat_output$bodytemp >= 38] <-1

#change the 9 into 0 for those symptoms
for (i in 4:22){
clinicdat_output[i][clinicdat_output[i] ==9] <- 0
}

#rearrange order
clinicdat_output <- as.data.frame(clinicdat_output[c(2,1,3,4,33:34,5:32)])
clinicdat_output<-clinicdat_output[order(clinicdat_output$scrID),]
#change male entry 1: male, 0:female
clinicdat_output$male <-ifelse(clinicdat_output$male =="m",1, ifelse(clinicdat_output$male =="f",0, NA) )

write.csv(clinicdat_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_clinicdat_h.csv", sep=""), row.names=FALSE)


# 2. get generaal household characteristics
housechar_name <- c("intervention", "q1_familysize","q2_familysize", "house_size",
"onsettime", "clinic_date","v1_date", "v2_date", "v3_date", "v4_date", 
"pets","dog","cat","turtle","bird","hamster","rabbit","otherpet","otherpet_no")

housechar_entry <- c("q2_arm", "q1_1", "q2_1", "q2_2", "q1_5",
"q1_date","q2_date","q3a_date","q3b_date","q4_date",
"q2_3", "q2_3_1","q2_3_2","q2_3_3","q2_3_4","q2_3_5","q2_3_6","q2_3_7_type","q2_3_7_no")

# store all demographic data into data frame "demographic_output
housechar_output <- getallhousehold_thing(hh, housechar_name, housechar_entry) 
#!!!housechar_output <-del_miss(housechar_output, "household") #delete the NA entries


# for substracting the date into days for homevisit
date_substract <- function(date1, date2){
  date1 <- dates(as.character(date1), format= "y-m-d")
  date2 <- dates(as.character(date2), format= "y-m-d")
  result<- date2 - date1
 
  return (result)
}

#add visit days since symptom onset

housechar_output$v1_day <- date_substract(housechar_output$clinic_date, housechar_output$v1_date)
housechar_output$v2_day <- date_substract(housechar_output$clinic_date, housechar_output$v2_date)
housechar_output$v3_day <- date_substract(housechar_output$clinic_date, housechar_output$v3_date)
housechar_output$v4_day <- date_substract(housechar_output$clinic_date, housechar_output$v4_date)

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
housechar_output$v4_day <- housechar_output$v4_day + housechar_output$onsetday

#change the 9 into 0 for those pets
for (i in 14:21){
housechar_output[i][is.na(housechar_output[i])] <- 0
}

#change onsettime 1 to 12, 2 to 24, 3 to 36, 4 to 48 and 5 to 60, 9 and NA to NA
housechar_output$onsettime[housechar_output$onsettime ==1 & !is.na(housechar_output$onsettime)]<- 12 
housechar_output$onsettime[housechar_output$onsettime ==2 & !is.na(housechar_output$onsettime)]<- 24
housechar_output$onsettime[housechar_output$onsettime ==3 & !is.na(housechar_output$onsettime)]<- 36
housechar_output$onsettime[housechar_output$onsettime ==4 & !is.na(housechar_output$onsettime)]<- 48
housechar_output$onsettime[housechar_output$onsettime ==5 & !is.na(housechar_output$onsettime)]<- 60
housechar_output$onsettime[housechar_output$onsettime ==9| 
is.na(housechar_output$onsettime)]<- NA

# delete onsettime and onsetday column
#num <- get_entry_num(housechar_output,"onsettime" )
#housechar_output <- housechar_output[-c(num)]
num <- get_entry_num(housechar_output,"onsetday" )
housechar_output <- housechar_output[-c(num)]

#rearrange order
housechar_output <- as.data.frame(housechar_output[c(1:21,26,22:25)])
housechar_output<-housechar_output[order(housechar_output$hhID),]

write.csv(housechar_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_housechar_h.csv", sep=""), row.names=FALSE)


# 3. get the objective hand hygiene compliance
handcom_name <- c("soap_given", "soap_remain", "alcohol_given","alcohol_remain")
handcom_entry <- c("q4_14_11", "q4_14_12","q4_14_21","q4_14_22")

handcom_output <- getallhousehold_thing(hand, handcom_name, handcom_entry) 
#!!!handcom_output <-del_miss(handcom_output, "household") #delete the NA entries

handcom_output<-handcom_output[order(handcom_output$hhID),]
write.csv(handcom_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_handcom_h.csv", sep=""), row.names=FALSE)
