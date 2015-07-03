#
# Extract information for household characteristics
#
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")


  cdcdat<-read.csv("P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_18_hku.csv")
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
clinicdat_output$dob <- dates(as.character(clinicdat_output$dob), format= "d/m/y")
#clinicdat_output$age <-round((dates("07-07-10", format= "y-m-d") - clinicdat_output$dob)/365,0)
clinicdat_output$age <-round((dates(as.character(clinicdat_output$date), format= "d/m/y") - clinicdat_output$dob)/365,0)
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

#date
 
year <- paste("20", substr(clinicdat_output$date,7,8), sep="")
month <- substr(clinicdat_output$date,4,5)
day <- substr(clinicdat_output$date,1,2)
clinicdat_output$date <- paste(year,month,day,sep="-")

write.csv(clinicdat_output, "P:\\GROUP\\HKUpilot\\data\\2007_12_18_clinicdat_h.csv", row.names=FALSE)