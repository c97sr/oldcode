
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
#

  cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""))
  hh <- cdcdat[cdcdat$q2_h_no!="      ",]

# 1. get the clinic data from index patient
clinicdat_entry <- c("q1_date","q1_p2","q1_p4","q1_4_1p","q1_4_2p","q1_4_3p","q1_4_4p","q1_4_5p","q1_4_6p","q1_4_7p",
   "q1_4_8p","q1_4_9p","q1_4_10p","q1_4_11p","q1_4_12p","q1_4_13p","q1_4_14p","q1_4_15p","q1_4_16p","q1_4_17p",
   "q1_4_18p","q1_4_19p","q1_5","q1_s1","q1_v1","q1_av")
clinicdat_name <- c("date","dob","sex","self_fever","headache","lnode","rnose","hvoice","sthroat","cough","phlegm","sbreath",
   "cpain","apain","vomit","diarrhoea","sjoint","pmuscle","dizzy","tired","chills","nsweat","onsettime", "QVres", "bodytemp","antiviral")
clinicdat_output <- getallhousehold_thing(cdcdat, clinicdat_name, clinicdat_entry) 

# get the age
clinicdat_output$dob <- dates(as.character(clinicdat_output$dob), format= "y-m-d")
clinicdat_output$age <-round((dates(as.character(clinicdat_output$date), format= "y-m-d") - clinicdat_output$dob)/365,0)
clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] <- 
                                       clinicdat_output$age[clinicdat_output$age<0 & !is.na(clinicdat_output$age)] + 100

# delete dob column
clinicdat_output <- clinicdat_output[,names(clinicdat_output)!="dob"]

#change the 9 into 0 for those symptoms
for (i in 5:23){
  clinicdat_output[i][clinicdat_output[i] ==9] <- 0
}
#change the 9 into NA for onsettime,QVres
clinicdat_output$onsettime[clinicdat_output$onsettime ==9] <- NA
clinicdat_output$QVres[clinicdat_output$QVres ==9] <- NA
#change the NA into 0 for antiviral
clinicdat_output$antiviral[is.na(clinicdat_output$antiviral)] <- 0
#rearrange order
clinicdat_output <- as.data.frame(clinicdat_output[c(2,1,4,28,5,26,6:25,27)])
clinicdat_output <- clinicdat_output[order(clinicdat_output$scrID),]

#write.csv(clinicdat_output, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\clinicdat_h.csv", sep="", row.names=FALSE)


# 2. get general household characteristics and hand wash data
housechar_name <- c("intervention","familysize", "house_size", "clinic_date","v1_date", "v2_date", "v3_date", "v4_date",
                    "soap_given", "soap_remain", "alcohol_given","alcohol_remain")
housechar_entry <- c("q2_arm", "q2_1", "q2_2","q1_date","q2_date","q3a_date","q3b_date","q4_date",
                     "q4_14_11", "q4_14_12","q4_14_21","q4_14_22")
housechar_output <- getallhousehold_thing(hh, housechar_name, housechar_entry) 

# for substracting the date into days for homevisit
date_substract <- function(date1, date2){
  date1 <- dates(as.character(date1), format= "y-m-d")
  date2 <- dates(as.character(date2), format= "y-m-d")
  result <- date2 - date1 
  return (result)
}

#add visit days since symptom onset

housechar_output$v1_day <- date_substract(housechar_output$clinic_date, housechar_output$v1_date)
housechar_output$v2_day <- date_substract(housechar_output$clinic_date, housechar_output$v2_date)
housechar_output$v3_day <- date_substract(housechar_output$clinic_date, housechar_output$v3_date)
housechar_output$v4_day <- date_substract(housechar_output$clinic_date, housechar_output$v4_date)

housechar_output$soap_usage <- housechar_output$soap_given-housechar_output$soap_remain
housechar_output$alcohol_usage <- housechar_output$alcohol_given-housechar_output$alcohol_remain

#rearrange order
housechar_output_1 <- as.data.frame(housechar_output[c(1,3:5,15:18)])
housechar_output_1 <- housechar_output_1[order(housechar_output_1$hhID),]

write.csv(housechar_output_1, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\hchar_h.csv", sep="", row.names=FALSE)

housechar_output_2 <- as.data.frame(housechar_output[housechar_output$intervention==3,c(1,19,20)])
housechar_output_2 <- housechar_output_2[order(housechar_output_2$hhID),]
write.csv(housechar_output_2, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\adherence_h.csv", sep="", row.names=FALSE)
