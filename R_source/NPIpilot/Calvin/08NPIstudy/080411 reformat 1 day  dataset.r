######################################################
#program begin
require(chron)
source ("D:/Work/SVNrepository/NPIpilot/Calvin/08NPIstudy/071008 reformat function.r")

#location of input and output file and today's date
###########################################
input_filename <- "080505 dat.csv"
input_dir <- "P:/GROUP/NPIstudy/data/raw/"
today <- "2008_05_05"
output_dir <-"P:/GROUP/NPIstudy/data/"
###########################################
#

cdcdat<-read.csv(paste(input_dir, input_filename, sep=""), )
hh <- cdcdat[!is.na(cdcdat$q2_h_no),]
  
 
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==3 | hh$q2_arm==4) & !is.na(hh$q2_arm),]
  control <- hh[hh$q2_arm==1,]  # for control group only
  mask <- hh[hh$q2_arm==3,] # for mask group only
  hand <- hh[hh$q2_arm==4,] # for hand hygiene group only


##########################################

#1. for getting the symptoms per day

symptomday_name <-c("bodytemp","headache","sthroat","cough","pmuscle","rnose","phlegm")
symptomday_entry <-c("temp","hdache","sthroat","cough","aches","nose","phlegm")

symptomday_output <- getallday_thing(hh, symptomday_name,symptomday_entry)

symptomday_output$hhID <- as.numeric(as.character(symptomday_output$hhID))
symptomday_output<-symptomday_output[order(symptomday_output$hhID),]
write.csv(symptomday_output, paste(output_dir, today,"_symptomday_d.csv", sep=""), row.names=FALSE)



##################################################
#old data of 2007, database modified in 2008, need to change

#1. for getting the daily hand hygiene compliance
handday_name <-c("soap_day","handrub_day","washhand_day") 
handday_entry <- c("21h","22h","23h")           #as all are in entry q3, only question number is needed
handday_output <- getallday_thing(hand, handday_name,handday_entry)
#!!!handday_output <- del_miss(handday_output, "day")

handday_output<-handday_output[order(handday_output$hhID),]
write.csv(handday_output, paste(output_dir, today,"_handday_d.csv", sep=""), row.names=FALSE)

#2. for getting the daily mask compliance
maskday_name <-c("mask_day","wheremask_day") 
maskday_entry <-c("21m","22m")
maskday_output <- getallday_thing(mask, maskday_name,maskday_entry)
#!!!maskday_output <- del_miss(maskday_output, "day")

# break the where to wear mask list into boolean entries
wheremask_var <- c(01,02,03,04,05,06,07,08,09)
wheremask_name <- c("m_crowded","m_workplace","m_home","m_transport","m_street",
                    "m_lift","m_school","m_clinic","m_others")

maskday_output <-addall_list_data(maskday_output, 
maskday_output$wheremask_day, wheremask_name, wheremask_var)

# delete wheremask_day column
num <- get_entry_num(maskday_output,"wheremask_day" )
maskday_output <- maskday_output[-c(num)]


maskday_output<-maskday_output[order(maskday_output$hhID),]
write.csv(maskday_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_maskday_d.csv", sep=""), row.names=FALSE)

#
