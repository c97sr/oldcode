######################################################
#program begin
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
today <- "2007_11_27"
#

cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""), )
  hh <- cdcdat[cdcdat$q2_h_no!="      ",]
  
 # household <- cdcdat[cdcdat$q2_h_no!="",]
  
  #names(household)[2] <- "ID"
  #hh<-household[2:dim(household)[1]]
  #hh <-household
  #names(hh)[2] <- "ID"

  #homevisit <-cdc[(cdc$household!="" & !is.na(cdc$ID)),]
  #hh <- merge(household, homevisit, by="ID", all=TRUE)
  
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==2 | hh$q2_arm==3) & !is.na(hh$q2_arm),]
  control <- hh[hh$q2_arm==1,]  # for control group only
  mask <- hh[hh$q2_arm==2,] # for mask group only
  hand <- hh[hh$q2_arm==3,] # for hand hygiene group only

#1. for getting the daily hand hygiene compliance
handday_name <-c("soap_day","handrub_day","washhand_day") 
handday_entry <- c("21h","22h","23h")           #as all are in entry q3, only question number is needed
handday_output <- getallday_thing(hand, handday_name,handday_entry)
#!!!handday_output <- del_miss(handday_output, "day")

handday_output<-handday_output[order(handday_output$hhID),]
write.csv(handday_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_handday_d.csv", sep=""), row.names=FALSE)

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

#3. for getting the symptoms per day
#symptomday_name <-c("bodytemp","symptom")
#symptomday_entry <-c("1","3to20")
symptomday_name <-c("bodytemp","symptom")
symptomday_entry <-c("1","3to20")

symptomday_output <- getallday_thing(hh, symptomday_name,symptomday_entry)
#!!!symptomday_output <- del_miss(symptomday_output, "day")

# break the symptoms list 3to20 into boolean entries
symptom_var <- c(03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20)
symptom_name <- c("headache", "lnode","rnose","hvoice","sthroat","cough",
   "phlegm","sbreath","cpain","apain","vomit",
   "diarrhoea","sjoint","pmuscle","dizzy","tired",
   "chill","nsweat")

symptomday_output_t <-addall_list_data(symptomday_output, 
symptomday_output$symptom , symptom_name, symptom_var)

# delete the symptom column
num <- get_entry_num(symptomday_output_t,"symptom" )
symptomday_output_t <- symptomday_output_t[-c(num)]

symptomday_output_t<-symptomday_output_t[order(symptomday_output_t$hhID),]
write.csv(symptomday_output_t, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_symptomday_d.csv", sep=""), row.names=FALSE)
