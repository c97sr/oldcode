
require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
#

cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""), )
cday <-read.csv("P:\\GROUP\\NPIpilot\\data\\cdc_raw\\hchar_h.csv", sep="", )
  hh <- cdcdat[cdcdat$q2_h_no!="      ",]

#1. for getting the daily hand hygiene and mask compliance
day_name <-c("soap_day","handrub_day","washhand_day","mask_day","wheremask_day") 
day_entry <- c("21h","22h","23h","21m","22m")           #as all are in entry q3, only question number is needed
day_output <- getallday_thing(hh, day_name,day_entry)
wheremask_var <- c(01,02,03,04,05,06,07,08,09)
wheremask_name <- c("m_crowded","m_workplace","m_home","m_transport","m_street","m_lift","m_school","m_clinic","m_others")
day_output <-addall_list_data(day_output, day_output$wheremask_day, wheremask_name, wheremask_var)

day_output <- day_output[,names(day_output)!="wheremask_day"]
arm <- data.frame(hhID = hh$q2_h_no)
arm$arm <- hh$q2_arm
output <- merge(arm,day_output,by="hhID")
output <- output[,c(1,3,2,4:18)]
output <- output[order(output$hhID),]
for(i in 1:nrow(output)){
    if(output$mask_day[i]==4&!is.na(output$mask_day[i])) output[i,c(10:18)] <- 0   # m_crowed - m_others
}
day <- rep(cday$v1_day,cday$familysize)
daycount <- rep(day,each=10)
output$daycount <- daycount
output$day <- as.numeric(as.character(output$day))+output$daycount
output_mask_hand <- output[output$arm!=1,c(1,4:9,12)]
write.csv(output_mask_hand, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\adherence_d.csv", sep="", row.names=FALSE)

#2. for getting the symptoms per day
symptomday_name <-c("bodytemp","symptom")
symptomday_entry <-c("1","3to20")
symptomday_output <- getallday_thing(hh, symptomday_name,symptomday_entry)

# break the symptoms list 3to20 into boolean entries
symptom_var <- c(03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20)
symptom_name <- c("headache", "lnode","rnose","hvoice","sthroat","cough",
   "phlegm","sbreath","cpain","apain","vomit",
   "diarrhoea","sjoint","pmuscle","dizzy","tired",
   "chills","nsweat")

symptomday_output_t <-addall_list_data(symptomday_output, 
symptomday_output$symptom , symptom_name, symptom_var)

symptomday_output_t <- symptomday_output_t[order(symptomday_output_t$hhID),]
symptomday_output_t$day <- output$day
symptomday_output_t <- symptomday_output_t[,-c(2,6)]
write.csv(symptomday_output_t, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\symptomday_d.csv", sep="", row.names=FALSE)
