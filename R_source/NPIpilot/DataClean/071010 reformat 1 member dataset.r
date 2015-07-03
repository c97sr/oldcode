require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
today <- "2007_11_27"
#

  cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""))
  hh<-cdcdat[cdcdat$q2_h_no!="      ",]
 
  
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==2 | hh$q2_arm==3) & !is.na(hh$q2_arm),]
  control <- hh[hh$q2_arm==1,]  # for control group only
  mask <- hh[hh$q2_arm==2,] # for mask group only
  hand <- hh[hh$q2_arm==3,] # for hand hygiene group only


#program begins

#1. For getting the demographic data 
demographic_name <- c( "dob", "male","occupation","education",
   "ever_smoke","current_smoke","when_smoke","dose_smoke", 
   "chor_disease","what_disease",
   "working","indoor_work","outdoor_work","work_place","work_people",
   "schooling", "school_type","school_people")

demographic_entry <- c( "q2_8m", "q2_9m", "q2_10m", "q2_11m","q2_12m","q2_12_1m","q2_12_2m",
"q2_13m","q2_14_1m","q2_14_2m","q2_18bm", "q2_20_1am", "q2_20_1bm", "q2_20_2m", "q2_20_3m", "q2_18am", 
"q2_19_1m", "q2_19_2m")

# store all demographic data into data frame "demographic_output
demographic_output <- getallmember_thing(hh, demographic_name, demographic_entry) 
#!!!demographic_output <-del_miss(demographic_output, "member") #delete the NA entries

#demographic_output <-demographic_output[!is.na(demographic_output$hhID),] 
#demographic_output <-del_missing(demographic_output)

# adding the chronic disease boolean entries
chronic_var <- c(01,02,03,11,12,21,22,23,31,32,33,41,42,43,51,52,53,61,62,71,72,73,81,82,83,84,90)
chronic_name <- c("c1","c2","c3","c11","c12","c21","c22","c23",
               "c31","c32","c33","c41","c42","c43","c51","c52","c53",
           "c61","c62","c71","c72","c73","c81","c82","c83","c84","c90")


demographic_output <-addall_list_data(demographic_output, 
demographic_output$what_disease, chronic_name, chronic_var)

#adding the age entry 
# for inputting age in all entries
#!!!demographic_output <-demographic_output[!is.na(demographic_output$dob),]  # delete those members not participated in the study
demographic_output$dob <- dates(as.character(demographic_output$dob), format= "y-m-d")
#demographic_output$age <-round((dates("07-07-01", format= "y-m-d") - dates(demographic_output$dob))/365,0)
hhsize <- hh$q2_1              # family size
visitdate <- hh$q2_date        # date of first home visit
demographic_output$visitdate <- rep(hh$q2_date,hh$q2_1)
demographic_output$visitdate <- dates(as.character(demographic_output$visitdate), format= "y-m-d")
demographic_output$age <-round((dates(demographic_output$visitdate) - dates(demographic_output$dob))/365,0)
demographic_output$age[demographic_output$age<0 & !is.na(demographic_output$age)] <- 
demographic_output$age[demographic_output$age<0 & !is.na(demographic_output$age)] + 100

#delete dob entry
num <- get_entry_num(demographic_output,"dob" )
demographic_output <- demographic_output[-c(num)]
#delete visitdate entry
num <- get_entry_num(demographic_output,"visitdate" )
demographic_output <- demographic_output[-c(num)]
# delete what_disease entry
#num <- get_entry_num(demographic_output,"what_disease" )
demographic_output <- demographic_output[-c(num)]

# change those 2 entries to 0
demographic_output$ever_smoke <- change2to0(demographic_output$ever_smoke)
demographic_output$current_smoke <- change2to0(demographic_output$current_smoke)
demographic_output$chor_disease <- change2to0(demographic_output$chor_disease)
demographic_output$indoor_work <- change2to0(demographic_output$indoor_work)
#demographic_output$indoor_work <- changeNAto0(demographic_output$indoor_work)
demographic_output$outdoor_work <- change2to0(demographic_output$outdoor_work)
#demographic_output$outdoor_work <- changeNAto0(demographic_output$outdoor_work)
demographic_output$working <- change2to0(demographic_output$working)
demographic_output$schooling <- change2to0(demographic_output$schooling)
#reorder the data
demographic_output <- demographic_output[c(1:3,47,4:10, 12:19,11,20:46)]
#sort by hhID
demographic_output<-demographic_output[order(demographic_output$hhID),]
# change sex to male, m = 1, f = 0
demographic_output$male <-ifelse(demographic_output$male =="m",1, ifelse(demographic_output$male =="f",0, NA) )

write.csv(demographic_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_demographic_m.csv", sep=""), row.names=FALSE)

write.csv(demographic_output,"d:/r/temp01.csv")


# 2. For getting the baseline data that related to flu studies
baseflu_name <- c("west_med", "health_supp", "chin_med", "vaccine07", 
"vaccine06", "vaccine05", "ILI_2week","ILI_1year")

baseflu_entry <-c("q2_15am","q2_15bm","q2_15cm","q2_16_1am","q2_16_1bm","q2_16_1cm",
"q2_16_2m","q2_16_3m")

# store all basline data that related to flu studies into "baseflu_output database 
baseflu_output <-getallmember_thing(hh, baseflu_name, baseflu_entry) 
#!!! baseflu_output <-del_miss(baseflu_output,"member") #delete the NA entries 

#change those 2 entries to 0
baseflu_output$west_med <- change2to0(baseflu_output$west_med)
baseflu_output$health_supp <- change2to0(baseflu_output$health_supp)
baseflu_output$chin_med <- change2to0(baseflu_output$chin_med)
baseflu_output$vaccine07 <- change2to0(baseflu_output$vaccine07)
baseflu_output$vaccine06 <- change2to0(baseflu_output$vaccine06)
baseflu_output$vaccine05 <- change2to0(baseflu_output$vaccine05)
baseflu_output$ILI_2week <- change2to0(baseflu_output$ILI_2week)

#sort by hhID
baseflu_output<-baseflu_output[order(baseflu_output$hhID),]

write.csv(baseflu_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_baseflu_m.csv", sep=""), row.names=FALSE)


#3. get follow up visits data (q3a, q3b, q4 first part)
followup_name <- c("medic_v2", "medic_v3", "medic_v4", "add_med_v2", 
"add_med_v3", "add_med_v4", "index_hour_v2", "index_hour_v3", "index_hour_v4")

followup_entry <-c("q3a_4m","q3b_4m","q4_4m","q3a_5m","q3b_5m","q4_5m",
"q3a_6m","q3b_6m","q4_6m")

# store all baseline data that related to flu studies into "baseflu_output" database 
followup_output <-getallmember_thing(hh, followup_name, followup_entry) 
#!!!followup_output <-del_miss(followup_output,"member") #delete the NA entries 

#change those 2 entries to 0
followup_output$medic_v2 <- change2to0(followup_output$medic_v2)
followup_output$medic_v3 <- change2to0(followup_output$medic_v3)
followup_output$medic_v4 <- change2to0(followup_output$medic_v4)
followup_output$add_med_v2 <- change2to0(followup_output$add_med_v2)
followup_output$add_med_v3 <- change2to0(followup_output$add_med_v3)
followup_output$add_med_v4 <- change2to0(followup_output$add_med_v4)

followup_output<-followup_output[order(followup_output$hhID),] #sort by hhID
write.csv(followup_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_followup_m.csv", sep=""), row.names=FALSE)


#4 get those data after all home visits (q4 2nd part)
aftervisit_name <- c("reduce_smoke", "avoid_2nd_smoke", "bal_diet", "rest","exercise", 
"ventilation", "soap", "handrub", "washhand", "mask", "where_mask")

aftervisit_entry <-c("q4_8_1m","q4_8_2m","q4_8_3m","q4_8_4m","q4_8_5m","q4_8_6m",
"q4_9_m","q4_10_m","q4_11_m","q4_12_m","q4_13_m")

aftervisit_output <-getallmember_thing(hh, aftervisit_name, aftervisit_entry) 
#!!!aftervisit_output <-del_miss(aftervisit_output, "member")

# adding where put on mask boolean entries
wheremask_var <- c(01,02,03,04,05,06,07,08,09)
wheremask_name <- c("m_crowded","m_workplace","m_home","m_transport","m_street",
                    "m_lift","m_school","m_clinic","m_others")

aftervisit_output <-addall_list_data(aftervisit_output, 
aftervisit_output$where_mask, wheremask_name, wheremask_var)

#delete wear mask entry
num <- get_entry_num(aftervisit_output,"where_mask" )
aftervisit_output <- aftervisit_output[-c(num)]

#change those 2 entries to 0
aftervisit_output$reduce_smoke <- change2to0(aftervisit_output$reduce_smoke)
aftervisit_output$avoid_2nd_smoke <- change2to0(aftervisit_output$avoid_2nd_smoke)
aftervisit_output$bal_diet <- change2to0(aftervisit_output$bal_diet)
aftervisit_output$rest <- change2to0(aftervisit_output$rest)
aftervisit_output$exercise <- change2to0(aftervisit_output$exercise)
aftervisit_output$ventilation <- change2to0(aftervisit_output$ventilation)

aftervisit_output<-aftervisit_output[order(aftervisit_output$hhID),]
write.csv(aftervisit_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_aftervisit_m.csv", sep=""), row.names=FALSE)

#5. get those mask compliance data

maskcom_name <- c("when_mask", "where_mask","given_mask","remain_mask")
maskcom_entry <-c("q4_12_m","q4_13_m", "q4_15_11m","q4_15_12m")

maskcom_output <-getallmember_thing(mask, maskcom_name, maskcom_entry) 
#!!!maskcom_output <-del_miss(maskcom_output,"member")


# adding where put on mask boolean entries
wheremask_var <- c(01,02,03,04,05,06,07,08,09)
wheremask_name <- c("m_crowded","m_workplace","m_home","m_transport","m_street",
                    "m_lift","m_school","m_clinic","m_others")

maskcom_output <-addall_list_data(maskcom_output, 
maskcom_output$where_mask, wheremask_name, wheremask_var)

#delete wear mask entry
num <- get_entry_num(maskcom_output,"where_mask" )
maskcom_output <- maskcom_output[-c(num)]

maskcom_output<-maskcom_output[order(maskcom_output$hhID),]
write.csv(maskcom_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_maskcom_m.csv", sep=""), row.names=FALSE)

#6. get those hand hygiene compliance data

handcom_name <- c( "soap", "handrub", "washhand", "smallgel_given", "smallgel_remain" )
handcom_entry <-c("q4_9_m","q4_10_m","q4_11_m", "q4_14_31m","q4_14_32m")

handcom_output <-getallmember_thing(hand, handcom_name, handcom_entry) 
#handcom_output <-del_miss(handcom_output,"member")

handcom_output<-handcom_output[order(handcom_output$hhID),]
write.csv(handcom_output, paste("P:\\GROUP\\NPIpilot\\data\\", today,"_handcom_m.csv", sep=""), row.names=FALSE)

#
#
#
