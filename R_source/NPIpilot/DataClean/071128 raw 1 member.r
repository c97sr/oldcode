require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "2007_11_27_cdc.csv"
#
  cdcdat<-read.csv(paste("P:\\GROUP\\NPIpilot\\data\\raw_data\\", input_filename, sep=""))
  hh<-cdcdat[cdcdat$q2_h_no!="      ",]

#1. For getting the demographic data 
demographic_name <- c("dob", "sex","vaccine07")
demographic_entry <- c("q2_8m", "q2_9m","q2_16_1am")
demographic_output <- getallmember_thing(hh, demographic_name, demographic_entry) 
demographic_output$vaccine07 <- change2to0(baseflu_output$vaccine07)

#adding the age entry 
demographic_output$dob <- dates(as.character(demographic_output$dob), format= "y-m-d")
hhsize <- hh$q2_1              # family size
visitdate <- hh$q2_date        # date of first home visit
demographic_output$visitdate <- rep(hh$q2_date,hh$q2_1)
demographic_output$visitdate <- dates(as.character(demographic_output$visitdate), format= "y-m-d")
demographic_output$age <- round((dates(demographic_output$visitdate) - dates(demographic_output$dob))/365,0)
demographic_output$age[demographic_output$age<0 & !is.na(demographic_output$age)] <- 
demographic_output$age[demographic_output$age<0 & !is.na(demographic_output$age)] + 100

demographic_output <- demographic_output[c(1,3,5,8,6)]
demographic_output <- demographic_output[order(demographic_output$hhID),]

#2. get those mask compliance data

maskcom_name <- c("when_mask", "where_mask","given_mask","remain_mask")
maskcom_entry <-c("q4_12_m","q4_13_m", "q4_15_11m","q4_15_12m")
maskcom_output <-getallmember_thing(hh, maskcom_name, maskcom_entry) 
maskcom_output$mask_usage <- as.numeric(as.character(maskcom_output$given_mask))-as.numeric(as.character(maskcom_output$remain_mask))

# adding where put on mask boolean entries
wheremask_var <- c(01,02,03,04,05,06,07,08,09)
wheremask_name <- c("m_crowded","m_workplace","m_home","m_transport","m_street",
                    "m_lift","m_school","m_clinic","m_others")
maskcom_output <-addall_list_data(maskcom_output,maskcom_output$where_mask, wheremask_name, wheremask_var)

for(i in 1:nrow(maskcom_output)){
    if(maskcom_output$when_mask[i]==4&!is.na(maskcom_output$when_mask[i])) maskcom_output[i,c(9:17)] <- 0   # m_crowed - m_others
}

maskcom_output <- maskcom_output[order(maskcom_output$hhID),c(1:4,8,11)]
output1 <- cbind(demographic_output,maskcom_output[,-c(1:3)])

#3. get those hand hygiene compliance data

handcom_name <- c( "soap", "handrub", "washhand", "smallgel_given", "smallgel_remain" )
handcom_entry <-c("q4_9_m","q4_10_m","q4_11_m", "q4_14_31m","q4_14_32m")
handcom_output <- getallmember_thing(hh, handcom_name, handcom_entry) 
handcom_output$smallgel_usage <- 
                   as.numeric(as.character(handcom_output$smallgel_given))-as.numeric(as.character(handcom_output$smallgel_remain))
handcom_output <- handcom_output[,names(handcom_output)!=c("smallgel_given","smallgel_remain")]
handcom_output <- handcom_output[order(handcom_output$hhID),]
output2 <- cbind(output1,handcom_output[,-c(1:3)])

write.csv(output2, "P:\\GROUP\\NPIpilot\\data\\cdc_raw\\adherence_m.csv", sep="", row.names=FALSE)

