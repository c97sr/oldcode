require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

# for file input and output, only change the date for today
input_filename <- "080327 dat.dat"
today <- "2008_3_27"
#
  #for 2008 data
  cdcdat<-read.delim(paste("P:/GROUP/CDC/NPI data/raw/", input_filename, sep=""))
  hh<-cdcdat[!is.na(cdcdat$q2_h_no),]
 
  hh <-  hh[(hh$q2_arm==1 | hh$q2_arm==3 | hh$q2_arm==4) & !is.na(hh$q2_arm),]
  
  #rename q2_hhsize to q2_1 to fit previous 07 program
  names(hh)[37] <- "q2_1"

# get base line characteristics about flu 
  baseflu_name <- c("west_med", "health_supp", "chin_med", "vaccine08", 
"vaccine07", "vaccine06", "ILI_2week","ILI_1year")

baseflu_entry <-c("q2_wmed_m","q2_hsupp_m","q2_cmed_m","q2_vac08_m","q2_vac07_m","q2_vac06_m",
"q2_ili2w_m","q2_ili1y_m")

# store all basline data that related to flu studies into "baseflu_output database 
baseflu_output <-getallmember_thing(hh, baseflu_name, baseflu_entry) 
#!!! baseflu_output <-del_miss(baseflu_output,"member") #delete the NA entries 

#change those 2 entries to 0
baseflu_output$west_med <- change2to0(baseflu_output$west_med)
baseflu_output$health_supp <- change2to0(baseflu_output$health_supp)
baseflu_output$chin_med <- change2to0(baseflu_output$chin_med)
baseflu_output$vaccine08 <- change2to0(baseflu_output$vaccine08)
baseflu_output$vaccine07 <- change2to0(baseflu_output$vaccine07)
baseflu_output$vaccine06 <- change2to0(baseflu_output$vaccine06)

baseflu_output$ILI_2week <- change2to0(baseflu_output$ILI_2week)

#sort by hhID
baseflu_output<-baseflu_output[order(baseflu_output$hhID),]
write.csv(baseflu_output, "D:/r/080327 antiviral vaccine/080327baseflu.csv", row.names=FALSE)

