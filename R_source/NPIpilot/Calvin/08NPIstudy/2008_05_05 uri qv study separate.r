require(chron)
source ("D:\\Work\\SVNrepository\\NPIpilot\\DataClean\\071008 reformat function.r")

#location of input and output file and today's date
###########################################
input_filename <- "080516 dat all.dat"
input_dir <- "P:/GROUP/NPIstudy/data/raw/"
today <- "080519"
output_dir <-"P:/GROUP/NPIstudy/data/"
###########################################

#separate NPI study and URI study
cdcdat<-read.delim(paste(input_dir, input_filename, sep=""))

########################
# for extracting personal info and delete it in the original dataset
personal <- cdcdat[c("q1_scr_no","q2_h_no","q1_namee", "q1_hkid",
"q2_cname_m0","q2_cname_m1","q2_cname_m2","q2_cname_m3", "q2_cname_m4",
"q2_cname_m5","q2_cname_m6","q2_cname_m7","q2_cname_m8")]
personal <- personal[order(personal$q1_scr_no),]
write.csv(personal, paste(input_dir, today," personal_info.csv", sep=""), row.names=FALSE) # personal info
########################
#delete personal info
cdcdat$q1_namee <-NA
cdcdat$q1_hkid <- NA
cdcdat$q2_cname_m0 <-NA
cdcdat$q2_cname_m1 <-NA
cdcdat$q2_cname_m2 <-NA
cdcdat$q2_cname_m3 <-NA
cdcdat$q2_cname_m4 <-NA
cdcdat$q2_cname_m5 <-NA
cdcdat$q2_cname_m6 <-NA
cdcdat$q2_cname_m7 <-NA
cdcdat$q2_cname_m8 <-NA

cdcdat$q1_scr_no <- tolower(cdcdat$q1_scr_no)
cdcuri <- cdcdat[
 #cdcdat$q1_recmethod ==1 & 
(substr(cdcdat$q1_scr_no,1,3) =="olp"|
 substr(cdcdat$q1_scr_no,1,3) =="tmh"|
 substr(cdcdat$q1_scr_no,1,3) =="ths"|
 substr(cdcdat$q1_scr_no,1,3) =="ucw"|
 substr(cdcdat$q1_scr_no,1,3) =="uyl"
)
 ,] # NPI quickvue study, NA treat as NPI study
cdcqv <- cdcdat[
!(substr(cdcdat$q1_scr_no,1,3) =="olp"|
 substr(cdcdat$q1_scr_no,1,3) =="tmh"|
 substr(cdcdat$q1_scr_no,1,3) =="ths"|
 substr(cdcdat$q1_scr_no,1,3) =="ucw"|
 substr(cdcdat$q1_scr_no,1,3) =="uyl"
),]
#cdcdat$q1_recmethod ==0 &!is.na(cdcdat$q1_recmethod),] #URI study

write.csv(cdcqv, paste(input_dir, today," dat.csv", sep=""), row.names=FALSE) # NPI quickvue study
write.csv(cdcuri, paste(input_dir, today," uri_dat.csv", sep=""), row.names=FALSE) # URI study


