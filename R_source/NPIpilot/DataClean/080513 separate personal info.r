# delete personal data in the raw data set

#CDC flu study pre-pilot data 
cdcdat<-read.csv("P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_18_hku.csv")

personal <- cdcdat[,c("Record_ID","q1_p1","q1_p3","q1_p5","q1_p6","q1_p7","q1_p8","q1_p9",
"q2_h_no","q2_5m0","q2_5m1","q2_5m2","q2_5m3","q2_5m4","q2_5m5","q2_5m6","q2_5m7","q2_5m8",
"q2_6m0","q2_6m1","q2_6m2","q2_6m3","q2_6m4","q2_6m5","q2_6m6","q2_6m7","q2_6m8")]

cdcdat[,c("q1_p1","q1_p3","q1_p5","q1_p6","q1_p7","q1_p8","q1_p9",
"q2_5m0","q2_5m1","q2_5m2","q2_5m3","q2_5m4","q2_5m5","q2_5m6","q2_5m7","q2_5m8",
"q2_6m0","q2_6m1","q2_6m2","q2_6m3","q2_6m4","q2_6m5","q2_6m6","q2_6m7","q2_6m8")] <- NA

write.csv(personal, "P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_18_hku_personal.csv", row.names=FALSE)
write.csv(cdcdat, "P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_18_hku_data.csv", row.names=FALSE)


#cdc flu study pilot 2007 data
cdcdat<-read.csv("P:/GROUP/NPIpilot/data/raw_data/2007_11_27_cdc.csv")

personal <- cdcdat[,c("Record_ID","q1_p1","q1_p3","q1_p5","q1_p6","q1_p7","q1_p8","q1_p9",
"q2_h_no","q2_5m0","q2_5m1","q2_5m2","q2_5m3","q2_5m4","q2_5m5","q2_5m6","q2_5m7","q2_5m8",
"q2_6m0","q2_6m1","q2_6m2","q2_6m3","q2_6m4","q2_6m5","q2_6m6","q2_6m7","q2_6m8")]

cdcdat[,c("q1_p1","q1_p3","q1_p5","q1_p6","q1_p7","q1_p8","q1_p9",
"q2_5m0","q2_5m1","q2_5m2","q2_5m3","q2_5m4","q2_5m5","q2_5m6","q2_5m7","q2_5m8",
"q2_6m0","q2_6m1","q2_6m2","q2_6m3","q2_6m4","q2_6m5","q2_6m6","q2_6m7","q2_6m8")] <- "info at 2007_11_27_cdc_personal.csv"

write.csv(personal, "P:/GROUP/NPIpilot/data/raw_data/2007_11_27_cdc_personal.csv", row.names=FALSE)
write.csv(cdcdat, "P:/GROUP/NPIpilot/data/raw_data/2007_11_27_cdc_data.csv", row.names=FALSE)

#cdc flu study main study 2008 data
cdcdat<-read.delim("P:/GROUP/NPIstudy/data/raw/080505 dat.dat")

personal <-cdcdat[,c("q1_namee","q1_hkid", "q1_scr_no", "q2_h_no","q2_cname_m0",
"q2_cname_m1","q2_cname_m2","q2_cname_m3","q2_cname_m4",
"q2_cname_m5","q2_cname_m6","q2_cname_m7","q2_cname_m8")]

cdcdat[,c("q1_namee","q1_hkid", "q1_scr_no", "q2_h_no","q2_cname_m0",
"q2_cname_m1","q2_cname_m2","q2_cname_m3","q2_cname_m4",
"q2_cname_m5","q2_cname_m6","q2_cname_m7","q2_cname_m8")] <-NA

write.csv(personal, "P:/GROUP/NPIstudy/data/raw/080505 dat_personal.csv", row.names=FALSE)
write.csv(cdcdat, "P:/GROUP/NPIstudy/data/raw/080505 dat_data.csv", row.names=FALSE)


