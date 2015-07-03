# separate this year's data from last year's data 

######################################################
#program begin
require(chron)

#location of input and output file and today's date
###########################################
input_filename <- "2008_05_05e_0308 delia.csv"
input_dir <- "D:/Calvin flu studies/School absenteeism study/Previous years data/edited/"
today <- "2008_05_08"
output_dir <-"D:/Calvin flu studies/School absenteeism study/Previous years data/edited/"
###########################################

dat<-read.csv(paste(input_dir, input_filename, sep=""), )
dat$date <-dates(as.character(dat$date), format ="y/m/d")

thisyr <- dat[dat$date >= dates("2008/2/18", format ="y/m/d"),]
prevyr <- dat[dat$date < dates("2008/2/18", format ="y/m/d"),]

write.csv(thisyr, paste(output_dir,today,"thisyr_delia.csv", sep=""), row.names=FALSE)
write.csv(prevyr, paste(output_dir,today,"prevyr_delia.csv", sep=""), row.names=FALSE)