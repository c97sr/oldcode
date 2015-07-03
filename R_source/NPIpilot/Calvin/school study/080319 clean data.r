# r program to combine daily school absenteeism data
# need to rename daily data to filename 1.csv, 2.csv...etc

require(chron)

dat <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613 renamed data/1.csv", header = TRUE)
dire = "D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613 renamed data/"

#combine all data 
no_of_file <- 150
for (i in 2:no_of_file){

 filename = paste (i, ".csv", sep="")
 temp <- read.csv( (paste (dire,filename, sep="")) , header = TRUE)
 dat <- rbind(temp, dat)
}

#rename the entries
names(dat) <-c("date", "school", "type", "level", "ili", "gi", "other_ill","not_ill","undefined","total_num")

dat[is.na(dat)] <-0

# treat half day absent as 1 day absent
#for (i in 5:10){
#  dat[i] <- ceiling(dat[i])
#}

#add day, month, year entries 
#for (i in 1:dim(dat)[1]){
#dat$day[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][3]
#dat$month[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][2]
#dat$year[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][1]
#}

#add total absentees entry 
#dat$total_abs <- dat$ili+ dat$gi+dat$other_ill+dat$not_ill+dat$undefined

#add short form school name 
schoolcode <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/080411 school code.csv", header = TRUE)
schoolcode <-as.data.frame(schoolcode)
schoolcode$name <- as.character(schoolcode$name)
schoolcode$short <- as.character(schoolcode$short)
dat$school <-as.character(dat$school)

for (i in 1:dim(schoolcode)[1]){
  dat$schoolcode[dat$school == schoolcode[i,2]  ] <-schoolcode[i,1]
}
#

dat$non_eclass <- 0 # data generate from e-class database

#combine old to newly received data
olddat <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/080602 data process/080602_clean_sch_data.csv", header = TRUE)
dat <- rbind(olddat, dat)

dat <- dat[order(dat$schoolcode),]

#check duplicate
dat$datesch <- paste(dat$date, dat$schoolcode, dat$type, dat$level, sep="_")

dup <- dat[duplicated(dat$datesch),] #extra duplicated rows
uni <- dat[!duplicated(dat$datesch),] #unique rows
uni <- uni[-13]

write.csv(uni, "D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613_clean_sch_data.csv", row.names = FALSE)

