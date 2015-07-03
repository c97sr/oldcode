###### for school TKPS reformat 
#for combining all TKPS files (need to rename those files to 1.csv, 2.csv...etc


dat <- read.csv("D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613 raw data/tkps/tkps1.csv", header=TRUE)
dire = "D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613 raw data/tkps/"

#combine all data
for (i in 2:99){

 filename = paste ("tkps",i, ".csv", sep="")
 temp <- read.csv( (paste (dire,filename, sep="")) , header = TRUE)
 dat <- rbind(temp, dat)
}

dat <- dat[dat[1] !="Date",]
dat <- dat[-c(11)]

#reformat the date
for (i in 1:dim(dat)[1]){
dat$day[i] <- strsplit(as.character(dat$Date[i]), "/", fixed = TRUE)[[1]][1]
dat$month[i] <- strsplit(as.character(dat$Date[i]), "/", fixed = TRUE)[[1]][2]
dat$year[i] <- strsplit(as.character(dat$Date[i]), "/", fixed = TRUE)[[1]][3]
}

dat$Date <- paste(dat$year,"/",dat$month,"/",dat$day, sep="")
dat <- dat[-c(11:13)]

write.csv( dat, "D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/080613tkps_reformat.csv", row.names=FALSE)

######
#end