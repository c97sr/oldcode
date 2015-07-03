#
# Script to read raw data on lab tests, clean and restructure format
#

#
# first sort out the home cultures ...
#
#home.culture07 <- read.csv("P:\\GROUP\\NPIpilot\\data\\raw_data\\2007_10_27_home_culture.csv", header=TRUE,
#  colClasses=rep(c("character", "integer", "character", "integer", "character"), c(1,1,1,1,11)))

#make 08 data become the format of home.culture07
home.culture <-read.delim("P:/GROUP/NPIstudy/data/raw/080505 household qpcr data.dat",header=TRUE,
  colClasses=rep(c("character"), c(14)))

home.culture$visit <- as.numeric(substr(home.culture$Record_ID,9,9))
home.culture$Record_ID <- substr(home.culture$Record_ID,1,7)

home.culture <-home.culture[c(1,15,3,4,2,5:14)]
names(home.culture) <- names(home.culture07)

#use previous 07 program


for(i in c(3,6:14)){
 home.culture[home.culture[,i]=="" & !is.na(home.culture[,i]),i] <- NA
 home.culture[home.culture[,i]=="NA " & !is.na(home.culture[,i]),i] <- NA
}

names(home.culture)[1] <- "hhID"
names(home.culture)[5] <- "scrID"

m0 <- cbind(home.culture[,c(1:5, 6)], member=0)
m1 <- cbind(home.culture[,c(1:5, 7)], member=1)
m2 <- cbind(home.culture[,c(1:5, 8)], member=2)
m3 <- cbind(home.culture[,c(1:5, 9)], member=3)
m4 <- cbind(home.culture[,c(1:5,10)], member=4)
m5 <- cbind(home.culture[,c(1:5,11)], member=5)
m6 <- cbind(home.culture[,c(1:5,12)], member=6)
m7 <- cbind(home.culture[,c(1:5,13)], member=7)
m8 <- cbind(home.culture[,c(1:5,14)], member=8)

names(m0)[6] <- names(m1)[6] <- names(m2)[6] <- names(m3)[6] <- names(m4)[6] <-
names(m5)[6] <- names(m6)[6] <- names(m7)[6] <- names(m8)[6] <- "qPCR"

home.culture.res <- rbind(m0, m1, m2, m3, m4, m5, m6, m7, m8)
home.culture.res <- home.culture.res[!is.na(home.culture.res$qPCR),]

home.culture.res <- home.culture.res[order(home.culture.res$hhID, home.culture.res$member, home.culture.res$visit),c(1,7,2:4,6,5)]

names(home.culture.res)[5] <- "flu type"

#check
home.culture.res[1:10,]
#table(home.culture.res$culture)
#sum(table(home.culture.res$culture))

write.csv(home.culture.res, "P:/GROUP/NPIstudy/data/2008_05_05_home_culture.csv", row.names=FALSE)

#
# the end
#
#
