#
#for prepilot data not including PCR. Script to read raw data on lab tests, clean and restructure format
#

#
# first sort out the home cultures ...
#

home.culture <- read.csv("P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_12_home_culture.csv", header=TRUE,
  colClasses=rep(c("character", "integer", "character", "integer", "character"), c(1,1,1,1,11)))

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
  names(m5)[6] <- names(m6)[6] <- names(m7)[6] <- names(m8)[6] <- "culture"

home.culture.res <- rbind(m0, m1, m2, m3, m4, m5, m6, m7, m8)
home.culture.res <- home.culture.res[!is.na(home.culture.res$culture),]

home.culture.res <- home.culture.res[order(home.culture.res$hhID, home.culture.res$member, home.culture.res$visit),c(1,7,2:4,6,5)]

home.culture.res[1:10,]
table(home.culture.res$culture)
sum(table(home.culture.res$culture))

#
# then sort out the clinic samples ...
#

clinic.culture <- read.csv("P:\\GROUP\\HKUpilot\\data\\raw_data\\2007_12_12_clinic_culture.csv", header=TRUE,
  colClasses=rep(c("character"), c(15)))

clinic.culture[1:10,]

max.n <- dim(clinic.culture)[1]
n.clinic <- dim(clinic.culture)[2] - 1

tmp <- rep(NA, max.n*n.clinic)
clinic.res <- data.frame(scrID=tmp, culture=tmp)

for(i in 1:n.clinic){
  for(j in 1:max.n){
    clinic.res$scrID[(i-1)*max.n + j] <- tolower(paste(names(clinic.culture)[i+1], clinic.culture[j,1], sep=""))
    clinic.res$culture[(i-1)*max.n + j] <- clinic.culture[j,i+1]
  }
}

clinic.res[1:20,]

dim(clinic.res)
clinic.res <- clinic.res[!(clinic.res$culture=="" | clinic.res$culture=="P"),]
dim(clinic.res)

clinic.res$culture[clinic.res$culture=="Pa2"] <- "0"     # what is Pa2 -- parainfluenza?

table(clinic.res$culture)
sum(table(clinic.res$culture))

#
# then check for double counting and remove duplicates from the _clinic_ file ...
#

home.clinic <- home.culture.res[home.culture.res$visit==0,c(7,6)]
names(home.clinic)[2] <- "culture.home"

table(merge(clinic.res, home.clinic)[,-1])

new.clinic <- merge(clinic.res, home.clinic, all.x=TRUE)
new.clinic <- new.clinic[is.na(new.clinic$culture.home),1:2]

# write finish dataset

dim(new.clinic)
dim(home.culture.res)

home.culture.res$AqPCR <-NA
home.culture.res$BqPCR <-NA
home.culture.res$AqPCR2<-NA
home.culture.res$AqPCR3 <-NA

write.csv(new.clinic, "P:\\GROUP\\HKUpilot\\data\\2007_12_17_clinic_culture.csv", row.names=FALSE)
write.csv(home.culture.res,"P:\\GROUP\\HKUpilot\\data\\2007_12_17_home_culture.csv", row.names=FALSE)

#
# the end
#
#
