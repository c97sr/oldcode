# convert websms data to our format

######################################################
#program begin
require(chron)

#location of input and output file and today's date
###########################################
input_filename <- "0308 delia.csv"
input_dir <- "D:/Calvin flu studies/School absenteeism study/Previous years data/dmsb/"
today <- "2008_05_05"
output_dir <-"D:/Calvin flu studies/School absenteeism study/Previous years data/edited/"
###########################################
#
alldat<-read.csv(paste(input_dir, input_filename, sep=""), )
alldat <- alldat[alldat$EN_REASON != "Other reason (waived record)",] #del wavied record

record_1 <- function(dat) {
r1 <- data.frame(date=0,school=0,type=0,level=0,ili=0,gi=0,other_ill=0,not_ill=0,undefined=0,total_num=0,schoolcode=0)

r1$date <- dat$NONATTDATE
r1$school <- "dmsb"
r1$type <- substr(dat$CLASSLVL,1,1)
r1$level <- substr(dat$CLASSLVL,2,2)
# Reasons
if (dat$EN_REASON =="Reason unknown"){
    r1$undefined <- 0.5 #unknown is underfine
}
if (dat$EN_REASON == "Sick leave" | dat$EN_REASON == "Sick (Early leave)" ) {
    r1$other_ill <-0.5 #treat sick leave as other illness at the time being
}
#else {r1$not_ill <-1}
#}
if (dat$EN_REASON =="Other reason" |
    dat$EN_REASON =="Parent's application"|
    dat$EN_REASON =="ss"|
    dat$EN_REASON =="Suspension"){
    r1$not_ill <-0.5
} 

#total number of students in each form
if (dat$CLASSLVL =="S1"){
    r1$total_num <- 160
} 
if (dat$CLASSLVL =="S2"){
    r1$total_num <- 160
} 
if (dat$CLASSLVL =="S3"){
    r1$total_num <- 160
} 
if (dat$CLASSLVL =="S4"){
    r1$total_num <- 160
} 
if (dat$CLASSLVL =="S5"){
    r1$total_num <- 160
} 
if (dat$CLASSLVL =="S6"){
    r1$total_num <- 30
} 
if (dat$CLASSLVL =="S7"){
    r1$total_num <- 30
} 
r1$schoolcode <-"dmsb"

return(r1)
}

allinone <- record_1(alldat[1,]) #process, change the format

for (i in 2:dim(alldat)[1]){
  tmp2 <- record_1(alldat[i,])
  allinone <- rbind(allinone,tmp2)
}

 allinone$schlv <- paste(allinone$date, allinone$level, sep="_")
 #combine same date and level to 1 row

comdat <- function(tmp){
  r1 <- tmp[1,]
  for (x in 5:9){
  r1[x] <- sum(tmp[x]) #add up those illness
 }
 return(r1)
}
 tmp <- allinone[allinone$schlv ==names(table(allinone$schlv))[1],]

finalone <- comdat(tmp)

for (i in 2:length(names(table(allinone$schlv))) ){
  tmp <- allinone[allinone$schlv ==names(table(allinone$schlv))[i],]

  finaltmp <- comdat(tmp)
  finalone <- rbind(finalone, finaltmp)
}
finalone <- finalone[-c(12)]
write.csv(finalone, paste(output_dir,today,"e_0308 delia.csv", sep=""), row.names=FALSE)