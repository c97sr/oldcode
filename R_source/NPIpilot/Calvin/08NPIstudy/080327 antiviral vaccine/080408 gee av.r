#NPI study analysis for vaccination and antiviral
###############

#alldat <- read.csv("P:/GROUP/CDC/NPI data/combined/080407_av_vac.csv", header=TRUE)
alldat <- read.csv("P:/GROUP/NPIstudy/data/combined/080417_av_vac.csv", header=TRUE)



#function for showing result
show_res <-function(inf.gee){
results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))

results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2)                 # Odds ratio and 95% CI
}

#function for inputting different paramaters in GEE and show it
diff_gee <- function (dat,score) {
  inf.gee <- gee(score~factor(av), id=factor(hhIDmem), data=dat )
  summary(inf.gee)
  show_res(inf.gee)
 }

#data format conversion
alldat$hhIDmem <- paste(alldat$hhID, alldat$member , sep="_")
#alldat$all_score <- alldat$all_score/7


alldat$av <- as.character(alldat$av)

alldat$av[alldat$av =="none"] <-0
alldat$av[alldat$av =="amantadine"] <-1
alldat$av[alldat$av =="tamiflu"] <-2
alldat$av[alldat$av =="relenza"] <-3
alldat$av[alldat$av =="ribavirin"] <-3
#exclude relenza and ribavirin
alldat <- alldat[alldat$av !=3,]

alldat$av <-as.factor(alldat$av)

# add dataset is adult and is child
isadult <- alldat[alldat$age >=16,]
ischild <- alldat[alldat$age <16,]

#linear model

library(gee) 
 #overall symptom 
 #overall age
 inf.gee <- gee(all_score~factor(av), id=factor(hhIDmem), data=alldat )
 #summary(inf.gee)
 show_res(inf.gee)
 #adult
 inf.gee <- gee(all_score~factor(av), id=factor(hhIDmem),data=alldat,  subset = (age>=16) )
  show_res(inf.gee)
 # children 
 inf.gee <- gee(all_score~factor(av), id=factor(hhIDmem),data=alldat,  subset = (age<16) )
  show_res(inf.gee)

###############################

#for overall symptom duration
diff_gee(alldat,alldat$all_duration)
diff_gee(isadult,isadult$all_duration)
diff_gee(ischild,ischild$all_duration) 

# for all members
#for overall symptom score
diff_gee(alldat,alldat$all_score)
diff_gee(isadult,isadult$all_score)
diff_gee(ischild,ischild$all_score) 

#for upper respiratory symptom score 
diff_gee(alldat,alldat$up_score)
diff_gee(isadult,isadult$up_score)
diff_gee(ischild,ischild$up_score)

#for lower respiratory symptom score 
diff_gee(alldat,alldat$low_score)
diff_gee(isadult,isadult$low_score)
diff_gee(ischild,ischild$low_score)

#for systemic respiratory symptom score 
diff_gee(alldat,alldat$sys_score)
diff_gee(isadult,isadult$sys_score)
diff_gee(ischild,ischild$sys_score)
#################################

#for all members who have symptoms
 symdat <-alldat[alldat$all_score !=0 & !is.na(alldat$all_score),]
 isadult <- symdat[symdat$isadult ==1,]
ischild <- symdat[symdat$isadult ==0,]
diff_gee(symdat,symdat$all_score)
diff_gee(isadult,isadult$all_score)
diff_gee(ischild,ischild$all_score) 

###############################
# for index only
indexdat <- alldat[alldat$member ==0,] #get index patient only
isadult <- indexdat[indexdat$isadult ==1,]
ischild <- indexdat[indexdat$isadult ==0,]

#for overall symptom score
diff_gee(indexdat,indexdat$all_score)
diff_gee(isadult,isadult$all_score)
diff_gee(ischild,ischild$all_score) 

#for upper respiratory symptom score 
diff_gee(indexdat,indexdat$up_score)
diff_gee(isadult,isadult$lup_score)
diff_gee(ischild,ischild$lup_score)

#for upper respiratory symptom score 
diff_gee(indexdat,indexdat$low_score)
diff_gee(isadult,isadult$low_score)
diff_gee(ischild,ischild$low_score)

#for symstemic respiratory symptom score 
diff_gee(indexdat,indexdat$sys_score)
diff_gee(isadult,isadult$sys_score)
diff_gee(ischild,ischild$sys_score)
#################################



#######################################
# 080417analysis by day

#symptom <-read.csv("P:/GROUP/CDC/NPI data/combined/2008_4_2_symptomday0708_d.csv" , header=TRUE)
symptom <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_17_symwfclin0708_d.csv" , header=TRUE)

symptom$bodytemp <-as.numeric(as.character(symptom$bodytemp))
symptom$scrID <-as.character(symptom$scrID)

symptom$fever <- 1*(symptom$bodytemp >=37.8) #create overall symptom score

#cut the last 3 days from 2007 data
symptom <- symptom[symptom$day !=7 & symptom$day !=8 & symptom$day !=9,] 


symptom$bodytemp <-as.numeric(as.character(symptom$bodytemp))
symptom$scrID <-as.character(symptom$scrID)

symptom$fever <- 1*(symptom$bodytemp >=37.8) #create overall symptom score
for ( i in 1:nrow(symptom)){
symptom$score[i] <- sum(symptom$fever[i],symptom$headache[i],symptom$sthroat[i],symptom$cough[i],
                       symptom$pmuscle[i],symptom$rnose[i],symptom$phlegm[i],na.rm=TRUE)

symptom$up_score[i] <- sum(symptom$sthroat[i], symptom$rnose[i], na.rm=TRUE)
symptom$lo_score[i] <- sum(symptom$cough[i], symptom$phlegm[i] ,na.rm=TRUE)
symptom$sy_score[i] <- sum(symptom$fever[i], symptom$headache[i], symptom$pmuscle[i] ,na.rm=TRUE)
}
#


dat <- symptom
#add antiviral
av <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_av0708_m.csv", header=TRUE)
dat <- merge(av, dat, by= c("hhID","member") ,  all.y=TRUE)
dat$av <- as.character(dat$av)
dat$av[is.na(dat$av)] <- "none"

#for demographics
demographic <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_demo0708_m.csv", header=TRUE)
demographic <-demographic[1:5] #get age and sex and other
dat <- merge( demographic, dat, by= c("hhID","member"),  all.y=TRUE)

#add isadult column
dat$isadult[dat$age >=16] <-1
dat$isadult[is.na(dat$isadult)] <-0


dat$hhmemID <- paste(dat$hhID, dat$member, sep="_")

dat$av <- as.character(dat$av)
dat$av[dat$av =="none"] <-0
dat$av[dat$av =="amantadine"] <-1
dat$av[dat$av =="tamiflu"] <-2
dat$av[dat$av =="relenza"] <-3
dat$av[dat$av =="ribavirin"] <-3

dat <- dat[dat$av != 3,] 

dat$av <-as.factor(dat$av)

#normalize the score
dat$score <- dat$score/7
dat$up_score <- dat$up_score/2
dat$lo_score <- dat$lo_score/2
dat$sy_score <- dat$sy_score/3

#check for logistic regression
#overall
library(gee) 
inf.gee <- gee(score~day*factor(av), id=factor(hhID), data=dat)
 show_res(inf.gee)

 inf.gee <- gee(up_score~day*factor(av), id=factor(hhID), data=dat)
 show_res(inf.gee)

 inf.gee <- gee(lo_score~day*factor(av), id=factor(hhID), data=dat)
 show_res(inf.gee)

 inf.gee <- gee(sy_score~day*factor(av), id=factor(hhID), data=dat)
 show_res(inf.gee)

inf.gee <- gee(fever~day*factor(av), id=factor(hhID), data=dat)
 show_res(inf.gee)

#adult
isadult <- dat[dat$isadult==1,]
inf.gee <- gee(score~day*factor(av), id=factor(hhmemID), data=isadult)
 show_res(inf.gee)

 inf.gee <- gee(up_score~day*factor(av), id=factor(hhmemID), data=isadult)
 show_res(inf.gee)

 inf.gee <- gee(lo_score~day*factor(av), id=factor(hhmemID), data=isadult)
 show_res(inf.gee)

 inf.gee <- gee(sy_score~day*factor(av), id=factor(hhmemID), data=isadult)
 show_res(inf.gee)

inf.gee <- gee(fever~day*factor(av), id=factor(hhmemID), data=isadult)
 show_res(inf.gee)


#children
ischild <- dat[dat$isadult==0,]
inf.gee <- gee(score~day*factor(av), id=factor(hhmemID), data=ischild)
 show_res(inf.gee)

 inf.gee <- gee(up_score~day*factor(av), id=factor(hhmemID), data=ischild)
 show_res(inf.gee)

 inf.gee <- gee(lo_score~day*factor(av), id=factor(hhmemID), data=ischild)
 show_res(inf.gee)

 inf.gee <- gee(sy_score~day*factor(av), id=factor(hhmemID), data=ischild)
 show_res(inf.gee)

inf.gee <- gee(fever~day*factor(av), id=factor(hhmemID), data=ischild)
 show_res(inf.gee)