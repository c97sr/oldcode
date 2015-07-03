


#NPI study analysis for vaccination 
###############

alldat <- read.csv("P:/GROUP/NPIstudy/data/combined/080414_av_vac.csv", header=TRUE)
dat <- alldat[alldat$member ==0,] #get index patient only

isadult <- dat[dat$isadult ==1,] #compare those child and adult with different antiviral taken
ischild <- dat[dat$isadult ==0,]
table(isadult$av)
table(ischild$av)

#antiviral precribed or not for index
in_noav <-dat[is.na(dat$av) &  dat$member==0,]
in_av <-dat[!is.na(dat$av) & dat$member==0,]

#antiviral perscribed 24hr before and after 24hrs
#dat24noav <- dat[dat$onsettime<24 & is.na(dat$av),] 
dat24av <- dat[dat$onsettime<24 & !is.na(dat$av),] 
#dat48noav <- dat[dat$onsettime>=24 & is.na(dat$av),] 
dat48av <- dat[dat$onsettime>=24 & !is.na(dat$av),] 


###############
#others
#vaccinated and non-vaccinated index
in_novac <-dat[dat$vac1==0 &  dat$member==0,]
in_vac <-dat[dat$vac1==1 & dat$member==0,]


#histogram for symptom score and duration
hist(in_novac$duration)
hist(in_vac$duration)
hist(in_novac$score)
hist(in_vac$score)


hist(in_noav$duration)
hist(in_av$duration)
hist(in_noav$score)
hist(in_av$score)

hist(dat24av$duration)
hist(dat48av$duration)
hist(dat24av$score)
hist(dat48av$score)


#Vaccination, for analysing symptom duration and symptom score started at 1st home visit.
#vaccination index only

mean(na.exclude(dat$duration[dat$vac1==0 & dat$duration != 0 & dat$member==0]))  #mean duration of symptoms of non-vaccinated index
length(na.exclude(dat$duration[dat$vac1==0 & dat$duration != 0 & dat$member==0]))

mean(na.exclude(dat$duration[dat$vac1==1 & dat$duration != 0 & dat$member==0]))  #mean duration of symptoms of vaccinated index
length(na.exclude(dat$duration[dat$vac1==1 & dat$duration != 0 & dat$member==0]))

mean(na.exclude(dat$score[dat$vac1==0 & dat$score != 0 & dat$member==0]))  #mean duration of symptoms of non-vaccinated members
length(na.exclude(dat$score[dat$vac1==0 & dat$score != 0 & dat$member==0]))

mean(na.exclude(dat$score[dat$vac1==1 & dat$score != 0 & dat$member==0])) 
length(na.exclude(dat$score[dat$vac1==1 & dat$score != 0 & dat$member==0]))

#vaccination index only 3 years
mean(na.exclude(dat$duration[dat$vac1==0 & dat$vac2==0 & dat$vac3==0  & dat$member==0 & dat$duration != 0]))  #mean duration of symptoms of non-vaccinated members
length(na.exclude(dat$duration[dat$vac1==0 & dat$vac2==0 & dat$vac3==0  & dat$member==0 & dat$duration != 0]))

mean(na.exclude(dat$duration[!(dat$vac1==0 & dat$vac2==0 & dat$vac3==0)  & dat$member==0 & dat$duration != 0]))  #mean duration of symptoms of vaccinated members
length(na.exclude(dat$duration[!(dat$vac1==0 & dat$vac2==0 & dat$vac3==0)  & dat$member==0 & dat$duration != 0]))

mean(na.exclude(dat$score[dat$vac1==0 & dat$vac2==0 & dat$vac3==0  & dat$member==0 & dat$score != 0]))  #mean symptom score of symptoms of non-vaccinated members
mean(na.exclude(dat$score[!(dat$vac1==0 & dat$vac2==0 & dat$vac3==0)  & dat$member==0 & dat$score != 0]))  #mean symptom score of symptoms of vaccinated members

# vaccination first year all members
mean(na.exclude(dat$duration[dat$vac1==0 & dat$duration != 0]))  #mean duration of symptoms of non-vaccinated members
mean(na.exclude(dat$duration[dat$vac1==1 & dat$duration != 0]))  #mean duration of symptoms of vaccinated members

mean(na.exclude(dat$score[dat$vac1==0 & dat$score != 0]))  #mean duration of symptoms of non-vaccinated members
mean(na.exclude(dat$score[dat$vac1==1 & dat$score != 0]))  #mean duration of symptoms of vaccinated members

#vaccination 3 years all members
mean(na.exclude(dat$duration[dat$vac1==0 & dat$vac2==0 & dat$vac3==0 & dat$duration != 0]))  #mean duration of symptoms of non-vaccinated members
mean(na.exclude(dat$duration[!(dat$vac1==0 & dat$vac2==0 & dat$vac3==0) & dat$duration != 0]))  #mean duration of symptoms of vaccinated members

mean(na.exclude(dat$score[dat$vac1==0 & dat$vac2==0 & dat$vac3==0 & dat$score != 0]))  #mean symptom score of symptoms of non-vaccinated members
mean(na.exclude(dat$score[!(dat$vac1==0 & dat$vac2==0 & dat$vac3==0) & dat$score != 0]))  #mean symptom score of symptoms of vaccinated members


