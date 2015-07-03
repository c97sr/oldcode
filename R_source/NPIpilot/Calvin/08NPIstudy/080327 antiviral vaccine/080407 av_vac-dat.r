#Combine this year and last year's symptom diary results
#for antiviral and vaccination study
#old
#symptom <-read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_symptomday0708_d.csv" , header=TRUE)

# 080417 new plots
symptom <-read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_17_symwfclin0708_d.csv" , header=TRUE)
symptom <- symptom[!is.na(symptom$hhID),]

#symptom <- symptom[symptom$member ==0,]

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

#for extracting the symptom duration to variable "dur"
#dur <- data.frame(scrID=0 , member=0, duration=0, score =0, hhID=0, 
#sc_day0=0,sc_day1=0,sc_day2=0,sc_day3=0,sc_day4=0,sc_day5=0,sc_day6=0)

#count <- 0
#cntscore <-0
#j <-1

#added 3103 for storing daily symptom score
# store_score <-NA

#for (i in 1:(nrow(symptom)-1)){
#  store_score <- cbind(store_score, symptom$score[i]) #add 3103
  
 # if ( symptom$score[i] >= 1 & !is.na(symptom$score[i]) ) {
 #   count <- count + 1	# symptom duration (1 or more ili)
 #   cntscore <- cntscore + symptom$score[i] # averge symptom score per day
 #   }
 
 # if (symptom$hhID[i] != symptom$hhID[i+1] | symptom$member[i] != symptom$member[i+1] ){
 #   # print (symptom$scrID[i])
  #   dur[j,1] <-symptom$scrID[i]
    # print (symptom$member[i])
  #   dur[j,2] <-symptom$member[i]
    # print (count)
  #   dur[j,3] <- count
    # print (cntscore)
   #  dur[j,4] <- cntscore/(symptom$day[i] +1) 
   #  dur[j,5] <- as.character(symptom$hhID[i])
    
      # daily average symptom score

    #   for (k in 1:7){
    #   dur[j,k+5] <- store_score[k+1]
     #  }

     #j <- j+1
     #count <- 0
     #cntscore <- 0
     #store_score <- NA

  #}

#}
#

#080409 add fever/lower/upper/systemic infection score

  #function for getting other scores
  get_score <- function(symptom, a_score, con){
   
    dur <- data.frame(hhID=0, scrID=0 , member=0, duration=0,   score =0,
    day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0)
    
    NA_count <-0
    count <- 0
    cntscore <-0
    j <-1

#added 3103 for storing daily symptom score
 store_score <-NA

for (i in 1:(nrow(symptom)-1)){
  store_score <- cbind(store_score, a_score[i]) 
  
  if ( a_score[i] >= 1 & !is.na(a_score[i]) ) {
    count <- count + 1	# symptom duration
    cntscore <- cntscore + a_score[i] # averge symptom score per day
    }
  else if (is.na(a_score[i])){  #NA entries count
    NA_count <-NA_count +1
  }
   
  if (symptom$hhID[i] != symptom$hhID[i+1] | symptom$member[i] != symptom$member[i+1] ){
     dur[j,1] <- as.character(symptom$hhID[i])
     dur[j,2] <-symptom$scrID[i]
     dur[j,3] <-symptom$member[i]
     dur[j,4] <- count
     dur[j,5] <- round(cntscore/(symptom$day[i] +1 -NA_count),2) 
   
    
      # daily average symptom score

       for (k in 1:7){
       dur[j,k+5] <- round(store_score[k+1],2)
       }

     j <- j+1
     NA_count <-0
     count <- 0
     cntscore <- 0
     store_score <- NA

  }

}
  names(dur)[4:12] <- paste( con , names(dur)[4:12], sep="_")
  return (dur)
} 
#
###########
 get_fever <- function(symptom, a_score, con){ #some problems of overall averaged temperature
 #somthing wrong with the score
   
    dur <- data.frame(hhID=0, scrID=0 , member=0, duration=0,   score =0,
    day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0)
    NA_count <-0 #counting the NAs
    count <- 0
    cntscore <-0
    j <-1

#added 3103 for storing daily symptom score
 store_score <-NA

for (i in 1:(nrow(symptom)-1)){
  store_score <- cbind(store_score, a_score[i]) 
   
   if (!is.na(a_score[i])){ 
    cntscore <- cntscore + a_score[i] # averge symptom score per day
  }
   
  if ( a_score[i] >= 37.8 & !is.na(a_score[i]) ) {
    count <- count + 1	# fever duration
    }
  else if (is.na(a_score[i]) ){ #NA in the entry
  NA_count <- NA_count + 1
  }
  
 
  if (symptom$hhID[i] != symptom$hhID[i+1] | symptom$member[i] != symptom$member[i+1] ){
     dur[j,1] <- as.character(symptom$hhID[i])
     dur[j,2] <-symptom$scrID[i]
     dur[j,3] <-symptom$member[i]
     dur[j,4] <- count
     dur[j,5] <- round(cntscore/(symptom$day[i] +1 - NA_count),2) 
   
      # daily average symptom score

       for (k in 1:7){
       dur[j,k+5] <- round(store_score[k+1],2)
       }

     j <- j+1
     count <- 0
     NA_count <-0
     cntscore <- 0
     store_score <- NA

  }

}
  names(dur)[4:12] <- paste( con , names(dur)[4:12], sep="_")
  return (dur)
} 
################

total <-get_score(symptom,symptom$score,"all")
fever <- get_fever(symptom,symptom$bodytemp,"f")
up_sym <- get_score(symptom,symptom$up_score,"up")
low_sym <- get_score(symptom,symptom$lo_score,"low")
sys_sym <- get_score(symptom,symptom$sy_score,"sys")

dur  <- merge(total, fever, by= c("hhID","scrID","member") )
dur <- merge(dur, up_sym, by= c("hhID","scrID","member") )
dur <- merge(dur, low_sym, by= c("hhID","scrID","member") )
dur <- merge(dur, sys_sym, by= c("hhID","scrID","member") )

#for vaccination status
baseflu <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_baseflu0708_m.csv", header=TRUE)
vac <- baseflu[,c(1,3,7:9)]
#compare with and without antiviral with symptom score and duration
#index case and household member case
dat <- merge(dur, vac, by= c("hhID","member") ,  all.x=TRUE)

####################
#for antivirals
av <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_av0708_m.csv", header=TRUE)
dat <- merge(av, dat, by= c("hhID","member") ,  all.y=TRUE)
dat$av <- as.character(dat$av)
dat$av[is.na(dat$av)] <- "none"

#add index antiviral type (for 2nd case analysis)
indexav <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_av0708_m.csv", header=TRUE)
indexav <- indexav[indexav$member==0,]
indexav <- indexav[c(1,4)]
names(indexav) <- c("hhID","indexav")

dat <- merge( dat,indexav,by= c("hhID") ,  all.x=TRUE)
dat$indexav <- as.character(dat$indexav)
dat$indexav[is.na(dat$indexav)] <- "none"
#######################################

#for demographics
demographic <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_demo0708_m.csv", header=TRUE)
demographic <-demographic[1:5] #get age and sex and other
dat <- merge( demographic, dat, by= c("hhID","member"),  all.y=TRUE)

#add isadult column
dat$isadult[dat$age >=16] <-1
dat$isadult[is.na(dat$isadult)] <-0

dat$dummy <- 1 #for getting all entries

# add onsettime
housechar <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_housechar0708_h.csv", header=TRUE)
onset <-housechar[c(1,2,6)]
dat <- merge( onset, dat, by= "hhID" ,  all.y=TRUE) #onsettime, if non-index, onsettime is follow index



#separate 2007 and 2008 data
#dat07 <- dat[grep("h07", dat$hhID) ,]
#dat08 <- dat[grep("h08", dat$hhID) ,]

dat <-dat[-c(5,8,14)] # delete duplicated scrID
#write.csv(dat, "P:/GROUP/NPIstudy/data/combined/080414_av_vac_allmem.csv", row.names=FALSE)

###################################################
# add more information base on index patient

#dat <- dat[dat$member ==0,] #get index patient only

# add clinicdat for index patient
#clinic <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_clinicdat0708_h.csv", header=TRUE)
#clinic <- clinic[clinic$hhID != "h080NA" & clinic$hhID != "      ",]

#for ( i in 1:nrow(clinic)){
#clinic$all_dayclin[i] <- sum(clinic$measure_fever[i],clinic$headache[i],clinic$sthroat[i],clinic$cough[i],
#                       clinic$aches[i],clinic$rnose[i],clinic$phlegm[i],na.rm=TRUE)

#clinic$f_dayclin[i] <- sum(clinic$bodytemp[i],na.rm=TRUE)
#clinic$f_dayclin[clinic$f_dayclin==0] <-NA
#clinic$up_dayclin[i] <- sum(clinic$sthroat[i],clinic$rnose[i],na.rm=TRUE)
#clinic$low_dayclin[i] <- sum(clinic$cough[i],clinic$phlegm[i],na.rm=TRUE)
#clinic$sys_dayclin[i] <- sum(clinic$measure_fever[i],clinic$headache[i],clinic$aches[i] ,na.rm=TRUE)

#} # aches = pmuscle, measure_fever = fever , clinicdat vs symptomday

#clinic <- clinic[c(1:2,19:23)] #add clinic symptom score
#dat <- merge(clinic, dat, by= c("hhID"),  all.y=TRUE)
#dat <- dat[-c(8)]
#dat <-dat[c(1,2,8:18,3,19:27,4,28:36,5,37:45,6,46:54,7,55:66)] #rearrange order


write.csv(dat, "P:/GROUP/NPIstudy/data/combined/080417_av_vac.csv", row.names=FALSE)