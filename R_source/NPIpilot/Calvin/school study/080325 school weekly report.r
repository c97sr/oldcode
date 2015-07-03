# Plot each school and compare all schools to chp

#parameters 
###############################
input_dir <- "D:/Calvin flu studies/School absenteeism study/Prospective data/080613 data process/"
file_name <- "080613_clean_sch_data.csv"
output_dir <- "D:/Calvin flu studies/School absenteeism study/Prospective data/"
today <- "080616"
##############################

require(chron)
dat <- read.csv( paste(input_dir, file_name, sep=""), header=TRUE)


#add day, month, year entries 
for (i in 1:dim(dat)[1]){
dat$day[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][3]
dat$month[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][2]
dat$year[i] <- strsplit(as.character(dat$date[i]), "/", fixed = TRUE)[[1]][1]
}
dat$year <- as.numeric(dat$year)
dat$month <- as.numeric(dat$month)
dat$day <- as.numeric(dat$day)
# add exclude data entry

dat$exclude <- 0
#dat$exclude[dat$non_eclass==1] <-1 # exclude those non e-class data 
dat$exclude[ dat$month ==2 & dat$day <18 ] <-1 #starting from 18/2
#dat$exclude[ dat$month ==3 & dat$day >=24 & dat$day <30 ] <-1 #exclude easter holiday

#dat$exclude[ dat$month ==4 & dat$day >26 ] <-1 #end at 26/4
 

#dat$total_abs <- dat$ili+ dat$gi+dat$other_ill+dat$not_ill+dat$undefined  # raw absenteeism rate
dat$total_abs <-0

#for non_eclass data, they only define illness (grouped in other illness entry), non-illness and undefine,
#therefore we can only exclude non-illness in the analysis 

for (i in 1:dim(dat)[1]){
  if (dat$non_eclass[i] == 1){
      dat$total_abs[i] <- dat$other_ill[i]+dat$undefined[i] 
  }
  else{ #for eclass data we can include only ili illness 
        dat$total_abs[i] <- dat$ili[i]+dat$undefined[i] 
  }#adjusted absenteeism rate
}

dat$exclude[(dat$level ==5 & dat$type=="S") | (dat$level==7 & dat$type=="S")] <-1 #exlcudeabsent rate for f.5 and f.7   

dat$exclude[ dat$total_abs / dat$total_num >0.5] <-1 #exlcude >50% absent rate for any class level due to special function
dat <- dat[dat$exclude ==0,]


#add total absent
#dat$total_abs <- dat$ili+ dat$gi+dat$other_ill+dat$not_ill+dat$undefined

#total ILI illness
#dat$total_abs <- dat$ili+dat$undefined

#get data from each school
dat$school_code <- as.numeric(dat$schoolcode)
dat$date <- dates(as.character(dat$date), format="y/m/d")

#define weeks, monday to friday
for (i in 1:dim(dat)[1]){

#if ( dat$date[i] >= dates("2008/2/11", format ="y/m/d")  &
#     dat$date[i] < dates("2008/2/18", format="y/m/d") 
#   )  {dat$week[i] <-7} 

 if ( dat$date[i] >= dates("2008/2/18", format ="y/m/d")  &
     dat$date[i] < dates("2008/2/25", format="y/m/d") 
   )  {dat$week[i] <-8} 

else if ( dat$date[i] >= dates("2008/2/25", format ="y/m/d")  &
     dat$date[i] < dates("2008/3/3", format="y/m/d") 
   ) { dat$week[i] <-9}
   
else if ( dat$date[i] >= dates("2008/3/3", format ="y/m/d")  &
     dat$date[i] < dates("2008/3/10" , format="y/m/d") 
   ) { dat$week[i] <-10}
   
else if ( dat$date[i] >= dates("2008/3/10", format ="y/m/d")  &
     dat$date[i] < dates("2008/3/17", format="y/m/d") 
   )  {dat$week[i] <-11}

else if ( dat$date[i] >= dates("2008/3/17", format ="y/m/d")  &
     dat$date[i] < dates("2008/3/23", format="y/m/d") 
   )  {dat$week[i] <-12}

else if ( dat$date[i] >= dates("2008/3/24", format ="y/m/d")  &
     dat$date[i] < dates("2008/3/30", format="y/m/d") 
   )  {dat$week[i] <-13}

else if ( dat$date[i] >= dates("2008/3/31", format ="y/m/d")  &
     dat$date[i] < dates("2008/4/6", format="y/m/d") 
   )  {dat$week[i] <-14}
else if ( dat$date[i] >= dates("2008/4/7", format ="y/m/d")  &
     dat$date[i] < dates("2008/4/13", format="y/m/d") 
   )  {dat$week[i] <-15}

else if ( dat$date[i] >= dates("2008/4/14", format ="y/m/d")  &
     dat$date[i] < dates("2008/4/20", format="y/m/d") 
   )  {dat$week[i] <-16}
else if ( dat$date[i] >= dates("2008/4/21", format ="y/m/d")  &
     dat$date[i] < dates("2008/4/27", format="y/m/d") 
   )  {dat$week[i] <-17}
else if ( dat$date[i] >= dates("2008/4/28", format ="y/m/d")  &
     dat$date[i] < dates("2008/5/4", format="y/m/d") 
   )  {dat$week[i] <-18}
else if ( dat$date[i] >= dates("2008/5/5", format ="y/m/d")  &
     dat$date[i] < dates("2008/5/11", format="y/m/d") 
   )  {dat$week[i] <-19}
else if ( dat$date[i] >= dates("2008/5/12", format ="y/m/d")  &
     dat$date[i] < dates("2008/5/18", format="y/m/d") 
   )  {dat$week[i] <-20}
else if ( dat$date[i] >= dates("2008/5/19", format ="y/m/d")  &
     dat$date[i] < dates("2008/5/25", format="y/m/d") 
   )  {dat$week[i] <-21}
else if ( dat$date[i] >= dates("2008/5/26", format ="y/m/d")  &
     dat$date[i] < dates("2008/6/1", format="y/m/d") 
   )  {dat$week[i] <-22}
else if ( dat$date[i] >= dates("2008/6/2", format ="y/m/d")  &
     dat$date[i] < dates("2008/6/8", format="y/m/d") 
   )  {dat$week[i] <-23}
else if ( dat$date[i] >= dates("2008/6/9", format ="y/m/d")  &
     dat$date[i] < dates("2008/6/15", format="y/m/d") 
   )  {dat$week[i] <-24}


 else dat$week[i] <-0  

}

#overall absence rate
overall <- vector(length = length(table(dat$week)))

for (i in 1:length(table(dat$week))){
 num <- sum(dat$total_num[dat$week ==(i+7)])
 abs <- sum (dat$total_abs[dat$week ==(i+7)])
 overall[i] <- abs/num 
}

names(overall) <- names(table(dat$week))
abs_dat <- overall

#by school
for (i in 1: length(table(dat$school_code))){
  
  tmp <- vector(length = length(table(dat$week)))

  for (j in 1:length(table(dat$week))){
    
    num <- sum(dat$total_num[dat$week ==(j+7) & dat$school_code == i])
    abs <- sum(dat$total_abs[dat$week ==(j+7) & dat$school_code == i])
    tmp[j] <- abs/num 
  
  }
   
  abs_dat <- rbind(abs_dat,tmp)
}

abs_dat <- t(abs_dat)
abs_dat <-as.data.frame(abs_dat)
abs_dat$month <- row.names(abs_dat)
names(abs_dat) <- c("overall", names(table(dat$schoolcode)), "week")

write.csv(abs_dat, paste(output_dir,today,"sch_abs_total.csv", sep=""), row.names=FALSE)
###############
