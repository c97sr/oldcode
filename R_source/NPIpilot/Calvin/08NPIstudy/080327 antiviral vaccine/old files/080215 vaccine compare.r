
#comparing with and without vaccination with PCR, symptoms, fever, culture ...etc
#for PCR
pcr <- read.csv("P:\\GROUP\\NPIpilot\\data\\2008_01_home_PCR.csv", header=TRUE)
lab <- read.csv("P:\\GROUP\\NPIpilot\\data\\lab2nd.csv", header=TRUE)
visitday <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_housechar_h.csv", header=TRUE)

day <- data.frame(hhID <- rep(visitday$hhID,each=5))
names(day) <- "hhID"

day$visit <- rep(0:4,128) 
j <- 1
for (j in 1:128){
    day$day[5*j-4] <- visitday$clinic_day[j]
    day$day[5*j-3] <- visitday$v1_day[j]
    day$day[5*j-2] <- visitday$v2_day[j]
    day$day[5*j-1] <- visitday$v3_day[j]
    day$day[5*j] <- visitday$v4_day[j]
}
pcr <- merge(pcr,day,by=c("hhID","visit"),all.x=TRUE)
pcr <- pcr[order(pcr$hhID,pcr$member,pcr$visit),]

pcrindex <- pcr[pcr$member ==0,] # for PCR 


lab <- read.csv("P:\\GROUP\\NPIpilot\\data\\lab2nd.csv", header=TRUE)

dat <- lab[lab$member == 0 & lab$V0 ==1 ,]

#extract vaccination information
baseflu <- read.csv("P:/GROUP/NPIpilot/data/2007_11_20_baseflu_m.csv", header=TRUE)
baseflu <- baseflu[,c(1,3,7:9)]
out2 <- merge (dat, baseflu, by=c("hhID","member"), all.x=TRUE)

out2<- out2[,c(1:2,13:15)]

abc <- merge(pcrindex, out2, by=c("hhID","member"), all.x=TRUE)
log(mean(na.exclude(abc$AqPCR[abc$vaccine07==0])),10)	#a positive group
log(mean(na.exclude(abc$AqPCR[abc$vaccine07==1])),10) 

mean(na.exclude(abc$BqPCR[abc$vaccine07==0]))# b positive group
mean(na.exclude(abc$BqPCR[abc$vaccine07==1]))

# for symptoms and fever
symptom <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_symptomday_d.csv", header=TRUE)
symptom$fever <- 1*(symptom$bodytemp>=37.8)
for ( i in 1:nrow(symptom)){
symptom$score[i] <- sum(symptom$fever[i],symptom$headache[i],symptom$rnose[i],symptom$hvoice[i],symptom$sthroat[i],
                        symptom$cough[i],symptom$sbreath[i],symptom$cpain[i],symptom$pmuscle[i],symptom$tired[i],na.rm=TRUE)
}

plot_sym <- data.frame(hhID=rep(lab2nd$hhID,each=10))
plot_sym$member <- rep(lab2nd$member,each=10)
plot_sym$day <- rep(0:9,101)
plot_fever <- merge(plot_sym,symptom[,c(1,3:5,25)],by=c("hhID","member","day"),all.x=TRUE)   # 25 is score
plot_fever$bodytemp <- (plot_fever$bodytemp-34.5)/2.0

# match symptom diary day (0-9) with index symptom onset day (as day 0)
fever_day <- pcr_age[pcr_age$visit==1,c(1,5)]  # hhID, day
# delete duplicated records (more than 1 seconday cases in the same household)
fever_day <- fever_day[order(fever_day$hhID),]
fever_day$mark <- 1
for (i in 2:99){
    if (fever_day$hhID[i]==fever_day$hhID[i-1]) fever_day$mark[i] <- 0
}
fever_day <- fever_day[fever_day$mark==1,-3]
# end
plot_fever <- merge(plot_fever,fever_day, by="hhID",all.x=T)
plot_fever$day <- plot_fever$day.x+plot_fever$day.y
plot_fever$pscore <- (plot_fever$score+2.5)/1.25

### Construct data frame for index bodytemp
plot_syma <- data.frame(hhID=rep(unique(pcr_indexw2$hhID),each=10))
plot_syma$day <- rep(0:9,length(unique(pcr_indexw2$hhID)))
plot_syma$member <- 0
plot_fever_a <- merge(plot_syma,symptom[,c(1,3:5,25)],by=c("hhID","member","day"),all.x=T)   # 25 is score
plot_fever_a$bodytemp <- (plot_fever_a$bodytemp-34.5)/2.0

# match symptom diary day (0-9) with index symptom onset day (as day 0)
plot_fever_a <- merge(plot_fever_a,fever_day, by="hhID",all.x=T)
plot_fever_a$day <- plot_fever_a$day.x + plot_fever_a$day.y

#symptom
abc <- merge(ss_index_plot, out2, by=c("hhID","member"), all.x=TRUE)
mean(na.exclude(abc$pscore[abc$vaccine07==0]))	#symptom score
mean(na.exclude(abc$pscore[abc$vaccine07==1])) 

#fever
abc <- merge(plot_fever_a, out2, by=c("hhID","member"), all.x=TRUE)
mean(na.exclude(abc$score[abc$vaccine07==0]))	 #fever score
mean(na.exclude(abc$score[abc$vaccine07==1])) 

mean(na.exclude(abc$bodytemp[abc$vaccine07==0])) #bodytemp
mean(na.exclude(abc$bodytemp[abc$vaccine07==1])) 




write.csv(out2, "d:/r/temp/2ndcaseID.csv",  row.names=FALSE)