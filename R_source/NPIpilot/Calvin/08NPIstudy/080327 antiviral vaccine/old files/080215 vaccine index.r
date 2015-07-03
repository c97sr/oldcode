
#
# Viral load for secondary cases
#

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

pcrindex <- pcr[pcr$member ==0,]


#lab2nd <- lab[lab$labsedcase==1,]
lab2nd <- lab[lab$member == 0 & lab$V0 ==1 ,] #get the index cases only

hhIDindex <- unique(lab2nd$hhID)
pcr_2nd <- pcr[as.character(pcr$hhID)==as.character(lab2nd$hhID)[1] & pcr$member==lab2nd$member[1],]
#    pcr_indexw2 <- pcr[as.character(pcr$hhID)==as.character(hhIDindex)[1] & pcr$member==0,]   # corresponding index

for (i in 2:nrow(lab2nd)){
   temp <- pcr[as.character(pcr$hhID)==as.character(lab2nd$hhID)[i] & pcr$member==lab2nd$member[i],]
   pcr_2nd <- rbind(pcr_2nd,temp)
}
#    for (i in 2:length(hhIDindex)){
#        temp <- pcr[as.character(pcr$hhID)==as.character(hhIDindex)[i] & pcr$member==0,]
#        pcr_indexw2 <- rbind(pcr_indexw2,temp)                                                 # corresponding index
#    }

for ( i in 1:nrow(pcr_2nd)){
   pcr_2nd$pcrplot[i] <- max(max(pcr_2nd$AqPCR[i],pcr_2nd$BqPCR[i],pcr_2nd$AqPCR2[i],na.rm=TRUE),900)
}
     for ( i in 1:nrow(pcr_indexw2)){                                                          # corresponding index
        pcr_indexw2$pcrplot[i] <- max(max(pcr_indexw2$AqPCR[i],pcr_indexw2$BqPCR[i],pcr_indexw2$AqPCR2[i],na.rm=TRUE),900)
     }


## plot the graph of 2nd case separately by age (show fever)

###### add age for pcr results
demog <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_27_demographic_m.csv", header=TRUE)
pcr_2nd_tmp <- pcr_2nd[,c(1:4,11,12)]
pcr_2nd_age <- merge(pcr_2nd_tmp,demog[,c(1,3,4)],by=c("hhID","member"),all.x=TRUE)
pcr_age <- pcr_2nd_age[order(pcr_2nd_age$age),]


## Get age and flu type of corresponding index
#hhID2nd_df <- data.frame(hhID=unique(pcr_2nd_age$hhID))
#age_indexa <- merge(hhID2nd_df,demog[demog$member==0,c(1,4)],by="hhID",all.x=T)
#names(age_indexa)[2] <- "indexage"
# add culture results
culture_0a <- pcr_indexw2[pcr_indexw2$visit==0,c(1,4)]  # hhID,culture
culture_1a <- pcr_indexw2[pcr_indexw2$visit==1,c(1,4)]  # hhID,culture
age_indexa$indexculture <- culture_0a$culture
for (i in 1:nrow(age_indexa)){
    if (age_indexa$indexculture[i]==0) age_indexa$indexculture[i] <- culture_1a$culture[i]
}

## Construct a data frame as reference: hhID,member,culture(A/B),age. Ordered by culture then by age.
# start from visit 2, since baseline culture should be negative!
culture_2 <- pcr_age[pcr_age$visit==2,c(1,2,4,7)]  # hhID, member, culture, age
culture_3 <- pcr_age[pcr_age$visit==3,c(1,2,4)]
culture_4 <- pcr_age[pcr_age$visit==4,c(1,2,4)]
culture_tmp23 <- merge(culture_2,culture_3, by=c("hhID","member"),all.x=TRUE)
culture_tmp234 <- merge(culture_tmp23,culture_4, by=c("hhID","member"),all.x=TRUE)
names(culture_tmp234)[6] <- "culture.z"

age_contact <- culture_tmp234[,c(1,2,4)] # hhID, member, age
age_contact$culture <- culture_tmp234$culture.x

# add by me
age_contact$culture [is.na(age_contact$culture)] <-0
culture_tmp234$culture.y[is.na(culture_tmp234$culture.y)] <-0
culture_tmp234$culture.z[is.na(culture_tmp234$culture.z)] <-0

for (i in 1:nrow(age_contact)){
     if(age_contact$culture[i]==0) age_contact$culture[i] <- culture_tmp234$culture.y[i]
     if(age_contact$culture[i]==0) age_contact$culture[i] <- culture_tmp234$culture.z[i]    
}

## There are 2 secondary cases with culture -ve for all visits but have +ve qPCR results, so need to replace it
id_cul0 <- age_contact[age_contact$culture==0,]
id_cul0_tmp <- merge(id_cul0,pcr_2nd,by=c("hhID","member"))
for (i in 1:nrow(id_cul0)){
    id_cul0$pcrA[i] <- sum(id_cul0_tmp$AqPCR[id_cul0_tmp$hhID==id_cul0$hhID[i] & id_cul0_tmp$member==id_cul0$member[i]],na.rm=T)
    id_cul0$pcrB[i] <- sum(id_cul0_tmp$BqPCR[id_cul0_tmp$hhID==id_cul0$hhID[i] & id_cul0_tmp$member==id_cul0$member[i]],na.rm=T)
}
for (i in 1:nrow(id_cul0)){
    if(id_cul0$pcrA[i]>0) id_cul0$culture[i] <- "A"
    if(id_cul0$pcrB[i]>0) id_cul0$culture[i] <- "B"
}

age_contact_tmp <- merge(age_contact,id_cul0[,c(1,2,4)],by=c("hhID","member"),all.x=TRUE)
for (i in 1:nrow(age_contact_tmp)){
     if(age_contact_tmp$culture.x[i]==0) age_contact_tmp$culture.x[i] <- age_contact_tmp$culture.y[i]
}
age_cul_cont <- age_contact_tmp[,-5]
names(age_cul_cont)[4] <- "culture"
age_cul_cont <- merge(age_cul_cont,age_indexa,by="hhID",all.x=TRUE)  # hhID,contact age & flu type, index age & flu type
age_cul_cont <- age_cul_cont[order(age_cul_cont$indexculture,age_cul_cont$indexage),]

###### construct info. for fever and symptom score
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


#plot_sym <- data.frame(hhID=rep(lab2nd$hhID,each=10))
#plot_sym$member <- rep(lab2nd$member,each=10)
#plot_sym$day <- rep(0:9,21)
#plot_fever <- merge(plot_sym,symptom[,c(1,3:5,25)],by=c("hhID","member","day"),all.x=TRUE)   # 25 is score
#
#plot_fever$fever <- 2.5*(plot_fever$bodytemp>=37.8)
#plot_fever$fever[is.na(plot_fever$fever)] <- 0
#
## match symptom diary day (0-9) with index symptom onset day (as day 0)
#fever_day <- pcr_age[pcr_age$visit==1,c(1,5)]  # hhID, day
## delete duplicated records (more than 1 seconday cases in the same household)
#fever_day <- fever_day[order(fever_day$hhID),]
#fever_day$mark <- 1
#for (i in 2:21){
#    if (fever_day$hhID[i]==fever_day$hhID[i-1]) fever_day$mark[i] <- 0
#}
#fever_day <- fever_day[fever_day$mark==1,-3]
## end
#plot_fever <- merge(plot_fever,fever_day, by="hhID",all.x=T)
#plot_fever$day <- plot_fever$day.x+plot_fever$day.y
#plot_fever$pscore <- (plot_fever$score+2.5)/1.25

###### construct info. for index symptom score
ss_index <- data.frame(hhID=rep(unique(pcr_2nd$hhID),each=10))
ss_index$member <- 0
ss_index$day <- rep(0:9,101)
ss_index_plot <- merge(ss_index,symptom[,c(1,3,4,25)],by=c("hhID","member","day"),all.x=TRUE)
ss_index_plot <- merge(ss_index_plot,fever_day,by="hhID",all.x=T)
ss_index_plot$day <- ss_index_plot$day.x + ss_index_plot$day.y
ss_index_plot$pscore <- (ss_index_plot$score+2.5)/1.25


pdf("D:\\r\\temp\\index_vload_symp.pdf", width=12, height=108, version="1.4")

pcr_age$pch_tmp <- NA
for ( i in 1:nrow(pcr_age)){
    if(pcr_age$culture[i] ==0) pcr_age$pch_tmp[i] <- "-"   
    else pcr_age$pch_tmp[i] <- "+"
}

layout(matrix(c(1,1,1,2:103), ncol=3, byrow=T),heights = c(0.5, rep(2.5,7)))
par(mar=c(0,4,2,4))
plot(0, type="n", axes=FALSE, xlab="", ylab="")
#mtext("Viral loads and symptom scores for 21 observed transmission pairs",cex=1.5,line=-1)

for (i in 1:101){
par(mar=c(4,6,0,4))
plot(0, type="n", axes=FALSE, xlim=c(0,14.8), ylim=c(1,10), xlab="", ylab="")

## symptom score for corresponding index
polygon(ss_index_plot$day[ss_index_plot$hhID==as.character(age_cul_cont$hhID[i]) & ss_index_plot$member==0][c(1,1:10,10)], 
        c(2,ss_index_plot$pscore[ss_index_plot$hhID==as.character(age_cul_cont$hhID[i]) & ss_index_plot$member==0],2),
        col=rgb(0,1,0,0.4),lty=0)

## symptom score for contact
#polygon(plot_fever$day[plot_fever$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever$member==age_cul_cont$member[i]][c(1,1:10,10)], 
#        c(2,plot_fever$pscore[plot_fever$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever$member==age_cul_cont$member[i]],2),
#        col=rgb(0,0,1,0.2),lty=0)

## Viral load for corresponding index
lines(type= "o", pcr_indexw2$day[pcr_indexw2$hhID==as.character(age_cul_cont$hhID[i])], 
                   log10( pcr_indexw2$pcrplot[pcr_indexw2$hhID==as.character(age_cul_cont$hhID[i])] ),
		   lty=2,col="purple",cex=1.5,pch=17)

## viral load for contact
#lines(type= "o", pcr_age$day[pcr_age$hhID==as.character(age_cul_cont$hhID[i]) & pcr_age$member==age_cul_cont$member[i]], 
#                   log10( pcr_age$pcrplot[pcr_age$hhID==as.character(age_cul_cont$hhID[i]) & pcr_age$member==age_cul_cont$member[i]] ),
#		   lty=1,col="blue",cex=1.8,pch=19)

## culture results (+/-)
lines( type= "p", pcr_age$day[pcr_age$hhID==as.character(age_cul_cont$hhID[i]) & pcr_age$member==age_cul_cont$member[i]]+0.5, 
                 log10(pcr_age$pcrplot[pcr_age$hhID==as.character(age_cul_cont$hhID[i]) & pcr_age$member==age_cul_cont$member[i]])+0.5,
		 lty=1,col="black",cex=2.0, 
		 pch=pcr_age$pch_tmp[pcr_age$hhID==as.character(age_cul_cont$hhID[i]) & pcr_age$member==age_cul_cont$member[i]] )

## line for fever definition (bodytemp>=37.8)
abline(h=(37.8-34.5)/2.0, col=gray(0.5),lty=2)
## body temperature for index
lines(type= "l", plot_fever_a$day[plot_fever_a$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever_a$member==0], 
                    plot_fever_a$bodytemp[plot_fever_a$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever_a$member==0] ,
		   lty=1,col="orange",lwd=2)
## body temperature for contact
#lines(type= "l", plot_fever$day[plot_fever$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever$member==age_cul_cont$member[i]], 
#                    plot_fever$bodytemp[plot_fever$hhID==as.character(age_cul_cont$hhID[i]) & plot_fever$member==age_cul_cont$member[i]] ,
#		   lty=1,col="red",lwd=2)

## bottom line
abline(h=log10(900), col=gray(0.5))

axis(1, at=2*0:7,labels=c(0,2,4,6,8,10,12,14),cex.axis=1.5)         
axis(2, at=2*1:5,labels=c(expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10) ), cex.axis=1.5,las=1)
axis(4, at=2*1:5,labels=c(2.5*(0:4)),col=gray(0.6),las=1,cex.axis=1.5)
legend(7,11,paste("Flu",age_cul_cont$culture[i],", age",age_cul_cont$age[i]),bty="n",cex=1.5,text.col="blue")
legend(-1,11,paste("Flu",age_cul_cont$indexculture[i],", age",age_cul_cont$indexage[i]),bty="n",cex=1.5,text.col="purple")
mtext("37.8",side=1,line=-2.5,at=13,cex=0.8)

#if (i>=19) mtext ("Day from index symptom onset", side=1,line=2.8,cex=1.2)
if ((i+2)/3==round((i+2)/3)) mtext("Viral load",side=2,line=3.5,cex=1.5)
}

dev.off()

