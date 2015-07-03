###### Include all index subjects (944) ######

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_14_home_culture.csv", header=TRUE)
cdat <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_24_clinicdat_h.csv", header=TRUE)
cc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_14_clinic_culture.csv", header=TRUE)

# Define culture result by combining culture, qPCR results together
hc0 <- hc[hc$member==0&hc$visit==0,-c(2:5)]
for (i in 1:nrow(hc0)){
      if( (!is.na(hc0$AqPCR[i])&as.character(hc0$AqPCR[i])>0) |
          (!is.na(hc0$AqPCR2[i])&as.character(hc0$AqPCR2[i])>0) | (!is.na(hc0$AqPCR3[i])&as.character(hc0$AqPCR3[i])>0))
	   hc0$culture[i]<- "A"
      if( !is.na(hc0$BqPCR[i])&as.character(hc0$BqPCR[i])>0)
	   hc0$culture[i]<- "B"
}
hc1 <- hc0[,-1]

# Combine clinic and home culture results 
hc1$scrID <- tolower(hc1$scrID)
cc$AqPCR <- cc$AqPCR2 <- cc$AqPCR3 <- cc$BqPCR <- NA
qvc.c <- rbind(cc,hc1)

# Merge QV and culture results
qvc <- data.frame(scrID = cdat$scrID)
qvc$hhID <- cdat$hhID
qvc$QVres <- cdat$QVres
qvc.all <- merge(qvc,qvc.c,by="scrID",all.x=T)

# Delete duplicated rows
qvc.all$mark <- rep(NA,nrow(qvc.all))
for (i in 2:nrow(qvc.all)){
  if(qvc.all$scrID[i]==qvc.all$scrID[i-1]) qvc.all$mark[i]<-1
}
qvc.all <- qvc.all[is.na(qvc.all$mark),-9]

# Assign A/B/0 for QuickVue results
qvc.new <- qvc.all
for (i in 1:nrow(qvc.new)){
     if(qvc.new$QVres[i]==1) qvc.new$QVres[i] <- "A"
     else if(qvc.new$QVres[i]==2) qvc.new$QVres[i] <- "B"
     else if(qvc.new$QVres[i]==3) qvc.new$QVres[i] <- 0
     else if(qvc.new$QVres[i]==4) qvc.new$QVres[i] <- "Inconclusive"
     else qvc.new$QVres[i] <- NA
}

# qvc.new

# random sample QV+ culture+

sample.size <- 20

qvc.pp <- qvc.new[!is.na(qvc.new$QVres) & (qvc.new$QVres=="A" | qvc.new$QVres=="B" ) 
  & !is.na(qvc.new$culture) & (qvc.new$culture=="A" | qvc.new$culture=="B") & qvc.new$QVres==qvc.new$culture,] 
dim(qvc.pp)
table(qvc.pp$QVres, qvc.pp$culture)

set.seed(10)
qvc.pp.selected <- qvc.pp[sample(1:nrow(qvc.pp), sample.size, replace=FALSE),]
#qvc.pp.selected <- qvc.pp.selected[order(qvc.pp.selected$scrID),]
qvc.pp.selected

# random sample QV- culture+

qvc.np <- qvc.new[!is.na(qvc.new$QVres) & qvc.new$QVres=="0"  
  & !is.na(qvc.new$culture) & (qvc.new$culture=="A" | qvc.new$culture=="B"),] 
dim(qvc.np)

set.seed(20)
qvc.np.selected <- qvc.np[sample(1:nrow(qvc.np), sample.size, replace=FALSE),]
#qvc.np.selected <- qvc.np.selected[order(qvc.np.selected$scrID),]
qvc.np.selected

# random sample QV- culture-

qvc.nn <- qvc.new[!is.na(qvc.new$QVres) & qvc.new$QVres=="0"  
  & !is.na(qvc.new$culture) & qvc.new$culture==0,] 
dim(qvc.nn)

set.seed(30)
qvc.nn.selected <- qvc.nn[sample(1:nrow(qvc.nn), sample.size, replace=FALSE),]
#qvc.nn.selected <- qvc.nn.selected[order(qvc.nn.selected$scrID),]
qvc.nn.selected

qvc.a <- rbind(qvc.pp.selected,qvc.np.selected,qvc.nn.selected)
qvc.sample <- qvc.a[is.na(qvc.a$AqPCR)&is.na(qvc.a$AqPCR2)&is.na(qvc.a$AqPCR3)&is.na(qvc.a$BqPCR),-c(2,5:8)]
qvc.sample <- qvc.sample[order(qvc.sample$scrID),]

write.csv(qvc.sample,"D:\\Work\\Influenza\\data\\2007_10_11\\2007_11_15_sample.csv",row.names=FALSE)

#
# the end
#
#