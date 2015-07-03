#
# Lab-confirmed non-secondary case: 
#    index has +ve culture at baseline (V0 or V1 based on culture and qPCR), contact has -ve culture for all 4 visits
#

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_16_home_culture.csv", header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_20_baseflu_m.csv", header=TRUE)

hculture <- data.frame(hhID = baseflu$hhID)
hculture$member <- 0
for(i in 2:nrow(hculture)){
  if(hculture$hhID[i]==hculture$hhID[i-1]) hculture$member[i] <- hculture$member[i-1]+1
}

j <- 1
i <- 1

## Read from dataset, one member per row
while(i<nrow(hculture)+1){  
  while(as.character(hculture$hhID[i])==as.character(hc$hhID[j])){
       if(hculture$member[i]==hc$member[j]){
          if((!is.na(hc$culture[j])&as.character(hc$culture[j])>0) | (!is.na(hc$AqPCR[j])&as.character(hc$AqPCR[j])>0) |
             (!is.na(hc$AqPCR2[j])&as.character(hc$AqPCR2[j])>0) | (!is.na(hc$AqPCR3[j])&as.character(hc$AqPCR3[j])>0) |
	     (!is.na(hc$BqPCR[j])&as.character(hc$BqPCR[j])>0) )
	     { hculture[i,hc$visit[j]+3] <- 1}
	     else { hculture[i,hc$visit[j]+3] <- 0}          
       j <- j+1
       } 
       else{i<-i+1}            
  }
  while(as.character(hculture$hhID[i])!=as.character(hc$hhID[j])){
       while(as.character(hculture$hhID[i])<as.character(hc$hhID[j])){i<-i+1}
       while(as.character(hculture$hhID[i])>as.character(hc$hhID[j])){j<-j+1}
  }
}

names(hculture) <- c("hhID","member","V0","V1","V2","V3","V4")

## Define excluded index and contact case
## index exclusion: none of V0/V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) 
                                  | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
				  | (is.na(hculture$V1[i]) & hculture$V0[i]==0)
				  | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
       {hculture$exd_index[i]<-1}
       else {hculture$exd_index[i]<-0}
}

## Define the household which should be excluded as long index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
dim(hculture)

## Define non-secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exd_index[i] == 0 & !is.na(hculture$V1[i]) & hculture$V1[i] ==0 &
           hculture$V2[i] ==0 & !is.na(hculture$V2[i]) & hculture$V3[i] ==0 & !is.na(hculture$V3[i])
               & hculture$V4[i] ==0 & !is.na(hculture$V4[i])  )
              {hculture$nsedcase[i] <- 1}
	 else {hculture$nsedcase[i] <- 0}
}

output <- hculture[hculture$nsedcase==1,-c(3,8,9)]

#######################################################
hculture$nsedcase[hculture$member==0] <- 1

i<- 5
while (i <1+nrow(hculture)){
      if (hculture$nsedcase[i]==0){ 
       while (hculture$hhID[i-1]==hculture$hhID[i]){
          i <- i-1
	  if (i==0) break
	  }
        hculture$nsedcase[i] <- 0
       while (hculture$hhID[i+1]==hculture$hhID[i]){
          i <- i+1
	  }
      }
      i<- i+1
}
hculture$nsedcase[1] <- 0
hculture.alln2nd <- hculture[hculture$member==0&hculture$nsedcase==1,]

set.seed(10)
hhID.sample <- hculture.alln2nd$hhID[sample(1:nrow(hculture.alln2nd),20,replace=FALSE)]
hhIDs <- rep(hhID.sample,each=5)

hc.sample <- data.frame( hhID = hhIDs)
hc.sample$member <- 0
hc.sample$visit <- 0
for(i in 2:nrow(hc.sample)){
  if(hc.sample$hhID[i]==hc.sample$hhID[i-1]) hc.sample$visit[i] <- hc.sample$visit[i-1]+1
}

hc.all <- merge(hc.sample,hc,all.x=T)
hc.alltest <- hc.all[is.na(hc.all$AqPCR)&is.na(hc.all$AqPCR2)&is.na(hc.all$AqPCR3)&is.na(hc.all$BqPCR),c(1,2,3,6)]

#write.csv(hc.alltest,"D:\\Work\\Influenza\\data\\2007_10_11\\2007_11_16_index_no2ndcase.csv",row.names=FALSE)