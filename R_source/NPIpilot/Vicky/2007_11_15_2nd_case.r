
#
# Exclude households (if index visit 0/1 culture is negative/NA&NA or if contact visit 1 culture is positive)
# Mark secondary cases (contacts +ve during followed-up days) and summarize their infected days
#

hc <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_11_14_home_culture.csv", header=TRUE)
baseflu <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_24_baseflu_m.csv", header=TRUE)

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
## Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) 
                                  | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
				  | (is.na(hculture$V1[i]) & hculture$V0[i]==0)
				  | (hculture$V0[i]==0 & hculture$V1[i]==0) ))     
          {hculture$exd_index[i]<-1}
       else {hculture$exd_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]==1) ))
          {hculture$exd_contact[i]=1}
	  else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
dim(hculture)

for (i in 1:nrow(hculture)){
    if ( hculture$exd_index[i]==1 | hculture$exd_contact[i] ==1)
       {hculture$exclude[i] <-1}
       else  {hculture$exclude[i] <-0}
}

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exclude[i] == 0 & !is.na(hculture$V1[i]) & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) | (hculture$V4[i] !=0 & !is.na(hculture$V4[i])) ) )
              {hculture$labsedcase[i] <- 1}
	 else {hculture$labsedcase[i] <- 0}
}

i<- 1
while (i <1+nrow(hculture)){
      if (hculture$labsedcase[i]==1){ 
      i<- i-1
       while (hculture$hhID[i]==hculture$hhID[i+1]){
          i <- i-1
	  if (i==0) break
	  }
	  i <- i+1
        hculture$labsedcase[i] <- 1
       i <- i+1
       while (hculture$hhID[i]==hculture$hhID[i-1]){
          i <- i+1
	  }
      }
      i<- i+1
}

#lab2nd <- unique(hculture$hhID[hculture$labsedcase==1])
#lab2nd <- as.character(lab2nd)
#
#hc.new <- hc[hc$hhID==lab2nd[1] | hc$hhID==lab2nd[2] | hc$hhID==lab2nd[3] | hc$hhID==lab2nd[4] | hc$hhID==lab2nd[5] | 
#             hc$hhID==lab2nd[6] | hc$hhID==lab2nd[7] | hc$hhID==lab2nd[8] | hc$hhID==lab2nd[9] | hc$hhID==lab2nd[10] |
#	     hc$hhID==lab2nd[11] | hc$hhID==lab2nd[12] | hc$hhID==lab2nd[13] | hc$hhID==lab2nd[14] | hc$hhID==lab2nd[15] |
#	     hc$hhID==lab2nd[16] | hc$hhID==lab2nd[17],]
#
#hc.new1 <- hc.new[is.na(hc.new$AqPCR),]

lab2ndID <- rep(hculture$hhID[hculture$labsedcase==1],each=4)
lab2ndm <- rep(hculture$member[hculture$labsedcase==1],each=4)
hc.new <- data.frame(hhID=lab2ndID)
hc.new$member <- lab2ndm
hc.new$visit <- 1
for(i in 2:nrow(hc.new)){
  if(hc.new$hhID[i]==hc.new$hhID[i-1]&hc.new$member[i]==hc.new$member[i-1]) hc.new$visit[i] <- hc.new$visit[i-1]+1
}

hc.index <- hculture[hculture$member==0&hculture$labsedcase==1,c(1,2)]
hc.index$visit <- 0

hc.list <- rbind(hc.new,hc.index)
hc.list <- hc.list[order(hc.list$hhID,hc.list$member,hc.list$visit),]

hc.all <- merge(hc.list,hc,all.x=T)
hc.alltest <- hc.all[is.na(hc.all$AqPCR)&is.na(hc.all$AqPCR2)&is.na(hc.all$AqPCR3)&is.na(hc.all$BqPCR),c(1,2,3,6)]

write.csv(hc.alltest,"D:\\Work\\Influenza\\data\\2007_10_11\\2007_11_16_2nd_case.csv",row.names=FALSE)
