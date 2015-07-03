
dat <- read.csv("P:/GROUP/NPIpilot/data/2007_11_27_housechar_h.csv", header=TRUE)


# Plot side by side 3 barcharts 
# 1. the delay of symptom onset to clinic visit
# 2. the delay of clinic visit to first home visit
# 3. the delay from symptom onset to the first home visit (in days)

#function for getting the delay days for plotting
delayday <- function(dataset){
  temp1 <- as.matrix(dataset[dataset==0])
  d0 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==1])
  d1 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==2])
  d2 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==3])
  d3 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==4])
  d4 <- dim(na.omit(temp1))[1]
  result <- c(d0,d1,d2,d3,d4)
  return (result) 
 }
#end

#plot side by side 
windows(width=8, height=12)
layout(matrix(1:3, nrow=3))
par(mar=c(5,5,2,2))

#horizonal barplot of 1st symptom appear to clinic of subject recruited by day

#get the delay days data

# plot the data

da <- delayday(dat$clinic_day)
names(da)<-c(0,1,2,3,4)

 barplot(rev(da), width = 1, horiz = TRUE, xlim=(c(0,80)), ylim=(c(0,5)), 
 xlab="", axisnames=FALSE, axes= FALSE,
 ylab="Delay (days)", cex.lab=1.5, cex.axis=1.5,
 space=0, las=1)
title("(a)",cex.main=1.5)
axis(1, at=0:8*10, 
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)


#horizonal barplot of subjects from clinic recruitment to 1st home visit 

between <-  dat$v1_day - dat$clinic_day
da <- delayday(between)
names(da)<- c(0,1,2,3,4)

barplot(rev(da), , width = 1, space=0, xlim=c(0,80), ylim =c(0,5), horiz = TRUE,
  xlab="", ylab="Delay (days)",las=1, axisnames=FALSE, axes= FALSE,
  cex.lab=1.5, cex.axis=1.5)
title("(b)",cex.main=1.5)
axis(1, at=0:8*10, 
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)

#horizonal barplot of subjects from symptom first appear to 1st home visit 

da <- delayday(dat$v1_day)

names(da)<- c(0,1,2,3,4)
barplot(rev(da), , width = 1, space=0, xlim=c(0,80), ylim =c(0,5), horiz = TRUE,
  axisnames=FALSE, axes= FALSE, 
  xlab="Number of index cases", cex.axis=1.5, ylab="Delay (days)",las=1, cex.lab=1.5 )
title("(c)",cex.main=1.5)
axis(1, at=0:8*10,
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)
