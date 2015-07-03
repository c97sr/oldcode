#080121 get those with PCR results 
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

cdcPCR <-cdc[!is.na(cdc$qPCR),]
cdcPCR$logqPCR <- log(cdcPCR$qPCR,10)
##############
# area under roc plot log(qPCR) vs QV 

library(ROCR)
dat <-cdcPCR[cdcPCR$culture=="A",]

p1 <- prediction(dat$logqPCR, dat$QVpos )
ss <- performance(p1, "sens", "spec")
round(ss@y.values[[1]][2],2)  # sensitivity
round(ss@x.values[[1]][2],2)  # specificity
p2 <- performance(p1, "auc")
round(p2@y.values[[1]],2)     # area under roc

plot(ss)

dat <-cdcPCR[cdcPCR$culture=="B",]
p1 <- prediction(dat$logqPCR, dat$QVpos )
ss <- performance(p1, "sens", "spec")
round(ss@y.values[[1]][2],2)  # sensitivity
round(ss@x.values[[1]][2],2)  # specificity
p2 <- performance(p1, "auc")
round(p2@y.values[[1]],2)     # area under roc

plot(ss)

#lines(ksmooth(ss$x.values,"box",2),col="BLACK")
