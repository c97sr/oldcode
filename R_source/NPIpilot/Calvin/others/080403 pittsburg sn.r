source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/070801 qv function.r")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

cdc <- read.csv("D:/Calvin flu studies/NPI study/080328 CDC all study summary/Pittsburgh/pittsburgh08data.csv")


############################
#2 by 2 table for sn and sp#
############################

#overall using culture as reference 
oritab <-table(cdc[c("quicktest","cultAB")])
maketab2(oritab)

#overall using pcr as reference
cdc$pcrpos[cdc$pcrAB !=0] <-1
cdc$pcrpos[is.na(cdc$pcrpos)] <-0
oritab <-table(cdc[c("quicktest","pcrpos")])
maketab2(oritab)

#separate by PCR intensity
pcrgp1 <- cdc[cdc$pcrAB ==1,]
oritab <-table(pcrgp1[c("quicktest","pcrAB")])
oritab

pcrgp2 <-cdc[cdc$pcrAB ==2,]
oritab <-table(pcrgp2[c("quicktest","pcrAB")])
oritab

pcrgp3 <-cdc[cdc$pcrAB ==3,]
oritab <-table(pcrgp3[c("quicktest","pcrAB")])
oritab

#separate by flu A and B
groupa <- cdc[cdc$pcrA > 0,]
oritab <-table(groupa[c("quicktest","pcrAB")])
oritab

groupb <- cdc[cdc$pcrB > 0,]
oritab <-table(groupb[c("quicktest","pcrAB")])
oritab