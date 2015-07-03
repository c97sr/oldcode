#sensitivity of each  

cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

oritab <-table(cdc[c("QVpos","cultpos")])
maketab2(oritab)
# 071123 for flu A and flu B
oritab <-table(cdc[c("QVposA","cultposA")])
maketab2(oritab)

oritab <-table(cdc[c("QVposB","cultposB")])
maketab2(oritab)


#extract only culture +ve A and culture +ve B

cultA <- cdc[cdc$cultposA ==1,]
cultB <- cdc[cdc$cultposB ==1,]

#function for getting the sensitivity for A +ve
snres <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]

print(table(tab$QVpos))
#print (length(tab$QVpos))
sn <- binom.test(x=table(tab$QVpos)[[2]], dim(tab)[1])


output <- data.frame(total = length(tab$QVpos),
  percent =  length(tab$QVpos)/dim(cdc)[1],
  Value = c(sn[[1]] / sn[[2]]),
  lowerCI = c(sn[[4]][1]),
  upperCI = c(sn[[4]][2]),
  row.names=c("Sensitivity"))
round(output, 2)
}

# for agegroup
a1 <- snres(cultA, cultA$agegroup, 1)
a2 <- snres(cultA, cultA$agegroup, 2)
a3 <- snres(cultA, cultA$agegroup, 3)
a4 <- snres(cultA, cultA$agegroup, 4)
a5 <- snres(cultA, cultA$agegroup, 5)
a6 <- snres(cultA, cultA$agegroup, 6)
out <- rbind (a1,a2,a3,a4,a5,a6)
row.names(out)<- c("age 0-5","5-10","10-15","16-30","31-50","50+")
out

a1 <- snres(cultB, cultB$agegroup, 1)
a2 <- snres(cultB, cultB$agegroup, 2)
a3 <- snres(cultB, cultB$agegroup, 3)
a4 <- snres(cultB, cultB$agegroup, 4)
a5 <- snres(cultB, cultB$agegroup, 5)
a6 <- snres(cultB, cultB$agegroup, 6)
out <- rbind (a1,a2,a3,a4,a5,a6)
row.names(out)<- c("age 0-5","5-10","10-15","16-30","31-50","50+")
out

#sex
a1 <- snres(cultA, cultA$male, 1)
a2 <- snres(cultA, cultA$male, 0)
out <- rbind (a1,a2)
row.names(out)<- c("male","female")
out

a1 <- snres(cultB, cultB$male, 1)
a2 <- snres(cultB, cultB$male, 0)
out <- rbind (a1,a2)
row.names(out)<- c("male","female")
out

#symptom onset
a1 <- snres(cultA, cultA$onsettime, 1)
a2 <- snres(cultA, cultA$onsettime, 2)
a3 <- snres(cultA, cultA$onsettime, 3)
a4 <- snres(cultA, cultA$onsettime, 4)
# over 48hrs too few cases
out <- rbind (a1,a2,a3,a4)
row.names(out)<- c("onset <12hrs","12-24hrs","24-36hrs","36-48hrs")
out

a1 <- snres(cultB, cultB$onsettime, 1)
a2 <- snres(cultB, cultB$onsettime, 2)
a3 <- snres(cultB, cultB$onsettime, 3)
a4 <- snres(cultB, cultB$onsettime, 4)
# over 48hrs too few cases

out <- rbind (a1,a2,a3,a4)
row.names(out)<- c("onset <12hrs","12-24hrs","24-36hrs","36-48hrs")
out

#clinic experience
a1 <- snres(cultA, cultA$qvdonegp, 1)
a2 <- snres(cultA, cultA$qvdonegp, 2)
a3 <- snres(cultA, cultA$qvdonegp, 3)
a4 <- snres(cultA, cultA$qvdonegp, 4)
out <- rbind (a1,a2,a3,a4)
row.names(out)<- c("QV <30","30-60","60-90",">90")
out

a1 <- snres(cultB, cultB$qvdonegp, 1)
a2 <- snres(cultB, cultB$qvdonegp, 2)
a3 <- snres(cultB, cultB$qvdonegp, 3)
a4 <- snres(cultB, cultB$qvdonegp, 4)

out <- rbind (a1,a2,a3,a4)
row.names(out)<- c("QV <30","30-60","60-90",">90")
out