source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/070801 qv function.r")
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")

# 080211 add the groups for the tables
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

#070824 get those culture result available data, show demographic result for presentation tables
ca <- cdc[(cdc$culture=="A"|cdc$culture=="B"|cdc$culture=="0")& !is.na(cdc$culture),] #get data from culture available subjects


 #change NA into 9
 ca$agegroup[is.na(ca$agegroup)] <-9
 ca$male[is.na(ca$male)] <-9
 ca$onsettime[is.na(ca$onsettime)] <-9


#define data into 4 groups QV, culture --,-+, +-, ++
#QV -ve, Culture -ve
qvncun <- ca[ca$QVpos==0 & ca$cultpos==0,]
qvncup <- ca[ca$QVpos==0 & ca$cultpos==1,]
qvpcun <- ca[ca$QVpos==1 & ca$cultpos==0,]
qvpcup <- ca[ca$QVpos==1 & ca$cultpos==1,]

# for getting one column of table 1
dat_tab <- function( totgp , retgp1, retgp2,retgp3,retgp4, num){
  re <- vector(length=8)
  #totnum <-  sum(1*(na.exclude(totgp == as.character (names(table(totgp))[num]))))
  retnum1 <- sum(1*(na.exclude(retgp1==as.character (names(table(totgp))[num]))))
  retnum2 <- sum(1*(na.exclude(retgp2==as.character (names(table(totgp))[num]))))
  retnum3 <- sum(1*(na.exclude(retgp3==as.character (names(table(totgp))[num]))))
  retnum4 <- sum(1*(na.exclude(retgp4==as.character (names(table(totgp))[num]))))
  re[1] <- retnum1
  re[2] <- retnum1/length(retgp1)
  re[3] <- retnum2
  re[4] <- retnum2/length(retgp2)
  re[5] <- retnum3
  re[6] <- retnum3/length(retgp3)
  re[7] <- retnum4
  re[8] <- retnum4/length(retgp4)

  return (re)
}


#for combining 1 group of table 1
com_tab <- function(cagp,qvncungp,qvpcungp, qvncupgp,qvpcupgp){
 out <- dat_tab(cagp, qvncungp, qvpcungp,qvncupgp,qvpcupgp,  1)
 
 for (i in 2: dim(table(cagp))) {
 temp <- dat_tab(cagp, qvncungp, qvpcungp,qvncupgp, qvpcupgp, i)
 out <- rbind(out, temp)
 }
 return (round(out,2))
}

#for different groups
com_tab(ca$agegroup, qvncun$agegroup,  qvpcun$agegroup,qvncup$agegroup, qvpcup$agegroup)
com_tab(ca$male, qvncun$male, qvpcun$male,qvncup$male,qvpcup$male)
com_tab(ca$onsettime, qvncun$onsettime, qvpcun$onsettime,qvncup$onsettime,qvpcup$onsettime)


#for symptoms
pos_dat <- function(gp0,gp1,gp2,gp3,gp4){
  
  total <- dim(gp0)[1]
  re <- vector(length=8)
  re[1] <- table(gp1)[[2]]
  re[2] <- table(gp1)[[2]]/ sum(table(gp1, exclude=NULL))
  re[3] <- table(gp2)[[2]]
  re[4] <- table(gp2)[[2]]/sum(table(gp2, exclude=NULL))
  re[5] <- table(gp3)[[2]]
  re[6] <- table(gp3)[[2]]/sum(table(gp3, exclude=NULL))
  re[7] <- table(gp4)[[2]]
  re[8] <- table(gp4)[[2]]/sum(table(gp4, exclude=NULL))

  return(re)
}

p_rnose <- pos_dat(ca, qvncun$rnose, qvpcun$rnose,  qvncup$rnose, qvpcup$rnose)
p_cough <- pos_dat(ca, qvncun$cough, qvpcun$cough,  qvncup$cough, qvpcup$cough)
p_sthroat <- pos_dat(ca, qvncun$sthroat, qvpcun$sthroat,  qvncup$sthroat, qvpcup$sthroat)
p_H_baseeadache <- pos_dat(ca, qvncun$headache, qvpcun$headache,  qvncup$headache, qvpcup$headache)
p_jointmuscle <- pos_dat(ca, qvncun$jointmuscle, qvpcun$jointmuscle,  qvncup$jointmuscle, qvpcup$jointmuscle)
p_phlegm <- pos_dat(ca, qvncun$phlegm, qvpcun$phlegm,  qvncup$phlegm, qvpcup$phlegm)
p_fever <- pos_dat(ca, qvncun$measure_fever, qvpcun$measure_fever,  qvncup$measure_fever, qvpcup$measure_fever)
p_fcs <- pos_dat(ca, qvncun$fcs, qvpcun$fcs,  qvncup$fcs, qvpcup$fcs)

round(rbind(p_rnose, p_cough, p_sthroat, p_H_baseeadache, p_jointmuscle, p_fever, p_fcs),2)

