source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/070801 qv function.r")


cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

############################
#2 by 2 table for sn and sp#
############################

#overall 
oritab <-table(cdc[c("QVpos","cultpos")])
maketab2(oritab)


#for adult and children
a1<- snspresult (cdc, cdc$agegroup,1) #
tmp1 <- cbind( a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$agegroup,2) #
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<- snspresult (cdc, cdc$agegroup,3) #
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$agegroup,4) #
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
a5<-snspresult (cdc, cdc$agegroup,5) #
tmp5 <- cbind(a5[1,1:3],a5[2,1:4])
a6<-snspresult (cdc, cdc$agegroup,6) #
tmp6 <- cbind(a6[1,1:3],a6[2,1:4])


out <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6) #combining
out <-out[c(7,1:6)] #rearranging order
names(out)<- c( "Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI") #renaming
row.names(out)<- c("age 0-5","6-10","11-15","16-30","31-50","50+") #renaming
out

# for sex  
a1<- snspresult (cdc, cdc$male,1) #male
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2 <- snspresult (cdc, cdc$male,0) #female
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
out <- rbind(tmp1,tmp2)
out <-out[c(7,1:6)]
#renaming
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("male","female")
out

#for symptom onset 12,24,36,48 hrs 
a1<- snspresult (cdc, cdc$onsettime,1) #symptom onet <12 hours
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$onsettime,2) #symptom onet <24 hours
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<-snspresult (cdc, cdc$onsettime,3) #symptom onet <36 hours
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$onsettime,4) #symptom onet <48 hours
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
a5<- snspresult (cdc, cdc$onsettime,5) 
tmp5 <- cbind(a5[1,1:3],a5[2,1:4])
out <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5)

out <-out[c(7,1:6)]
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("onset <12hrs","12-24hrs","24-36hrs","36-48hrs",">48hrs")
out

#071204 for experienced in qv done
a1<- snspresult (cdc, cdc$qvdonegp,1)
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$qvdonegp,2)
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<- snspresult (cdc, cdc$qvdonegp,3)
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$qvdonegp,4)
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
out <- rbind(tmp1,tmp2,tmp3,tmp4)
out <-out[c(7,1:6)]
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("QV <30","30-60","60-90",">90")
out



#add 080502 p value (fisher's test, chi-sq test)
# age sn
chisq.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$agegroup[cdc$cultpos!=0]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$agegroup[cdc$cultpos!=0]))

#sex
chisq.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$male[cdc$cultpos!=0]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$male[cdc$cultpos!=0]))

#onsettime
chisq.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$onsettime[cdc$cultpos!=0]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$onsettime[cdc$cultpos!=0]))

#qvdonegp
chisq.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$qvdonegp[cdc$cultpos!=0]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=0], cdc$qvdonegp[cdc$cultpos!=0]))



# age sp
chisq.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$agegroup[cdc$cultpos!=1]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$agegroup[cdc$cultpos!=1]))

#sex
chisq.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$male[cdc$cultpos!=1]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$male[cdc$cultpos!=1]))

#onsettime
chisq.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$onsettime[cdc$cultpos!=1]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$onsettime[cdc$cultpos!=1]))

#qvdonegp
chisq.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$qvdonegp[cdc$cultpos!=1]))
fisher.test(table(cdc$QVpos[cdc$cultpos!=1], cdc$qvdonegp[cdc$cultpos!=1]))