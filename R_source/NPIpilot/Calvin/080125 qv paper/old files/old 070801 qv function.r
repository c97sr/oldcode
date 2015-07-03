
# useful functions

# generate random number
rand <-function(x){
set.seed(100)
sort(sample(1:200, x))
}
#end

#
#function for display results when put in different parameters 
#

snspresult <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]
  tab <- table(tab[c("QVpos","cultpos")])
  maketab2(tab)
  #070927 add more test stat such as ppv, npv lr+ lr-...etc 
  #teststat2by2(maketab2(tab))
}
#end

snspresultA <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]
  tab <- table(tab[c("QVposA","cultposA")])
  maketab2(tab)
  #070927 add more test stat such as ppv, npv lr+ lr-...etc 
  #teststat2by2(maketab2(tab))
}
#end

snspresultB <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]
  tab <- table(tab[c("QVposB","cultposB")])
  maketab2(tab)
  #070927 add more test stat such as ppv, npv lr+ lr-...etc 
  #teststat2by2(maketab2(tab))
}
#end


#
#function for making a 2 by 2 table "maketab"
#

maketab <- function(oritab){
  print(oritab)
  tab <-matrix(NA, ncol=2, nrow=2)
  tab[1,1]<-oritab[2,2]+ oritab[2,3]+ oritab[3,2]+ oritab[3,3]
  tab[1,2]<-oritab[2,1]+ oritab[3,1]
  tab[2,1]<-oritab[4,2]+ oritab[4,3]
  tab[2,2]<-oritab[4,1]
  print(tab)
  print(sum(tab))
  snsp(tab)
}
# end

#function for rearrange the 2 by 2 table to ++ +- -+ -- and display the result
maketab2 <-function(oritab){
  tab <-matrix(NA, ncol=2, nrow=2)
  tab[1,1]<-oritab[2,2]
  tab[1,2]<-oritab[2,1]
  tab[2,1]<-oritab[1,2]
  tab[2,2]<-oritab[1,1]
  print(tab)
  print(sum(tab))
  
  snsp(tab)
  
  #return (tab)
}

#
# function for calculating sens, spec and exact 95% CIs "snsp"
#

snsp <- function(tab){
cdc <-tab
total <- sum(tab)
# 95% CI for sensitivity
sn <- binom.test(x=cdc[1,1], n=cdc[1,1] + cdc[2,1])

# 95% CI for specificity
sp <- binom.test(x=cdc[2,2], n=cdc[1,2] + cdc[2,2])

#95% CI for PPV
ppv <- binom.test(x=cdc[1,1], n=cdc[1,1] + cdc[1,2])
#95% CI for NPV
npv <- binom.test(x=cdc[2,2], n=cdc[2,1] + cdc[2,2])


# Reformat and present results
output <- data.frame(
  Value = c(sn[[1]] / sn[[2]], sp[[1]] / sp[[2]] ,  
  ppv[[1]] / ppv[[2]], npv[[1]] / npv[[2]]),
 lowerCI = c(sn[[4]][1], sp[[4]][1], ppv[[4]][1], npv[[4]][1]),
  upperCI = c(sn[[4]][2], sp[[4]][2], ppv[[4]][2], npv[[4]][2]),
  Total = c(total, total, total,total),
  row.names=c("Sensitivity", "Specificity","PPV","NPV"))

round(output, 2)
}

# end


# function for printing the discriptive demographics data
show_data <- function(dataset, group){

print( table(group))
print ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

groupnum <- dim(table(group))

totalnum <- length(na.exclude(group))	#total number of culture available subjects
table(group, exclude=NULL)		#show age group

for (i in 1:groupnum){					
  
  which_group <- as.character (names(table(group))[i]) # get the entry name of each group
  groupdat <- dataset[group == which_group,]  # get the data from that group
  grouptotal <- sum(1*(na.exclude(group==which_group)))

  percentage <- sum(1*(na.exclude(group== which_group))) / totalnum  
  percentage <- round(percentage, digits = 2) *100
  print (paste ("group", which_group ," % :", percentage, "%")) #show % of each group
 
  print("culture result :")
  print (table(dataset$cultpos[group == which_group], exclude=NULL) )
  
  cultnegnum <- sum(1*(groupdat$cultpos ==0 & !is.na(groupdat$cultpos))) # number of culture positives
  cultposnum <- sum(1*(groupdat$cultpos ==1 & !is.na(groupdat$cultpos))) # number of culture negatives

  print (paste ("% of -ve culture: ",  round(cultnegnum / grouptotal, digits =2)*100, "%") )
  print (paste ("% of +ve culture: ",  round(cultposnum / grouptotal, digits =2)*100, "%" ) )

  # for quickvue result
  print("quickvue result :")
  print (table(dataset$QVpos[group == which_group], exclude=NULL) )
  
  QVnegnum <- sum(1*(groupdat$QVpos ==0 & !is.na(groupdat$QVpos))) # number of QV positives
  QVposnum <- sum(1*(groupdat$QVpos ==1 & !is.na(groupdat$QVpos))) # number of QV negatives

  print (paste ("% of -ve QV: ",  round(QVnegnum / grouptotal, digits =2)*100, "%") )
  print (paste ("% of +ve QV: ",  round(QVposnum / grouptotal, digits =2)*100, "%" ) )
 
 
  print (" ************************************** ")

}

}

#end


#
#function for display all the test statistics
#
teststat2by2 <- function(tab){
output <- vector(length=12)
names(output) <- c("sn","sp","lrpos","lrneg","ppv","npv", "prevalence", 
"pretest odd", "post test odd pos", "post pro pos", "post test odd neg", "post test neg") 

  cdc <-matrix(NA, ncol=2, nrow=2)
  cdc[1,1]<-tab[2,2]
  cdc[1,2]<-tab[2,1]
  cdc[2,1]<-tab[1,2]
  cdc[2,2]<-tab[1,1]

# sensitivity
sn <- cdc[1,1] / (cdc[1,1] + cdc[2,1])
# specificity
sp <- cdc[2,2] / (cdc[1,2] + cdc[2,2])
# likelihood ratio for a positive test result
lrpos <- sn/(1-sp)
# likelihood ratio for a negative test result
lrneg <- (1-sn)/sp
#ppv
ppv <- cdc[1,1] / (cdc[1,1] + cdc[1,2])
#npv
npv <- cdc[2,2] / (cdc[2,1] + cdc[2,2])

#prevalence 
prevalence <- (cdc[1,1] + cdc[2,1]) / (cdc[1,1] + cdc[1,2] + cdc[2,1] + cdc[2,2])

#pre-test odds
pretest_odd <- prevalence/(1-prevalence)
#post-test odds for positive test
post_odd_pos <- pretest_odd * lrpos
#post-test odds for positive test
post_pro_pos <- post_odd_pos / (post_odd_pos +1)

#post-test odds for negative test
post_odd_neg <- pretest_odd * lrneg
#post-test odds for negative test
post_pro_neg <- post_odd_neg / (post_odd_neg +1)

output[1] <- sn
output[2] <- sp 
output[3] <- lrpos
output[4] <- lrneg
output[5] <- ppv
output[6] <- npv
output[7] <- prevalence
output[8] <- pretest_odd
output[9] <- post_odd_pos 
output[10] <- post_pro_pos
output[11] <- post_odd_neg
output[12] <- post_pro_neg

return(round(output,3))
}


teststat <- function(sn, sp, prevalence){

# likelihood ratio for a positive test result
lrpos <- sn/(1-sp)
# likelihood ratio for a negative test result
lrneg <- (1-sn)/sp
#pre-test odds
pretest_odd <- prevalence/(1-prevalence)
#post-test odds for positive test
post_odd_pos <- pretest_odd * lrpos
#post-test odds for positive test
post_pro_pos <- post_odd_pos / (post_odd_pos +1)

#post-test odds for negative test
post_odd_neg <- pretest_odd * lrneg
#post-test odds for negative test
post_pro_neg <- post_odd_neg / (post_odd_neg +1)

print("post test probablility for +ve test: ")
print(post_pro_pos)

print("post test probablility for -ve test: ")
print(post_pro_neg)

}
  

#totalnum <- length(na.exclude(ca$agegroup))	#total number of culture available subjects
#table(ca$agegroup, exclude=NULL)		#show age group

#for (i in 1:4){					
#  groupdat <- ca[ca$agegroup ==i,]		# get the data from that group
#  grouptotal <- sum(1*(na.exclude(ca$agegroup==i)))

#  percentage <- sum(1*(na.exclude(ca$agegroup==i))) / totalnum  
#  percentage <- round(percentage, digits = 2) *100
#  print (paste ("group", i ," % :", percentage, "%")) #show % of each group
 
#  print("culture result :")
# print (table(ca$cultpos[ca$agegroup ==i], exclude=NULL) )
  
   # cultposnum <- sum(1*(groupdat$cultpos ==0 & !is.na(groupdat$cultpos))) # number of culture positives
  #cultnegnum <- sum(1*(groupdat$cultpos ==1 & !is.na(groupdat$cultpos))) # number of culture negatives

  #print (paste ("% of -ve culture: ",  round(cultposnum / grouptotal, digits =2)*100, "%") )
  #print (paste ("% of +ve culture: ",  round(cultnegnum / grouptotal, digits =2)*100, "%" ) )
 

  #print (" ************************************** ")

#}



# 071121 function for return the discriptive demographics data
return_data <- function(dataset, group){

print( table(group))
print ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

dat <- data.frame()

groupnum <- dim(table(group))

totalnum <- length(na.exclude(group))	#total number of culture available subjects
table(group, exclude=NULL)		#show age group

for (i in 1:groupnum){					
  
  which_group <- as.character (names(table(group))[i]) # get the entry name of each group
  dat[i,1] <- which_group

  groupdat <- dataset[group == which_group,]  # get the data from that group
  
  grouptotal <- sum(1*(na.exclude(group==which_group)))
  dat[i,2] <-grouptotal

  percentage <- sum(1*(na.exclude(group== which_group))) / totalnum  
  percentage <- round(percentage, digits = 2) *100
  dat[i,3] <- percentage
  
  print (paste ("group", which_group ," % :", percentage, "%")) #show % of each group
 
  print("culture result :")
  print (table(dataset$cultpos[group == which_group], exclude=NULL) )
  
  cultnegnum <- sum(1*(groupdat$cultpos ==0 & !is.na(groupdat$cultpos))) # number of culture positives
  dat[i,4] <-cultnegnum
  
  cultposnum <- sum(1*(groupdat$cultpos ==1 & !is.na(groupdat$cultpos))) # number of culture negatives
  dat[i,5] <-cultposnum

  print (paste ("% of -ve culture: ",  round(cultnegnum / grouptotal, digits =2)*100, "%") )
  print (paste ("% of +ve culture: ",  round(cultposnum / grouptotal, digits =2)*100, "%" ) )
  dat[i,6] <- round(cultnegnum / grouptotal, digits =2)*100
  dat[i,7] <- round(cultposnum / grouptotal, digits =2)*100

  # for quickvue result
  print("quickvue result :")
  print (table(dataset$QVpos[group == which_group], exclude=NULL) )
  
  QVnegnum <- sum(1*(groupdat$QVpos ==0 & !is.na(groupdat$QVpos))) # number of QV positives
  dat[i,8] <- QVnegnum
  
  QVposnum <- sum(1*(groupdat$QVpos ==1 & !is.na(groupdat$QVpos))) # number of QV negatives
  dat[i,9] <- QVposnum

  print (paste ("% of -ve QV: ",  round(QVnegnum / grouptotal, digits =2)*100, "%") )
  print (paste ("% of +ve QV: ",  round(QVposnum / grouptotal, digits =2)*100, "%" ) )
  dat[i,10] <- round(QVnegnum / grouptotal, digits =2)*100
  dat[i,11] <- round(QVposnum / grouptotal, digits =2)*100
 
  print (" ************************************** ")

}
names(dat) <- c("group,","subject num,","subject %,",
"cult-ve no,","cult+ve no,","cult-ve %,","cult+ve %,",
"QV-ve no,","QV+ve no,","QV-ve %,","QV+ve %,")
return (dat[c(1:3,5,7,9,11)])
}
