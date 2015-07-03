# Plot attack rates by age, education levels, and family size
rm(list=ls(all=TRUE))
options(error=NULL)

#source("funcs_hksero.r")
source("scrap.R")

base_raw <- read.csv("../anon_data/base_ind.csv")
labs1_raw <- read.csv("../anon_data/lab_results1.csv")
labs2_raw <- read.csv("../anon_data/lab_results2.csv")
symp_raw <- read.csv("../anon_data/phone_symp.csv")
diary_raw <- read.csv("../anon_data/diary.csv")
quest_raw <- read.csv("../anon_data/questionnaire_symptoms.csv")
recruit_raw <- read.csv("../anon_data/subj_recruit.csv")

pp_data <- post_proc_sero(base_raw, labs1_raw, labs2_raw, symp_raw, diary_raw, quest_raw, recruit_raw)

base <- pp_data$base
labs1 <- pp_data$labs1
labs2 <- pp_data$labs2
symp <- pp_data$symp
diary <- pp_data$diary
quest <- pp_data$quest
recruit <- pp_data$recruit

# names(base)
# [1] "ind_id"          "hh_id"           "hh_index"        "hh_size"         "base_attendance" "base_blood"     
# [7] "fu_attendance"   "fu_blood"        "sex"             "relation"        "base_dob"        "base_hkid"      
#[13] "base_ind_date"   "fu_ind_date"     "district"        "education"       "profession"      "profession_text"
#[19] "indoor"          "bldg"            "bldg_text"       "coworker_num"    "school"          "classmates_num" 
#[25] "eversmoked"      "smoking"         "quit_y"          "quit_m"          "smoke_num"       "chronic"        
#[31] "adultdis"        "adultdis_text"   "childdis"        "childdis_text"   "westmed"         "suppl"          
#[37] "chinmed"         "med_text"        "vac0809"         "vac0708"         "vac0607"         "Coldnum"        
#[43] "Cold_date1"      "Cold_date2"      "Cold_date3"      "cold1_symptoms"  "cold2_symptoms"  "cold3_symptoms" 
#[49] "remark"          "lab_tested"      "lab_positive"    "ag"              "hhsize"          "hhsgt2"         
#[55] "tested"          "fourfold"        "b20"             "f40"             "b40"             "f20"            
#[61] "age"             "ag1"             "sy_ph_crude"     "sy_ph_tpu"       "child_hh"        "blood.pair"     
#[67] "symp.reported"   "symp.fever.rp"   "symp.ILI.rp"     "symp.ARI.rp"     "diary.reported"  "diary.fever.rp" 
#[73] "diary.ILI.rp"    "diary.ARI.rp"    "quest.reported"  "quest.fever"     "quest.ILI"       "quest.ARI" 


# Make subselections
t_base <- base[base[,"fourfold"] > -1,]

#Initalize variables
norowtbase <- dim(t_base)[1]
batchsize <- 30
batchnum <- floor(norowtbase/batchsize) 
batchnum.age <- 4 #Number of age groups
age.matrix <- matrix(999, nrow=batchsize, ncol=batchnum)
meanage <- array(999, c(1,batchnum))
attackrate.age2 <- array(999, c(1,batchnum))

batchnum.edu <- 7 #Number of education groups
edu.matrix <- matrix(999, nrow=batchsize, ncol=batchnum)
meanedu <- array(999, c(1,batchnum))
attackrate.edu2 <- array(999, c(1,batchnum))

batchnum.famsize <- 6 #Number of family size groups
famsize.matrix <- matrix(999, nrow=batchsize, ncol=batchnum)
meanfamsize <- array(999, c(1,batchnum))
attackrate.famsize2 <- array(999, c(1,batchnum))


# Compute overall attack rates
print('Individuals who have been infected (fourfold) versus those who have not:')
print(table(t_base$fourfold))
attack.rate <- as.vector(table(t_base$fourfold))[2]/norowtbase
print(noquote(paste("Overall attack rate (%): ", attack.rate*100)))

## -- Tables --##
print('No. of people with four-fold rise by age groups:')
print(table(t_base$ag1, t_base$fourfold))
print('No. of people with four-fold rise by education levels:')
print(table(t_base$education, t_base$fourfold))
print('No. of people with four-fold rise by family size:')
print(table(t_base$hh_size, t_base$fourfold))
print('No. of people with four-fold rise by vaccination status:')
print(table(t_base$vac0809, t_base$fourfold))
#Among the young adults (19-64 yrs)
ya_base <- t_base[t_base[,"ag1"] == 2,]
print(table(ya_base$vac0809, ya_base$fourfold))

#Compute the main effect means
t_base$ag1 <- relevel(factor(t_base$ag1),ref="1")
print('Main effect means by age groups:')
tapply(t_base$fourfold, t_base$ag1, mean)
print('Main effect means by education levels:')
tapply(t_base$fourfold, t_base$education, mean)
print('Main effect means by family size:')
tapply(t_base$fourfold, t_base$hh_size, mean)
print('Main effect means by vaccination status:')
tapply(t_base$fourfold, t_base$vac0809, mean)


## --  Attack rate by ages -- ##
t_base.sortage <- t_base[order(t_base$age),] #Sort dataset by age

# Compute attack rates by age
for (i in 1:batchnum) {
    age.matrix[,i] <- t_base.sortage$age[(batchsize*(i-1)+1): (batchsize*i)] #Divide dataset by ages
    meanage[i] <- mean(age.matrix[,i])
    tmp <- table(t_base.sortage$fourfold[(batchsize*(i-1)+1): (batchsize*i)])[2] #Compute attack rates of different ages
    if (!is.na(tmp)) attackrate.age2[i] <- tmp/ batchsize
    else attackrate.age2[i] <- 0}
    
    
## -- Attack rates by education levels -- ##
t_base.sortedu <- t_base[order(t_base$education),] #Sort dataset by education

# Compute attack rates by education levels
for (i in 1:batchnum) {
    edu.matrix[,i] <- t_base.sortedu$education[(batchsize*(i-1)+1): (batchsize*i)]
    meanedu[i] <- mean(edu.matrix[,i])
    tmp <- table(t_base.sortedu$fourfold[(batchsize*(i-1)+1): (batchsize*i)])[2]
    if (!is.na(tmp)) attackrate.edu2[i] <- tmp/ batchsize
    else attackrate.edu2[i] <- 0
}


## -- Attack rates by family size -- ##
t_base.sortfamsize <- t_base[order(t_base$hh_size),] #Sort dataset by family size

# Compute attack rates by family size
for (i in 1:batchnum) {
    famsize.matrix[,i] <- t_base.sortfamsize$hh_size[(batchsize*(i-1)+1): (batchsize*i)] #Divide dataset by ages
    meanfamsize[i] <- mean(famsize.matrix[,i])
    tmp <- table(t_base.sortfamsize$fourfold[(batchsize*(i-1)+1): (batchsize*i)])[2] #Compute attack rates of different ages
    if (!is.na(tmp)) attackrate.famsize2[i] <- tmp/ batchsize
    else attackrate.famsize2[i] <- 0}
    

## -- Plots -- ##
#Variables for plots
windows()
op <- par(mfrow = c(1,2))

#Plot attack rates by age
title = "Four-fold rise by Age"
xlabel = "Age (years)"
ylabel = "Number of people with four-fold rise"
barplot(table(t_base$ag1), main=title, xlab=xlabel, ylab=ylabel, col="grey")
ylabel = "Proposition of people with four-fold rise"
plot(meanage, attackrate.age2, ylim=c(0,0.55), main=title, xlab=xlabel, ylab=ylabel, type="b")
grid(20, 10, lwd=2) # 20 grids in x-direction, 10 in y-direction

# Plot attack rate by education levels
windows()
op <- par(mfrow = c(1,2))
title = "Four-fold rise by Education levels"
xlabel = "Education levels"
ylabel = "Number of people with four-fold rise"
barplot(table(t_base$education), main=title, xlab=xlabel, ylab=ylabel, col="grey")
ylabel = "Proposition of people with four-fold rise"
plot(meanedu, attackrate.edu2, ylim=c(0,0.55), main=title, xlab=xlabel, ylab=ylabel, type="b")
grid(24, 14, lwd=2) # 24 grids in x-direction, 14 in y-direction

# Plot attack rate by family size
windows()
op <- par(mfrow = c(1,2))
title = "Four-fold rise by Family Size"
xlabel = "Family Size"
ylabel = "Number of people with four-fold rise"
barplot(table(t_base$hh_size), main=title, xlab=xlabel, ylab=ylabel, col="grey")
ylabel = "Proposition of people with four-fold rise"
plot(meanfamsize, attackrate.famsize2, ylim=c(0,0.55), main=title, xlab=xlabel, ylab=ylabel, type="b")
grid(24, 14, lwd=2) # 24 grids in x-direction, 14 in y-direction
