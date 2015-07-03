## -- This file assign new categorical variables to create Table 1 of serosurvey -- ##
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
# [1] "ind_id"          "hh_id"           "hh_index"        "hh_size"         "base_attendance"
# [6] "base_blood"      "fu_attendance"   "fu_blood"        "sex"             "relation"       
#[11] "base_dob"        "base_hkid"       "base_ind_date"   "fu_ind_date"     "district"       
#[16] "education"       "profession"      "profession_text" "indoor"          "bldg"           
#[21] "bldg_text"       "coworker_num"    "school"          "classmates_num"  "eversmoked"     
#[26] "smoking"         "quit_y"          "quit_m"          "smoke_num"       "chronic"        
#[31] "adultdis"        "adultdis_text"   "childdis"        "childdis_text"   "westmed"        
#[36] "suppl"           "chinmed"         "med_text"        "vac0809"         "vac0708"        
#[41] "vac0607"         "Coldnum"         "Cold_date1"      "Cold_date2"      "Cold_date3"     
#[46] "cold1_symptoms"  "cold2_symptoms"  "cold3_symptoms"  "remark"          "lab_tested"     
#[51] "lab_positive"    "ag"              "hhsize"          "hhsgt2"          "pop_recruit"    
#[56] "recruit_source"  "hh_recruit"      "hh_telsource"    "tested1"         "tested2"        
#[61] "final_tested"    "hh_tested"       "test1_fourfold"  "test2_fourfold"  "final_fourfold" 
#[66] "b20"             "f40"             "b40"             "f20"             "age"            
#[71] "ag1"             "sy_ph_crude"     "sy_ph_tpu"       "child_hh"        "blood.pair"     
#[76] "symp.reported"   "symp.fever.rp"   "symp.ILI.rp"     "symp.ARI.rp"     "diary.reported" 
#[81] "diary.fever.rp"  "diary.ILI.rp"    "diary.ARI.rp"    "quest.reported"  "quest.fever"    
#[86] "quest.ILI"       "quest.ARI"       "any.fever"       "any.ILI"         "any.ARI"        
#[91] "district5" 


# Variables to be tested
# ag1[Y], phone.ili[N], phone.ari[N], diary.ili[N], diary.ari[N], symp.ili[N], symp.ari[N]

#Initalize variables
norowbaseline <- dim(base)[1]
agegroup <- array(888,c(1,norowbaseline))
age.hh <- array(888,c(1,norowbaseline))
age.pair <- array(888,c(1,norowbaseline))
age.pos <- array(888,c(1,norowbaseline))
sexgroup <- as.factor(base$sex)
edugroup <- array(888, c(1,norowbaseline)) 
edu.hh <- array(888,c(1,norowbaseline))
edu.pair <- array(888,c(1,norowbaseline))
edu.pos <- array(888,c(1,norowbaseline))
jobgroup <- array(888,c(1,norowbaseline)) 
job.hh <- array(888,c(1,norowbaseline))
job.pair <- array(888,c(1,norowbaseline))
job.pos <- array(888,c(1,norowbaseline)) 
famgroup <- array(888,c(1,norowbaseline)) 
fam.hh <- array(888,c(1,norowbaseline))
fam.pair <- array(888,c(1,norowbaseline))
fam.pos <- array(888,c(1,norowbaseline))
base <- cbind(base,district5=rep(-1))
district.hh <- array(888,c(1,norowbaseline))
district.pair <- array(888,c(1,norowbaseline))
district.pos <- array(888,c(1,norowbaseline))


for (i in 1:norowbaseline) {
        # Define age classes
        tmp <- base$age[i]
        if (is.na(tmp) || tmp < 0) agegroup[i] <- "NK" else
        if (tmp < 20) agegroup[i] <- 1 else
        if (tmp < 40) agegroup[i] <- 2 else
        if (tmp < 60) agegroup[i] <- 3 else
        if (tmp < 110) agegroup[i] <- 4 else
        stop("Problem with age allocation")
        
        age.hh <- ifelse (base$hh_tested==TRUE, agegroup, "NA") #Age of individuals whose at least one member in household provided blood samples
        age.pair <- ifelse (base$final_tested==TRUE, agegroup, "NA") #Age of individuals based on whether or not they had provided blood samples
        age.pos <- ifelse (base$final_fourfold==1, agegroup, "NA") #Age of individuals based on whether or not they had provided blood samples
        
        
        # Define sex groups
        sex.hh <- ifelse(base$hh_tested==TRUE, base$sex, 999)
        sex.pair <- ifelse(base$final_tested==TRUE, base$sex, 999)
        sex.pos <- ifelse(base$final_fourfold==1, base$sex,999)
        
        
        # Define education groups
        tmp <- base$education[i]
        if (is.na(tmp) || tmp == 0) edugroup[i] <- "NK" else
        if (tmp == 1) edugroup[i] <- 1 else
        if (tmp == 2) edugroup[i] <- 2 else
        if (tmp == 3 || tmp == 4 || tmp == 5) edugroup[i] <- 3 else
        if (tmp == 6) edugroup[i] <- 4 else
        stop("Problem with education allocation")
        
        edu.hh <- ifelse(base$hh_tested==TRUE, edugroup, 999)
        edu.pair <- ifelse (base$final_tested==TRUE, edugroup, 999) 
        edu.pos <- ifelse (base$final_fourfold==1, edugroup, 999)     
        
        
        #Define occupation class
        tmp <- base$profession[i]
        if (is.na(tmp) || tmp == 0) jobgroup[i] <- "NK" else
        if (tmp == 1 || tmp == 2) jobgroup[i] <- 1 else
        if (tmp == 3 || tmp == 4) jobgroup[i] <- 2 else
        if (tmp == 5 || tmp == 6) jobgroup[i] <- 3 else
        if (tmp == 7 || tmp == 8 || tmp == 9) jobgroup[i] <- 4 else 
        if (tmp == 10) jobgroup[i] <- 5 else #reassign housewives (10) to group (5)
        if (tmp == 11) jobgroup[i] <- 5 else #reassign retirees (11) to group (5)
        if (tmp == 12) jobgroup[i] <- 6 else #reassign students (12) to group (6)
        if (tmp == 13) jobgroup[i] <- 5 else #reassign unemployed (13) to group (5)
        if (tmp == 14) jobgroup[i] <- 5 else #reassign housekeepers (14) to group (5)
        stop("Problem with occupation allocation")
                
        job.hh <- ifelse (base$hh_tested==TRUE, jobgroup, 999) 
        job.pair <- ifelse (base$final_tested==TRUE, jobgroup, 999)         
        job.pos <- ifelse (base$final_fourfold==1, jobgroup, 999) 

        
        #Define family size class
        tmp <- base$hh_size[i]
        if (is.na(tmp) | tmp == 0) famgroup[i] <- "NA" else
        if (tmp == 1) famgroup[i] <- 1 else
        if (tmp == 2) famgroup[i] <- 2 else
        if (tmp == 3) famgroup[i] <- 3 else
        if (tmp == 4) famgroup[i] <- 4 else
        if (tmp == 5) famgroup[i] <- 5 else
        if (tmp == 6) famgroup[i] <- 6 else
        if (tmp == 7) famgroup[i] <- 6 else
        if (tmp == 8) famgroup[i] <- 6 else
        if (tmp == 9) famgroup[i] <- 6 else
        stop("Problem with family size allocation")
        
        fam.hh <- ifelse (base$hh_tested==TRUE, famgroup, 999) 
        fam.pair <- ifelse (base$final_tested==TRUE, famgroup, 999) 
        fam.pos <- ifelse (base$final_fourfold==1, famgroup, 999) 
        
        #Recategorize residential districts
        tmp <- as.character(base$district[i])
        if (is.na(tmp)) base$district5[i] <- "NK" else
        if (tmp=="#N/A") base$district5[i] <- "NK" else
        if (tmp=="ABERDEEN") base$district5[i] <- "HK Island" else
        if (tmp=="CENTRAL AND WESTERN" || tmp=="CENTRAL AND WESTERN ") base$district5[i] <- "HK Island" else
        if (tmp=="EASTERN") base$district5[i] <- "HK Island" else
        if (tmp=="ISLAND") base$district5[i] <- "NT West" else
        if (tmp=="KOWLOON CITY" || tmp=="KOWLOON CITY ") base$district5[i] <- "KLN West" else
        if (tmp=="KWAI TSING" || tmp=="KWAI TSING ") base$district5[i] <- "NT West" else
        if (tmp=="KWUN TONG") base$district5[i] <- "KLN East" else
        if (tmp=="NORTH") base$district5[i] <- "NT East" else
        if (tmp=="SAI KUNG") base$district5[i] <- "NT East" else
        if (tmp=="EASTERN / SAI KUNG") base$district5[i] <- "NT East" else
        if (tmp=="SHA TIN" || tmp=="SHA TIN ") base$district5[i] <- "NT East" else
        if (tmp=="SHAM SHUI PO") base$district5[i] <- "KLN West" else
        if (tmp=="SOUTHERN") base$district5[i] <- "HK Island" else
        if (tmp=="TAI PO") base$district5[i] <- "NT East" else
        if (tmp=="TUEN MUN") base$district5[i] <- "NT West" else
        if (tmp=="TSUEN WAN") base$district5[i] <- "NT West" else
        if (tmp=="WAN CHAI") base$district5[i] <- "HK Island" else
        if (tmp=="WONG TAI SIN") base$district5[i] <- "KLN East" else
        if (tmp=="YAU TSIM MONG") base$district5[i] <- "KLN West" else
        if (tmp=="YUEN LONG") base$district5[i] <- "NT West" else
        stop("Problem with residential district allocation")
        
        district.hh <- ifelse (base$hh_tested==TRUE, base$district5, "NA")
        district.pair <- ifelse (base$final_tested==TRUE, base$district5, "NA") 
        district.pos <- ifelse (base$final_fourfold==1, base$district5, "NA") 
}        

agegroup <- as.factor(agegroup)
levels(agegroup)[1] <- "0-19yrs"
levels(agegroup)[2] <- "20-39yrs"
levels(agegroup)[3] <- "40-59yrs"
levels(agegroup)[4] <- "60+ yrs"

age.hh <- as.factor(age.hh)
levels(age.hh)[1] <- "0-19yrs"
levels(age.hh)[2] <- "20-39yrs"
levels(age.hh)[3] <- "40-59yrs"
levels(age.hh)[4] <- "60+ yrs"

age.pair <- as.factor(age.pair)
levels(age.pair)[1] <- "0-19yrs"
levels(age.pair)[2] <- "20-39yrs"
levels(age.pair)[3] <- "40-59yrs"
levels(age.pair)[4] <- "60+ yrs"

age.pos <- as.factor(age.pos)
levels(age.pos)[1] <- "0-19yrs"
levels(age.pos)[2] <- "20-39yrs"
levels(age.pos)[3] <- "40-59yrs"
levels(age.pos)[4] <- "60+ yrs"

sex.hh <- as.factor(sex.hh)
sex.pair <- as.factor(sex.pair)
sex.pos <- as.factor(sex.pos)

edugroup <- as.factor(edugroup)
levels(edugroup)[1] <- "Kindergarten"
levels(edugroup)[2] <- "Primary"
levels(edugroup)[3] <- "Secondary"
levels(edugroup)[4] <- "Degrees"

edu.hh <- as.factor(edu.hh)
levels(edu.hh)[1] <- "Kindergarten"
levels(edu.hh)[2] <- "Primary"
levels(edu.hh)[3] <- "Secondary"
levels(edu.hh)[4] <- "Degrees"

edu.pair <- as.factor(edu.pair)
levels(edu.pair)[1] <- "Kindergarten"
levels(edu.pair)[2] <- "Primary"
levels(edu.pair)[3] <- "Secondary"
levels(edu.pair)[4] <- "Degrees"

edu.pos <- as.factor(edu.pos)
levels(edu.pos)[1] <- "Kindergarten"
levels(edu.pos)[2] <- "Primary"
levels(edu.pos)[3] <- "Secondary"
levels(edu.pos)[4] <- "Degrees"

jobgroup <- as.factor(jobgroup)
#levels(jobgroup)[1] <- "Managers, administators, professionals" #jobgroup = 1
#levels(jobgroup)[2] <- "Associate professionals, clerks" #jobgroup = 2
#levels(jobgroup)[3] <- "Service workers & shop sales workers/ Craft & related workers" #jobgroup = 3
#levels(jobgroup)[4] <- "Plant & machine operators & assemblers/Elementary occupations/Skilled agriculural & fishery workers" #jobgroup = 4
#levels(jobgroup)[5] <- "Retirees/ Unemployed/ home-makers" #jobgroup = 5
#levels(jobgroup)[6] <- "Students" #jobgroup = 6

job.hh <- as.factor(job.hh)
job.pair <- as.factor(job.pair)
job.pos <- as.factor(job.pos)

famgroup <- as.factor(famgroup)
levels(famgroup)[1] <- "1"
levels(famgroup)[2] <- "2"
levels(famgroup)[3] <- "3"
levels(famgroup)[4] <- "4"
levels(famgroup)[5] <- "5"
levels(famgroup)[6] <- "6 or above"

fam.hh <- as.factor(fam.hh)
levels(fam.hh)[1] <- "1"
levels(fam.hh)[2] <- "2"
levels(fam.hh)[3] <- "3"
levels(fam.hh)[4] <- "4"
levels(fam.hh)[5] <- "5"
levels(fam.hh)[6] <- "6 or above"

fam.pair <- as.factor(fam.pair)
levels(fam.pair)[1] <- "1"
levels(fam.pair)[2] <- "2"
levels(fam.pair)[3] <- "3"
levels(fam.pair)[4] <- "4"
levels(fam.pair)[5] <- "5"
levels(fam.pair)[6] <- "6 or above"

fam.pos <- as.factor(fam.pos)

district.hh <- as.factor(district.hh)
district.pair <- as.factor(district.pair)
district.pos <- as.factor(district.pos)

print('Age of Cohort:')
print(table(agegroup))
print('Age of those whose family members provided paired samples:')
print(ftable(age.hh,base$hh_telsource))
print('Age of those with paired samples:')
print(table(age.pair))
print('Age of those who tested positive:')
print(table(age.pos))
print('Sex of Cohort:')
print(table(sexgroup))
print('Sex of those whose family members provided paired samples:')
print(ftable(sex.hh,base$hh_telsource))
print('Sex of those with paired samples:')
print(table(sex.pair))
print('Sex of those who tested positive:')
print(table(sex.pos))
print('Education level of Cohort:')
print(table(edugroup))
print('Education level of those family members have provided paired samples:')
print(ftable(edu.hh,base$hh_telsource))
print('Education level of those who have provided paired samples:')
print(table(edu.pair))
print('Education level of those who are sero-converted:')
print(table(edu.pos))
print('Occupation of Cohort:')
print(table(jobgroup))
print('Occupation of those whose family members have provided paired samples:')
print(ftable(job.hh,base$hh_telsource))
print('Occupation of those who have provided paired samples:')
print(table(job.pair))
print('Occupation of those who are sero-converted:')
print(table(job.pos))
print('Family size of Cohort:')
print(table(famgroup))
print('Family size of those whose family members have provided paired samples:')
print(ftable(fam.hh,base$hh_telsource))
print('Family size of those who have provided paired samples:')
print(table(fam.pair))
print('Family size of those who are sero-converted:')
print(table(fam.pos))
print('Residential district of Cohort:')
print(table(base$district5))
print('Residentail district of those whose family members have provided paired samples:')
print(ftable(district.hh,base$hh_telsource))
print('Residential district of those who have provided paired samples:')
print(table(district.pair))
print('Residential district of those who are sero-converted:')
print(table(district.pos))

# Save subdatasets
write.csv(base, file="M:\\post_proc_data\\base.csv")
t_base <- base[base[,"final_fourfold"] > -1,]
write.csv(t_base, file="M:\\post_proc_data\\t_base.csv")
