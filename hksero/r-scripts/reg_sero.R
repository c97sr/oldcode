rm(list=ls(all=TRUE))
options(error=NULL)

# To dos
# - make a variable for household contains child

# Change working directory if required
setwd("./hksero/")

source("./hksero/r-scripts/scrap.r")

#<<<<<<< .mine
base_raw <- read.csv("./hksero/anon_data/base_ind.csv")
labs_raw_old <- read.csv("./hksero/anon_data/lab_results1.csv")
labs_raw_new <- read.csv("./hksero/anon_data/lab_results2.csv")
symp_raw <- read.csv("./hksero/anon_data/phone_symp.csv")
diar_raw <-read.csv("./hksero/anon_data/diary.csv")
ques_raw <-read.csv( "./hksero/anon_data/questionnaire_symptoms.csv")
rec_raw  <-read.csv("./hksero/anon_data/subj_recruit.csv")
#=======
#base_raw <- read.csv("../anon_data/base_ind.csv")
#labs_raw_old <- read.csv("../anon_data/lab_results1.csv")
#labs_raw_new <- read.csv("../anon_data/lab_results2.csv")
#symp_raw <- read.csv("../anon_data/phone_symp.csv")
#diar_raw <-read.csv("../anon_data/diary.csv")
#ques_raw <-read.csv( "../anon_data/questionnaire_symptoms.csv")
#>>>>>>> .r69

# Comment out below as required
pp_data <- post_proc_sero(base_raw, labs_raw_old,labs_raw_new,symp_raw,diar_raw,ques_raw,rec_raw)
#tmpfile <- "tmp_reg_sero_pp_data.Rdata"
#save(pp_data,file=tmpfile)
#load("./tmp_reg_sero_pp_data.Rdata")

base <- pp_data$base
labs1 <- pp_data$labs1
labs2 <- pp_data$labs2
symp <- pp_data$symp
diary<- pp_data$diary
quest<- pp_data$quest
recruit <- pp_data$recruit

# names(base)
# [1] "ind_id"          "hh_id"           "hh_index"        "hh_size"        
# [5] "base_attendance" "base_blood"      "fu_attendance"   "fu_blood"       
# [9] "sex"             "relation"        "base_dob"        "base_hkid"      
#[13] "base_ind_date"   "fu_ind_date"     "district"        "education"      
#[17] "profession"      "profession_text" "indoor"          "bldg"           
#[21] "bldg_text"       "coworker_num"    "school"          "classmates_num" 
#[25] "eversmoked"      "smoking"         "quit_y"          "quit_m"         
#[29] "smoke_num"       "chronic"         "adultdis"        "adultdis_text"  
#[33] "childdis"        "childdis_text"   "westmed"         "suppl"          
#[37] "chinmed"         "med_text"        "vac0809"         "vac0708"        
#[41] "vac0607"         "Coldnum"         "Cold_date1"      "Cold_date2"     
#[45] "Cold_date3"      "cold1_symptoms"  "cold2_symptoms"  "cold3_symptoms" 
#[49] "remark"          "lab_tested"      "lab_positive"    "ag"             
#[53] "hhsize"          "hhsgt2"          "pop_recruit"     "recruit_source" 
#[57] "hh_recruit"      "hh_telsource"    "tested1"         "tested2"        
#[61] "final_tested"    "hh_tested"       "test1_fourfold"  "test2_fourfold" 
#[65] "final_fourfold"  "b20"             "f40"             "b40"            
#[69] "f20"             "age"             "ag1"             "sy_ph_crude"    
#[73] "sy_ph_tpu"       "child_hh"        "blood.pair"      "symp.reported"  
#[77] "symp.fever.rp"   "symp.ILI.rp"     "symp.ARI.rp"     "diary.reported" 
#[81] "diary.fever.rp"  "diary.ILI.rp"    "diary.ARI.rp"    "quest.reported" 
#[85] "quest.fever"     "quest.ILI"       "quest.ARI"       "any.fever"      
#[89] "any.ILI"         "any.ARI"  
# Variables to be tested
# ag1[Y], hh_size[N], sex[N], education[N], vaccination0809[N],smoking[N], presence_children

# Make subselections
# tested

# regroup district
base  <-cbind(base,district5=rep(-1,dim(base)[1]))

for (i in 1:dim(base)[1]){

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

}




#source("./hksero/r-scripts/subj_characteristics.r")

# exclude three subject, two missing education level, one missing vaccination status
t_base_1 <- base[base[,"final_fourfold" ] > -1,]
t_base_2 <- t_base_1[t_base_1[,"education"]!=0,]
t_base   <-t_base_2[t_base_2[,"ind_id"]!="S090189-1",]


t_base_adults <- t_base[t_base[,"ag20"]!=1,]
t_base_children <- t_base[t_base[,"ag20"]==1,]
t_base_adult2 <- t_base[t_base[,"ag20"]==2,]
t_base_adult3 <- t_base[t_base[,"ag20"]==3,]
t_base_adult4 <- t_base[t_base[,"ag20"]==4,]

#t_base_source <- t_base[t_base[,"hh_telsource" ] > -1,]

t_base_yadults <- t_base[t_base[,"ag20"] %in% c(2,3),]



# not tested
#n_base <- base[base[,"blood.pair"] == 1,]
#m_base <- n_base[n_base[,"tested"]==FALSE,]
#m_base_adults <-m_base[m_base[,"ag1"]!=1,]
#m_base_yadults <- m_base_adults[m_base_adults[,"ag1"] %in% c(2,3),]
#ftable(m_base_adults$child_hh,m_base_adults$ag1)



# Make some preliminary comparisons
#table(t_base$sy_ph_tpu,t_base$final_fourfold)

# First, age groups
table_age <-table(t_base$ag20,t_base$final_fourfold)
cbind(table_age,percent=table_age[,2]/(table_age[,2]+table_age[,1]))



t_base$ag20 <- relevel(factor(t_base$ag20),ref="1")
model1 <- glm(formula = t_base$final_fourfold ~ factor(t_base$ag20), family = binomial(logit))
summary(model1)


# OR
round(exp(model1$coefficients),digits=4)

# 95% CI and p Values
coef(summary(model1))[2,][[4]]

round(exp(confint.default(model1)),digits=4)



# sex

table_sex <-table(t_base$sex,t_base$final_fourfold)
table_sex
cbind(table_sex,percent=table_sex[,2]/(table_sex[,2]+table_sex[,1]))

model_sex <- glm(formula = t_base$final_fourfold ~ factor(t_base$sex), family = binomial(logit))
summary(model_sex)

# OR
round(exp(model_sex$coefficients),digits=4)

# 95% CI and p Values
coef(summary(model_sex))[2,][[4]]

round(exp(confint.default(model_sex)),digits=4)



# presence of children in househould for all adults 

table_child<- ftable(t_base_adults$child_hh,t_base_adults$final_fourfold)
table_child
cbind(table_child,percent=table_child[,2]/(table_child[,2]+table_child[,1]))

#ftable(t_base_adults$child_hh,t_base_adults$ag1)


# refernce group is adults without children in the household, 2 stands for No children
t_base_adults$child_hh <- relevel(factor(t_base_adults$child_hh),ref="2")

model_child <- glm(formula = t_base_adults$final_fourfold ~ factor(t_base_adults$child_hh), family = binomial(logit))
summary(model_child)


# OR
round(exp(model_child$coefficients),digits=4)


# 95%CI and p Values
coef(summary(model_child))[2,][[4]]

round(exp(confint.default(model_child)),digits=4)




# smoking status
# use the following command to complete the data set, blank cells are "No",
# weak association(p<0.2)


for (i in 1:length(t_base$smoking)){
	if (t_base$smoking[i]!="Yes"){t_base$smoking[i] <-"No"}
}	

table_smoking <-table(t_base$smoking,t_base$final_fourfold)
table_smoking
cbind(table_smoking,percent=table_smoking[,2]/(table_smoking[,2]+table_smoking[,1]))


model_smoking <- glm(formula = t_base$final_fourfold ~ factor(t_base$smoking), family = binomial(logit))
summary(model_smoking)


# OR
round(exp(model_smoking$coefficients),digits=4)


# 95%CI and p Values
coef(summary(model_smoking))[2,][[4]]

round(exp(confint.default(model_smoking)),digits=4)




# household size
# reset houldhold groups
t_base  <-cbind(t_base,hh_size2=rep(-1,dim(t_base)[1]))

for (i in 1:dim(t_base)[1]){
        if (is.na(t_base[i,"hh_size"])) t_base$hh_size2[i] <- -1
        
        else if (t_base[i,"hh_size"] ==1) t_base[i,"hh_size2"] <- 1
        else if (t_base[i,"hh_size"] ==2) t_base[i,"hh_size2"] <- 2
        else if (t_base[i,"hh_size"] ==3) t_base[i,"hh_size2"] <- 3
        else if(t_base[i,"hh_size"] ==4) t_base[i,"hh_size2"] <- 4
        else if(t_base[i,"hh_size"] ==5) t_base[i,"hh_size2"] <- 5
        else if(t_base[i,"hh_size"] ==6) t_base[i,"hh_size2"] <- 6
        else if(t_base[i,"hh_size"] ==7) t_base[i,"hh_size2"] <- 6
        else if(t_base[i,"hh_size"] ==8) t_base[i,"hh_size2"] <- 6
        else if(t_base[i,"hh_size"] ==9) t_base[i,"hh_size2"] <- 6

        else stop("Problem with the age groups")
 }       
    
table_hh <-table(t_base$hh_size2,t_base$final_fourfold)
table_hh
cbind(table_hh,percent=table_hh[,2]/(table_hh[,2]+table_hh[,1]))


t_base$hh_size2 <- relevel(factor(t_base$hh_size2),ref="2")

model_hh <- glm(formula = t_base$final_fourfold ~ factor(t_base$hh_size2), family = binomial(logit))
summary(model_hh)


# OR
round(exp(model_hh$coefficients),digits=4)

# 95%CI and p Values
coef(summary(model_hh))[2,][[4]]

round(exp(confint.default(model_hh)),digits=4)


# education and profession are not associated with infection at significant level 0.2
# education level
# two missing education level

# regroup education level
t_base  <-cbind(t_base,education2=rep(-1,dim(t_base)[1]))

for (i in 1:dim(t_base)[1]){
        #if (is.na(t_base[i,"education"])) t_base$education2[i] <- -1
        if (t_base[i,"education"] ==0) t_base[i,"education2"] <- 0
        else if (t_base[i,"education"] ==1) t_base[i,"education2"] <- 1
        else if (t_base[i,"education"] ==2) t_base[i,"education2"] <- 2
        else if (t_base[i,"education"] ==3) t_base[i,"education2"] <- 3
        else if(t_base[i,"education"] ==4) t_base[i,"education2"] <- 3
        else if(t_base[i,"education"] ==5) t_base[i,"education2"] <- 3
        else if(t_base[i,"education"] ==6) t_base[i,"education2"] <- 4
        
        #else stop("Problem with the grouping educaiton level")
 } 
 
 
 
 
table_edu <-table(t_base$education2,t_base$final_fourfold)
table_edu
cbind(table_edu,percent=table_edu[,2]/(table_edu[,2]+table_edu[,1]))


t_base$education2 <- relevel(factor(t_base$education2),ref="1")

model_edu <- glm(formula = t_base$final_fourfold ~ factor(t_base$education2), family = binomial(logit))
summary(model_edu)

# OR
round(exp(model_edu$coefficients),digits=4)

# 95%CI and p Values
coef(summary(model_edu))[2,][[4]]

round(exp(confint.default(model_edu)),digits=4)





# profession

# regroup professon categories
t_base  <-cbind(t_base,profession2=rep(-1,dim(t_base)[1]))

for (i in 1:dim(t_base)[1]){
       
        if (t_base[i,"profession"] %in% c(1,2)) t_base[i,"profession2"] <- 1
        else if (t_base[i,"profession"] %in% c(3,4)) t_base[i,"profession2"] <- 2
        else if (t_base[i,"profession"] %in% c(5,6)) t_base[i,"profession2"] <- 3
        else if (t_base[i,"profession"] %in% c(7,8,9)) t_base[i,"profession2"] <-4
        else if(t_base[i,"profession"] %in% c(10,11,13,14)) t_base[i,"profession2"] <- 5
        else if(t_base[i,"profession"] %in% c(12)) t_base[i,"profession2"] <- 6
                
 } 
 
table_prof <-table(t_base$profession2,t_base$final_fourfold)
table_prof
cbind(table_prof,percent=table_prof[,2]/(table_prof[,2]+table_prof[,1]))

t_base$profession2 <- relevel(factor(t_base$profession2),ref="1")

model_prof <- glm(formula = t_base$final_fourfold ~ factor(t_base$profession2), family = binomial(logit))
summary(model_prof)

# OR
round(exp(model_prof$coefficients),digits=4)

# 95%CI and p Values
coef(summary(model_prof))[2,][[4]]

round(exp(confint.default(model_prof)),digits=4)

# vaccinated with seasonal influenza vaccine during year 2008 and 2009 for all subjects, not adjusted for other variables
for (i in 1:length(t_base$vac0809)){
	if (t_base$vac0809[i]=="yes"){t_base$vac0809[i] <-"Yes"}
	if (t_base$vac0809[i]=="no"){t_base$vac0809[i] <-"No"}
}	
 
table_vac0809 <-ftable(t_base$vac0809,t_base$final_fourfold)
cbind(table_vac0809,percent=table_vac0809[,2]/(table_vac0809[,2]+table_vac0809[,1]))

t_base$vac0809 <- relevel(factor(t_base$vac0809),ref="No")
model_vac0809 <- glm(formula = t_base$final_fourfold ~ t_base$vac0809, family = binomial(logit))
summary(model_vac0809)

# OR
round(exp(model_vac0809$coefficients),digits=2)[2]

# 95% CI and p Values
coef(summary(model_vac0809))[2,][[4]]
round(exp(confint.default(model_vac0809)),digits=4)

# vaccinated with seasonal flu vaccine during year 2008 and 2009 for children(age group 1)
table_vac0809 <-ftable(t_base$vac0809,t_base$final_fourfold,t_base$ag20)
cbind(table_vac0809,percent=table_vac0809[,2]/(table_vac0809[,2]+table_vac0809[,1]))

t_base_children$vac0809 <- relevel(factor(t_base_children$vac0809),ref="No")
model_children <- glm(formula = t_base_children$final_fourfold ~ t_base_children$vac0809, family = binomial(logit))
summary(model_children)

# OR
round(exp(model_children$coefficients),digits=4)

# 95%CI and p Values
P_vac <-coef(summary(model_children))[3,][[4]]
CI_vac<-cbind(exp(coef(model_children)[[3]]-1.96*coef(summary(model_children))[3,][[2]]),exp(coef(model_children)[[3]]+1.96*coef(summary(model_children))[3,][[2]]),P_vac)
colnames(CI_vac) <-c("Lower CI","Upper CI","p value")
round(CI_vac,4)

# vaccinated with seasonal flu vaccine during year 2008 and 2009 for age group 2
t_base_adult2$vac0809 <- relevel(factor(t_base_adult2$vac0809),ref="No")
model4 <- glm(formula = t_base_adult2$final_fourfold ~ factor(t_base_adult2$vac0809), family = binomial(logit))
summary(model4)

# OR
round(exp(model4$coefficients),digits=4)

# 95% CI and p Values
coef(summary(model4))[2,][[4]]

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)
round(exp(confint.default(model4)),digits=4)


# vaccinated with seasonal flu vaccine during year 2008 and 2009 for age group 3
t_base_adult3$vac0809 <- relevel(factor(t_base_adult3$vac0809),ref="No")
model4 <- glm(formula = t_base_adult3$final_fourfold ~ factor(t_base_adult3$vac0809), family = binomial(logit))
summary(model4)


# OR
round(exp(model4$coefficients),digits=4)

# p Values
round(coef(summary(model4))[3,][[4]],digits=4)

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)

round(exp(confint.default(model4)),digits=4)


# vaccinated with seasonal flu vaccine during year 2008 and 2009 for age group 4
t_base_adult4$vac0809 <- relevel(factor(t_base_adult4$vac0809),ref="No")
model4 <- glm(formula = t_base_adult4$final_fourfold ~ factor(t_base_adult4$vac0809), family = binomial(logit))
summary(model4)


# OR
round(exp(model4$coefficients),digits=4)

# P Values
round(coef(summary(model4))[2,][[4]],digits=4)

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)

round(exp(confint.default(model4)),digits=4)







# all the varaible are not significant at 0.05 level and household size is significant with infection risk at 0.2 significance level
#thus,household size will be adjusted for further analyses performed as follows



tab1 <- ftable(t_base$ag20,t_base$vac0809,t_base$final_fourfold,row.vars=c(1,2))
tab1_dash <- cbind(tab1,percent=tab1[,2]/(tab1[,2]+tab1[,1]))
table(t_base$final_fourfold,t_base$sy_ph_crude)
table(t_base$ag20,t_base$final_fourfold)
table(t_base$ag20,t_base$b40)
table(t_base$ag20,t_base$f40)

table(t_base_children$vac0809,t_base_children$final_fourfold)
table(t_base_adult2$vac0809,t_base_adult2$final_fourfold)
table(t_base_adult3$vac0809,t_base_adult3$final_fourfold)
table(t_base_adult4$vac0809,t_base_adult4$final_fourfold)
table(t_base$vac0809,t_base$final_fourfold)

table(t_base$sex,t_base$final_fourfold)




# Newly added by Danny, three regresson model: vaccination, presence of children, age structure,smoking status and household clustering.
# variabls that have slightly association with risk of infection(p<0.2) will be adjusted in multivariate regression analysis.
# group household size


#  source of recruitment

# univariate analysis of comparision between ramdom dialing and swine flu cohort
tab_source <- ftable(t_base$hh_telsource,t_base$final_fourfold)
cbind(tab_source,percent=tab_source[,2]/(tab_source[,2]+tab_source[,1]))


t_base$hh_telsource <- relevel(factor(t_base$hh_telsource),ref="1")
model_source <- glm(formula = t_base$final_fourfold ~ factor(t_base$hh_telsource), family = binomial(logit))
summary(model_source)


# OR
round(exp(model_source$coefficients),digits=4)

# P Values
round(coef(summary(model_source))[2,][[4]],digits=4)

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)

round(exp(confint.default(model_source)),digits=4)






# date of clinic visit, use the median date as the cutoff date 


t_base  <-cbind(t_base,date_visit=rep(-1,dim(t_base)[1]))

 
 # cute off for every two weeks

for (i in 1:dim(t_base)[1]){
       
        if(t_base[i,"base_ind_date"] < as.date("4jul9")+14) t_base[i,"date_visit"] <- 1
        else if(t_base[i,"base_ind_date"] < as.date("18jul9")+14) t_base[i,"date_visit"] <-2
		else if(t_base[i,"base_ind_date"] < as.date("1Aug9")+14) t_base[i,"date_visit"] <-3
		else if(t_base[i,"base_ind_date"] < as.date("14Aug9")+14) t_base[i,"date_visit"] <-4
		else if(t_base[i,"base_ind_date"] < as.date("28Aug9")+28) t_base[i,"date_visit"] <-5
		#else if(t_base[i,"base_ind_date"] < as.date("11Sep2009")+14) t_base[i,"date_visit"] <-6



 } 
 
 
#t_base$date_visit


tab_visit <- ftable(t_base$date_visit,t_base$final_fourfold)
cbind(tab_visit,percent=tab_visit[,2]/(tab_visit[,2]+tab_visit[,1]))


# timing of clinic visit
t_base$date_visit <- relevel(factor(t_base$date_visit),ref="1")
model_source <- glm(formula = t_base$final_fourfold ~ factor(t_base$date_visit), family = binomial(logit))
summary(model_source)


# OR
round(exp(model_source$coefficients),digits=4)

# P Values
round(coef(summary(model_source))[2,][[4]],digits=4)

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)

round(exp(confint.default(model_source)),digits=4)



# univiate for district 

tab_dist <- ftable(t_base$district5,t_base$final_fourfold)
cbind(tab_dist,percent=tab_dist[,2]/(tab_dist[,2]+tab_dist[,1]))


# timing of clinic visit
t_base$district5 <- relevel(factor(t_base$district5),ref="HK Island")
model_dist <- glm(formula = t_base$final_fourfold ~ factor(t_base$district5), family = binomial(logit))
summary(model_dist)


# OR
round(exp(model_dist$coefficients),digits=4)

# P Values
round(coef(summary(model_dist))[2,][[4]],digits=4)

#CI_vac_y<-cbind(exp(coef(model4)[[2]]-1.96*coef(summary(model4))[2,][[2]]),exp(coef(model4)[[2]]+1.96*coef(summary(model4))[2,][[2]]),P_vac_y)

#colnames(CI_vac_y) <-c("Lower CI","Upper CI","p value")
#round(CI_vac_y,2)

round(exp(confint.default(model_dist)),digits=4)




# Multivarate regression starts here
################################################################  
library(gee)
t_base$vac0809 <- relevel(factor(t_base$vac0809),ref="No")
t_base$smoking <- relevel(factor(t_base$smoking),ref="No")
t_base$ag20 <- relevel(factor(t_base$ag20),ref="1")
t_base$child_hh <- relevel(factor(t_base$child_hh),ref="2")
t_base$hh_size2 <- relevel(factor(t_base$hh_size2),ref="2")
t_base$education2 <- relevel(factor(t_base$education2),ref="1")
t_base$hh_telsource <- relevel(factor(t_base$hh_telsource),ref="1")
t_base$date_visit <- relevel(factor(t_base$date_visit),ref="1")


model_gee <- gee(formula = t_base$final_fourfold ~ factor(t_base$ag20)+factor(t_base$hh_size2)+factor(t_base$child_hh)+factor(t_base$education2)+factor(t_base$hh_telsource)+factor(t_base$profession2)+factor(t_base$district5),id=factor(t_base$hh_id),family = binomial(logit))



summary(model_gee)




results <- data.frame(beta=model_gee$coef, se=sqrt(diag(model_gee[[20]])),
row.names=names(model_gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results <- round(results, 4)
results



#regression analysis ends
#############################################################




























# binomial confidence intrvals 
#######################################
# ORs

round(exp(model9$coefficients),digits=2)



# 95% CI AND p values
# p Values
coef(summary(model9))[2,][[4]]
coef(summary(model9))[3,][[4]]
coef(summary(model9))[4,][[4]]


#alternative way to compute 95%CI
round(exp(confint.default(model9)),digits=2)


# Code to compute Binomial 95% Confidence Interval

SR_BinOptim <- function(p,P,N,n) {
guess <- pbinom(n,N,p,lower.tail=TRUE)
rtn <- P-guess
rtn
}



binCI <- function(vecN,vecn,P,min=0) {
noObs <- length(vecN)
rtn <- array(NA,c(noObs))
for (i in 1:noObs) {
if (vecN[i] > 0) {
sol <- uniroot(SR_BinOptim,c(0,1),P=P,N=vecN[i],n=vecn[i])
if (sol$root > min) rtn[i] <- sol$root
else rtn[i] <- min
} else rtn[i] <- min
}
rtn
}

# 95%CI for age groups

round(binCI(c(98,294,293,82,767),c(41,27,17,1,86),0.975),digits=3)
round(binCI(c(98,294,293,82,767),c(41,27,17,1,86),0.025),digits=3)

# 95%CI for sex
round(binCI(c(462,305),c(47,39),0.975),digits=3)
round(binCI(c(462,305),c(47,39),0.025),digits=3)

# smoking

round(binCI(c(722,45),c(84,2),0.975),digits=3)
round(binCI(c(722,45),c(84,2),0.025),digits=3)

# vaccine
round(binCI(c(628,139),c(67,19),0.975),digits=3)
round(binCI(c(628,139),c(67,19),0.025),digits=3)



# 95%CI for presence of children
round(binCI(c(411,258,669),c(19,26,45),0.975),digits=3)
round(binCI(c(411,258,669),c(19,26,45),0.025),digits=3)






# 95%CI for vaccination status of children
round(binCI(c(80,18),c(33,8),0.975),digits=3)
round(binCI(c(80,18),c(33,8),0.025),digits=3)

# 95%CI for vaccination status of age group2
round(binCI(c(254,41),c(22,5),0.975),digits=3)
round(binCI(c(254,41),c(22,5),0.025),digits=3)

# 95%CI for vaccination status of age group3
round(binCI(c(254,40),c(12,5),0.975),digits=3)
round(binCI(c(254,40),c(12,5),0.025),digits=3)

# 95%CI for vaccination status of age group4
round(binCI(c(42,40),c(0,1),0.975),digits=3)
round(binCI(c(42,40),c(0,1),0.025),digits=3)

# family size
round(binCI(c(139,32,207,254,103,32),c(3,3,19,37,19,5),0.975),digits=3)
round(binCI(c(139,32,207,254,103,32),c(3,3,19,37,19,5),0.025),digits=3)



# education
round(binCI(c(16,128,478,145),c(6,17,52,11),0.975),digits=3)
round(binCI(c(16,128,478,145),c(6,17,52,11),0.025),digits=3)


round(binCI(769,86,0.975),digits=3)

round(binCI(769,86,0.025),digits=3)

