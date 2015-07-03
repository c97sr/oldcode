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


# Variables to be tested
# ag1[Y], 
# symp.fever.rp[N], symp.ILI.rp[N], symp.ARI.rp[N], 
# diary.fever.rp[N], diary.ILI.rp[N], diary.ARI.rp[N], 
# quest.fever[N], quest.ILI[N], quest.ARI[N]


# Make subselections
t_base <- base[base[,"final_fourfold"] > -1,]
t1_base <- t_base[t_base[,"any.fever"] > -1,]
t2_base <- t_base[t_base[,"any.ILI"] > -1,]
t3_base <- t_base[t_base[,"any.ARI"] > -1,]
s0_base <- t_base[t_base[,"symp.fever.rp"] > -1,]
s1_base <- t_base[t_base[,"symp.ILI.rp"] > -1,]
s2_base <- t_base[t_base[,"symp.ARI.rp"] > -1,]
d0_base <- t_base[t_base[,"diary.fever.rp"] > -1,]
d1_base <- t_base[t_base[,"diary.ILI.rp"] > -1,]
d2_base <- t_base[t_base[,"diary.ARI.rp"] > -1,]
q0_base <- t_base[t_base[,"quest.fever"] > -1,]
q1_base <- t_base[t_base[,"quest.ILI"] > -1,]
q2_base <- t_base[t_base[,"quest.ARI"] > -1,]

# Output tables
print('Laboratory tested:')
print(ftable(t_base$ag1, t_base$final_fourfold))
print('Symp.fever:')
print(ftable(s0_base$ag1, s1_base$symp.fever.rp, s1_base$final_fourfold))
print('Symp.ILI:')
print(ftable(s1_base$ag1, s1_base$symp.ILI.rp, s1_base$final_fourfold))
print('Symp.ARI:')
print(ftable(s2_base$ag1, s2_base$symp.ARI.rp, s2_base$final_fourfold))
print('Diary.fever:')
print(ftable(d0_base$ag1, d0_base$diary.fever.rp, d0_base$final_fourfold))
print('Diary.ILI:')
print(ftable(d1_base$ag1, d1_base$diary.ILI.rp, d1_base$final_fourfold))
print('Diary.ARI:')
print(ftable(d2_base$ag1, d2_base$diary.ARI.rp, d2_base$final_fourfold))
print('Questionnaire.fever:')
print(ftable(q0_base$ag1, q0_base$quest.fever, q0_base$final_fourfold))
print('Questionniare.ILI:')
print(ftable(q1_base$ag1, q1_base$quest.ILI, q1_base$final_fourfold))
print('Questionnaire.ARI:')
print(ftable(q2_base$ag1, q2_base$quest.ARI, q2_base$final_fourfold))
print('Reported fever in any source:')
print(ftable(t1_base$ag1, t1_base$any.fever, t1_base$final_fourfold))
print('Reported ILI in any source:')
print(ftable(t1_base$ag1, t1_base$any.ILI, t1_base$final_fourfold))
print('Reported ARI in any source:')
print(ftable(t1_base$ag1, t1_base$any.ARI, t1_base$final_fourfold))

# To compute the p-value
pos.ili.phone <- matrix(c(41,7,45,1),nrow=2)
pos.ili.diary <- matrix(c(41,9,45,3),nrow=2)
pos.ili.quest <- matrix(c(41,19,45,11),nrow=2)
pos.ili.any <- matrix(c(41,24,45,11),nrow=2)

pos.ari.phone <- matrix(c(41,9,45,2),nrow=2)
pos.ari.diary <- matrix(c(41,18,45,11),nrow=2)
pos.ari.quest <- matrix(c(41,24,45,26),nrow=2)
pos.ari.any <- matrix(c(41,31,45,26),nrow=2)

pos.fever.phone <- matrix(c(41,10,45,1),nrow=2)
pos.fever.diary <- matrix(c(41,9,45,5),nrow=2)
pos.fever.quest <- matrix(c(41,21,45,13),nrow=2)
pos.fever.any <- matrix(c(41,27,45,14),nrow=2)

neg.ili.phone <- matrix(c(57,2,627,12),nrow=2)
neg.ili.diary <- matrix(c(57,3,627,10),nrow=2)
neg.ili.quest <- matrix(c(57,7,627,44),nrow=2)
neg.ili.any <- matrix(c(57,11,627,57),nrow=2)

neg.ari.phone <- matrix(c(57,4,627,27),nrow=2)
neg.ari.diary <- matrix(c(57,9,627,80),nrow=2)
neg.ari.quest <- matrix(c(57,26,627,192),nrow=2)
neg.ari.any <- matrix(c(57,29,627,204),nrow=2)

neg.fever.phone <- matrix(c(57,2,627,14),nrow=2)
neg.fever.diary <- matrix(c(57,5,627,20),nrow=2)
neg.fever.quest <- matrix(c(57,8,627,53),nrow=2)
neg.fever.any <- matrix(c(57,13,627,67),nrow=2)

pos.ili.phone
chisq.test(pos.ili.phone)
fisher.test(pos.ili.phone)

# Save subdatasets
write.csv(base, file="M:\\post_proc_data\\base.csv")
write.csv(t_base, file="M:\\post_proc_data\\t_base.csv")
write.csv(labs1, file="M:\\post_proc_data\\labs1.csv")
write.csv(labs2, file="M:\\post_proc_data\\labs2.csv")
write.csv(symp, file="M:\\post_proc_data\\symp.csv")
write.csv(diary, file="M:\\post_proc_data\\diary.csv")
write.csv(quest, file="M:\\post_proc_data\\quest.csv")
