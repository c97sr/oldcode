rm(list=ls(all=TRUE))
options(error=NULL)

# To dos
# - make a variable for household contains child
# - 

source("scrap.R")

base_raw <- read.csv("../anon_data/base_ind.csv")
labs_raw <- read.csv("../anon_data/lab_results.csv")
symp_raw <- read.csv("../anon_data/phone_symp.csv")

pp_data <- post_proc_sero(base_raw, labs_raw, symp_raw)

base <- pp_data$base
labs <- pp_data$labs
symp <- pp_data$symp

# names(base)

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
#[67] "symp.reported"   "symp.fever"      "symp.cough"      "symp.sputum"     "symp.sore_t"     "symp.running_n" 
#[73] "symp.myalgia"    "symp.ILI"        "symp.ARI"   

# Variables to be tested
# ag1[Y], phone.ili[N], phone.ari[N], diary.ili[N], diary.ari[N], symp.ili[N], symp.ari[N]

norowbaseline <- dim(base)[1]
agneg.symp.ILI <- array(-1,c(1,norowbaseline))
agneg.symp.ARI <- array(-1,c(1,norowbaseline))
agpos.symp.ILI <- array(-1,c(1,norowbaseline))
agpos.symp.ARI <- array(-1,c(1,norowbaseline))


# Put in the tested field
for (i in 1:norowbaseline) {

        age.labneg <- ifelse (base$lab_positive==0, base$age, 999)
        tmp <- age.labneg[i]
        
        if (is.na(tmp) || tmp < 0 || base$symp.ILI==0) agneg.symp.ILI[i] <- "NK"
        else if (tmp == 999) agneg.symp.ILI[i] <- "NA"
        else if (tmp < 19) agneg.symp.ILI[i] <- 1
        else if (tmp < 50) agneg.symp.ILI[i] <- 2
        else if (tmp < 65) agneg.symp.ILI[i] <- 3
        else if (tmp < 110) agneg.symp.ILI[i] <- 4
        else stop("Problem with age allocation for agneg.symp.ILI")
        
        if (is.na(tmp) || tmp < 0 || base$symp.ARI==0) agneg.symp.ARI[i] <- "NK"
        else if (tmp == 999) agneg.symp.ARI[i] <- "NA"
        else if (tmp < 19) agneg.symp.ARI[i] <- 1
        else if (tmp < 50) agneg.symp.ARI[i] <- 2
        else if (tmp < 65) agneg.symp.ARI[i] <- 3
        else if (tmp < 110) agneg.symp.ARI[i] <- 4
        else stop("Problem with age allocation for agneg.symp.ARI")
        
        
        age.labpos <- ifelse (base$lab_positive==1, base$age, 999) 
        tmp <- age.labpos[i]
        
        if (is.na(tmp) || tmp < 0 || base$symp.ILI==0) agpos.symp.ILI[i] <- "NK"
        else if (tmp == 999) agpos.symp.ILI[i] <- "NA"
        else if (tmp < 19) agpos.symp.ILI[i] <- 1
        else if (tmp < 50) agpos.symp.ILI[i] <- 2
        else if (tmp < 65) agpos.symp.ILI[i] <- 3
        else if (tmp < 110) agpos.symp.ILI[i] <- 4
        else stop("Problem with age allocation for agpos.symp.ILI")
        
        if (is.na(tmp) || tmp < 0 || base$symp.ARI==0) agpos.symp.ARI[i] <- "NK"
        else if (tmp == 999) agpos.symp.ARI[i] <- "NA"
        else if (tmp < 19) agpos.symp.ARI[i] <- 1
        else if (tmp < 50) agpos.symp.ARI[i] <- 2
        else if (tmp < 65) agpos.symp.ARI[i] <- 3
        else if (tmp < 110) agpos.symp.ARI[i] <- 4
        else stop("Problem with age allocation for agpos.symp.ARI")


}

print(table(agneg.symp.ILI))
print(table(agneg.symp.ARI))
print(table(agpos.symp.ILI))
print(table(agpos.symp.ARI))

# Save data
#age.frame <- as.data.frame(cbind(base$ind_id, base$age, base$ag1))
#colnames(age.frame, do.NULL = TRUE)
#colnames(age.frame) <- c("ind_id", "age", "ag1")    
#write.csv(age.frame, file="M:\\post_proc_data\\age.csv")
