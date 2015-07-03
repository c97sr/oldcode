# A single library file for finctions used for the hk sero survey analysis
require("date")

multimatch <- function(x,vec) {
	rtn <- NULL
	size <- length(vec)
	start <- 1	
	while (start <= size && !is.na(match(x,vec[start:size]))) {
		current <- match(x,vec[start:size]) + start - 1
		rtn <- c(rtn,current)
		start <- current + 1
	}
	rtn
}

post_proc_sero <- function(baseline,laboratory1,laboratory2,symptoms,sympdiary,questionnaire,recruitment,base_date=as.date("04Jul2009")) {
	
	# This function takes the baseline individual data, lab data and symptom reports as input
	# It does all the required post-processing, and adds the data onto the baseline return 
	# - generating a tested / not tested field
	# - generating an age group field
	# - generating a reported symptoms/ not reported symptoms field...	
	
	# Put in the tested field
	norowbaseline <- dim(baseline)[1]
	baseline <- cbind(baseline,pop_recruit=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,recruit_source=rep(-1,norowbaseline))
	baseline <- cbind(baseline,hh_recruit=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,hh_telsource=rep(-1,norowbaseline))
	baseline <- cbind(baseline,tested1=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,tested2=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,final_tested=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,hh_tested=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,test1_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,test2_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,final_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,b20=rep(NA,norowbaseline))
	baseline <- cbind(baseline,f40=rep(NA,norowbaseline))
	baseline <- cbind(baseline,b40=rep(NA,norowbaseline))
	baseline <- cbind(baseline,f20=rep(NA,norowbaseline))
	baseline <- cbind(baseline,pre_titre=rep(-1,norowbaseline))
	baseline <- cbind(baseline,pre_titre_agg=rep(-1,norowbaseline))
	baseline <- cbind(baseline,post_titre=rep(-1,norowbaseline))
	baseline <- cbind(baseline,post_titre_agg=rep(-1,norowbaseline))
	baseline <- cbind(baseline,age=rep(-1,norowbaseline))
	baseline <- cbind(baseline,ag1=rep(0,norowbaseline))
	baseline <- cbind(baseline,ag10=rep(0,norowbaseline))
	baseline <- cbind(baseline,ag20=rep(0,norowbaseline))
	baseline <- cbind(baseline,sy_ph_crude=rep(-1,norowbaseline))
	baseline <- cbind(baseline,sy_ph_tpu=rep(-1,norowbaseline))
	baseline <- cbind(baseline,child_hh=rep(-1,norowbaseline))
	
	# New SR edits April 30
	# Will do these manipulations out of the setup file now
	
	# baseline <- cbind(baseline,bl_wk=rep(-1,norowbaseline))
	# baseline <- cbind(baseline,fu_wk=rep(-1,norowbaseline))
	
	baseline <- cbind(baseline,blood.pair=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,symp.fever.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.ILI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.ARI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,diary.fever.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.ILI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.ARI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,quest.fever=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.ILI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.ARI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.fever=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.ILI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.ARI=rep(-1,norowbaseline))
	
	# Some rows added by SR to consolidate Danny's edits for the regression model
	baseline  <-cbind(baseline,profession2=rep(-1,norowbaseline))
	baseline  <-cbind(baseline,education2=rep(-1,norowbaseline))
			
	norowsymptoms <- dim(symptoms)[1]
	symptoms <- cbind(symptoms,id=rep(NA,norowsymptoms))
	symptoms <- cbind(symptoms,fever.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,cough.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,sputum.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,sore_t.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,running_n.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,shortness_b.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,chills.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,headache.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,vomit.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,diarrehea.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,myalgia.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ILI=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ILI.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ARI=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ARI.rp=rep(-1,norowsymptoms))
	symptoms$id <- paste(symptoms$houseid,symptoms$sub_no,sep="-")
	
	sympdiary <- sympdiary[,1:13] #remove unused columns
	sympdiary <- sympdiary[1:1962,] #remove unused rows
	norowdiary <- dim(sympdiary)[1]
	sympdiary <- cbind(sympdiary,temperature.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,temperature_up.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,temperature_low.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,fever=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,fever.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,chills=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,chills.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,headache=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,headache.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sore_t=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sore_t.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,cough=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,cough.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sputum=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sputum.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,stuffy_n=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,stuffy_n.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,running_n=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,running_n.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,myalgia=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,myalgia.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,vomit=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,vomit.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,diarrehea=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,diarrehea.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ILI=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ILI.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ARI=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ARI.rp=rep(-1,norowdiary))
	
	questionnaire <- questionnaire[1:861,] #remove unused rows
	norowquest <- dim(questionnaire)[1]
	questionnaire <- cbind(questionnaire,subid=rep(NA,norowquest))
	questionnaire <- cbind(questionnaire,ILI=rep(-1,norowquest))
	questionnaire <- cbind(questionnaire,ARI=rep(-1,norowquest))
	
	# Correct some formats
	baseline[,"base_ind_date"] <- as.date(as.character(baseline[,"base_ind_date"]),order="dmy")
	baseline[,"fu_ind_date"] <- as.date(as.character(baseline[,"fu_ind_date"]), order="dmy")
	baseline[,"base_dob"] <- as.date(as.character(baseline[,"base_dob"]),order="dmy")
	baseline$age <- (base_date - baseline$base_dob)/ 365.25
	
	baseline$blood.pair <- ifelse (baseline$base_blood == "Yes" & baseline$fu_blood == "Yes", 1, 0)
	
	str1<-substr(questionnaire$id,1,5)
	str2<-substr(questionnaire$id,6,6)
	str0 <- noquote(paste("S0",str1))
	str1.new <- gsub( "[^[:alnum:]]", "", str0)
	questionnaire$subid <- paste(str1.new,str2,sep="-")
	
	temp.tmp <- as.numeric(levels(sympdiary$temp)[sympdiary$temp])
	
	# Determine prevalence of ILI and ARI on phone symptoms
	for (m in 1:norowsymptoms) {
		if (symptoms$fever[m]==1 && (symptoms$cough[m]==1 || symptoms$sore_t[m]==1)) 
			symptoms$ILI[m] <- 1
		else symptoms$ILI[m] <- 0
		
		symp.sum <- symptoms$fever[m]+symptoms$cough[m]+symptoms$sputum[m]+symptoms$sore_t[m]+symptoms$running_n[m]+symptoms$myalgia[m]
		if (symp.sum>=2) symptoms$ARI[m] <- 1
		else symptoms$ARI[m] <- 0
		
	}
	
	
	# Set symptoms on phone reporting during reporting period
	for (m in 1:norowsymptoms) {
		if (!is.na(symptoms$record[m]) && symptoms$record[m] == 1) {
			record_start <- m
			record_finish <- m -1 + symptoms$record.num[m]
			for (k in record_start: record_finish) {
				if (symptoms$fever[k]==1) symptoms$fever.rp[m] <- 1
				else symptoms$fever.rp[m] <- 0
				if (symptoms$cough[k]==1) symptoms$cough.rp[m] <- 1
				else symptoms$cough.rp[m] <- 0
				if (symptoms$sputum[k]==1) symptoms$sputum.rp[m] <- 1
				else symptoms$sputum.rp[m] <- 0
				if (symptoms$sore_t[k]==1) symptoms$sore_t.rp[m] <- 1
				else symptoms$sore_t.rp[m] <- 0
				if (symptoms$running_n[k]==1) symptoms$running_n.rp[m] <- 1
				else symptoms$running_n.rp[m] <- 0
				if (symptoms$shortness_b[k]==1) symptoms$shortness_b.rp[m] <- 1
				else symptoms$shortness_b.rp[m] <- 0
				if (symptoms$chills[k]==1) symptoms$chills.rp[m] <- 1
				else symptoms$chills.rp[m] <- 0
				if (symptoms$headache[k]==1) symptoms$headache.rp[m] <- 1
				else symptoms$headache.rp[m] <- 0
				if (symptoms$vomit[k]==1) symptoms$vomit.rp[m] <- 1
				else symptoms$vomit.rp[m] <- 0
				if (symptoms$diarrehea[k]==1) symptoms$diarrehea.rp[m] <- 1
				else symptoms$diarrehea.rp[m] <- 0
				if (symptoms$myalgia[k]==1) symptoms$myalgia.rp[m] <- 1                   
				else symptoms$myalgia.rp[m] <- 0
				if (symptoms$ILI[k]==1) symptoms$ILI.rp[m] <- 1
				else symptoms$ILI.rp[m] <- 0
				if (symptoms$ARI[k]==1) symptoms$ARI.rp[m] <- 1
				else symptoms$ARI.rp[m]}
		}
	}                  
	
	
	# Perform computation on symptoms diary
	temp.tmp <- as.numeric(levels(sympdiary$temp)[sympdiary$temp])
	#sympdiary$temp_up and sympdiary$temp_low are numerics, instead of factor
	
	for (j in 1:norowdiary) {    
		
		# Convert fahrenheit to celsius
		if (is.na(temp.tmp[j])) sympdiary$temperature.clean[j] <- -1
		else {if (temp.tmp[j] >= 45) 
				sympdiary$temperature.clean[j] <- (temp.tmp[j] - 32) *5/9
			else sympdiary$temperature.clean[j] <- temp.tmp[j]}
		
		if (is.na(sympdiary$temp_up[j])) sympdiary$temperature_up.clean[j] <- -1
		else {if (sympdiary$temp_up[j] >= 45) 
				sympdiary$temperature_up.clean[j] <- (sympdiary$temp_up[j] - 32) *5/9
			else sympdiary$temperature_up.clean[j] <- sympdiary$temp_up[j]}
		
		if (is.na(sympdiary$temp_low[j])) sympdiary$temperature_low.clean[j] <- -1
		else {if (sympdiary$temp_low[j] >= 45) 
				sympdiary$temperature_low.clean[j] <- (sympdiary$temp_low[j] - 32) *5/9
			else sympdiary$temperature_low.clean[j] <- sympdiary$temp_low[j]}
		
		# Assign symptoms on diary
		if (sympdiary$temperature.clean[j] >= 37.5) sympdiary$fever[j] <- 1        
		else {if (sympdiary$temperature.clean[j] == -1)
			{if (sympdiary$temperature_up.clean[j] == -1 && sympdiary$temperature_low.clean[j] == -1) 
					sympdiary$fever[j] <- -1
				else {if (sympdiary$temperature_low.clean[j] >= 37.5) sympdiary$fever[j] <- 1
					else sympdiary$fever[j] <- 0}
			}
			else sympdiary$fever[j] <- 0
		}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$chills[j] <- -1 
		else {tmp <- regexpr("1",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$chills[j] <- 1
			else if (tmp[1] == -1) sympdiary$chills[j] <- 0
			else sympdiary$chills[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$headache[j] <- -1
		else {tmp <- regexpr("2",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$headache[j] <- 1
			else if (tmp[1] == -1) sympdiary$headache[j] <- 0        
			else sympdiary$headache[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$sore_t[j] <- -1
		else {tmp <- regexpr("3", sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$sore_t[j] <- 1
			else if (tmp[1] == -1) sympdiary$sore_t[j] <- 0
			else sympdiary$sore_t[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$cough[j] <- -1
		else {tmp <- regexpr("4",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$cough[j] <- 1
			else if (tmp[1] == -1) sympdiary$cough[j] <- 0
			else sympdiary$cough[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$sputum[j] <- -1
		else {tmp <- regexpr("5",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$sputum[j] <- 1
			else if (tmp[1] == -1) sympdiary$sputum[j] <- 0
			else sympdiary$sputum[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$stuffy_n[j] <- -1
		else {tmp <- regexpr("6",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$stuffy_n[j] <- 1
			else if (tmp[1] == -1) sympdiary$stuffy_n[j] <- 0
			else sympdiary$stuffy_n[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$running_n[j] <- -1
		else {tmp <- regexpr("7",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$running_n[j] <- 1
			else if (tmp[1] == -1) sympdiary$running_n[j] <- 0
			else sympdiary$runnning_n[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$myalgia[j] <- -1
		else {tmp <- regexpr("8",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$myalgia[j] <- 1
			else if (tmp[1] == -1) sympdiary$myalgia[j] <- 0
			else sympdiary$myalgia[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$vomit[j] <- -1
		else {tmp <- regexpr("9",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$vomit[j] <- 1
			else if (tmp[1] == -1) sympdiary$vomit[j] <- 0
			else sympdiary$vomit[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$diarrehea[j] <- -1
		else {tmp <- regexpr("10",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$diarrehea[j] <- 1
			else if (tmp[1] == -1) sympdiary$diarrehea[j] <- 0
			else sympdiary$diarrehea[j] <- -1}
		
		
		# Determine prevalence of ILI and ARI on symptoms diary
		if (sympdiary$fever[j]==1 && (sympdiary$cough[j]==1 || sympdiary$sore_t[j]==1)) 
			sympdiary$ILI[j] <- 1
		else sympdiary$ILI[j] <- 0
		
		diary.sum <- sympdiary$fever[j]+sympdiary$cough[j]+sympdiary$sputum[j]+sympdiary$sore_t[j]+sympdiary$running_n[j]+sympdiary$myalgia[j]
		if (diary.sum>=2) sympdiary$ARI[j] <- 1
		else sympdiary$ARI[j] <- 0
		
	}
	
	# Set symptoms on symptoms diary during reporting period
	for (j in 1:norowdiary) {
		if (!is.na(sympdiary$record.clean[j]) && sympdiary$record.clean[j] == 1) {
			record_start <- j
			record_finish <- j -1 + sympdiary$record.num[j]
			for (k in record_start: record_finish) {
				if (sympdiary$fever[k]==1) sympdiary$fever.rp[j] <- 1
				else sympdiary$fever.rp[j] <- 0
				if (sympdiary$chills[k]==1) sympdiary$chills.rp[j] <- 1
				else sympdiary$chills.rp[j] <- 0
				if (sympdiary$headache[k]==1) sympdiary$headache.rp[j] <- 1
				else sympdiary$headache.rp[j] <- 0
				if (sympdiary$sore_t[k]==1) sympdiary$sore_t.rp[j] <- 1
				else sympdiary$sore_t.rp[j] <- 0
				if (sympdiary$cough[k]==1) sympdiary$cough.rp[j] <- 1
				else sympdiary$cough.rp[j] <- 0
				if (sympdiary$sputum[k]==1) sympdiary$sputum.rp[j] <- 1
				else sympdiary$sputum.rp[j] <- 0
				if (sympdiary$stuffy_n[k]==1) sympdiary$stuffy_n.rp[j] <- 1
				else sympdiary$stuffy_n.rp[j] <- 0
				if (sympdiary$running_n[k]==1) sympdiary$running_n.rp[j] <- 1
				else sympdiary$running_n.rp[j] <- 0
				if (sympdiary$myalgia[k]==1) sympdiary$myalgia.rp[j] <- 1
				else sympdiary$myalgia.rp[j] <- 0
				if (sympdiary$vomit[k]==1) sympdiary$vomit.rp[j] <- 1
				else sympdiary$vomit.rp[j] <- 0
				if (sympdiary$diarrehea[k]==1) sympdiary$diarrehea.rp[j] <- 1
				else sympdiary$diarrehea.rp[j] <- 0
				if (sympdiary$ILI[k]==1) sympdiary$ILI.rp[j] <- 1
				else sympdiary$ILI.rp[j] <- 0
				if (sympdiary$ARI[k]==1) sympdiary$ARI.rp[j] <- 1
				else sympdiary$ARI.rp[j] <- 0}
		}
	}          
	
	for (r in 1:norowquest) {
		
		#Correct some typos from inconsistencies
		if (questionnaire$id[r] == 900902) questionnaire$no_doctor[r] <- 1
		
		# Determine prevalence of ILI and ARI on contact questionnaire 
		if (questionnaire$fever.clean[r]==1 && (questionnaire$cough.clean[r]==1 || questionnaire$sore_t.clean[r]==1))
			questionnaire$ILI[r] <- 1
		else questionnaire$ILI[r] <- 0
		
		quest.sum <- questionnaire$fever.clean[r]+questionnaire$cough.clean[r]+questionnaire$sputum.clean[r]+questionnaire$sore_t.clean[r]+questionnaire$running_n.clean[r]+questionnaire$myalgia.clean[r]
		if (quest.sum>=2) questionnaire$ARI[r] <- 1
		else questionnaire$ARI[r] <- 0
		
	}
	
	## -- Kendra's editing ends -- ##
	for (i in 1:norowbaseline) {
		
		# Some FU interview dates are in year 2000
		# KW added on 4 June 2010
		basedate1 <- as.date("01Jan1900"); basedate1 <- as.numeric(basedate1)
		basedate2 <- as.date("31Dec1900"); basedate2 <- as.numeric(basedate2)
		
		if (!is.na(baseline$fu_ind_date[i])) {
			if (baseline$fu_ind_date[i] <= basedate2 && baseline$fu_ind_date[i] >= basedate1) {
				tmp <- as.numeric(baseline$fu_ind_date[i])
				tmp <- baseline$fu_ind_date[i] + 40177 #convert those in 1900 to 2010
				baseline$fu_ind_date[i] <- as.date(as.character(tmp, order="dmy"))
			}
		}
		
		
		# Tested or not during test1
		if (baseline$ind_id[i] %in% laboratory1$Our_ref) baseline$tested1[i] <- TRUE
		
		# Result of the testing
		if (baseline$tested1[i]) {           
			lab1row <- match(baseline$ind_id[i],laboratory1$Our_ref)
			if (is.na(lab1row)) stop("There shouldn't be an NA here.")
			baseline$b20[i] <- laboratory1$Base20[lab1row]
			baseline$f40[i] <- laboratory1$Post40[lab1row]
			baseline$b40[i] <- laboratory1$Base40[lab1row]
			baseline$f20[i] <- laboratory1$Post20[lab1row]
			if (baseline$b20[i]==0 || baseline$b40[i]==0) baseline$pre_titre[i] <- 0
			if (baseline$b20[i]==1) baseline$pre_titre[i] <- 20
			if (baseline$b40[i]==1) baseline$pre_titre[i] <- 40
			if (baseline$f20[i]==0 || baseline$f40[i]==0) baseline$post_titre[i] <- 0
			if (baseline$f20[i]==1) baseline$post_titre[i] <- 20
			if (baseline$f40[i]==1) baseline$post_titre[i] <- 40
			if (baseline$b20[i]==0 && baseline$f40[i]==1) baseline$test1_fourfold[i] <- 1
			else baseline$test1_fourfold[i] <- 0
			
		}
		
		## -- Kendra edited: -- ##
		
		# Tested or not
		if (baseline$ind_id[i] %in% laboratory2$ind_id) baseline$tested2[i] <- TRUE
		
		# Result of the testing
		if (baseline$tested2[i]) {           
			lab2row <- match(baseline$ind_id[i],laboratory2$ind_id)
			if (is.na(lab2row)) stop("There shouldn't be an NA here.")
			baseline$test2_fourfold[i] <- laboratory2$Kpos[lab2row]
			baseline$pre_titre[i] <- laboratory2$Base_titre[lab2row]
			baseline$post_titre[i] <- laboratory2$Post_titre[lab2row]
		}
		
		baseline$pre_titre_agg[i] <- baseline$pre_titre[i]
		baseline$post_titre_agg[i] <- baseline$post_titre[i]
		
		if (baseline$pre_titre[i]==0 || baseline$pre_titre[i]==5) baseline$pre_titre_agg[i] <- 10
		if (baseline$post_titre[i]==0 || baseline$post_titre[i]==5) baseline$post_titre_agg[i] <-10
		
		if (baseline$tested1[i]==TRUE || baseline$tested2[i]==TRUE) baseline$final_tested[i] <- TRUE
		
		if (baseline$tested1[i]==TRUE && baseline$tested2[i]==FALSE) baseline$final_fourfold[i] <- baseline$test1_fourfold[i]
		if (baseline$tested1[i]==FALSE && baseline$tested2[i]==TRUE) baseline$final_fourfold[i] <- baseline$test2_fourfold[i] 
		if (baseline$tested1[i]==TRUE && baseline$tested2[i]==TRUE) baseline$final_fourfold[i] <- baseline$test2_fourfold[i]
		
		#Correct some formats
		if (baseline$ind_id[i] == "S090433-0") baseline$base_ind_date[i] <- as.date("27Aug2009")
		
		if (baseline$ind_id[i] == "S090002-0") {baseline$base_dob[i] <- as.date("04Oct2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090008-0") {baseline$base_dob[i] <- as.date("11Feb2004")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090033-0") {baseline$base_dob[i] <- as.date("10Apr2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090048-2") {baseline$base_dob[i] <- as.date("10Aug2001")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090056-0") {baseline$base_dob[i] <- as.date("01Apr2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090061-4") {baseline$base_dob[i] <- as.date("10Oct1918")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090066-3") {baseline$base_dob[i] <- as.date("01Jul1914")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090076-0") {baseline$base_dob[i] <- as.date("12Jun2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090093-0") {baseline$base_dob[i] <- as.date("12Feb2009")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090093-1") {baseline$base_dob[i] <- as.date("12Jan2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090094-2") {baseline$base_dob[i] <- as.date("01Jan1928")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090095-0") {baseline$base_dob[i] <- as.date("09Sep2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090097-0") {baseline$base_dob[i] <- as.date("01Apr2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090097-1") {baseline$base_dob[i] <- as.date("01Apr2006")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090102-3") {baseline$base_dob[i] <- as.date("01Jan1924")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-0") {baseline$base_dob[i] <- as.date("01Feb2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-1") {baseline$base_dob[i] <- as.date("04Feb2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-2") {baseline$base_dob[i] <- as.date("05Jul2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090116-0") {baseline$base_dob[i] <- as.date("09Jun2000")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090117-4") {baseline$base_dob[i] <- as.date("06Apr1927")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090117-5") {baseline$base_dob[i] <- as.date("26Nov1921")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090125-0") {baseline$base_dob[i] <- as.date("06Jul2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-0") {baseline$base_dob[i] <- as.date("01Nov2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-1") {baseline$base_dob[i] <- as.date("01May2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-2") {baseline$base_dob[i] <- as.date("01Dec2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090143-0") {baseline$base_dob[i] <- as.date("04Dec2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090144-3") {baseline$base_dob[i] <- as.date("01Jan1924")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090152-0") {baseline$base_dob[i] <- as.date("07Jan2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090157-0") {baseline$base_dob[i] <- as.date("07Nov2007")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090168-0") {baseline$base_dob[i] <- as.date("06Dec2007")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090173-3") {baseline$base_dob[i] <- as.date("06Feb1921")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090178-0") {baseline$base_dob[i] <- as.date("01Jan2005")   
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090194-1") {baseline$base_dob[i] <- as.date("04Aug1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090202-2") {baseline$base_dob[i] <- as.date("01Jan1911") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090203-1") {baseline$base_dob[i] <- as.date("06Aug2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090204-2") {baseline$base_dob[i] <- as.date("29Jul1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090207-3") {baseline$base_dob[i] <- as.date("01Jan1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090213-0") {baseline$base_dob[i] <- as.date("07Feb2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090222-1") {baseline$base_dob[i] <- as.date("07Dec2000")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090245-0") {baseline$base_dob[i] <- as.date("09Apr2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090249-1") {baseline$base_dob[i] <- as.date("01Sep1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090266-1") {baseline$base_dob[i] <- as.date("03Oct2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090267-4") {baseline$base_dob[i] <- as.date("01Jan1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090275-0") {baseline$base_dob[i] <- as.date("06Apr2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090283-3") {baseline$base_dob[i] <- as.date("01Jan1929")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090285-0") {baseline$base_dob[i] <- as.date("10Feb2009") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090293-0") {baseline$base_dob[i] <- as.date("29Aug2007") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090300-1") {baseline$base_dob[i] <- as.date("11Feb2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090332-1") {baseline$base_dob[i] <- as.date("10Aug1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090332-2") {baseline$base_dob[i] <- as.date("15Apr1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090339-2") {baseline$base_dob[i] <- as.date("01Mar1924") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090334-1") {baseline$base_dob[i] <- as.date("01Jan1926") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090363-0") {baseline$base_dob[i] <- as.date("01Jan1908") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090371-4") {baseline$base_dob[i] <- as.date("01Jan1920")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090372-0") {baseline$base_dob[i] <- as.date("10Jul2000") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090373-2") {baseline$base_dob[i] <- as.date("01Jan1922")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090379-0") {baseline$base_dob[i] <- as.date("01May2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090379-1") {baseline$base_dob[i] <- as.date("01May2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090401-1") {baseline$base_dob[i] <- as.date("01Jan1920") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090412-1") {baseline$base_dob[i] <- as.date("04May2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090413-0") {baseline$base_dob[i] <- as.date("01Jan2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090415-0") {baseline$base_dob[i] <- as.date("30Aug2003") #updated on 7 June 2010
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090418-0") {baseline$base_dob[i] <- as.date("01Jan2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090418-1") {baseline$base_dob[i] <- as.date("01Jan2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090430-0") {baseline$base_dob[i] <- as.date("06Apr2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}    
		if (baseline$ind_id[i] == "S090430-1") {baseline$base_dob[i] <- as.date("09Apr2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090454-0") {baseline$base_dob[i] <- as.date("06Apr2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090464-1") {baseline$base_dob[i] <- as.date("01Jan1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090470-4") {baseline$base_dob[i] <- as.date("01Sep1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090476-4") {baseline$base_dob[i] <- as.date("01Jan1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090486-0") {baseline$base_dob[i] <- as.date("01Jan2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090513-1") {baseline$base_dob[i] <- as.date("05Jan2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090518-3") {baseline$base_dob[i] <- as.date("01Jan1913") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090525-0") {baseline$base_dob[i] <- as.date("11Jun2001")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090525-3") {baseline$base_dob[i] <- as.date("01Jan1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090531-4") {baseline$base_dob[i] <- as.date("01Jan1917") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090550-1") {baseline$base_dob[i] <- as.date("08May2000") #updated on 7 June 2010
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090556-0") {baseline$base_dob[i] <- as.date("04Dec2008")  
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090583-1") {baseline$base_dob[i] <- as.date("06Nov2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090584-2") {baseline$base_dob[i] <- as.date("29May1926") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090587-0") {baseline$base_dob[i] <- as.date("04Aug2003") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		
		if (baseline$ind_id[i] == "S090324-2") {
			baseline$fu_ind_date[i] <- as.date("29Dec2009") #added on 06 July 2010
			baseline$fu_attendance[i] <- "Yes"
			baseline$fu_blood[i] <- "Yes"}
		
		# Add by KMW on 1 June 2010
		if (baseline$ind_id[i] == "S090008-0") baseline$profession[i] <- 12
		if (baseline$ind_id[i] == "S090008-1") baseline$profession[i] <- 12
		if (baseline$ind_id[i] == "S090358-0") baseline$profession[i] <- 3
		
		## -- Kendra's editing ends -- ##
		
		# Set age and age groups
		if (is.na(baseline[i,"age"])) baseline$ag1[i] <- -1
		else if (baseline[i,"age"] < 19) baseline[i,"ag1"] <- 1
		else if (baseline[i,"age"] <= 48) baseline[i,"ag1"] <- 2
		else if (baseline[i,"age"] < 65) baseline[i,"ag1"] <- 3
		else if (baseline[i,"age"] < 110) baseline[i,"ag1"] <- 4
		else stop("Problem with the age group ag1 allocation")
		
		#Add by KMW on 8 June 2010
		if (is.na(baseline[i,"age"])) baseline$ag10[i] <- -1
		else if (baseline[i,"age"] < 10) baseline$ag10[i] <- 1
		else if (baseline[i,"age"] < 20) baseline$ag10[i] <-  2
		else if (baseline[i,"age"] < 30) baseline$ag10[i] <- 3
		else if (baseline[i,"age"] < 40) baseline$ag10[i] <- 4
		else if (baseline[i,"age"] < 50) baseline$ag10[i] <- 5
		else if (baseline[i,"age"] < 60) baseline$ag10[i] <- 6
		else if (baseline[i,"age"] < 70) baseline$ag10[i] <- 7
		else if (baseline[i,"age"] < 80)  baseline$ag10[i] <- 8
		else if (baseline[i,"age"] < 110)  baseline$ag10[i] <- 9
		
		if (is.na(baseline[i,"age"])) baseline$ag20[i] <- -1
		else if (baseline[i,"age"] < 20) baseline[i,"ag20"] <- 1
		else if (baseline[i,"age"] < 40) baseline[i,"ag20"] <- 2
		else if (baseline[i,"age"] < 60) baseline[i,"ag20"] <- 3
		else if (baseline[i,"age"] < 110) baseline[i,"ag20"] <- 4
		else stop("Problem with ag20 allocation")
		
		## -- Kendra edited: -- ##
		# Reported symptoms or not through the phone
		if (baseline$ind_id[i] %in% symptoms$id) baseline$symp.reported[i] <- TRUE
		
		# Result of the reporting
		sympindex <- match(baseline$ind_id[i],symptoms$id)
		if (baseline$symp.reported[i]) {
			if (is.na(sympindex)) stop("There shouldn't be an NA in sympindex among reported")
			baseline$symp.fever.rp[i] <- symptoms$fever.rp[sympindex]
			baseline$symp.ILI.rp[i] <- symptoms$ILI.rp[sympindex]
			baseline$symp.ARI.rp[i] <- symptoms$ARI.rp[sympindex]          
		}
		
		# Reported symptoms or not on symptoms diary 
		if (baseline$ind_id[i] %in% sympdiary$id) baseline$diary.reported[i] <- TRUE
		
		# Result of the reporting
		diaryindex <- match(baseline$ind_id[i], sympdiary$id)
		if (baseline$diary.reported[i]) {
			if (is.na(diaryindex)) stop("There shouldn't be an NA in diaryindex among reported")
			baseline$diary.fever.rp[i] <- sympdiary$fever.rp[diaryindex]
			baseline$diary.ILI.rp[i] <- sympdiary$ILI.rp[diaryindex]
			baseline$diary.ARI.rp[i] <- sympdiary$ARI.rp[diaryindex]
		}
		
		
		# Reported symptoms or not on contact questionniare
		if (baseline$ind_id[i] %in% questionnaire$subid) baseline$quest.reported[i] <- TRUE
		
		# Result of the reporting
		questindex <- match(baseline$ind_id[i], questionnaire$subid)
		if (baseline$quest.reported[i]) {
			if (is.na(questindex)) stop("There shouldn't be an NA in questindex among reported")
			baseline$quest.fever[i] <- questionnaire$fever[questindex]
			baseline$quest.ILI[i] <- questionnaire$ILI[questindex]
			baseline$quest.ARI[i] <- questionnaire$ARI[questindex]
		}
		
		if (baseline$symp.fever.rp[i]==1 || baseline$diary.fever.rp[i]==1 || baseline$quest.fever[i]==1) 
			baseline$any.fever[i] <- 1
		else baseline$any.fever[i] <- 0
		
		if (baseline$symp.ILI.rp[i]==1 || baseline$diary.ILI.rp[i]==1 || baseline$quest.ILI[i]==1)
			baseline$any.ILI[i] <- 1
		else baseline$any.ILI[i] <- 0
		
		if (baseline$symp.ARI.rp[i]==1 || baseline$diary.ARI.rp[i]==1 || baseline$quest.ARI[i]==1) 
			baseline$any.ARI[i] <- 1
		else baseline$any.ARI[i] <- 0
		
		
		#Where did POP recruit this household
		if (baseline$ind_id[i] %in% recruitment$ind_id) baseline$pop_recruit[i] <- recruitment$Match[i]
		recruitindex <- match(baseline$ind_id[i], recruitment$ind_id)
		if (baseline$pop_recruit[i]) baseline$recruit_source[i] <- recruitment$Source[recruitindex]
		
		
		## -- Kendra's editing ends -- ##
		
		
		# Most crude phone symptoms
		if (!is.na(sympindex)) baseline[i,"sy_ph_crude"] <- 1
		else baseline[i,"sy_ph_crude"] <- 0
		
		if (!is.na(sympindex) && symptoms[sympindex,"self_r"] == 1) baseline[i,"sy_ph_tpu"] <- 1
		else baseline[i,"sy_ph_tpu"] <- 0
		
		# Copying in Danny's reclassifications
		if (baseline[i,"profession"] %in% c(1,2)) baseline[i,"profession2"] <- 1
		else if (baseline[i,"profession"] %in% c(3,4)) baseline[i,"profession2"] <- 2
		else if (baseline[i,"profession"] %in% c(5,6)) baseline[i,"profession2"] <- 3
		else if (baseline[i,"profession"] %in% c(7,8,9)) baseline[i,"profession2"] <-4
		else if(baseline[i,"profession"] %in% c(10,11,13,14)) baseline[i,"profession2"] <- 5
		else if(baseline[i,"profession"] %in% c(12)) baseline[i,"profession2"] <- 6
		
#		# date of clinic visit, use the median date as the cutoff date 	
#		baseline  <-cbind(baseline,date_visit=rep(-1,norowbaseline))
#	
#		# cute off for every two weeks
#
#		if (baseline[i,"base_ind_date"] < as.date("4jul9")+14) baseline[i,"date_visit"] <- 1
#		else if (baseline[i,"base_ind_date"] < as.date("18jul9")+14) baseline[i,"date_visit"] <-2
#		else if (baseline[i,"base_ind_date"] < as.date("1Aug9")+14) baseline[i,"date_visit"] <-3
#		else if (baseline[i,"base_ind_date"] < as.date("14Aug9")+14) baseline[i,"date_visit"] <-4
#		else if (baseline[i,"base_ind_date"] < as.date("28Aug9")+28) baseline[i,"date_visit"] <-5
#		#else if(baseline[i,"base_ind_date"] < as.date("11Sep2009")+14) baseline[i,"date_visit"] <-6

		#if (is.na(baseline[i,"education"])) baseline$education2[i] <- -1
		if (baseline[i,"education"] ==0) baseline[i,"education2"] <- 0
		else if (baseline[i,"education"] ==1) baseline[i,"education2"] <- 1
		else if (baseline[i,"education"] ==2) baseline[i,"education2"] <- 2
		else if (baseline[i,"education"] ==3) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==4) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==5) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==6) baseline[i,"education2"] <- 4
		else stop("Problem with the grouping education level")
		
	}          
	
	## More editing from Kendra
	# Household stuff
	for (i in 1:norowbaseline) {
		
		#Check whether someone in the household has had their paired samples tested
		if (baseline[i,"hh_index"]==0) {
			hh_index_start <- i
			hh_index_finish <- i -1 + baseline[i,"hh_size"]
		}
		
		# Set presence child_household
		if (1 %in% baseline[hh_index_start:hh_index_finish,"ag1"]) baseline[i,"child_hh"] <- 1
		else if (-1 %in% baseline[hh_index_start:hh_index_finish,"ag1"]) {
			if (1 %in% baseline[hh_index_start:hh_index_finish,"clean_is_child"]) 
				baseline[i,"child_hh"] <- 1
			else 
				baseline[i,"child_hh"] <- 2
		}
		else baseline[i,"child_hh"] <- 2
		
		if (TRUE %in% baseline[hh_index_start: hh_index_finish,"final_tested"]) baseline$hh_tested[i] <- TRUE
		
		if (baseline$ind_id[i] == "S090546-x") baseline$hh_tested[i] <- FALSE
		
		# Where did the household recruit from
		if (TRUE %in% baseline[hh_index_start: hh_index_finish, "pop_recruit"]) baseline$hh_recruit[i] <- TRUE
		if (1 %in% baseline[hh_index_start: hh_index_finish, "recruit_source"]) baseline$hh_telsource[i] <- 1
		else if (2 %in% baseline[hh_index_start: hh_index_finish, "recruit_source"]) baseline$hh_telsource[i] <- 2
		
		
	}
	
	# Set up year week of recruitment SR editing here
	
	#baseline$ag1 <- relevel(factor(baseline$ag1),ref="2")
	
	list(base=baseline,labs1=laboratory1,labs2=laboratory2,symp=symptoms,diary=sympdiary,quest=questionnaire,recruit=recruitment)
	
}


sev_g <- function(k,bi,fi,q1,q2,q3) {
	if (identical(bi-k,2)) 			rtn <- 1-q2
	else if (identical(bi-k,1)) 	rtn <- 1-q1
	else if (k >= bi && k < fi-2) 	rtn <- 1
	else if (identical(fi-k,2)) 	rtn <- q2
	else if (identical(fi-k,1)) 	rtn <- q1
	else 							rtn <- 0
	as.numeric(rtn)
}

sev_r <- function(p,index,ya,b,f,N,q1,q2,minw=1,maxw=length(ya)) {
	rtn <- 0
	bi <- b[index]
	fi <- f[index]
	for (i in minw:maxw) rtn <- rtn + sev_g(i,bi,fi,q1,q2)*ya[i]
	rtn <- rtn / p / N
	as.numeric(rtn)
}

sev_like <- function(p,q1,q2,xa,ba,fa,ya,N,noobs=length(xa),verbose=FALSE) {
	rtn <- 0
	for (i in 1:noobs) {
		prob <- sev_r(p,i,ya,ba,fa,N,q1,q2)
		if (prob > 1) rtn <- rtn - 1e10
		else {
			if (xa[i]==1) rtn <- rtn + log(prob)
			else rtn <- rtn + log(1-prob)
		}
	}
	if (verbose) {
		cat("p ",p," val ",rtn,"\n")
		flush.console()
	}
	rtn
}

sev_like_ci <- function(p,target,q1,q2,xa,ba,fa,ya,N,noobs=length(xa)) {
	raw <- sev_like(p,q1,q2,xa,ba,fa,ya,N)
	abs(raw-(target-1.96))
}

sim_study <- function(fnSimData,vecPs,vecAgs,matWkSev,vecWkBins,q1,q2) {
	
	# Function to simulate a sero study
	# Takes a study protocol, a vector of age specific probabilities, 
	# Not sure that this is working at all XXXX
	
	if (is.character(fnSimData)) {
	
		rtn <- read.csv(fnSimData)
		rtn$base_ind_date <- as.date(as.character(rtn$base_ind_date))
		rtn$fu_ind_date <- as.date(as.character(rtn$fu_ind_date))
		
	} else {
		
		rtn <- fnSimData
		
	}
	
	noSubs <- dim(rtn)[1]
	noWks <- length(vecWkBins)
	if (noWks != dim(matWkSev)[2]+1) stop("problem matching wk bins with severity mtrix")
	noAgs <- length(vecPs)
	if (noAgs != length(vecAgs)) stop("problem matching age group prob vector and ag size vector")
	stFirstWk <- vecWkBins[1]
	if (rtn$base_ind_date[1] < stFirstWk) stop("likely a problem with the date format in a study protocol file")
	
	# Cycle through each invidiual and each week
	for (i in 1:noSubs) {
		base_wk <- floor((rtn$base_ind_date[i]-stFirstWk)/7)+1
		fu_wk <- floor((rtn$fu_ind_date[i]-stFirstWk)/7)+1
		ag <- rtn$ag20[i]
		infs <- 0
		if (base_wk - 2 > 0) infs <- infs + matWkSev[ag,base_wk-2]* (1 - q2) / vecPs[ag]
		if (base_wk - 1 > 0) infs <- infs + matWkSev[ag,base_wk-1]* (1 - q1) / vecPs[ag]
		for (j in base_wk:(fu_wk-3)) infs <- infs + matWkSev[ag,j] / vecPs[ag]
		infs <- infs + matWkSev[ag,fu_wk-2] / vecPs[ag] * q2
		infs <- infs + matWkSev[ag,fu_wk-1] / vecPs[ag] * q1
		ind_p <- infs / vecAgs[ag]
		if (ind_p > 1) stop("problem with simulation")
		if (runif(1) < ind_p) rtn$final_fourfold[i] <- 1
		else rtn$final_fourfold[i] <- 0
	}
	
	rtn
	
}

summarize_models <- function(
						lstModels,
						outfunc=exp,
						writetab=TRUE,
						file="modsum.csv",
						sigdigits=3,
						transpose=FALSE) {
	
	# Figure out the number of models
	nomods <- length(lstModels) 
	
	# Make a vector of all coefficients
	allCoeffs <- c()
	for (i in 1:nomods) {
		
		# Select current model results
		mod <- lstModels[[i]]
		
		# Get a list of variables
		vars <- names(mod$coefficients)
		novars <- length(vars)
		
		# Go through each variabel and add it if its not already in the list
		for (j in 1:novars) {
			
			# Get the variable name
			curname <- vars[j]
			
			# Test for the presence of the variable in the master list
			var_present <- (curname %in% allCoeffs)
			
			# If not in the list add it
			if (!(var_present)) allCoeffs <- c(allCoeffs,curname)
			
		# Close the for loop for j
		}
		
	# Close the for loop for i	
	}
	
	# Define the data structures used to extract the information from the models 	
	noCoeffs <- length(allCoeffs)
	matPointEst <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	matLB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	matUB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	vecAIC <- vector(mode="numeric",length=nomods)
	
	# Loop back though the models and the coeffciients to populate the data structures
	for (i in 1:nomods) {
	
		# Select current model results
		mod <- lstModels[[i]]
		cis <- confint.default(mod)
	
		# Get a list of variables
		vars <- names(mod$coefficients)
		novars <- length(vars)

		# Record the AIC
		vecAIC[i] <- mod$aic
		
		# Go through each variabel and add it if its not already in the list
		for (j in 1:novars) {
		
			# Get the variable name
			curname <- vars[j]
					
			# Extract the point estimate and confidence intervals for the parameters
			matPointEst[i,curname] 	<- mod$coefficients[curname]  
			matLB[i,curname] 		<- cis[curname,1]  
			matUB[i,curname] 		<- cis[curname,2]  
						
		# Close the for loop for j
		}
	
	# Close the for loop for i	
	}
	
	# If selected, write a nicely formatted csv table for the parameters and models
	if (writetab) {
		
		if (transpose) {
			
			# Declare the output string
			strTable <- ""
			
			# Put in the first header row
			strTable <- paste(strTable,"Parameter",sep="")
			for (i in 1:noCoeffs) strTable <- paste(strTable,",",allCoeffs[i],",",allCoeffs[i],sep="")
			strTable <- paste(strTable,",AIC\n",sep="")
			
			# Put in the second header row
			strTable <- paste(strTable,"Model",sep="")
			for (i in 1:noCoeffs) strTable <- paste(strTable,",PE,CI",sep="")
			strTable <- paste(strTable,",AIC\n",sep="")
			
			# Output individual model lines, starting with coefficient loop
			for (i in 1:nomods) {
				
				# Pull the name of the current coefficient
				# curname <- allCoeffs[i]
				
				# Put in the name of the coefficient
				strTable <- paste(strTable,i,sep="")
				
				# Cycle through the tables looking at the different models
				for (j in 1:noCoeffs) {
					
					# Itentify the current coefficient
					curname <- allCoeffs[j]
					
					# Put in the point estimates and confidence intervals for each parameter / model combination
					curPE <- signif(outfunc(matPointEst[i,curname]),digits=sigdigits)
					curLB <- signif(outfunc(matLB[i,curname]),digits=sigdigits)
					curUB <- signif(outfunc(matUB[i,curname]),digits=sigdigits)
					
					# Paste in the parameter values and the confidence intervals
					if (is.na(curPE)) {
						
						# Put in the entry for NA results
						strTable <- paste(strTable,",","-",",","-",sep="")
						
					} else {
						
						# Put in the entry for non NA results
						strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
						
					}
					
				# End j loop for coefficients
				}
				
				# Add the AIC at the end of the line, with a return
				mod <- lstModels[[i]]
				curAIC <- round(mod$aic,digits=1)
				strTable <- paste(strTable,",",curAIC,"\n",sep="")
				
			# End the i for loop for models
			}
		
		# End the if clause for transpose
		} else {
			
			# Declare the output string
			strTable <- ""
			
			# Put in the first header row
			strTable <- paste(strTable,",Model 1",sep="")
			if (nomods>1) for (i in 2:nomods) strTable <- paste(strTable,",,Model ",i,sep="")
			strTable <- paste(strTable,"\n",sep="")
			
			# Put in the second header row
			if (nomods>1) for (i in 1:nomods) strTable <- paste(strTable,",Estimate,(95% CI)",sep="")
			strTable <- paste(strTable,"\n",sep="")
			
			# Output individual coefficient lines, starting with coefficient loop
			for (i in 1:noCoeffs) {
				
				# Pull the name of the current coefficient
				curname <- allCoeffs[i]
				
				# Put in the name of the coefficient
				strTable <- paste(strTable,curname,sep="")
				
				# Cycle through the tables looking at the different models
				for (j in 1:nomods) {
					
					# Put in the point estimates and confidence intervals for each parameter / model combination
					curPE <- signif(outfunc(matPointEst[j,curname]),digits=sigdigits)
					curLB <- signif(outfunc(matLB[j,curname]),digits=sigdigits)
					curUB <- signif(outfunc(matUB[j,curname]),digits=sigdigits)
					
					# Paste in the parameter values and the confidence intervals
					if (is.na(curPE)) {
						
						# Put in the entry for NA results
						strTable <- paste(strTable,",","-",",","-",sep="")
						
					} else {
						
						# Put in the entry for non NA results
						strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
						
					}
					
					# End model for loop
				}
				
				# Return at the end of the line
				strTable <- paste(strTable,"\n",sep="")
				
				# End for for coeffs	
			}
			
			
			# Write the row name for the AICs	
			strTable <- paste(strTable,"AIC",sep="")
			
			# Start the for loop for the AICs
			for (i in 1:nomods) {
				
				# Get the current model
				mod <- lstModels[[i]]
				
				# Format the AIC for the current model
				curAIC <- round(mod$aic,digits=1)
				
				# Write the value and the space
				strTable <- paste(strTable,",",curAIC,",",sep="")
			}
			
			# Return at the end of the AIC line
			strTable <- paste(strTable,"\n",sep="")
		
		# End else statement for transpose of 
		}
		
		# Write the string to the selected file
		cat(strTable,file=file)
		
	}
	
	# Return 
	list(pe=matPointEst,lb=matLB,ub=matUB,aic=vecAIC)
	
}

# A function to make a rolling average of p ff based on age of width width
# Doesn't assume that x is sorte, but does assume it can be ordered
# Uses the median value for the y axis if median is true, otherwise uses mean
windowInc <- function(x,p,width=100,use_median=TRUE) {
	
	# Set up some auxilliary variables
	noMeasures <- length(x)
	noRolling <- noMeasures - width + 1
	if (noRolling < 1) stop("not enough data points for that rolling average")
	
	# Setup return values
	vecX <- vector(mode="numeric",length=noRolling)
	vecY <- vector(mode="numeric",length=noRolling)
	
	# Sort x and p
	x_dash <- x[order(x)]
	p_dash <- p[order(x)]
	
	# Populate the rolling averagesim_study
	
	for (i in width:noMeasures) {
		
		# Setup the loop variabels
		start_window <- i - width + 1
		end_window <- i
		
		# Calculate the rolling average
		vecY[i] <- sum(p_dash[start_window:end_window]) / width
		if (use_median) vecX[i] <- median(x_dash[start_window:end_window])
		else vecX[i] <- mean(x_dash[start_window:end_window])
		
	}
	
	list(x=vecX,y=vecY)
	
}


# An aux routine to compute Binomial 95% Confidence Interval
SR_BinOptim <- function(p,P,N,n) {
	guess <- pbinom(n,N,p,lower.tail=TRUE)
	rtn <- P-guess
	rtn
}

# Computes Binomial 95% Confidence Interval
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

# Auxilliary function for the full set of results
srPieFunc <- function(x,y,vec,rscale=1/10,base=4) {
	agcounts <- vec
	totalcount <- sum(vec)
	cumulativecount <- 0
	for (ag in (1:length(vec))) {
		if (totalcount > 0) circle(x=x,y=y,r= rscale*(log(totalcount,base=base)+1),theta=c(cumulativecount/totalcount*360,(cumulativecount+agcounts[ag])/totalcount*360),col=vecColors[ag])
		cumulativecount <- cumulativecount + agcounts[ag] 
	}	
}

# Make the symptom plot
# Data frame for counts, total number of positives, filename, ylim for chart
plot.symptoms <- function(df, N, filename, ylim_in = c(0,1.0), ms=0.3, ...) {
	
	# Define the number of ...
	fw <- 8.3/cm(1)
	fh <- 8.3/cm(1)
	
	
	# Open the pdf file
	pdf(filename,height=fh,width=fw)
	
	# Set some standard parameter options
	par(mai=(c(0.15*fh,0.175*fw,0,0)))
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	par(cex = 10/12)
	
	# Set up the vectors
	vecLev1 <- names(df)
	noLev1 <- length(vecLev1)
	offvec1 <- (1:noLev1 - 0.5)
	
	vecLev2 <- row.names(df)
	noLev2 <- length(vecLev2)
	offvec2 <- ((0:(noLev2-1))/(noLev2-1)*ms - ms/2)
	
	colvec <- c("red","green","blue","cyan","magenta")
	
	
	# Setup the axes
	plot(1:2,type="n", xlim=c(0,noLev1), ylim=ylim_in, axes=FALSE)
	axis(2,las=1)
	veclabs <- rep(" ",noLev1*2+1)
	for (i in 1:noLev1) veclabs[i*2] <- vecLev1[i]
	axis(1,at=(0:(noLev1*2))/2,labels=veclabs)
	
	# Set up a loop for the main offset
	for (i in 1:noLev1) {
		
		# Set up a loop for the secondary offset
		for (j in 1:noLev2) {
			
			# Plot the points and the confidence intervals
			xoff <- offvec1[i] + offvec2[j]
			n <- df[j,i]
			pest <- n/N
			lb <- binCI(c(N),c(n),0.975)
			ub <- binCI(c(N),c(n),0.025)
			points(c(xoff),c(pest),col=colvec[j],pch=22,bg=colvec[j])
			points(c(xoff,xoff),c(lb,ub),type="l",col=colvec[j])
			
		}
		
	}
	
	legend(0,max(ylim_in),legend=vecLev2,col=colvec,pt.bg=colvec,pch=22,bty="n")
	
	# Close the pdf device
	dev.off()
	
}

# A function to take a character string of a typical serological result and convert it to a raw titre
# Different lookup-tables need to be coded in different ways
# Initially, assumes a full dilution schedule starting at 1:10 
convertRawTitre <- function(charTitre,serotable="full_start_1_10") {
	
	# Define the table
	if (serotable=="full_start_1_10") {
		luTable <- data.frame(	raw=c("<1:20","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560"),
								lu=c(1,2,3,4,5,6,7,8,9)
								)
	}
	
	# Define return vector
	tmp <- luTable[luTable$raw==charTitre,"lu"]
	if (length(tmp) != 1) stop("Wrong value in convertRawTitre")
	
	# Return the correct value
	tmp
	
}

FigPersistenceCompare <- function(solutionA, solutionB, solutionN,
		file="fig1.pdf") {
	
	
	# Plot the illustartive scenarios
	pdf(file,height=10/cm(1),width=10/cm(1))
	
	# standard figure parameter adjustments
	par(	cex=0.8,
			mgp=c(2.5,1,0),
			mai=(c(0,0,0,0)),
			cex.axis=0.9,
			fig=c(0.2,1.0,0.1,1.0),
			las=1)
	
	# Set up values for the axes
	yaxisvals <- 10^(0:5)
	yaxislabs <- c("1","10","100","1,000","10,000","100,000")
	
	# Plot the x axis
	plot(	1:2,
			#log="y",
			xlab="Time (years)",
			ylab="Incidence (per day)",
			type="n",
			xlim=c(0,ill_params["noyears"]),
			ylim=c(1,20),
			axes=FALSE)
	axis(1)
	axis(2)
	#axis(2,at=yaxisvals,labels=yaxislabs)
	points(solutionA[,"time"]/364,(solutionA[,"dS"]+1),type="l",col="red")
	points(solutionB[,"time"]/364,(solutionB[,"dS"]+1),type="l",col="green")
	points(solutionN[,"time"]/364,(solutionN[,"dS"]+1),type="l",col="black")
	dev.off()
	
}