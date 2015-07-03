# To dos: 
# - Next : tidy up the data code for plotting the confidence bounds then estimate the correct adjustment for the 
# - growth rate based on the daily data
# - next, setup growth rate

rm(list=ls(all=TRUE))
options(error=NULL)

source("./funcs.r")

ind_data  <- LoadAndCleanHKSeroInd("/Volumes/NO NAME/data/influenza/hk_serosurvey/fu_attendence_sr_noident.csv")
symp_data <- LoadAndCleanHKSeroSymp("/Volumes/NO NAME/data/influenza/hk_serosurvey/fu_symptoms_sr_noident.csv")

# write.csv(symp_data,"~/Dropbox/tmp/tmp.csv")

# Set up some basic parameters
no_people <- dim(ind_data)[1]
no_reports <- dim(symp_data)[1]

# Set up the time scale
timesc.fday 		<- as.date("06Jul2009")
timesc.stwkno		<- 28
timesc.nowks		<- 30
timesc.vecstart 	<- 0:(timesc.nowks-1)*7+timesc.fday
number_age_classes	<- 5

# The lines below generate a matrix of age group and week of report for the symptom data
raw_inc <- matrix(data=0,nrow=timesc.nowks,ncol=number_age_classes)
for (i in 1:no_reports) {
	dayrep 		<- symp_data$onset_d[i]
	hhrep 		<- symp_data$household_id[i]
	indrep 		<- symp_data$subject_id[i]
	index_ind	<- LookUpMain(hhrep,indrep,ind_data)
	agrep 		<- ind_data$AG[index_ind]
	weekrep		<- round((dayrep-timesc.fday)/7)
	raw_inc[weekrep,indrep] <- raw_inc[weekrep,agrep] + 1
}

# Generate the cohort size matrix

# Generate some runs form a Hong Kong - like model
