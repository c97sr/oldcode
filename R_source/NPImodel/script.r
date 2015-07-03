# Script to run the MCMC sampler for NPI study
# - make up 100 houses of 4 as 2 example sets. One with complete data and one with partial data
# - check that you can recover q_base, q_int_1 and q_int_2 from sim data
# - find out why simulated data set falls over occasionally with the debug check every 50

rm(list=ls())
options(error=NULL)
vecArgs <- commandArgs(trailingOnly=TRUE)

# Set broad parameters for machine, data and simulation mode, and seed
if (length(vecArgs) != 10) {
	runtype <- "local" # local | server
 	dirData <- "20_h_of_2" # 2008_full | 20_h_of_2 | 50_h_of_2 | 1_h_of_10 | 20_h_of_2_poor_comp | 2008_no_null_index
	simulate <- TRUE # TRUE | FALSE
	set.seed(-38273910)
	strOutstem <- "sample"
	strRps <- "rps.txt"
	strBps <- "bps.txt"
	strBpsRange <- "bps_ranges.txt"
	strSnapshotName <- "NULL"
	stop(paste(" ",vecArgs))
} else {
	runtype <- vecArgs[2]
	dirData <- vecArgs[3]
	if (vecArgs[4]=="TRUE") simulate <- TRUE
	else simulate <- FALSE
	set.seed(-1*as.numeric(vecArgs[5]))
	strOutstem <- vecArgs[6]
	strRps <- vecArgs[7]
	strBps <- vecArgs[8]
	strBpsRange <- vecArgs[9]
	strSnapshotName <- vecArgs[10]
}

# Set up directories and read in source, data param files
if (runtype=="local") {
	dirRCodeRoot <- "./"
	dirWorking <- "./NPImodel/working/"
} else if (runtype=="server") {
	dirRCodeRoot <- "../../r-project/"
	dirWorking <- "./"
} else stop("Unknown runtype.")
source(paste(dirRCodeRoot,"SR_misc.r",sep=""))
source(paste(dirRCodeRoot,"NPImodel/time_since_functions.r",sep=""))
fnhc 		<- paste(dirRCodeRoot,"NPImodel/data/",dirData,"/tests.csv",sep="")
fnhh 		<- paste(dirRCodeRoot,"NPImodel/data/",dirData,"/hchar.csv",sep="")
fnic 		<- paste(dirRCodeRoot,"NPImodel/data/",dirData,"/ichar.csv",sep="")
fnbiops 	<- paste(dirWorking,strBps,sep="")
fnrunps 	<- paste(dirWorking,strRps,sep="")
fnptab 		<- paste(dirWorking,strBpsRange,sep="")
fnMcmcTrace <- paste(dirWorking,strOutstem,sep="")
studyData 	<- read_study_data(fnhc,fnhh,fnic)
dataUsed 	<- studyData$indDb
runPs 		<- fnReadFileAsParamVector(fnrunps)
bioPs 		<- fnReadFileAsParamVector(fnbiops)
tabParams 	<- read.table(fnptab)
hhlu 		<- fnMakeHHLookup(studyData$hh)

# Simulate if needed and set initial infection times and infection states
if (simulate) dataUsed <- simStudy(bioPs,runPs,studyData$hhDb,studyData$indDb)	
dataUsed1 <- fnAssignInitialTimeValues(dataUsed,runPs)
dataUsed2 <- fnAssignInitialInfectionStates(dataUsed1)

write.csv(studyData$indDb,file="tmp.csv")

fnRunMcmc(	bioPs,runPs,tabParams,dataUsed2,hhlu,
			fnOutputStem=fnMcmcTrace, fnSnapshot=strSnapshotName,
			db_like=TRUE, db_freq=1000)
