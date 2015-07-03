# Script to load, view and then fit the timeseries of cases
# Next 	- enter the presanis ratios and look at a few fits
#		- change the output formatting to be consistent with tracer
#		- run MCMC fits to obtain estimates of key parameters
#  		- add parameter to turn off eigen vector "fitting"

rm(list=ls(all=TRUE))
source("../funcs.R")
require("date")
options(error=NULL)

# Load the raw data file and process comment the three lines below as needed
# hosp.data <- procEFlu("~sriley/Dropbox/projects/influenza/age_peak_attack/anon_data/totalcases.coded(2010-02-08).csv")
# save(hosp.data,file="tmp_hosp_data.Rdata")
load("tmp_hosp_data.Rdata")

# Some useful aux quantities
ind_sev <- hosp.data$sev == TRUE

# Some preliminary analyses
cat("Date of first admission",date.ddmmmyy(as.date(min(hosp.data$adm))),"\n")
cat("Date of last admission",date.ddmmmyy(as.date(max(hosp.data$adm))),"\n")

# Generate age specific vectors of frequency for full period
stFirstWeek <- as.date("27Apr2009")		# monday of week of first admission
dtLastAdmission <- as.date("08Feb2010")	# day of last admission
noWeeks <- ceiling((dtLastAdmission - stFirstWeek)/7)
wkBins <- stFirstWeek + (0:(noWeeks+1))*7

# Make a matrix for the data for the full period
agcats <- names(table(hosp.data$ag1))
noagegroups <- length(agcats)
agdata <- matrix(0,nrow = noagegroups,ncol=noWeeks+1)
row.names(agdata) <- agcats
for (cat in agcats) {
	catmask <- hosp.data[,"ag1"]==cat
	agdata[cat,] <- (hist(hosp.data[catmask,"adm"],breaks=wkBins,plot=FALSE))$counts
}

# Extract mid-section for fitting
noWeeksSkipStart <- 3
noWeeksSkipEnd <- 3
noFittedWeeks <- noWeeks - noWeeksSkipStart - noWeeksSkipEnd
agdata_fit <- matrix(0,nrow = noagegroups,ncol=noWeeks+1-noWeeksSkipStart - noWeeksSkipEnd)
row.names(agdata_fit) <- agcats
wkBins_fit <- stFirstWeek + (noWeeksSkipStart:(noWeeks-noWeeksSkipEnd+1))*7
agdata_fit[,] <- agdata[,(noWeeksSkipStart+1):(noWeeks-noWeeksSkipEnd+1)]	

save(agdata_fit,wkBins_fit,agcats,file="data_bins_cats_fit.Rdata")
