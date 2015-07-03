bioPs$q_ptp <- 0.3
dataUsed <- simStudy(bioPs,runPs,studyData$hhDb,studyData$indDb)	

vals <- 0.8+1:100/100*0.3
likes <- fnMargLikeIndexInfs(bioPs,runPs,dataUsed,vals=vals,
		t_factor=1,v_factor=1,s_factor=1)
plot(vals,likes,type="l")

plot(vals,fnMargLike("q_ptp",vals,bioPs,runPs,dataUsed),type="l")

# Current objective - generate simulated data in exactly the same format as the real data
runPs$totSamples <- 20
biopstmp <- bioPs
tmpseed <- .Random.seed
set.seed(-928340294)
fnRunMcmc(bioPs,runPs,tabParams,dataUsed,fnOutputStem = fnMcmcTrace)
runPs$totSamples <- 10
bioPs <- biopstmp
set.seed(-928340294)
fnRunMcmc(bioPs,runPs,tabParams,dataUsed,fnOutputStem = fnMcmcTraceA)
runPs$totSamples <- 10
fnRunMcmc(bioPs,runPs,tabParams,dataUsed,fnOutputStem = fnMcmcTraceB, 
		fnSnapshot = "./NPImodel/working/mcmc_out_10a_snapshot.out")


rm(list=ls(all=TRUE))
options(error=recover)
# options(error=NULL)

setwd("C:\\eclipse\\R Source")
fnOutput <- "c:\\tmp\\mcmc.out"

source("NPImodel\\time_since_functions.r")
source("SR_misc.r")

fnhc <- "NPImodel\\data\\pilot_home_culture.csv"
fnhh <- "NPImodel\\data\\hchar_h.csv"

source("NPImodel\\local_copy_runPs.r")
source("NPImodel\\local_copy_bioPs.r")

x <- simStudy(runPs,bioPs,fnhh,fnhc)

# Checking runtime efficiency and profiling
system.time(simStudy(runPs,bioPs,studyData$hhDb,studyData$indDb))
system.time(fnRunMcmc(bioPs,runPs,tabParams,x,fnOutputIn = fnMcmcTrace))
proffile <- "c:\\tmp\\profout.txt"
Rprof(filename=proffile)
fnRunMcmc(bioPs,runPs)
Rprof()
summaryRprof(proffile)

x <- readProfileData(proffile)
flatProfile(x)
printProfileCallGraph(x)
summaryRprof(proffile)

# Do repeated simulations and calcuate the attack rate (hh=2 only here)
set.seed(-4567)
notests <- 100
testrates <- array(0,c(notests))
for (i in 1:notests) {	
	dataUsed <- simStudy(bioPs,runPs,studyData$hhDb,studyData$indDb)
	testrates[i] <- fnCalcAttackRate(dataUsed)
}
hist(testrates,breaks= -0.025 + 0.05*(0:20))

fnDebugLike(bioPs,runPs,dataUsed2,fnMakeHHLookup(studyData$hh),t_factor=1,v_factor=1,s_factor=1)

# Section to plot marginal likelihoods
pname <- "p1_symp"
tmin <- 0.1
tmax <- 0.9
nsamp <- 20
tf <- 1
vf <- 1
sf <- 1
xvals <- tmin+0:nsamp/nsamp*(tmax-tmin)
yvals <- fnMargLike(pname,xvals,bioPs,runPs,dataUsed,fnMakeHHLookup(studyData$hh),t_factor=tf,v_factor=vf,s_factor=sf)
yvals <- yvals - max(yvals)
plot(xvals,yvals,type="l")

for (i in 1:5) {
	if (simulate) dataUsed <- simStudy(bioPs,runPs,studyData$hhDb,studyData$indDb)
	yvals <- fnMargLike(pname,xvals,bioPs,runPs,dataUsed,fnMakeHHLookup(studyData$hh),t_factor=tf,v_factor=vf,s_factor=sf)
	yvals <- yvals - max(yvals)
	points(xvals,yvals,type="l")
}



proffile <- "profout.txt"
Rprof(filename=proffile)
for (i in 1:100000) runPs[1,"t_min"] <- runif(1)
Rprof()
summaryRprof(proffile)

