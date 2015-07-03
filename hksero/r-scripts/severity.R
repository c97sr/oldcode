# Script to integrate study data with HK severity data
# You would usually run this one section at a time
# All outputs go to ../outputs/
# If needed: setwd("~/Dropbox/svneclipse/hksero/r-scripts")
# or a similar command for your machine, the script is designed to 
# be run either interactively or from via "source" from the 
# directory in which it is placed

# Next
# - Remove the _dummy from the ICU matrices

rm(list=ls(all=TRUE))
# source("../../H1N1pdm/funcs.R")
# source("funcs_hksero.R")
source("../../idsource/R/stevensRfunctions.R")
require("date")
require("Design")
require("Hmisc")
options(warn=0)
options(error=NULL)

vecArgs <- commandArgs(trailingOnly=TRUE)

if (length(vecArgs) == 1 && vecArgs[1] == "MAKEALL") {
	makeall = TRUE
} else {makeall = FALSE}

remakestudy =		FALSE
remakeeflu =		FALSE
remakefigure =		FALSE
remakegrid =		FALSE
remakesevtab =		FALSE
plotsimprotocol =	FALSE
remakemultiplot =	FALSE
remakeserotable = 	FALSE
remakeregression = 	TRUE
remakespline =		FALSE
makeplottiming = 	FALSE
makesymptomtable = 	FALSE
datashares = 		TRUE

# Load the study data (as needed)
if (makeall || remakestudy) {
	
	# Routine to read in the raw data files from the lab and field work and then make the objects used for the analysis
	base_raw 	<- read.csv("../anon_data/base_ind.csv")
	labs_raw1	<- read.csv("../anon_data/lab_results1.csv")
	labs_raw2 	<- read.csv("../anon_data/lab_results2.csv")
	symp_raw 	<- read.csv("../anon_data/phone_symp.csv")
	diar_raw 	<- read.csv("../anon_data/diary.csv")
	ques_raw 	<- read.csv("../anon_data/questionnaire_symptoms.csv")
	recruit_raw <- read.csv("../anon_data/subj_recruit.csv")
	pp_data 	<- post_proc_sero(base_raw, labs_raw1, labs_raw2, symp_raw,diar_raw,ques_raw,recruit_raw)
	save(pp_data,file="../anon_data/tmp_pp_data.Rdata")
	write.csv(pp_data$base,file="pp_data_base.csv")
	
} else {
	
	# If the objects are not to be remade, load from file
	load("../anon_data/tmp_pp_data.Rdata")
	
}

# Load the eflu severity data
if (makeall || remakeeflu) {
	
	# Comb through the eflu data and make the subset needed for the analysis
	hosp.data_new <- procEFlu("/Volumes/NO NAME/data/influenza/eflu/totalcases.coded.201007-08.csv")
	hosp.data_old <- procEFlu("/Volumes/NO NAME/data/influenza/eflu/totalcases.coded(2010-02-08).unique.csv")
	save(hosp.data_new,hosp.data_old,file="../anon_data/tmp_hosp_data.Rdata")

} else {
	
	# If the objects are not to be remade, laod from file
	load("../anon_data/tmp_hosp_data.Rdata")

}

# Make the extra object combining info from the two sources
# Generate age specific vectors of frequency for full period

# Number of simulated studies 10 for the initial paper
noSims <- 10

# Use 20 year age gaps
agstudy <- "ag20"

# Values for q1 and q2 from Miller et al
q1 <- 0.403871
q2 <- 0.9419355

# Scaling factors for units
stdfactors <- c(100,100000,100000)
hosp.data <- hosp.data_new

# Define the weeks of the study by their first and last days, and other aux variablers
stFirstWeek <- as.date("27Apr2009")							# monday of week of first onset
dtLastAdmission <- as.date(max(hosp.data$first.onset))
noWeeks <- ceiling((dtLastAdmission - stFirstWeek)/7)
wkBins <- stFirstWeek + (0:(noWeeks+1))*7
agcats <- names(table(hosp.data[,agstudy]))
noagegroups <- length(agcats)

# Define matrices for the severe events with age group rows and week columns
matHOS <- matrix(0,nrow = noagegroups,ncol=noWeeks+1)
matICU <- matrix(0,nrow = noagegroups,ncol=noWeeks+1)
matDTH <- matrix(0,nrow = noagegroups,ncol=noWeeks+1)
row.names(matHOS) <- agcats
row.names(matICU) <- agcats
row.names(matDTH) <- agcats

# Populate the matrices using age group masks and the R histogram function
strForTiming <- "first.onset"
offsetForTiming <- 0
# strForTiming <- "adm"
# offsetForTiming <- 3
for (cat in agcats) {
	catmask <- hosp.data[, agstudy]==cat
	matHOS[cat,] <- (hist(hosp.data[catmask,strForTiming]-offsetForTiming,breaks=wkBins,plot=FALSE))$counts
	ICUmask <- hosp.data[, agstudy]==cat & hosp.data$sev 
	matICU[cat,] <- (hist(hosp.data[ICUmask,strForTiming]-offsetForTiming,breaks=wkBins,plot=FALSE))$counts
	DTHmask <- hosp.data[, agstudy]==cat & hosp.data$dth 
	matDTH[cat,] <- (hist(hosp.data[DTHmask,strForTiming]-offsetForTiming,breaks=wkBins,plot=FALSE))$counts
}

# Below used occasionally for debugging
matICU_dummy <- matICU
matICU_dummy[,] <- 0
matICU_dummy[,25] <- rowSums(matICU)

# Getting the +ve s out of the principal database plus some other housekeeping stuff
base 		<- pp_data$base
basetested 	<- base[base$final_fourfold %in% c(0,1),]
notest 		<- sum(table(base$final_fourfold)[c("0","1")])
vecColors 	<- c("red","green","blue","magenta")	
agvector 	<- as.vector(c(1119*1000, 2105*1000, 2367*1000, 1260*1000))

# Remake the simulated study
if (makeall || remakesevtab) {
		
	# Make a simulted version of the data using the study protocol implied by sim_study.csv
	lsSims <- list()
	for (i in 1:noSims) {
		simData <- sim_study("../anon_data/sim_study_revised.csv",
			c(7.983119605, 7.687006357, 14.70931798, 22.9364505)/100000,
			agvector,
			matICU, 
			wkBins, 
			q1=q1, 
			q2=q2)
		lsSims <- c(lsSims,list(simData))
	}
	
	sum(simData$final_fourfold[simData$ag20==1])
	
	# table((lsSims[[1]])$ag20,(lsSims[[1]])$final_fourfold)
	save(lsSims,file="../anon_data/tmp_sim_data.Rdata")
	
} else {
	
	load("../anon_data/tmp_sim_data.Rdata")

}

if (makeall || remakefigure) {
	
	# Make the figure for the 

	# Set up the figure
	xaxeA <- c(as.date("1May2009"),as.date("31Jan2010"))
	ymaxA <- 800
	
	# Setup a multipart figure to show the study data and the HK hosptialization data
	# Prepare a different version of this for what _could_ have been done
	# Keep these in landscape for slide presentations# Figure not good enough. Can use illustrator from here.
	
	height 	<- 10
	width 	<- 7.6
	relh	<- 2.5
	xg1 <- 1.0 / width
	xg2 <- 0.1 / width
	yg1 <- 1.5 / height
	yg2 <- 0.0 / height
	yg3 <- 0.1 / height
	
	xf <-   1 - xg1 - xg2
	yf <- ( 1 - yg1 - yg2 - yg3 ) / ( relh + 1 )
	
	posA <- c(xg1,1-xg2,yg1+relh*yf+yg2,1-yg3)
	posB <- c(xg1,1-xg2,yg1,yg1+relh*yf)
	
	# head(base)
	# noparts <- dim(base)[1]
	
	# setup some auxilliary things for part B of the figure
	xdates <- as.date(c("1May2009","1Jun2009","1Jul2009","1Aug2009","1Sep2009","1Oct2009","1Nov2009","1Dec2009","1Jan2010","1Feb2010"))
	xdatelab <- c("1 May","","1 Jul","","1 Sep","","1 Nov","","1 Jan","")
	line_counter <- 1

	if (plotsimprotocol) base_tmp <- simData
	else base_tmp <- basetested	
	
	nobase <- dim(base_tmp)[1]
	
	# Set up the pdf file
	pdf("../outputs/model_severity.pdf",height=height/cm(1),width=width/cm(1))
	
	# Standard margin zeroing and setting axis label formats
	par(mai=(c(0,0,0,0)))
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	
	# Make the initial plot
	par(fig=posB)
	plot(1:2,type="n",axes=FALSE,xlim=xaxeA,ylim=c(0,ymaxA),xlab="Time",ylab="Incidence of hospitalized cases by week of onset")
	axis(1,at=xdates,labels=xdatelab,las=2,cex=10/12)
	axis(2,cex=10/12,las=1)
	plotDiscInc(matHOS,wkBins,vecColors)
	
	par(fig=posA,new=TRUE)
	plot(1:2,type="n",axes=FALSE,xlim=xaxeA,ylim=c(0,1),ylab="")
	for (ag in c(1,2,3,4)) {
		for (i in 1:nobase) {
			if (base_tmp[i,agstudy]==ag) {
					points(
						c(base_tmp[i,"base_ind_date"],base_tmp[i,"fu_ind_date"]),
						c(line_counter/notest,line_counter/notest),
						type="l",col=vecColors[ag],lwd=0.1)
						line_counter <- line_counter+1

			}
		}
	}
	
	# Close the pdf file
	dev.off()

}

# Make the scatter plot of all the different results
if (makeall || remakegrid) {
	
	# aux stuff for the chart
	axticks <- 0:8
	axlabs <- c("","10","20","40","80","160","320","640","")
	height 	<- 15
	width 	<- 15
	rscale 	<- 1/6
	base 	<- 10
	
	# Dummy data
	pdf("../outputs/full_results.pdf",height=height/cm(1),width=width/cm(1))
	par(mgp=c(2,0.4,0))
	par(mai=(c(0,0,0,0)))
	par(fig=c(0.1,0.95,0.1,0.95))

	# Set up the chart
	plot(1:2,type="n",axes=FALSE,xlim=c(0,8),ylim=c(0,8),xlab="Baseline",ylab="Follow-up")
	axis(1,at=axticks,labels=axlabs)
	axis(2,at=axticks,labels=axlabs)
	abline(h=1:8,v=1:8,col="grey",lwd=0.5)
	abline(a=1.5,b=1,col="grey",lwd=2)
	
	# Plot the pies
	tabage <- table(basetested$pre_titre_agg,basetested$post_titre_agg,basetested[,agstudy])
	for (bt in 1:7) for (ft in 1:6) srPieFunc(bt,ft,tabage[bt,ft,],base=base,rscale=rscale)

	# Add cicles for legend
	polygon(c(5,7,7,5),c(0.5,0.5,3.5,3.5),col="white",border=NULL)
	srPieFunc(6,3,c(1,0,0,0),base=base,rscale=rscale)
	srPieFunc(6,2,c(3,3,2,2),base=base,rscale=rscale)
	srPieFunc(6,1,c(25,25,25,25),base=base,rscale=rscale)

	dev.off()

}

# Estimate the rates based on the actual data and the simulated data and then outputs
# a severity table
# Need to check that it works well on the
if (makeall || remakesevtab) {
	
	# The indexes are 1:4 for the age groups and then -1 for the total 
	agindexes <- c(1,2,3,4,-1)
	# agindexes <- c(1)
	# agindexes <- c(-1)

	# Estimation results for point estimate/ub/lb, age group and simulation number (-1 for actual data)
	simRes <- array(-1,dim=c(3,length(agindexes),noSims),dimnames=c("est_lb_ub","age_group","sim"))
	
	# Which datasets should be run? -1 for the actual data and then an index for reach possible simulation
	seqToRun <- c(-1,1:noSims)
	# seqToRun <- c(-1)
	# seqToRun <- c(1:noSims)

	# Cycle thourgh all the selected datasets
	for (r in seqToRun) {
		
		# If simulated data then set up a boolean
		if (r < 0) {	
			fitsimulation <- FALSE
		} else {
			fitsimulation <- TRUE
		}
		
		# Use check for use of simulated data
		if (fitsimulation) {
			
			# Assign the study data, the outcomes, the multiplier factors and the number of outcomes
			studData <- lsSims[[r]]
			poss_outcomes <- list(matICU)
			multfactors <- c(100000)
			noouts <- length(poss_outcomes)
		
		# Else use the actual study data 
		} else {
			
			studData <- basetested
			poss_outcomes <- list(matHOS,matICU,matDTH)
			# poss_outcomes <- list(matICU)
			multfactors <- stdfactors
			noouts <- length(poss_outcomes)
			nullcol <- rep(-999,noouts*length(agindexes))
			sevtable <- as.data.frame(cbind(outcome=nullcol,agegroup=nullcol,
							rateestimate=nullcol,ratelb=nullcol,rateub=nullcol,
							extestimate=nullcol,extlb=nullcol,extub=nullcol
					))
			
		# End else for the use of actual study data	
		}
	
		# This relies on numeric age cats with increasing age for increasing number
		# Numbers of individuals in each of the age categoriea 
		studyvector <- as.vector(table(studData[,agstudy]))
	
		# Normalize the population vector
		popsize <- sum(agvector)
		studysize <- sum(studyvector)
		# Try this the other way around
		# Do we need to adjust for the number of positives in each age group?
		# For the average value to be correct, we would need to have sampled the population in the 
		# same ratio as the number of infections, not the same ratio as the size of the population
		# Construct a posvec at this point
		# posvec <- c(44/(68+44),13/(13+133),20/381,1/131)
		# browser()
		# posvec <- (table(studData$final_fourfold,studData$ag20))[2,] / studyvector
		adjvec <- (studyvector / studysize) / (agvector / sum(agvector))
		
		# Start for loop for the number of outcomes
		for (o in 1:noouts) {
		
			# Assign the first outcome
			outcome <- poss_outcomes[[o]]
			multfactor <- multfactors[o]
		
			# Calculate the total number of individuals in each age group in the data
			# Adjust the results of the study as if the same means had been obtained with the 
			# same overall size, but the study had been perfectly representative
			# This is needed to estimate overall attack rates
			outcome_adj <- outcome
			outcome_adj[,] <- -1
			for (i in 1:4) outcome_adj[i,] <- outcome[i,] * adjvec[i]
			# outcome_adj <- outcome_adj * sum(outcome) / sum(outcome_adj)
		
			# Loop around each possible age class
			for (ai in 1:(length(agindexes))) {
				
				# Get the correct age class for the 
				ageclass <- agindexes[ai]
			
				# For actual age classes, i.e. not for overall estimates
				if (ageclass > 0) {
					
					# Extract the vectors required for the likelihood estimation
					# Vector of cases per week
					ya <- as.vector(outcome[ageclass,])
					
					# Study data for the correct agegroup
					base_at <- studData[studData[,agstudy] == ageclass,]
					
					# The size of the population at large in this age group
					size_group <- agvector[ageclass]
					
					# The totsl number of severe events (used to doublecheck the calculation on dummy data)
					totalevents <- sum(outcome[ageclass,])
					vec_b_aux <- studData[studData[,agstudy]==ai,"base_ind_date"]
					vec_f_aux <- studData[studData[,agstudy]==ai,"fu_ind_date"]
				
				
				} else {
					
					# Make up an overall list of the number of events
					ya <- as.vector(colSums(outcome_adj))
					
					# Set the sero data to be all the age group data
					# This could well be where the error is
					base_at <- studData
					
					# Define the size of the population under consideration
					size_group <- sum(agvector)
					
					# Add up the total number of severe events
					totalevents <- sum(outcome_adj)	
					
					# Make the vector of baseline and followup sample weeks for the study 
					# participants
					vec_b_aux <- studData[,"base_ind_date"]
					vec_f_aux <- studData[,"fu_ind_date"]
					
				}
				
				vec_b <- ceiling((vec_b_aux - stFirstWeek)/7)
				vec_f <- ceiling((vec_f_aux - stFirstWeek)/7)
				xa <- as.vector(base_at$final_fourfold)
				
				# Some debug lines
				# f <- function(x){sev_like(x,q1=q1,q2=q2,xa,vec_b,vec_f,ya,agvector[ageclass],verbose=FALSE)}
				# xtmp <- as.vector(10^((-15:-1)/3))
				# plot(xtmp,-vapply(xtmp,f,c(1)),log="xy")

				# Start thinking about debugging here
				# You might want to evaluate sum(ya)/(sum(xa)/length(xa)* size_group)
				
				# Sewup the linear interpolation
				lb <- 0.00000001
				opt_pest <- optimize(
						f=sev_like,interval=c(lb,1.0-lb),q1=q1,q2=q2,xa=xa,
						ba=vec_b,fa=vec_f,ya=ya,size_group,tol=1e-10,maximum=TRUE,verbose=FALSE)
				p_est <- opt_pest$maximum
				p_est_lnlike <- opt_pest$objective
				
				# sev_like_ci(p_est,-5.4,0.5,0.5,xa,vec_b,vec_f,ya,7000000/4)
				opt_lb <- optimize(sev_like_ci,c(lb,p_est),target=p_est_lnlike,q1=q1,q2=q2,xa=xa,
						ba=vec_b,fa=vec_f,ya=ya,size_group,maximum=FALSE,tol=1e-10)
				opt_ub <- optimize(sev_like_ci,c(p_est,1.0),target=p_est_lnlike,q1=q1,q2=q2,xa=xa,
						ba=vec_b,fa=vec_f,ya=ya,size_group,maximum=FALSE,tol=1e-10)
				
				if (fitsimulation) {
					simRes[1,ai,r] <- p_est
					simRes[2,ai,r] <- opt_lb$minimum
					simRes[3,ai,r] <- opt_ub$minimum
				} else {					
					sevtable[(o-1)*length(agindexes)+ai,1:5] <- c(o,ageclass,multfactor*p_est,multfactor*opt_lb$minimum,multfactor*opt_ub$minimum)
					sevtable[(o-1)*length(agindexes)+ai,6:8] <- c(totalevents/p_est/size_group,totalevents/opt_ub$minimum/size_group,totalevents/opt_lb$minimum/size_group)
				}
				
				cat("r ",r," o ",o," g ",ageclass," ",multfactor*p_est," (",multfactor*opt_lb$minimum,",",multfactor*opt_ub$minimum,")\n")

				# if (ageclass < 0) browser()
				# posvec <- c(44/(68+44),13/(13+133),20/381,1/131)
				# ratevec <- c(0.01137379,0.01120409,0.02173004,0.1143651)
				# sum(ratevec*(posvec*agvector/sum(posvec*agvector)))
				# sum(outcome[4,])/((1/131)*agvector[4])*100
				
				flush.console()
				
			}
			
		}
			
	}
	
	# Write output
	if (seqToRun[1]< 0) write.csv(sevtable,file="../outputs/tab_severity.csv",row.names=FALSE)
	save(simRes,file="tmp_simRes.Rdata")
	
	# Report the number of severe cases used
	cat("HOS",sum(matHOS),"ICU",sum(matICU),"DTH",sum(matDTH),"\n")
	
} else {
	load("tmp_simRes.Rdata")
}

if (makeall || remakemultiplot) {

	# Routine to put the three line plots one above the other on the same page

	# Load up the required data
	tab1 <- read.csv("../outputs/tab_severity.csv")
	
	# Define some plotting variables
	sevindexoff <- c(0,5,10)
	sevplotoff <- c(-0.15,0,0.15)
	sevpch <- c(22,24,19) 
	groupindexoff <- c(5,1,2,3,4)
	nosevs <- length(sevindexoff)
	nogrps <- length(groupindexoff)
	simoff <- (0:noSims)/noSims*0.8
	simoff <- simoff - simoff[noSims %/% 2 + 1]

	# Set up the geometry
	height 	<- 16
	width 	<- 14-2.6
	xg1 <- 1.5 / width
	xg2 <- 0.1 / width
	yg1 <- 1.0 / height
	yg2 <- 0.2 / height
	yg3 <- 0.2 / height
	yg4 <- 0.2 / height
	yr1	<- 1
	yr2 <- 1
	yr3 <- 1
	xf <-   1 - xg1 - xg2
	yf <- (1 - yg1 - yg2 - yg3 - yg4 ) / ( yr1 + yr2 + yr3 )
	posA <- c(xg1,1-xg2,yg1+yf*yr1+yg2+yf*yr2+yg3,1-yg4)
	posB <- c(xg1,1-xg2,yg1+yf*yr1+yg2,yg1+yf*yr1+yg2+yf*yr2)
	posC <- c(xg1,1-xg2,yg1,yg1+yf*yr1)
	
	# Axes ranges and labels
	xrange 	<- c(0.5,5.5)
	xat1 	<- 1:5
	xlabs1	<- c("Overall","3-19","20-39","40-59","60 and over")
	xat2 	<- 1:6 - 0.5
	xlabs2	<- c("","","","","","")
	xat3 <- 2:5 - 0.5
	
	y1vals <- c(1e-6,1e-4,1e-2,1)
	y1labs <- c(expression(10^{-6}),expression(10^{-4}),0.01,1)
	y2vals <- (0:10)/10
	y2labs <- c("0%","","20%","","40%","","60%","","80%","","100%")
	y3vals <- c(1e-6,1e-4,1e-2,1)
	y3labs <- c(expression(10^{-6}),expression(10^{-4}),0.01,1)

	y1range <- c(min(y1vals),max(y1vals))
	y2range <- c(min(y2vals),max(y2vals))
	y3range <- c(min(y3vals),max(y3vals))
	
	pdf("../outputs/rates_real_and_sim.pdf",height=height/cm(1),width=width/cm(1))
	
	# Standard parameter calls
	par(mai=(c(0,0,0,0)))
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	
	# Make the initial plot
	par(fig=posA)
	plot(1:2,axes=FALSE,xlim=xrange,ylim=y1range,log="y",type="n")
	
	axis(2,at=y1vals,labels=y1labs,las=1)
	text(xrange[1],y1range[2],"A",adj=c(1,1))
	abline(v=xat3,lty=2)

	# Put the data opn the chart
	for (i in 1:nosevs) {

		for (j in 1:nogrps) {

			points(j+sevplotoff[i],tab1$rateestimate[groupindexoff[j]+(i-1)*5]/stdfactors[i],pch=sevpch[i],col="black")
			points(c(j+sevplotoff[i],j+sevplotoff[i]),c(tab1$ratelb[groupindexoff[j]+(i-1)*5]/stdfactors[i],tab1$rateub[groupindexoff[j]+(i-1)*5]/stdfactors[i]),type="l",col="black")

		}

	}

	par(fig=posB,new=TRUE)
	plot(1:2,axes=FALSE,xlim=xrange,ylim=y2range,type="n")
	axis(2,at=y2vals,labels=y2labs,las=1)
	text(xrange[1],y1range[2],"B",adj=c(1,1))
	abline(v=xat3,lty=2)

	# Put the data opn the chart
	for (i in 1:nosevs) {

		for (j in 1:nogrps) {

			points(j+sevplotoff[i],tab1$extestimate[groupindexoff[j]+(i-1)*5],pch=sevpch[i],col="black")
			points(c(j+sevplotoff[i],j+sevplotoff[i]),c(tab1$extlb[groupindexoff[j]+(i-1)*5],tab1$extub[groupindexoff[j]+(i-1)*5]),type="l",col="black")

		}

	}
	
	par(fig=posC,new=TRUE)
	plot(1:2,axes=FALSE,xlim=xrange,ylim=y1range,log="y",type="n")
	axis(2,at=y3vals,labels=y3labs,las=1)
	axis(1,at=xat1,labels=xlabs1,lty=0,line=0)
	axis(1,at=xat2,labels=xlabs2,lty=1)
	text(xrange[1],y3range[2],"C",adj=c(1,1))
	abline(v=xat3,lty=2)

	# Put the data actual data the chart
	for (i in c(2)) {
		for (j in 1:nogrps) {
			points(j+simoff[1],tab1$rateestimate[groupindexoff[j]+(i-1)*5]/stdfactors[i],pch=sevpch[i],col="black")
			points(c(j+simoff[1],j+simoff[1]),c(tab1$ratelb[groupindexoff[j]+(i-1)*5]/stdfactors[i],tab1$rateub[groupindexoff[j]+(i-1)*5]/stdfactors[i]),type="l",col="black")
		}
	}

	# Put the data opn the chart
	for (i in 2:length(simoff)) {
		for (j in 1:nogrps) {
			points(j+simoff[i],simRes[1,groupindexoff[j],i-1],pch=sevpch[2],col="grey")
			points(c(j+simoff[i],j+simoff[i]),c(simRes[2,groupindexoff[j],i-1],simRes[3,groupindexoff[j],i-1]),type="l",col="grey")
		}
	}
		
	dev.off()	

}

if (makeall || remakeserotable) {
	x <- ftable(basetested$pre_titre_agg>20,basetested$final_fourfold,basetested$ag20,col.vars=c(2),row.vars=c(3,1))
	write.csv(x,"../outputs/table_base_titres_age_seroconvert.csv")
}

if (makeall || remakeregression) {
	
	# Needs a pretty diagram to show the estimates and confidence intervals for the different
	# age-specific models
	
	# Need a figure with the spline fit and the rolling average, nicely formatted to go into the main text
	# just the raw data and the spline fit with confidence bounds
	
	# Regroup district
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
	
	for (i in 1:length(base$smoking)){
		if (base$smoking[i]!="Yes"){base$smoking[i] <-"No"}
	}	
		
	# reset houldhold groups
	base  <- cbind(base,hh_size2=rep(-1,dim(base)[1]))
	
	for (i in 1:dim(base)[1]){
		if (is.na(base[i,"hh_size"])) base$hh_size2[i] <- -1
		else if (base[i,"hh_size"] ==0) base[i,"hh_size2"] <- -1
		else if (base[i,"hh_size"] ==1) base[i,"hh_size2"] <- 1
		else if (base[i,"hh_size"] ==2) base[i,"hh_size2"] <- 2
		else if (base[i,"hh_size"] ==3) base[i,"hh_size2"] <- 3
		else if(base[i,"hh_size"] ==4) base[i,"hh_size2"] <- 4
		else if(base[i,"hh_size"] ==5) base[i,"hh_size2"] <- 5
		else if(base[i,"hh_size"] ==6) base[i,"hh_size2"] <- 6
		else if(base[i,"hh_size"] ==7) base[i,"hh_size2"] <- 6
		else if(base[i,"hh_size"] ==8) base[i,"hh_size2"] <- 6
		else if(base[i,"hh_size"] ==9) base[i,"hh_size2"] <- 6
		else stop("Problem with the hh sizes groups")
	}       
	
	base  <- cbind(base,vaccine.0809=rep(-1,dim(base)[1]))

	# Vaccinated with seasonal influenza vaccine during year 2008 and 2009 for all subjects, not adjusted for other variables
	for (i in 1:length(base$vaccine.0809)) {
		if (base$vac0809[i]=="yes" || base$vac0809[i]=="Yes")	{base$vaccine.0809[i] <- "Yes"} else
		if (base$vac0809[i]=="no" || base$vac0809[i]=="No")	{base$vaccine.0809[i] <- "No"} else {base$vaccine.0809[i] <- "NK"}
	}
		
	t_base_1 <- base[base[,"final_fourfold" ] > -1,]

	# Establish an "exposure" variable for low-baseline titre
	t_base_1$low_pre_titre <- (t_base_1$pre_titre_agg < 40)	
	t_base_1$exp.data.comp <- (t_base_1[,"education"]!=0) & (t_base_1[,"ind_id"] != "S090189-1")
	t_base <- t_base_1[t_base_1[,"exp.data.comp"]==1,]
			
	# Automate model fitting to be able to test all possible models
	# setup a vector of model terms
	vecTerms <- c(
		"factor(t_base$sex)",
		"factor(t_base$education2)",
		"factor(t_base$profession2)",
		"factor(t_base$hh_size2)",
		"factor(t_base$district5)",
		# "factor(t_base$low_pre_titre)",
		"factor(t_base$vaccine.0809)"
	)
		
	# To generate the set of all models, define some aux variables
	noModels <- 0
	AICBaseline <- 400
	AICThreshold <- 1000
	noTerms <- length(vecTerms)
	sysListMods <- list()
	
	# Open loop for all model sizes
	for (i in 1:noTerms) {
		
		# Make the list of parameter indices
		possMods <- combn(noTerms,i)
		noMods <- dim(possMods)[2]
		
		# Open loop for all model combinations of this size
		for (j in 1:noMods) {
			
			# Solve the regression model
			cur_mod <- glm(
				family = binomial(logit),
				as.formula(
					paste("t_base$final_fourfold ~ t_base$age + t_base$child_hh + ", 
						paste(vecTerms[possMods[,j]],collapse=" + "),
					collapse="")
				)
			)
			
			# Test conditions for AIC
			if ((cur_mod$aic - AICBaseline) < AICThreshold) sysListMods <- c(sysListMods,list(cur_mod))
			noModels <- noModels + 1		
		}
		
	}	
	
	tmp <- srg.summarise.glm(sysListMods,file="../outputs/sysmodsum.csv",transpose=TRUE)

	# Models to present or mention in the main text
	# Run regressions based on variables of interest
	tablemodels <- list(
			
			# Univariate models for Table
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							t_base$age
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$child_hh)
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$sex)
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$district5)
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$vaccine.0809)
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$hh_telsource)
			),
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							factor(t_base$low_pre_titre)
			),
			
			# Best model
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							t_base$age +
							factor(t_base$low_pre_titre) +
							factor(t_base$child_hh) 
			),
			
			# Complete "good" model
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							t_base$age +
							factor(t_base$child_hh) +
							factor(t_base$district5) +
							factor(t_base$sex) + 
							factor(t_base$vaccine.0809)
			),

			# Complete "good" model with source of recruitment as well
			glm(family = binomial(logit), 
					formula = t_base$final_fourfold ~ 
							t_base$age +
							factor(t_base$district5) +
							factor(t_base$child_hh) +
							factor(t_base$sex) + 
							factor(t_base$vaccine.0809) +
							factor(t_base$hh_telsource)				
			)
	)
	
	srg.summarise.glm(tablemodels,file="../outputs/tablemodels_2sf.csv",sigdigits=2)
	srg.summarise.glm(tablemodels,file="../outputs/tablemodels_3sf.csv",sigdigits=3)
	
	# Make the small tabels for the table, if needed
	table(t_base$sex)

	# For 5 knots the outer quantiles are 0.05 and 0.95 and the inner knots are equally spaced (I assume)
	# Based on example for rcspline.restate
	# Probably interesting to look at the spline fit after adjusting for the other variables
	# Not interpreting the intercept correctly at the moment
	lineartrans <- exp(0.26946 - 0.06918*t_base$age)
	age20coeffs_tmp <- c(-0.4353,-1.8901,-2.4523,-4.4245)
	age20coeffs <- exp(c(0,age20coeffs_tmp[1]+age20coeffs_tmp[c(2,3,4)]))
	age20trans <- age20coeffs[t_base$ag20]
	
	# XXXX up to here. Need to add in the coefficients for each individual age
	# plot(t_base$age,coeffs[1] + xtrans(t_base$age),type="n")
	# points(t_base$age,log(lineartrans),col="red")
	# points(t_base$age,log(age20trans),col="blue")
	
	# Plot with confidence intervals

	# rcspline.plot(t_base$age,t_base$final_fourfold,
	#	model=c("logistic", "cox", "ols"),xrange,event,nk=5,knots=NULL,
	#	show=c("xbeta","prob"),adj=NULL,xlab,ylab,ylim,plim=c(0,1),plotcl=TRUE,
	#	showknots=TRUE,add=FALSE,subset,lty=1,noprint=FALSE,m,smooth=FALSE,bass=1,
	#	main="auto",statloc)
	
	# Make a function for the average rolling incidence with age and then  super-impose it on the plot
	# Add an additional figure and then alter the regression table
	# Then move on to the custom inference method
	
	# Write the output to a nicely formatted file	

}

# Figure to make a nicely formatted pdf of the fitted spline function and the rolling average data
# Should be a two part figure with a probability scale at the top chart and a log-odds scale on the lower chart
# Do all the ptrettying up in illustrator
if (makeall || remakespline) {
	
	# Make the spline outputs
	plot(1:2)
	unadjsplineprob <- rcspline.plot(t_base$age,t_base$final_fourfold,show="prob",add=TRUE)
	unadjsplinebeta <- rcspline.plot(t_base$age,t_base$final_fourfold,show="xbeta",add=TRUE)
	rollave <- windowInc(t_base$age,t_base$final_fourfold,width=100)
	
	# Set up the knot locations for 5 knot rcs on t_base$age
	ageknts <- quantile(t_base$age,c(0.05,0.9/4+0.05,2*0.9/4+0.05,3*0.9/4+0.05,0.95))
	
	# Make the same spline using rcs and glm
	glmrcsbasic <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~ 
					rcs(t_base$age))

	# Make the same spline using rcs and glm
	glmrcsadj1 <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~ 
					rcs(t_base$age) + 
					factor(t_base$child_hh))	
	
	# Make the same spline using rcs and glm
	glmrcsadj2 <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~ 
					rcs(t_base$age)+ 
					factor(t_base$child_hh) +
					factor(t_base$district5))	
	
	# Make the same spline using rcs and glm
	glmrcsadj3 <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~ 
					rcs(t_base$age) +
					factor(t_base$district5) +
					factor(t_base$child_hh) +
					factor(t_base$sex) + 
					factor(t_base$vaccine.0809))
	
	# Make the same spline using rcs and glm
	glmrcsadj4 <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~  
					rcs(t_base$age) +
					factor(t_base$district5) +
					factor(t_base$child_hh) +
					factor(t_base$sex) + 
					factor(t_base$vaccine.0809) +
					factor(t_base$hh_telsource)	)	
	
	# Make the same spline using rcs and glm
	glmlinearadj <- glm(family = binomial(logit), 
			formula = t_base$final_fourfold ~  
					t_base$age +
					factor(t_base$child_hh))	
	
		
	# Regenerate the spline function
	wbasic <- rcsplineFunction(ageknts,glmrcsbasic$coefficients[1:5])
	wadj1 <- rcsplineFunction(ageknts,glmrcsadj1$coefficients[1:5])
	wadj2 <- rcsplineFunction(ageknts,glmrcsadj2$coefficients[2:5])
	wadj3 <- rcsplineFunction(ageknts,glmrcsadj3$coefficients[2:5])	
	
	# Figure height and width, one-col plos med
	fh <- 8.3
	fw <- 8.3
	
	# Set up the single position to use multiple y-axis
	pos <- c(0.25,1.0,0.2,1.0)
	pos2 <- c(0.35,0.7,0.3,0.7)
	
	# Open the file to write the pdf
	pdf("../outputs/spline.pdf",height=fh/cm(1),width=fw/cm(1))
	
		# Set the standard graphics parameters
		par(mai=(c(0,0,0,0)))
		par(mgp=c(2,0.4,0))
		par(tcl=-0.25)
	
		# Make the initial plot
		par(fig=pos)
		xat <- 10+0:6*10
		yat <- c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0)
		subspline <- unadjsplineprob$x < max(xat) & unadjsplineprob$x > min(xat)
		plot(	unadjsplineprob$x[subspline],
				unadjsplineprob$xbeta[subspline],
				type="n",log="y",
				ylim=c(min(yat),max(yat)),
				xlim=c(min(xat),max(xat)),
				axes=FALSE)
		polygon(	c(unadjsplineprob$x[subspline],rev(unadjsplineprob$x[subspline])),
					c(unadjsplineprob$upper[subspline],rev(unadjsplineprob$lower[subspline])),
					col="grey",border=NA)
		points(rollave)
		points(unadjsplineprob$x[subspline],unadjsplineprob$xbeta[subspline],type="l",col="red")
		axis(1,at=xat,las=1)
		axis(2,at=yat,las=1)
		
		# Plot an inset chart
		par(fig=pos2,new=TRUE)
		chtxran <- 10:70
		plot(chtxran,wbasic(chtxran),type="l",col="red",axes=FALSE,ylim=c(-5,0))
		points(chtxran,wadj1(chtxran),type="l",col="blue")
		points(chtxran,glmlinearadj$coefficients[1] + glmlinearadj$coefficients[2]*chtxran,type="l",col="green")
		# 	points(chtxran,wadj3(chtxran),type="l",col="cyan")	
		axis(1,las=1)
		axis(2,las=1)
						
	# Close the pdf
	dev.off()	
		
}

# Setup routine to make the scatter plot of recruitment and followup
if (makeall || makeplottiming) {
	
	# Setup the output file
	# Figure height and width, one-col plos med
	fh <- 15
	fw <- 15
	
	# Set up the axes	
	x_axis_str <- c("04 Jul 2009","01 Aug 2009","01 Sep 2009","19 Sep 2009")
	y_axis_str <- c("11 Nov 2009","01 Dec 2009","01 Jan 2010","01 Feb 2010","06 Feb 2010")
	x_axis_num <- as.numeric(as.date(x_axis_str))
	y_axis_num <- as.numeric(as.date(y_axis_str))
	
	# Open the file to write the pdf
	pdf("../outputs/scatter_timing.pdf",height=fh/cm(1),width=fw/cm(1))
	
		# Set the position in the frame, with no "default spacing"
		par(mai=(c(0,0,0,0)))
		par(fig=c(0.3,0.95,0.3,0.95))
	
		# Set up the chart
		plot(	1:2,
				xlim=c(min(x_axis_num),max(x_axis_num)),
				ylim=c(min(y_axis_num),max(y_axis_num)),
				axes=FALSE,
				type="n",
				xlab="",
				ylab=""
		)
		
		# Plot the points and the axes 
		points(	jitter(as.numeric(basetested$base_ind_date)),
				jitter(as.numeric(basetested$fu_ind_date)),
				cex=0.5)
		axis(	1,
				labels=x_axis_str,
				at=x_axis_num)
		axis(	2,
				at=y_axis_num,
				labels=y_axis_str,
				las=2)
		
		# Close the pdf device
		dev.off()
		
}

# Make the symptom chart
if (makeall || makesymptomtable) {
	
	# Define the data table for the chart
	dfSymp <- data.frame(	row.names=	c("ILI","ARI","Fever"),
							phone=		c(8,11,11),
							diary=		c(12,29,14),
							followup=	c(30,50,34),
							any=		c(35,57,41))
					
	# Plot the chart
	plot.symptoms(dfSymp,86,ylim_in=c(0,0.8),"../outputs/symptoms.pdf")
	
}

# Routine to construct a well named and exhaustive data set of the serological study data
if (makeall || datashares) {
	
	# Define the output data using the study identifiers
	seroshare <- data.frame(	id=t_base_1$ind_id,										# Individual identifier
								gender=t_base_1$sex,									# Male or female
								age=t_base_1$age,										# Age of particiopant
								date.of.baseline=as.numeric(t_base_1$base_ind_date),	# Date of first sample
								date.of.followup=as.numeric(t_base_1$fu_ind_date),		# Date of followup sample
								child.present=as.numeric(t_base_1$child_hh),			# Presence of a child in the household
								district=t_base_1$district5,							# District
								vaccine.0809=t_base_1$vaccine.0809,
								recruit.source=t_base_1$hh_telsource,
								low.pre.titre=t_base_1$low_pre_titre,
								fourfold.rise=t_base_1$final_fourfold
	)
	
	# Write the file to the output directory
	write.csv(seroshare,file="../outputs/sero_data_riley_and_others_plos_med_2011.csv",row.names=FALSE)

	# Define the output data using the matrices for the different outcomes
	severeshare <- data.frame (	monday.of.week=stFirstWeek+(0:(noWeeks-1))*7,
			hosp.3.19=matHOS["1_3t19",1:(noWeeks)],
			hosp.20.39=matHOS["2_20t39",1:(noWeeks)],
			hosp.40.59=matHOS["3_40t59",1:(noWeeks)],
			hosp.60.plus=matHOS["4_60plus",1:(noWeeks)],
			icu.3.19=matICU["1_3t19",1:(noWeeks)],
			icu.20.39=matICU["2_20t39",1:(noWeeks)],
			icu.40.59=matICU["3_40t59",1:(noWeeks)],
			icu.60.plus=matICU["4_60plus",1:(noWeeks)],
			death.3.19=matDTH["1_3t19",1:(noWeeks)],
			death.20.39=matDTH["2_20t39",1:(noWeeks)],
			death.40.59=matDTH["3_40t59",1:(noWeeks)],
			death.60.plus=matDTH["4_60plus",1:(noWeeks)]
	)
	
	# Write the file to the output directory
	write.csv(severeshare,file="../outputs/severe_data_riley_and_others_plos_med_2011.csv",row.names=FALSE)

}

