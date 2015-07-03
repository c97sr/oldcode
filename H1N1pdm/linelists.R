 # To dos: 
# - Next : tidy up the data code for plotting the confidence bounds then estimate the correct adjustment for the 
# - growth rate based on the daily data
# - next, setup growth rate

rm(list=ls(all=TRUE))
options(error=NULL)
source("./funcs.r")

# Load up the final database and the version immediately before that one
# Use "table(tabDataX$CaseStatus) to check the categories"
tabData3 	<- read.csv("/Volumes/NO NAME/data/influenza/us_h1n1/Hospitalized_052609.csv")
td3 		<- postProcLineList(tabData3,order="dmy")
tabData2 	<- read.csv("/Volumes/NO NAME/data/influenza/us_h1n1/linelists13may_new.csv")
td2 		<- postProcLineList(tabData2,order="dmy")
tabData1 	<- read.csv("/Volumes/NO NAME/data/influenza/us_h1n1/linelists08may.csv")
td1 		<- postProcLineList(tabData1,order="mdy")

sd 			<- as.date("1Apr2009")
ed 			<- as.date("10May2009")
freqdates 	<- sd:ed

# XXXX edit here to use srDailyInc

freq1 <- data.frame(lb=freqdates,ub=freqdates+1,f=srDailyInc(td1$cleanOnset,sd=sd,ed=ed))
freq2 <- data.frame(lb=freqdates,ub=freqdates+1,f=srDailyInc(td2$cleanOnset,sd=sd,ed=ed))
td2_dash <- td2[td2$cleanMexTrav != "Y",]
freq3 <- data.frame(lb=freqdates,ub=freqdates+1,f=srDailyInc(td2_dash$cleanOnset,sd=sd,ed=ed))
# freq1 <- sr1DHist(td1$cleanOnset,min=sd,max=ed,bins=ed-sd)
# freq2 <- sr1DHist(td2$cleanOnset,min=sd,max=ed,bins=ed-sd)

pvec <- as.vector(mode="numeric",data.frame(1,0.5))
names(pvec) <- c("A","r")
fb <- 13
lb <- 27

mle1 <- optim(pvec,likeInc,control=list(trace=0,fnscale=-1,maxit=1000),freqtab=freq1,firstbin=fb,lastbin=lb)
mle2 <- optim(pvec,likeInc,control=list(trace=0,fnscale=-1,maxit=1000),freqtab=freq2,firstbin=fb,lastbin=lb)
mle3 <- optim(pvec,likeInc,control=list(trace=0,fnscale=-1,maxit=1000),freqtab=freq3,firstbin=fb,lastbin=lb)

freq1 <- cbind(freq1,model=rep(NA,ed-sd+1))
freq2 <- cbind(freq2,model=rep(NA,ed-sd+1))

for (i in 1:(ed-sd)) {
	if (i >= fb && i <= lb) { 
		freq1$model[i] <- mle1$par["A"]*exp(mle1$par["r"]*(freq1$lb[i]-freq1$lb[fb-1]))
		freq2$model[i] <- mle2$par["A"]*exp(mle2$par["r"]*(freq2$lb[i]-freq2$lb[fb-1]))
	}
}

# Figure out the proportion of reports with onsets
tabOnsetDates			<- table(td2$cleanOnset)
tabOnsetDates_nmt		<- table(td2_dash$cleanOnset)
vecOnsetDates			<- as.numeric(names(tabOnsetDates))
noOnsetDates			<- length(vecOnsetDates)
tabOnsetDatesWk			<- table(td2$cleanOnsetWk)
vecOnsetDatesWk			<- as.numeric(names(tabOnsetDatesWk))
noOnsetDatesWk			<- length(vecOnsetDatesWk)

#Plot hospitalization status and ICU by week of onset
tabHospByOnsetWk 	<- table(td2$cleanOnsetWk,td2$cleanHosp)
tabICUByOnsetWk 	<- table(td2$cleanOnsetWk,td2$cleanICU)
tabHospByOnsetDay	<- table(td2$cleanOnset,td2$cleanHosp)
lbrates 			<- 0.01
pltHospWeek <- CIsForPlot(tabHospByOnsetWk[,"Y"]+tabHospByOnsetWk[,"N"],tabHospByOnsetWk[,"Y"],ratelb=lbrates)
pltICUWeek 	<- CIsForPlot(tabHospByOnsetWk[,"Y"],tabICUByOnsetWk[,"Y"],ratelb=lbrates)
pltHospDay 	<- CIsForPlot(tabOnsetDates[],tabHospByOnsetDay[,"Y"],ratelb=lbrates)

# Estimate overall and without onset values for hospitalization
hospital_all_n 	<- sum(tabHospByOnsetWk[,"Y"]) - tabHospByOnsetWk[1,"Y"]
hospital_all_N 	<- hospital_all_n + sum(tabHospByOnsetWk[,"N"]) - tabHospByOnsetWk[1,"N"]
hospital_all_pt	<- hospital_all_n / hospital_all_N
hospital_all_ub	<- binCI(c(hospital_all_N),c(hospital_all_n),0.025,min=0.000001)
if (hospital_all_n > 0) hospital_all_lb <- binCI(c(hospital_all_N),c(hospital_all_n),0.975,min=0.000001)
if (hospital_all_n == 0) hospital_all_lb <- lbrates

hosp_no_onset_n 	<- tabHospByOnsetWk[1,"Y"]
hosp_no_onset_N 	<- hosp_no_onset_n + tabHospByOnsetWk[1,"N"]
hosp_no_onset_pt	<- hosp_no_onset_n / hosp_no_onset_N
hosp_no_onset_ub	<- binCI(c(hosp_no_onset_N),c(hosp_no_onset_n),0.025,min=0.000001)
if (hosp_no_onset_n > 0) hosp_no_onset_lb <- binCI(c(hosp_no_onset_N),c(hosp_no_onset_n),0.975,min=0.000001)
if (hosp_no_onset_n == 0) hosp_no_onset_lb <- lbrates

# Estimate overall and without onset values for ICU
ICU_all_n 	<- sum(tabICUByOnsetWk[,"Y"]) - tabICUByOnsetWk[1,"Y"]
ICU_all_N 	<- hospital_all_n
ICU_all_pt	<- ICU_all_n / ICU_all_N
ICU_all_ub	<- binCI(c(ICU_all_N),c(ICU_all_n),0.025,min=0.000001)
if (ICU_all_n > 0) ICU_all_lb <- binCI(c(ICU_all_N),c(ICU_all_n),0.975,min=0.000001)
if (ICU_all_n == 0) ICU_all_lb <- lbrates

ICU_no_onset_n 	<- tabICUByOnsetWk[1,"Y"]
ICU_no_onset_N 	<- hosp_no_onset_n
ICU_no_onset_pt	<- ICU_no_onset_n / ICU_no_onset_N
ICU_no_onset_ub	<- binCI(c(ICU_no_onset_N),c(ICU_no_onset_n),0.025,min=0.000001)
if (ICU_no_onset_n > 0) ICU_no_onset_lb <- binCI(c(ICU_no_onset_N),c(ICU_no_onset_n),0.975,min=0.000001)
if (ICU_no_onset_n == 0) ICU_no_onset_lb <- lbrates

# Numbers needed for text
sum(table(td2$cleanOnset))-table(td2$cleanOnset)[1]
# The 8 May version of the US CDC line list data for the S-OIV outbreak contained...
sum(tabOnsetDates)
sum(tabOnsetDates) - tabOnsetDates[1]
# XX of these occurred between the XX and XX of ...
indexstart 	<- match(as.numeric(as.date("13Apr2009")),row.names(tabOnsetDates))
indexend	<- match(as.numeric(as.date("27Apr2009")),row.names(tabOnsetDates))
sum(tabOnsetDates[indexstart:indexend])
indexstart 	<- match(as.numeric(as.date("13Apr2009")),row.names(tabOnsetDates_nmt))
indexend	<- match(as.numeric(as.date("27Apr2009")),row.names(tabOnsetDates_nmt))
sum(tabOnsetDates_nmt[indexstart:indexend])
# We estimated the daily growth rate
mle2
mle1
# Hospitalization was reported for... 
sum(tabHospByOnsetWk) - sum(tabHospByOnsetWk[,"NK"])
sum(tabHospByOnsetWk) - sum(tabHospByOnsetWk[,"NK"]) - sum(tabHospByOnsetWk[1,"Y"]) - sum(tabHospByOnsetWk[1,"N"]) 
# Of those for whom ages were available
table(td2$cleanAG,td2$cleanHosp)
table(td2$cleanAG2)

table(td2$cleanContact,td2$cleanHosp)
binCI(127,4,0.025,min=0.00001)

# Do necessary regressions, following example from ben.
x <- data.frame(day			= as.numeric(row.names(tabHospByOnsetDay))[fb:lb],
		admitted	= tabHospByOnsetDay[fb:lb,"Y"],
		total		= tabHospByOnsetDay[fb:lb,"Y"] + tabHospByOnsetDay[fb:lb,"N"])

# plot(x$week, x$admitted/x$total)

# reshape data to "long" format i.e. one row per person
tmp <- c()
for(i in 1:nrow(x)) tmp <- c(tmp,rep(c(1,0),c(x$admitted[i],x$total[i]-x$admitted[i])))
x.reshaped <- data.frame(day=rep(x$day,x$total),admitted=tmp)

#check the reshaping has worked
x
# table(x.reshaped$day, x.reshaped$admitted)

x.glm <- glm(admitted~day, family=poisson, data=x.reshaped)
#summary(x.glm)
