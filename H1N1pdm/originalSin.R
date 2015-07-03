rm(list=ls(all=TRUE))
options(error=NULL)
require("date")
require("seqinr")
source("SR_misc.r")
source("swineflu/swine_flu_funcs.r")

# Load up the final database and the version immediately before that one
# Use "table(tabDataX$CaseStatus) to check the categories"
tabData2 	<- read.csv("../../projects/influenza/swineflu/data/cdc/linelists13may_new.csv")
td2 		<- postProcLineList(tabData2,order="dmy")

# Declare the required variables
sy			<- 1910
ey			<- 2009
possYears 	<- sy:ey
noYears 	<- length(possYears)
mainTab	 	<- data.frame(years = possYears, N = rep(0,noYears), n = rep(0,noYears), rounc = rep(0,noYears))
dummyTab	<- data.frame(years = possYears, N = rep(0,noYears), n = rep(0,noYears), rounc = rep(0,noYears))
noData 		<- dim(td2)[1]

# Set up actual data
for (i in 1:noData) {
	year <- round(ey-td2$cleanAge[i])
	if (year <= ey) {
		mainTab$N[year-sy+1] <- mainTab$N[year-sy+1] + 1
		if (td2$cleanHosp[i] == "Y") mainTab$n[year-sy+1] <- mainTab$n[year-sy+1] + 1
	}
}

# Set up dummy data
for (y in possYears) {
	dummyTab$N[y-sy+1] <- 100
	if (y >= 1973) dummyTab$n[y-sy+1] <- 10
	else dummyTab$n[y-sy+1] <- 5
}

# Likelihood function
test_function <- function(y,vecp) {
	if (y < vecp[1]) rtn <- 1
	else rtn <- vecp[2]
	rtn
}

oasLike <- function(vecp,allyears,func,data) {
	total_risk <- 0
	for (y in allyears) total_risk <- total_risk + func(y,vecp) * data$N[y-sy+1]
	total_risk
	noYears <- dim(data)[1]
	probVec <- rep(0,noYears)
	for (y in allyears) probVec[y-sy+1] <- 	func(y,vecp) * data$N[y-sy+1] / total_risk 
	# browser()
	like <- dmultinom(data$n,prob=probVec,log=TRUE)
	like
}

likeWrapper <- function(vecf,fittedp,vecp,allyears,func,data) {
	for (i in length(fittedp)) vecp[fittedp[i]] <- vecf[i]
	rtn <- oasLike(vecp,allyears,func,data)
	rtn
}

tmpLikeWrapper <- function(x,fittedindex,vecp,allyears,func,data) {
	vecp[fittedindex] <- x
	rtn <- oasLike(vecp,allyears,func,data)
	rtn
}

# Test the liklihood
oasLike(c(1987,2.3),sy:ey,test_function,dummyTab)
likeWrapper(c(2.0),c(2),c(1973,2),sy:ey,test_function,dummyTab)

# Optimize to find best parameter values
optim(	c(2.2),likeWrapper,fittedp=c(2),vecp=c(1988,2.0),allyears=sy:ey,func=test_function,data=dummyTab,
		list(trace=0,fnscale=-1,maxit=1000))
optimize(	f=tmpLikeWrapper,interval=c(0.1,4),maximum=TRUE,
			fittedindex=2,vecp=c(1973,2.0),allyears=sy:ey,func=test_function,data=dummyTab)
	
tmp <- mainTab$N/max(mainTab$N)
col.seq <- rgb( 1-tmp, 1-tmp, 1-tmp  )
plot(mainTab$year, mainTab$n/sum(mainTab$n), pch=16)

for (i in 1:noYears) mainTab$rounc[i] <- mainTab$years[i] %/% 3

mainTab