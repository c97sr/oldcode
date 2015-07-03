# Remove all objects 
rm(list=ls(all=TRUE))
require("odesolve")
setwd("C:\\eclipse\\R Source\\gz_syphilis")

tabData <- read.csv("data\\syphilis_gz_1994_to_2007.csv",header=TRUE)

# Command to make figure 2.1
windows(width=5.75,height=4)
plot(	tabData$year,
		tabData$cases.male+tabData$cases.female,
		ylab="Incidence",
		xlab="Year",
		type="l")

# Command to make figure 2.2
windows(width=5.75,height=4)
plot(	tabData$year,
		tabData$congenital,
		ylab="Incidence",
		xlab="Year",
		type="l")

sirsModel <- function(t,y,p) {
	rtn=array(0,c(4))
	names(rtn) <- c("S","I","R","dS")
	names(y) <- c("S","I","R","dS")
	lambda <- p["beta"]*y["I"]/(y["S"]+y["I"]+y["R"])
	rtn["S"] <- -lambda*y["S"] + p["omega"]*y["R"]+(1-p["p_i"])*p["gamma"]*y["I"]
	rtn["I"] <- lambda*y["S"] - p["gamma"]*y["I"]
	rtn["R"] <- p["p_i"]*p["gamma"]*y["I"] - p["omega"]*y["R"]
	rtn["dS"] <- lambda*y["S"]
	list(rtn)
}

sirsParams <- c(
		popsize = 1000,
		seed = 10,
		omega = 1/0.5,
		gamma = 1/2,
		p_i = 0.2,
		beta = 1
		)

initial_conditions <- c(sirsParams["popsize"]-sirsParams["seed"],sirsParams["seed"],0,0)
names(initial_conditions) <- c("S","I","R","dS")

no_years <- 14
time_points <- 0:no_years
timestep <- 1

solution <- lsoda(initial_conditions,time_points,sirsModel,sirsParams,atol=1e-80)

# Command to make figure 3.1
windows(width=5.75,height=4)
plot(	tabData$year,
		tabData$cases.male+tabData$cases.female,
		ylab="Incidence",
		xlab="Year",
		type="l",
		col="red",
		lwd=2)
points(	tabData$year,
		solution[2:(no_years+1),"dS"],
		type="l",
		col="green",
		lwd=2)


