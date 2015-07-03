rm(list=ls(all=TRUE))

require("odesolve")

seirModel <- function(t,y,p) {
	rtn=array(0,c(4))
	names(rtn) <- c("S","E","I","dS")
	names(y) <- c("S","E","I","dS")
	lambda <- p["beta"]*(1+p["a"]*cos((t+182)*pi/182))*y["I"]/p["N"]
	rtn["S"] <- p["omega"]*(p["N"]-y["S"]-y["E"]-y["I"])/364 - lambda*y["S"]
	rtn["E"] <- lambda*y["S"] - p["alpha"]*y["E"]
	rtn["I"] <- p["alpha"]*y["E"] - p["gamma"]*y["I"]
	rtn["dS"] <- lambda*y["S"]
	list(rtn)
}

# plot((1+0.05*cos((0:1000+182)*pi/182)))

seirParams <- c(
		N = 6800000,
		seed = 1,
		R0 = 1.8,
		omega = 1/2.4,
		alpha = 1/1.5,
		a = 0.0625,
		gamma = 1/0.5
		)

seirParams["S0"] = seirParams["N"] - seirParams["seed"]
seirParams["I0"] = seirParams["seed"]
seirParams["beta"] = seirParams["gamma"]*seirParams["R0"]

initial_conditions <- c(seirParams["S0"],0,seirParams["I0"],0)
names(initial_conditions) <- c("S","E","I","dS")

no_invals <- 200
no_years <- 40
time_points <- ((0:no_invals)*364*no_years)/no_invals
max_day <- max(time_points)
timestep <- time_points[2]-time_points[1]
last_index <- length(time_points+1)
inc_period <- 364/12.0

# Insert 
solution <- lsoda(initial_conditions,time_points,seirModel,seirParams,atol=1e-80)
solution <- cbind(solution, R = seirParams["N"] - solution[,"S"] - solution[,"E"] - solution[,"I"])
solution <- cbind(solution, incInf = rep(NA,length(time_points)))
solution[1,"incInf"] <- 0
for (i in 2:last_index) solution[i,"incInf"] <- (solution[i,"dS"] - solution[i-1,"dS"]) / timestep * inc_period 

plot(	0:1,
		xlim=c(0,max_day/364),
		ylim=c(0,0.75),
		type="n",
		xlab="Years",
		ylab="Number"
		)

pvars <- c("R","incInf")
pcols <- c("blue","black")

for (i in 1:(length(pvars))) {points(time_points[]/364,solution[,pvars[i]]/seirParams["N"],type="l",col=pcols[i],lwd=2)}

#legend(x=200,y=6000000,vecTimeSeries,lwd=2,col=vecColours)
