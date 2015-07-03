# Make a very simple version of the model that runs and gives output
# Next things to do
# - calculate quick likelihood
# - have a burnin period
# - implement an annual spline for seeding

rm(list=ls(all=TRUE))
options(error=NULL)
require("odesolve")

# Force of infection 
# Takes the time, a vector of current state variables and a vector of parameters
# Returns the hazard per susceptible individuals

foi <- function(t,y,p) {
	pop <- y["S"] + y["E"] + y["I"] + y["R"]
	beta_d <- p["R0"]/p["D_I"]
	if (p["direct_seas"]) rtn <- beta_d * y["I"] / pop * season(t,p["off_s"],p["d_s"],p["A"]) + p["beta_w"] * y["W"] 
	else rtn <- beta_d * y["I"] / pop + p["beta_w"] * y["W"] 
	as.numeric(rtn) 
}

# Seasonality
# Takes the current time, the offset time, the duration of the season and the amplitude of seasonality
# Returns either 1 or (1 - amplitude)
season <- function(t,off,d,a){
	rtn <- 1
	if ((t-off) %% 364 > d) rtn <- rtn * (1 - a)
	rtn
}

# Differential equation model
# Takes the current time, a vector of state variables and the parameters
# Returns the vector of time derivatives for the model
birdMod <- function(t,y,p) {
	rtn <- y
	lambda <- foi(t,y,p)
	N <- sum(y[c("S","E","I","R")])
	p_it <- splinefun(c(0,p["tmax"]/2,p["tmax"]),c(1,p["rp_i2"],p["rp_i3"]))
	rtn["S"]  <-  (1-p["p_i"]*p_it(t))*1/p["D_M"]*N - 1 / p["D_M"]*y["S"] - lambda*y["S"]
	rtn["E"]  <-  + lambda * y["S"] - 1 / p["D_M"]*y["E"] - 1 / p["D_E"] * y["E"]
	rtn["I"]  <-  {+ p["p_i"]*p_it(t)*1 / p["D_M"]*N - 1 / p["D_M"]*y["I"] +
					1 / p["D_E"] * y["E"] - 1 / p["D_I"] * y["I"]}
	rtn["R"]  <-  - 1 / p["D_M"]*y["R"] + 1 / p["D_I"] * y["I"] 
	rtn["W"]  <-  + p["gamma"]*y["I"] - p["mu_w"] * season(t,p["off_s"],p["d_s"],p["A"])*y["W"]
	rtn["dS"] <-  + lambda*y["S"]
	
	list(rtn)
}

# Model parameters
birdPs <- c(
	R0 = 5.5,
	D_M = 3,
	D_E = 1.0,
	D_I = 7.0,
	d_s = 364/4,
	A = 0.3,
	gamma = 0,
	mu_w = 0,
	beta_w = 0,
	p_i = 0.00001,
	rp_i2 = 1, 
	rp_i3 = 1, 
	p_I0 = 0.01,
	N0 = 10000,
	off_s = -120,
	direct_seas = TRUE,
	tmax = 364*5
	)

# Load in the surveillance data and make a single strain dataset
hosts <- c("Chicken","Waterfowl","Quail","Other","X_other")
strain <- "H5"
dat <- read.csv("./tab1.csv")[1:73,]
dat_times <- (0:72)*28
allN <- rowSums(dat[,paste(hosts,"N",sep=".")])
alln <- rowSums(dat[,paste(hosts,strain,sep=".")])

# Get a simple model result
mod_times <- 28*(0:72)
y0 <- c(S=as.numeric(birdPs["N0"]*(1-birdPs["p_I0"])),E=0,I=as.numeric(birdPs["p_I0"]*birdPs["N0"]),R=0,W=0,dS=0)
solution <- lsoda(y0,mod_times,birdMod,birdPs)

pdf("Model_and_data.pdf",width=10/cm(1),height=10/cm(1))
par(mgp=c(1.8,0.6,0))
plot(1:2,type="n",xlim=c(0,birdPs["tmax"]),ylim=c(0,0.1),axes=FALSE,xlab="Years",ylab="Prevalence Infectious")
axis(1,at=(0:20/4)*364,labels=c("0","","","","1","","","","2","","","","3","","","","4","","","","5"))
axis(2)
points(dat_times,alln/allN,type="l",col="red")
points(mod_times,solution[,"I"]/rowSums(solution[,c("S","E","I","R")]),type="l",col="blue")
points(spline(c(0,max(mod_times)/2,max(mod_times)),c(0.1,0.2,0.1)),type="l",col="green")

legend(birdPs["tmax"]*4/8,y=0.07,c("Data","Model"),
		col=c("red","blue"),
		lty=1,lwd=2,y.intersp=1.2,horiz=FALSE,yjust=0,bty="n")
dev.off()

# Play with a spline for changes in annual beta
x <- splinefun(c(0,36,72),c(1.03,1.01,1.01))

