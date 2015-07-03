require("date")
source("../H1N1pdm/funcs.R")

# The commands below seem to work pretty well, but the charts aren't the ones now needed.

# x <- read.csv("../../files/projects/influenza/swineflu/peak_demand/results/091015_mac/test_mcmc.csv")
x <- read.csv("~/tmpresults/091021/test_mcmc_tmp.csv")

# names(x)
# "sample"   "lnlike"   "R0"       "a_total"  "a_yc"     "a_oc"    
# "a_ya"     "a_oa"     "phi_yc"   "phi_oc"   "phi_ya"   "phi_oa"  
# "r"        "t0"       "p_H_base" "gamma_h" 

sizex <- dim(x)[1]
rburn <- 1:round(sizex/3)
rmid <- round(sizex/3+1):round(2*sizex/3)
rend <- round(2*sizex/3+1):sizex

param <- "r"
s1 <- x[rmid,param]
s2 <- x[rend,param]
plot(s1,type="l",col="red",ylim=c(min(c(s1,s2)),max(c(s1,s2))),main=param)
points(s2,type="l",col="blue")

outputTable	<- genMeanMedianCIs(x[rend,],c(
				"R0","p_H_base","a_total",
				"a_yc","a_oc","a_ya","a_oa",
				"gamma_h","t0"				
			))

OzHosp <- read.csv("../H1N1pdm/working/Pub_Aus_Data.csv")
OzHosp$Date <- as.date(as.character(OzHosp$Date),order="mdy")
source("../H1N1pdm/working/default_QLD_params.R")

pdf(file="../H1N1pdm/figs/modelfit.pdf",width=20/cm(1),height=10/cm(1),version="1.4")
plot(OzHosp$Date,OzHosp$Qld_Hosp,xlim=c(as.date("05May2009"),as.date("31Oct2009")),ylim=c(0,120),col="black")
# Test the basic setup of the model
DbParamNames 	<- c("r",		"t0",		"p_H_base",	"gamma_h")
sampstaken <- 2500+0:20*10  
# DbParamValues 	<- c(0.0435637341188034,18016.0708335772,0.0163009724031923,0.333355533974246)
for (i in sampstaken) {
	DbParamValues <- c()
	for (name in DbParamNames) {
			DbParamValues <- c(DbParamValues,x[i,name])
	}
	cat(i,DbParamValues,"\n")
	
	tmp 			<- SirModelFourAgeClasses(pname=DbParamNames,pvals=DbParamValues,
							casemix=c(7,81,118,15),vp=FourAgeParams(),NGM=MakeNgmFour)		 
	stnFour 		<- tmp$sol
	
	# plot(as.date(stnFour[,"time"]),(stnFour[,"Ih"]+stnFour[,"Iv"])*4380000,col="black",type="l")
	points(as.date(stnFour[,"time"]),(stnFour[,"Ih"]+stnFour[,"Iv"])*4380000,col="black",type="l")
	# points(OzHosp$Date,OzHosp$Qld_Hosp)
}

points(OzHosp$Date,OzHosp$Qld_Hosp,col="red")
dev.off()
