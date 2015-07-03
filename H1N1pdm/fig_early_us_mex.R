rm(list=ls(all=TRUE))
options(recover=NULL)
require(date)

source("../R Source/SR_misc.r")
source("./funcs.r")

sdate <- as.numeric(as.date("11Mar2009"))
edate <- as.numeric(as.date("31Jul2009"))

timepoints 	<- sdate + 0:((edate-sdate)*10)/10
notimepoints <- length(timepoints)-1
noReals		<- 10

varnames 	<- c(	"S_MM",	"I_MM", "R_MM", "S_UU","I_UU","R_UU",
					"S_MU",	"I_MU", "R_MU",	"S_UM","I_UM","R_UM")
			
nocols 		<- length(varnames)
vars		<- vector(mode="numeric",length=nocols)
names(vars) <- varnames

p <- fnReadFileAsParamVector("./params_base.txt")

output <- array(dim=c(noReals,notimepoints+1,3))
# Last index 	= 1 for infections in mexico
#				= 2 for infections in returning travellers
#				= 3 for infections within the united states

output[,1,] <- 0

# Outer realization lool
for (k in 1:noReals) {
	
	# Set up initial conditions for the variables
	vars[] 		 <- 0
	vars["S_MU"] <- round(p["N_T"]*(1-p["P_U"])*p["D_S"])
	vars["S_MM"] <- p["N_M"] - p["seed"] - vars["S_MU"]
	vars["I_MM"] <- p["seed"]
	vars["S_UM"] <- round(p["N_T"]*p["P_U"]*p["D_S"])
	vars["S_UU"] <- p["N_U"] - vars["S_MU"]
	
	# Set up initial conditions for run variables
	for (i in 1:notimepoints) {
		
		# Set up dt
		t <- timepoints[i]
		dt <- timepoints[i+1]-timepoints[i]
		
		# Calculate infection probabilities
		pInf_U 	<- 1-exp(-dt*lambda_U(t,p,vars))
		pInf_M 	<- 1-exp(-dt*lambda_M(t,p,vars))
		pRec	<- 1-exp(-dt/p["D_I"])
				
		# Generate infection events
		evInf_UU 	<- rbinom(1,vars["S_UU"],pInf_U)
		evInf_UM 	<- rbinom(1,vars["S_UM"],pInf_M)
		evInf_MU 	<- rbinom(1,vars["S_MU"],pInf_U)
		evInf_MM 	<- rbinom(1,vars["S_MM"],pInf_M)
		evRec_UU 	<- rbinom(1,vars["I_UU"],pRec)
		evRec_UM 	<- rbinom(1,vars["I_UM"],pRec)
		evRec_MU 	<- rbinom(1,vars["I_MU"],pRec)
		evRec_MM 	<- rbinom(1,vars["I_MM"],pRec)

		# Update state variables with infection events
		vars["S_UU"]	<- vars["S_UU"] - evInf_UU
		vars["S_UM"]	<- vars["S_UM"] - evInf_UM
		vars["S_MU"]	<- vars["S_MU"] - evInf_MU
		vars["S_MM"]	<- vars["S_MM"] - evInf_MM
		vars["I_UU"]	<- vars["I_UU"] + evInf_UU - evRec_UU
		vars["I_UM"]	<- vars["I_UM"] + evInf_UM - evRec_UM
		vars["I_MU"]	<- vars["I_MU"] + evInf_MU - evRec_MU
		vars["I_MM"]	<- vars["I_MM"] + evInf_MM - evRec_MM
		vars["R_UU"]	<- vars["R_UU"] + evRec_UU
		vars["R_UM"]	<- vars["R_UM"] + evRec_UM
		vars["R_MU"]	<- vars["R_MU"] + evRec_MU
		vars["R_MM"]	<- vars["R_MM"] + evRec_MM
		
		# Calculate probability of travel		
		pTrv	<- 1-exp(-dt*1/p["D_S"])				
		
		# Draw random travel events
		evTrav_SU <- rbinom(1,vars["S_UM"],pTrv)
		evTrav_IU <- rbinom(1,vars["I_UM"],pTrv)
		evTrav_RU <- rbinom(1,vars["R_UM"],pTrv)
		totTrav_U <- evTrav_SU + evTrav_IU + evTrav_RU
		vecRtnTrav_U <- rmultinom(1,totTrav_U,c(vars["S_UU"],vars["I_UU"],vars["R_UU"])) 
				
		evTrav_SM <- rbinom(1,vars["S_MU"],pTrv)
		evTrav_IM <- rbinom(1,vars["I_MU"],pTrv)
		evTrav_RM <- rbinom(1,vars["R_MU"],pTrv)
		totTrav_M <- evTrav_SM + evTrav_IM + evTrav_RM
		vecRtnTrav_M <- rmultinom(1,totTrav_M,c(vars["S_MM"],vars["I_MM"],vars["R_MM"])) 
		
		# Update state variables for travel events
		vars["S_UU"]	<- vars["S_UU"] - vecRtnTrav_U[1] + evTrav_SU
		vars["S_UM"]	<- vars["S_UM"] + vecRtnTrav_U[1] - evTrav_SU
		vars["S_MU"]	<- vars["S_MU"] + vecRtnTrav_M[1] - evTrav_SM
		vars["S_MM"]	<- vars["S_MM"] - vecRtnTrav_M[1] + evTrav_SM
		vars["I_UU"]	<- vars["I_UU"] - vecRtnTrav_U[2] + evTrav_IU
		vars["I_UM"]	<- vars["I_UM"] + vecRtnTrav_U[2] - evTrav_IU
		vars["I_MU"]	<- vars["I_MU"] + vecRtnTrav_M[2] - evTrav_IM
		vars["I_MM"]	<- vars["I_MM"] - vecRtnTrav_M[2] + evTrav_IM
		vars["R_UU"]	<- vars["R_UU"] - vecRtnTrav_U[3] + evTrav_RU
		vars["R_UM"]	<- vars["R_UM"] + vecRtnTrav_U[3] - evTrav_RU
		vars["R_MU"]	<- vars["R_MU"] + vecRtnTrav_M[3] - evTrav_RM
		vars["R_MM"]	<- vars["R_MM"] - vecRtnTrav_M[3] + evTrav_RM
		
		output[k,i+1,1] <- evInf_MM
		output[k,i+1,2] <- evTrav_IU
		output[k,i+1,3] <- evInf_UU
		# output[k,i+1,4:(3+nocols)] <- vars[]
		
	}
	
}

pdf(file="./figs/Simple_2_Pop.pdf",width=10,height=6)
plot(as.date(timepoints),output[1,,1]+1,type="n",log="y",xlab="Day",ylab="Incidence + 1 (log scale)")
for (k in 1:noReals) {
	points(timepoints,output[k,,1]+1,type="l",col="red")
	points(timepoints,output[k,,2]+1,type="l",col="blue")
	points(timepoints,output[k,,3]+1,type="l",col="yellow")
}

dev.off()
