rm(list=ls(all=TRUE))

#setwd("D:\\eclipse\\R Source")
setwd("D:\\Eric\\JE\\Source")

#source("je\\je_revised_funcs.r")
#source("je\\je_params.r")

source("je_revised_funcs.r")
source("je_params.r")

#pigData = read.table("je\\pig.in",header=TRUE)
pigData = read.table("pig.in",header=TRUE)

pigData$Day <- pigData$Day + 360
pigData$Prev <- pigData$n / pigData$N

p_mos_fitted_names <- c("v_s","N_site_TK","N_site_LK","N_site_PM","tw_g","tw_h","tw_f","tw_l")
p_mos_fitted_values <- c(68.517,191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)

# try other most fitted value: larger peak to trough ratio (1st par)
# not needed
# p_mos_fitted_values <- c(68.517*5,191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)


jeParams <- updatePVector(p_mos_fitted_values,p_mos_fitted_names,jeParams)

theta_names = c("N_m","omega","R0Max")
theta_in = c(10000,0.1,10)
jeDynLBounds = c(1,0.00001,5)
jeDynUBounds = c(2000000,1,10000)

# jeEst <- optim(theta_in,jeLike,method = "L-BFGS-B",lowe =jeDynLBounds,upper=jeDynUBounds,
#            names_theta=theta_names,allps = jeParams,x = pigData,
#            control=list(trace=0,fnscale=-1,maxit=10000,parscale=theta_in))


# theta_in <- c(30,0.1,40)


jeEst <- optim(theta_in,jeLike,method = "L-BFGS-B",lowe =jeDynLBounds,upper=jeDynUBounds,
            names_theta=theta_names,allps = jeParams,x = pigData,
            control=list(trace=0,fnscale=-1,maxit=10000,parscale=theta_in))

jeLike(theta_out,theta_names,pigData,jeParams,saturated=FALSE)


# updated estimate
theta_out <- c(39.0502051611532,0.097324661805472,93.1205756327272)

theta_out <- c(39.0502051611532*5,0.097324661805472,93.1205756327272)

# result
# c(39.0502051611532*5,0.097324661805472,93.1205756327272) jeLike -177.7535
# c(17.8164025659383,0.701461390992609,190.675306875091) jeLike -176.2515

##------
## add one more parameter to estimate: mos peak to trough ratio
p_mos_fitted_names <- c("N_site_TK","N_site_LK","N_site_PM","tw_g","tw_h","tw_f","tw_l")
p_mos_fitted_values <- c(191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)

jeParams <- updatePVector(p_mos_fitted_values,p_mos_fitted_names,jeParams)

theta_names = c("v_s","N_m","omega","R0Max")
theta_in = c(200,10000,0.1,10)
jeDynLBounds = c(100,1,0.00001,5)
jeDynUBounds = c(100000,2000000,1,10000)

# lower bound for v_s, should be higher, 10 or more?

# is there possible identifiability problem: N_m prop (R0Max)^2 ?


jeEst <- optim(theta_in,jeLike,method = "L-BFGS-B",lowe =jeDynLBounds,upper=jeDynUBounds,
            names_theta=theta_names,allps = jeParams,x = pigData,
            control=list(trace=0,fnscale=-1,maxit=10000,parscale=theta_in))

jeLike(theta_out,theta_names,pigData,jeParams,saturated=FALSE)


#  result
# -177.977545024687 c(66.6830265796294,42.1797003654916,0.0995112409505348,86.3747797716872)


#### plot result
p_mos_fitted_names <- c("v_s","N_site_TK","N_site_LK","N_site_PM","tw_g","tw_h","tw_f","tw_l")
p_mos_fitted_values <- c(68.517,191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)
jeParams <- updatePVector(p_mos_fitted_values,p_mos_fitted_names,jeParams)

theta_out <- c(66.6830265796294*10,42.1797003654916*10,0.0995112409505348,86.3747797716872)
params_out <- updatePVector(theta_out,theta_names,jeParams)

jeLike(theta_out,theta_names,pigData,jeParams,saturated=FALSE)

theta_names = c("start","N_m","omega","R0Max")
theta_out <- c(0,39.0502051611532,0.097324661805472,93.1205756327272)
params_out <- updatePVector(theta_out,theta_names,jeParams)

no_points <- 100
max_days <- 2000
ics <- jeMakeICS(params_out)
solution <- lsoda(ics,0:no_points/no_points*max_days,jeModel,params_out,hmax=1)

mospop <- rep(0,2000)
for (i in 0:2000)
{
mospop[i] <- jeMosquitoPop(i,params_out)
}

ps <- getPopSizes(solution,jeParams)
pp <- getPigPrevFinalAgeClass(solution,jeParams)
pi <- getPigImmFinalAgeClass(solution,jeParams)
ii <- getMosPrevFinalAgeClass(solution,jeParams)
mp <- getMosPop(solution,jeParams)
windows()
plot(pp,type="l",ylim=c(0,1),lwd=2)
lines(pi[,1],1-pi[,2], col="blue")
points(pigData$Day,pigData$Prev)
lines(ps$t, ps$M/max(ps$M), col="red", lty=1)
lines(ii[,1],ii[,2]/max(ii[,2]), col="brown", lty=2)

windows()
plot(pp,type="l",ylim=c(0,1),lwd=2)
lines(pi[,1],1-pi[,2], col="blue")
points(pigData$Day,pigData$Prev)
lines(ps$t+theta_out[1], ps$M/max(ps$M), col="red", lty=1)
lines(ii[,1],ii[,2]/max(ii[,2]), col="brown", lty=2)





##------

## allow the timing of mos cycle as par
p_mos_fitted_names <- c("v_s","N_site_TK","N_site_LK","N_site_PM","tw_g","tw_h","tw_f","tw_l")
p_mos_fitted_values <- c(68.517,191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)

jeParams <- updatePVector(p_mos_fitted_values,p_mos_fitted_names,jeParams)

theta_names = c("start","N_m","omega","R0Max")
theta_in = c(1,20000,0.5,10)
jeDynLBounds = c(0,100,0.00001,5)
jeDynUBounds = c(360,100000,1,1000)

jeEst <- optim(theta_in,jeLike,method = "L-BFGS-B",lowe =jeDynLBounds,upper=jeDynUBounds,
            names_theta=theta_names,allps = jeParams,x = pigData,
            control=list(trace=0,fnscale=-1,maxit=10000,parscale=theta_in))

# some estimate:
# (allowing large R0 (up to 1000) results in poor fit: likelihood ~ -100
# -101.341500194030 c(58.0262688681185,100,0.099950278530821,1121.40906368216)
# -69.859294159986 c(98.1119871858753,99587.6572584663,0.100179277211274,100) **
# -69.511273211595 c(91.3711010594994,41518.4198922334,0.875464362205162,100) 
# -69.8639087419294 c(98.1460927489414,99999.4717946982,0.499998823970087,99.99)
# -65.115786939943 c(107.318368052297,1e+05,0.950992761474657,140.720348476100)

theta_out <- c(98.1119871858753,99587.6572584663,0.100179277211274,100)
theta_out <- c(91.3711010594994,41518.4198922334,0.875464362205162,100)
theta_out <- c(107.318368052297,1e+05,0.950992761474657,140.720348476100)


# some estimate (with gamma_m = 1)
# -65.6524638001752 c(90.1797081048126,28595.0314997641,0.346321162428643,82.1830526863394)


# some estimate (with gamma_m = 1, sigma=1/42)
# -68.1864363638851 c(84.8512965296314,1e+05,0.518271479748789,84.2745280539394)

##------


theta_out <- c(62.3257916613551,0.115487082574281,11.4927939612781)
params_out <- updatePVector(theta_out,theta_names,jeParams)

# increase mos population
theta_out <- c(62.3257916613551*20000,0.115487082574281*8,11.4927939612781*2)
params_out <- updatePVector(theta_out,theta_names,jeParams)

no_points <- 100
max_days <- 2000
ics <- jeMakeICS(params_out)
solution <- lsoda(ics,0:no_points/no_points*max_days,jeModel,params_out,hmax=1)

ps <- getPopSizes(solution,jeParams)
pp <- getPigPrevFinalAgeClass(solution,jeParams)
pi <- getPigImmFinalAgeClass(solution,jeParams)
windows()
plot(pp,type="l",ylim=c(0,1))
lines(pi[,1],1-pi[,2], col="blue")
points(pigData$Day,pigData$Prev)


# Plot Prevelance at other age
a <- 2
ppa <- getPigPrevAnyAgeClass(solution,a,jeParams)
pia <- getPigImmAnyAgeClass(solution,a,jeParams)
windows()
plot(ppa,type="l",ylim=c(0,1))
lines(pia[,1],1-pia[,2], col="blue")
points(pigData$Day,pigData$Prev)


# Plot Prevalence at each age
windows()
par(mfrow=c(3,2))
for (a in 1:6)
{
ppa <- getPigPrevAnyAgeClass(solution,a,jeParams)
pia <- getPigImmAnyAgeClass(solution,a,jeParams)
plot(ppa,type="l",ylim=c(0,1), main=paste("Age =",a), las=1)
lines(pia[,1],1-pia[,2], col="blue")
points(pigData$Day,pigData$Prev)
}


# Plot Population
windows(height=12, width=10)
par(mfrow=c(2,1))
plot(ps$t,ps$M,type="l",main="mosquito")
plot(ps$t,ps$P,type="l",main="pigs",ylim=c(0,400000))
#plot(ps$t,ps$P,type="l",main="pigs")

windows()
plot(getInc(solution,"M",1:1,jeParams),type="l")
windows()
plot(getInc(solution,"P",1:jeParams["n_a"],jeParams),type="l")

jeLike(theta_out,theta_names,pigData,jeParams,saturated=FALSE)

# Likelihood with other parameters
# "original" Likelihood -330.4205 (with warning message)
# increase mos pop in "solution" to 10 times: -146.7745 (with warning message)
# increase mos pop in "solution" to 100 times: -111.6246 (with warning message)
# increase mos pop in "solution" to 1000 times: -123.2385 (no error message)
# increase mos pop in "solution" to 10000 times: -144.5840 (no error message)
# 100x mos pop, 2x omega: 
