rm(list=ls(all=TRUE))

# setwd("D:\\eclipse\\R Source")
setwd("D:\\Eric\\JE\\Source")

source("je_revised_funcs.r")
source("je_params.r")

jeParams <- jeUpdateAuxParams(jeParams)
culexData = read.table("CulexTri_sr_for_r.txt",header=TRUE)
pstofit = c("v_s","N_site_TK","N_site_LK","N_site_PM","tw_g","tw_h","tw_f","tw_l")
pstofitvals = c(10,1,1,1,987.77,1277.7,1290.54,172.57)
jeLbounds = c(1,1,1,1,0.00001,0.00001,0.00001,0.00001)
jeUbounds = c(1000,100000,100000,100000,100000,100000,100000,100000)

# jeEst <- optim(pstofitvals,jeMosquitoLikelihood,method = "L-BFGS-B",lowe =jeLbounds,upper=jeUbounds,
#            control=list(trace=0,fnscale=-1,maxit=100000),fitpn=pstofit,allps=jeParams,df=culexData)

#Other inital values
#pstofitvals2 = c(100,1,1,1,987.77,1277.7,1290.54,172.57)
#pstofitvals3 = c(10,1,1,1,8,8,8,1)
#pstofitvals4 = c(100,1,1,1,8,8,8,1)
#pstofitvals = c(50,200,200,6,5,5,5,1)


# new estimate a2
# Initial c(10,1,1,1,987.77,1277.7,1290.54,172.57)
# currentFita2 <- c(64.5876990,186.8703418,199.3023820,5.9730636,1467.9847867,1271.9154733,880.6454923,0.4177638)
# Likelihood: -3250.899

# new estimate a3
# Initial c(100,1,1,1,987.77,1277.7,1290.54,172.57)
# currentFita3 <- c(68.239878,190.161165,215.141342,5.819128,1491.006601,1354.584311,930.498574,0.000010)
# Likelihood: -3316.432

# new estimate a4
# Initial c(100,100,100,100,987.77,1277.7,1290.54,172.57)
# currentFita4 <- c(90.783910,207.178505,232.305937,9.967395,1490.142221,1217.173797,966.443209,0.000010)
# Likelihood: -3294.011

# new estimate a5
# Initial c(50,200,200,6,987.77,1277.7,1290.54,172.57)
# currentFita5 <- c(82.508844,216.178958,235.364620,7.010295,1468.436962,1221.131716,932.706499,0.000010)
# Likelihood: -3272.032

# new estimate a6
# Initial c(50,200,200,6,5,5,5,1)
# currentFita6 <- c(49.93068810,200.03401229,200.07302588,6.11006589,6.09440381,6.35589309,3.99135089,0.05198705)
# Likelihood: -3203.634

# b: even longer summer
# new estimate b1
# Initial c(10,1,1,1,8,8,8,1)
# currentFitb1 <- c(15.565993,181.113298,159.088082,4.461355,68.706787,134.994962,41.691973,12.464268)
# Likelihood :  -4270.181



currentFit <- c(68.517,191.19,199.01,10.790,1468.5,1260.8,893.74,1e-05)
p_out <- updatePVector(currentFit,pstofit,jeParams)
jeMosquitoLikelihood(currentFit,pstofit,jeParams,culexData)

culexData <- jeNormWithMax(culexData,updatePVector(currentFit,pstofit,jeParams))
plot(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",log="y")
points(culexData$Days360_01_11_04,culexData$NormMosPerTrap,type="p")

# Check that dynamic model has same properties
max_days <- 720
no_points <- 100
ics <- jeMakeICS(p_out)
solution <- lsoda(ics,0:no_points/no_points*max_days,jeModel,p_out,hmax=10)
solution
plot(solution,type="l")
points(0:720 ,sapply(0:720,jeMosquitoPop,p=p_out))


# ------ Other outputs

p_out <- updatePVector(currentFit,pstofit,jeParams)
plot(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l")
points(culexData$Days360_01_11_04,culexData$NormMosPerTrap/max(culexData$NormMosPerTrap)*1000,type="p")

p_out <- updatePVector(currentFita2,pstofit,jeParams)
lines(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",col="blue")

p_out <- updatePVector(currentFita3,pstofit,jeParams)
lines(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",col="green")

p_out <- updatePVector(currentFita4,pstofit,jeParams)
lines(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",col="brown")

p_out <- updatePVector(currentFita5,pstofit,jeParams)
lines(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",col="yellow")

p_out <- updatePVector(currentFita6,pstofit,jeParams)
lines(300:720,sapply(300:720,jeMosquitoPop,p=p_out),type="l",col="red")
