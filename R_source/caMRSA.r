# Quick model to fit caMRSA data for hong kong
rm(list=ls(all=TRUE))
require(odesolve)

# Lets check this time sync stuff by opening the same eclipse session form both machines

updVec <- function(values,names,pvec) {                                                                        
    no_updates <- length(values)
    for (i in 1:no_updates) pvec[names[i]]=values[i]
    pvec["beta"] <- pvec["R0"]/(1/pvec["kappa"]+ pvec["gamma"]/(pvec["kappa"]+pvec["gamma"])*pvec["epsilonZ"]/pvec["tau"])
    pvec
}

#Load data
hospData = read.table("D:\\files\\projects\\caMRSA\\data\\061011\\caMRSA_97-05.txt",header=TRUE)

ps <- c(
    R0 = 1,         
    epsilonZ = 1,   #Not important at this time
    tau = 1/10,     #Assumed ot be 10 days
    kappa = 1/370,  #Scanvic, CID, 32, 1393
    gamma = 0.005952,
    n = 2000000,    #Approximate population served by QMH
    s = 1,          #Initial seed size
    timestart = 365*5
)

caMRSA <- function(t,s,p) {

    # set up the variables
    rtn <- array(0,3)
    names(rtn) <- c("x","y","cumInc")
    names(s) <- c("x","y","cumInc")

    rtn["x"] <- -p["beta"]*s["x"]*(s["y"] + p["epsilonZ"]*(1-s["x"]-s["y"])) + p["kappa"]*s["y"] + p["tau"]*(1-s["x"]-s["y"])
    rtn["y"] <- +p["beta"]*s["x"]*(s["y"] + p["epsilonZ"]*(1-s["x"]-s["y"])) - p["kappa"]*s["y"] - p["gamma"]*s["y"]
    rtn["cumInc"] <- p["gamma"]*s["y"]*p["n"]

    list(rtn)

}

incFromSol <- function(sol) {
    n <- dim(sol)[1] 
    rtn <- array(0,c(n-2,2))
    for (i in 3:n) {
        deltaT = sol[i,"time"]-sol[i-1,"time"]
        rtn[i-2,1]=sol[i-1,"time"]+deltaT/2
        rtn[i-2,2]=sol[i,"cumInc"]-sol[i-1,"cumInc"]
    }
    rtn
}

debug(likeMRSA)
likeMRSA <- function(fitted_ps,fitted_names,all_ps,hdata) {
    if (min(fitted_ps) < 0) lnlike <- -1e10 else {
        noData = dim(hdata)[1]
        lnlike <- 0 
        p=updVec(fitted_ps,fitted_names,all_ps)
        ics <- c(1-p["s"]/p["n"],p["s"]/p["n"],0)
        names(ics) <- c("x","y","cumInc")
        sol <- lsoda(ics,c(-p["timestart"],0:(noData)*365),caMRSA,p)
        modInc <- incFromSol(sol)[,2]
        if (min(modInc) < 0.0000001) lnlike <- -1.1e10 else {
            for (i in 1:noData) lnlike = lnlike + dpois(hdata$Incidence[i],modInc[i],log=TRUE)
            }
            plot(hdata$Year,hdata$Incidence,ylim=c(0,max(c(hdata$Incidence,modInc))))
            points(hdata$Year,modInc,type="l")
        }
    cat(lnlike,fitted_ps,"\n")
    flush.console()
    lnlike
}



# XXXX Next, calculate the R0 for this model and recast the beta in terms of R0
# Then look at the bounds for R0 - must be greater than one
psToFit <- c("R0","timestart")
startVals <- c(2,365*5)
estMRSA <- optim(startVals,likeMRSA,control=list(trace=0,fnscale=-1,maxit=10000),
        fitted_names=psToFit,all_ps=ps,hdata=hospData)
       

#testPs = c(1.19,1.409761e-4,1.008782e+2)
testPs = c(1.19,1.409761e-4,1.008782e+2)

likeMRSA(estMRSA$par,psToFit,ps,hospData)

fitPs <- updVec(estMRSA$par,psToFit,ps)
ics <- c(1-ps["s"]/ps["n"],ps["s"]/ps["n"],0)
names(ics) <- c("x","y","cumInc")
sol = lsoda(ics,0:10*365,caMRSA,fitPs)
inc <- incFromSol(sol)
plot(sol[,"time"],sol[,"y"]*fitPs["n"],type="l",ylim=c(0,max(sol[,"y"]*fitPs["n"])))
