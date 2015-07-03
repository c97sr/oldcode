require("odesolve")
require("lattice")

posCh <- function(t,v) {
    if (t==0) rtn <- 0
    else rtn <- v
    rtn
}

updatePVector <- function(values,names,pvec) {                                                                        
    no_updates <- length(values)
    for (i in 1:no_updates) pvec[names[i]]=values[i]
    pvec <- jeUpdateAuxParams(pvec)
    pvec
}

jeModel <- function(t,y,p) {

    # set up the variables
    rtn <- array(0,c(p["nVar"]))
    N_m <- sum(y[1:3])
    N_p <- sum(y[jePI("M",1,1,p):jePI("R",p["n_a"],p["n_i"],p)])
    NI_p <- 0 
    for (a in 1:p["n_a"]) NI_p <- NI_p + y[jePI("I",a,1,p)]
    nu_m <- fnNu_m(t,p)

    # Forces of infection
    lambda_m = ( p["b"] * p["omega"] ) * (NI_p / N_p) * p["p_m"]
    lambda_p = ( p["b"] * p["omega"] / N_p) * y[jeMI("I")] * p["p_p"]
    
    # Model definition
    rtn[jeII("M",1,p)] <- lambda_m * y[jeMI("S")]
    rtn[jeMI("S")] <- nu_m * N_m - ( lambda_m + p["mu_m"] ) * y[jeMI("S")]
    rtn[jeMI("E")] <- lambda_m * y[jeMI("S")] - ( p["gamma_m"] + p["mu_m"] ) * y[jeMI("E")]
    rtn[jeMI("I")] <- p["gamma_m"] * y[jeMI("E")] - p["mu_m"] * y[jeMI("I")]
    for (a in 1:p["n_a"]) {
        for (i in 1:p["n_i"]) {
            rtn[jePI("M",a,i,p)] <- fnNu_p(a,i,p) * N_p + p["mu_age_p"] * (posCh(a-1,y[jePI("M",a-1,i,p)]) 
                - y[jePI("M",a,i,p)]) + p["mu_imm_p"] * (posCh(i-1,y[jePI("M",a,i-1,p)]) - y[jePI("M",a,i,p)]) 
        }
        rtn[jeII("P",a,p)] <- lambda_p * y[jePI("S",a,1,p)] 
        rtn[jePI("S",a,1,p)] <- p["mu_imm_p"]*y[jePI("M",a,p["n_i"],p)] -lambda_p * y[jePI("S",a,1,p)] + p["mu_age_p"]*(posCh(a-1,y[jePI("S",a-1,1,p)])-y[jePI("S",a,1,p)])
        rtn[jePI("I",a,1,p)] <- lambda_p * y[jePI("S",a,1,p)] - p["sigma_p"]  * y[jePI("I",a,1,p)] + p["mu_age_p"]*(posCh(a-1,y[jePI("I",a-1,1,p)])-y[jePI("I",a,1,p)])
        rtn[jePI("R",a,1,p)] <- p["sigma_p"] * y[jePI("I",a,1,p)] + p["mu_age_p"]*(posCh(a-1,y[jePI("R",a-1,1,p)])-y[jePI("R",a,1,p)])
    }
    
    #if (y[jeMI("S")] < 0) browser()

    list(rtn)

}

jeLike <- function(theta,names_theta,x,allps,saturated=FALSE) {
    p <- updatePVector(theta,names_theta,allps)
    ics <- jeMakeICS(p) 
    if (!saturated) {
        sol <- lsoda(ics,c(0,x$Day),jeModel,p,hmax=10)
        pp <- getPigPrevFinalAgeClass(sol,p)
    }
    sizeX <- dim(x)[1]
    lnlike <- 0
    for (i in 2:sizeX) {
        if (saturated) lnlike <- lnlike + dbinom(x$n[i-1],x$N[i-1],x$n[i-1]/x$N[i-1],log=TRUE)
        else lnlike <- lnlike + dbinom(x$n[i-1],x$N[i-1],pp$prev[i],log=TRUE)
    }
    fwdth=15
    no_fitted=length(theta)
    report <- format(lnlike,digits=fwdth)
    report <- paste(report,"c(")
    for (i in 1:no_fitted) {
        report <- paste(report,format(theta[i],digits=fwdth),sep="")
        if (i != no_fitted) report <- paste(report,",",sep="")
    }
    report <- paste(report,")",sep="")
    for (i in 1:no_fitted) report <- paste(report,names_theta[i])
    cat("jeLike",report,"\n")
    flush.console()
    plot(pp,type="l",ylim=c(0,1),main=theta)
    points(pigData$Day,pigData$Prev)
    lnlike
}

jeMakeICS <- function(p) {
    rtn <- array(0,p["nVar"])
    rtn[jeMI("S")] <- p["N_m"]/p["v_s"] - p["E_m0"]
    rtn[jeMI("E")] <- p["E_m0"]
    for (a in 1:p["n_a"]) rtn[jePI("S",a,1,p)] <- (p["N_p"] - p["I_p0"])/p["n_a"] 
    for (a in 1:p["n_a"]) rtn[jePI("I",a,1,p)] <- p["I_p0"]/p["n_a"] 
    rtn
}

jePI <- function(k,a,i,p) {
    n_states = 4
    if (k=="M") rtn <- 4 + (a-1)*p["n_i"]*n_states+(i-1)*n_states+0 else 
    if (k=="S") rtn <- 4 + (a-1)*p["n_i"]*n_states+(i-1)*n_states+1 else
    if (k=="I") rtn <- 4 + (a-1)*p["n_i"]*n_states+(i-1)*n_states+2 else
    if (k=="R") rtn <- 4 + (a-1)*p["n_i"]*n_states+(i-1)*n_states+3 else
    stop("Only states M,S,I, and R are allowed.")
    rtn
}

jeMI <- function(k,p) {
    if (k=="S") rtn <- 1 else
    if (k=="E") rtn <- 2 else
    if (k=="I") rtn <- 3 else
    stop("Problem here")
    rtn
}

jeII <- function(sp,a,p) {
    if (sp=="M") rtn <- 3+p["n_a"]*p["n_i"]*4 + 1
    if (sp=="P") rtn <- 3+p["n_a"]*p["n_i"]*4 + 1 + a
    rtn
}

jeUpdateAuxParams <- function(p) {
    sigma_tw = (1 + p["tw_g"]+ p["tw_h"] + p["tw_f"] + p["tw_l"])/360
    p["t_gs"] = 1/sigma_tw
    p["t_gf"] = p["t_gs"] + p["tw_g"]/sigma_tw
    p["t_fs"] = p["t_gf"] + p["tw_h"]/sigma_tw
    p["t_ff"] = p["t_fs"] + p["tw_f"]/sigma_tw
    D_mp1 <- p["t_gf"]-p["t_gs"]
    D_mp2 <- p["t_ff"]-p["t_fs"]
    p["D_mp1"] <- D_mp1
    p["D_mp2"] <- D_mp2
    p["alpha_1"] = p["mu_m"]+log(p["v_s"])/D_mp1
    if (p["mu_m"] > 1/D_mp2*log(1/p["v_s"])) p["alpha_2"] <- p["mu_m"] + 1/D_mp2*log(1/p["v_s"]) else stop("Implied alpha_2 value is not reasonable.")
    p["N_m0"] <- p["N_m"]/p["v_s"]
    p["mu_age_p"] <- p["n_a"] / p["a_max"] 
    p["mu_imm_p"] <- 1/p["ave_imm_p"]/p["n_i"]
    p["nVar"] <- 3+p["n_a"]*p["n_i"]*4 + 1 + p["n_a"]
    p["D_p"] <- 1/(p["sigma_p"]+1/p["a_max"])
    p["D_m"] <- p["gamma_m"]/p["mu_m"]/(p["gamma_m"]+p["mu_m"])
#    p["b"] <- p["R0Max"]/(p["omega"]*sqrt(p["D_p"]*p["D_m"]*p["N_m"]/p["N_p"]*p["p_m"]*p["p_p"])) # should be sqrt(R0Max)? according to the expression for R0
    p["b"] <- sqrt(p["R0Max"])/(p["omega"]*sqrt(p["D_p"]*p["D_m"]*p["N_m"]/p["N_p"]*p["p_m"]*p["p_p"])) # should be sqrt(R0Max)? according to the expression for R0
    p
}

# fnNu_m: New function with one more parameter (start of the cycle, assuming same shape)
#fnNu_m <- function(t,p) {
#    t <- t %% 360
fnNu_m <- function(t,p) {
    t <- (t+p["start"]) %% 360
    if (t >= p["t_gs"] && t < p["t_gf"]) {rtn <- p["alpha_1"]} else {if (t >= p["t_fs"] && t < p["t_ff"]) {rtn <- p["alpha_2"]} else {rtn <- p["mu_m"]}}
    rtn
}


fnNu_p <- function(a,i,p) {
    if (a==1 && i==1) rtn <- 1/p["a_max"]
    else rtn <- 0
    rtn
}

jeMosquitoPop <- function(t,p) {
    t <- t %% 360
    if (t >= p["t_gs"] && t < p["t_gf"]) {
        t <- t-p["t_gs"]
        rtn <- p["N_m0"]*exp(t*(p["alpha_1"]-p["mu_m"]))
    } else if (t >= p["t_fs"] && t < p["t_ff"]) {
        t <- t-p["t_fs"]
        rtn <- p["N_m"]*exp(t*(p["alpha_2"]-p["mu_m"]))
    } else if (t >= p["t_gf"] && t < p["t_fs"]) {
        rtn <- p["N_m"]
    } else {
        rtn <- p["N_m0"]
    }
    rtn
}

jeMosquitoLikelihood <- function(fitps,fitpn,allps,df,saturated=FALSE) {
    fwdth=5
    noData <- dim(df)[1]
    no_fitted <- length(fitps)
    lnlike <- 0
    satlike <- 0
    ps <- updatePVector(fitps,fitpn,allps)
    for (i in 1:noData) {
#        N_norm <- jeMosquitoPop(df$Days360_01_11_04[i],ps)
        N_norm <- jeMosquitoPop(df$Days360_01_11_04[i],ps)/ps["N_m"]   # standardized 
        if (df$Site[i]=="Luk_Keng") N_norm <- N_norm*ps["N_site_LK"]
        else if (df$Site[i]=="Tam_Kon_chau") N_norm <- N_norm*ps["N_site_TK"]
        else if (df$Site[i]=="Pak_Mong") N_norm <- N_norm*ps["N_site_PM"]
        else stop("Problem with the site cases")
        N_norm <- N_norm*df$No_Traps[i]
        lnlike <- lnlike+dpois(df$Mosquitoes_Collected[i],N_norm,log=TRUE)
        if (saturated) satlike <- satlike + dpois(df$Mosquitoes_Collected[i],df$Mosquitoes_Collected[i],log=TRUE)
    }
    if (ps["alpha_1"]*ps["D_mp1"] > ps["mu_m"]*(ps["D_mp1"]+ps["D_mp2"])) lnlike <- -10000000
    report <- format(lnlike,digits=fwdth)
    report <- paste(report,"c(")
    for (i in 1:no_fitted) {
        report <- paste(report,format(fitps[i],digits=fwdth),sep="")
        if (i != no_fitted) report <- paste(report,",",sep="")
    }
    report <- paste(report,")",sep="")
    for (i in 1:no_fitted) report <- paste(report,fitpn[i])
    cat("jeMosLike",report,"\n")
    flush.console()
    if (saturated) lnlike <- satlike
    lnlike
}

jeNormWithMax <- function(df,ps) {
    noData <- dim(df)[1]
    for (i in 1:noData) {
        N_norm = df$Mosquitoes_Collected[i]/df$No_Traps[i] 
        if (df$Site[i]=="Luk_Keng") N_norm = N_norm/ps["N_site_LK"]
        else if (df$Site[i]=="Tam_Kon_chau") N_norm = N_norm/ps["N_site_TK"]
        else if (df$Site[i]=="Pak_Mong") N_norm = N_norm/ps["N_site_PM"]
        df$NormMosPerTrap[i]=N_norm
    }
    df
}

getIncAux <- function(sol,cs,cf) {
    noPoints <- dim(sol)[1]-1
    noColumns <- cf - cs + 2
    rtn <- array(0,c(noPoints,noColumns))
    for (i in 1:noPoints) {
        dt <- (sol[i+1,1]-sol[i,1])/2
        rtn[i,1] = sol[i,1]+dt
        for (j in cs:cf) {
            rtn[i,2+j-cs] <- (sol[i+1,j] - sol[i,j]) / dt
        }
    }
    rtn
}

getInc <- function(sol,sp,ind,p) {
    noTPs = dim(sol)[1]
    
    if (sp=="M") rtn <- getIncAux(sol,jeII(sp,1,p)+1,jeII(sp,1,p)+1)
    else if (sp=="P") {
        rtn <- array(0,c(noTPs-1,2))
        rtn[,1] <- getIncAux(sol,jeII(sp,1,p)+1,jeII(sp,1,p)+1)[,1]
        for (a in ind) rtn[,2] <- rtn[,2] + getIncAux(sol,jeII(sp,a,p)+1,jeII(sp,a,p)+1)[,2] 
    }   
    rtn
}

getPopSizes <- function(sol,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,3)))
    names(rtn) <- c("t","M","P")
    rtn$t <- sol[,1]
    for (i in 1:noTPs) {
        rtn$M[i] <- sum(sol[i,2:4])
        for (a in 1:p["n_a"]) rtn$P[i] <- sum(solution[i,(jePI("M",1,1,p)+1):(jePI("R",p["n_a"],p["n_i"],p)+1)])
    }
    rtn
}

getPigPrevFinalAgeClass <- function(sol,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","prev")
    rtn$t <- sol[,1]
    for (j in 1:noTPs) {
        ptot <- 0
        pinf <- 0
        a <- p["n_a"]
        for (i in 1:p["n_i"]) {
            ptot <- ptot + sol[j,jePI("M",a,i,p)+1] + sol[j,jePI("S",a,i,p)+1]  + 
                sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
            pinf <- pinf + sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
        }
        rtn$prev[j] <- pinf / ptot
    } 
    rtn
}

getPigImmFinalAgeClass <- function(sol,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","prev")
    rtn$t <- sol[,1]
    for (j in 1:noTPs) {
        ptot <- 0
        pinf <- 0
        a <- p["n_a"]
        for (i in 1:p["n_i"]) {
            ptot <- ptot + sol[j,jePI("M",a,i,p)+1] + sol[j,jePI("S",a,i,p)+1]  + 
                sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
            pinf <- pinf + sol[j,jePI("M",a,i,p)+1]
        }
        rtn$prev[j] <- pinf / ptot
    } 
    rtn
}


getMosPrevFinalAgeClass <- function(sol,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","PrevM")
    rtn$t <- sol[,1]
    for (i in 1:noTPs) {
        rtn$PrevM[i] <- sol[i,4]/sum(sol[i,2:4])
    }
    rtn
}


getPigPrevAnyAgeClass <- function(sol,a,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","prev")
    rtn$t <- sol[,1]
    for (j in 1:noTPs) {
        ptot <- 0
        pinf <- 0
#        a <- p["n_a"]
        for (i in 1:p["n_i"]) {
            ptot <- ptot + sol[j,jePI("M",a,i,p)+1] + sol[j,jePI("S",a,i,p)+1]  + 
                sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
            pinf <- pinf + sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
        }
        rtn$prev[j] <- pinf / ptot
    } 
    rtn
}


getPigImmAnyAgeClass <- function(sol,a,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","prev")
    rtn$t <- sol[,1]
    for (j in 1:noTPs) {
        ptot <- 0
        pinf <- 0
#        a <- p["n_a"]
        for (i in 1:p["n_i"]) {
            ptot <- ptot + sol[j,jePI("M",a,i,p)+1] + sol[j,jePI("S",a,i,p)+1]  + 
                sol[j,jePI("I",a,i,p)+1] + sol[j,jePI("R",a,i,p)+1]
            pinf <- pinf + sol[j,jePI("M",a,i,p)+1]
        }
        rtn$prev[j] <- pinf / ptot
    } 
    rtn
}


getMosPop <- function(sol,p) {
    noTPs <- dim(sol)[1]
    rtn = as.data.frame(array(0,c(noTPs,2)))
    names(rtn) <- c("t","pop")
    rtn$t <- sol[,1]
    for (j in 1:noTPs) {
        rtn$pop[j] <- sol[j,jeII("M",1,p)]
    } 
    rtn
}
