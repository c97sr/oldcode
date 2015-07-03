require(odesolve)
require(Hmisc)

nameSpecies <- c("Humans","Cats","Dogs","Pigs","WaterB","Rats")
noSpecies <- length(nameSpecies)
glbTabVecComps <- array(NA,c(noSpecies))

vill_assign <- function(pms,df,index) {                                                                               
    pms["N_H"] = df$N_H[index]
    pms["N_D"] = df$N_D[index]
    pms["N_C"] = df$N_C[index]
    pms["N_P"] = df$N_P[index]
    pms["N_W"] = df$N_W[index]
    pms["N_S"] = df$N_S[index]
    pms
}

FOI_to_snails <- function(y,p) {
    rtn <-  p["beta_MS"] * (   	p["epsilon_L"]*y["I_L"] +
               					p["epsilon_H"]*y["I_H"] + 
               					p["epsilon_D"]*y["I_D"] +
               					p["epsilon_C"]*y["I_C"] +
               					p["epsilon_P"]*y["I_P"] +
               					p["epsilon_W"]*y["I_W"] +
               					p["epsilon_R"]*y["I_R"] ) 
    rtn
}

schistAdjMam <- function(t,y,p) {                                                                                     
    rtn 			<- 	c(0,0,0,0,0,0,0,0,0) 
    names(rtn) 		<- 	c("I_L","I_H","E_S","I_S","I_D","I_C","I_P","I_W","I_R") 
    names(y) 		<- 	c("I_L","I_H","E_S","I_S","I_D","I_C","I_P","I_W","I_R")
    aux 			<- 	c(0,0,0,0)
    names(aux) 		<- 	c("y_H_L","y_H_M","y_S","y_D")
 
    rtn["I_L"]  	<-  p["beta_SM"]*p["omega_L"]*y["I_S"]*(p["N_H"]-y["I_L"]-y["I_H"]) 
                            - p["beta_SM"]*p["omega_H"]*y["I_L"]*y["I_S"]
                            - p["gamma_L"]*y["I_L"] + p["gamma_H"]*y["I_H"]
    rtn["I_H"]  	<-  p["beta_SM"]*p["omega_H"]*y["I_L"]*y["I_S"] - p["gamma_H"]*y["I_H"]
    rtn["E_S"]  	<-  FOI_to_snails(y,p)*(p["N_S"]-y["E_S"]-y["I_S"]) - p["sigma_S"]*y["E_S"]
    rtn["I_S"]  	<-  p["sigma_S"]*y["E_S"] - p["alpha_S"]*y["I_S"]
    rtn["I_D"]  	<-  p["beta_SM"]*p["omega_D"]*y["I_S"]*(p["N_D"]-y["I_D"]) - p["gamma_D"]*y["I_D"]
    rtn["I_C"]  	<-  p["beta_SM"]*p["omega_C"]*y["I_S"]*(p["N_C"]-y["I_C"]) - p["gamma_C"]*y["I_C"] 
    rtn["I_P"]  	<-  p["beta_SM"]*p["omega_P"]*y["I_S"]*(p["N_P"]-y["I_P"]) - p["gamma_P"]*y["I_P"] 
    rtn["I_W"]  	<-  p["beta_SM"]*p["omega_W"]*y["I_S"]*(p["N_W"]-y["I_W"]) - p["gamma_W"]*y["I_W"] 
    rtn["I_R"]  	<-  p["beta_SM"]*p["omega_R"]*y["I_S"]*(p["N_R"]-y["I_R"]) - p["gamma_R"]*y["I_R"] 

    aux["y_H_L"] 	<- 	y["I_L"]/p["N_H"]
    aux["y_H_M"] 	<- 	y["I_H"]/p["N_H"]
    aux["y_D"] 		<- 	y["I_D"]/p["N_D"] 
    aux["y_S"] 		<- 	y["I_S"]/p["N_S"]

    list(rtn,aux)

}

fnRestGivenIs <- function(I_S,p) { # From multisite_simple.nb
    rtn <- c(0,0,0,0,0,0,0,0,0) 
    names(rtn) <- c("I_L","I_H","E_S","I_S","I_D","I_C","I_P","I_W","I_R") 
    rtn["I_L"] =  I_S*p["N_H"]*p["beta_SM"]*p["omega_L"]*p["gamma_H"] / 
                (   I_S*I_S*p["beta_SM"]*p["beta_SM"]*p["omega_L"]*p["omega_H"] +
                    I_S*p["beta_SM"]*p["omega_L"]*p["gamma_H"] +
                    p["gamma_L"]*p["gamma_H"])
    rtn["I_H"] =  I_S*I_S*p["N_H"]*p["beta_SM"]*p["beta_SM"]*p["omega_L"]*p["omega_H"] /
                (   I_S*I_S*p["beta_SM"]*p["beta_SM"]*p["omega_L"]*p["omega_H"] +
                    I_S*p["beta_SM"]*p["omega_L"]*p["gamma_H"] +
                    p["gamma_L"]*p["gamma_H"])
    rtn["E_S"] = I_S*p["alpha_S"]/p["sigma_S"]
    rtn["I_S"] = I_S
    rtn["I_D"] = I_S*p["N_D"]*p["omega_D"]*p["beta_SM"]/(I_S*p["omega_D"]*p["beta_SM"]+p["gamma_D"])
    rtn["I_C"] = I_S*p["N_C"]*p["omega_C"]*p["beta_SM"]/(I_S*p["omega_C"]*p["beta_SM"]+p["gamma_C"])
    rtn["I_P"] = I_S*p["N_P"]*p["omega_P"]*p["beta_SM"]/(I_S*p["omega_P"]*p["beta_SM"]+p["gamma_P"])
    rtn["I_W"] = I_S*p["N_W"]*p["omega_W"]*p["beta_SM"]/(I_S*p["omega_W"]*p["beta_SM"]+p["gamma_W"])
    rtn["I_R"] = I_S*p["N_R"]*p["omega_R"]*p["beta_SM"]/(I_S*p["omega_R"]*p["beta_SM"]+p["gamma_R"])
    rtn
}

fnSchistIsAnalytical <- function(I_S,p) {                                                                             
    rtn <- fnRestGivenIs(I_S,p)
    delta_current = abs(    FOI_to_snails(rtn,p)*(p["N_S"]-rtn["E_S"]-rtn["I_S"]) - 
                                p["sigma_S"]*rtn["E_S"] )
    list(delta_current,rtn)
}

fnSchistFullAnalytical <- function(p) {
    I_S <- fnIsGivenP(p)
    rtn <- fnRestGivenIs(I_S,p)
    rtn
}

fnSchistIsAnalyticalIsOnly <- function(I_S,p) {                                                                       
    x <- fnSchistIsAnalytical(I_S,p)
    x[[1]]
}

fnSchistAnalyticalSteady <- function(pars) {
  testmin = optimise(fnSchistIsAnalyticalIsOnly,c(0,pars["N_S"]),tol=0.001,p=pars)
}

fnSteadyState <- function(ics,params,model,tol,dt,max_it) {                                                           
    current_tol <- tol+1
    current_it <- 0
    no_vars <- length(ics)
    current_state <- ics
    names(current_state) <- c("I_L","I_H","E_S","I_S","I_D","I_C","I_P","I_W","I_R")
    while ((current_tol > tol) && (current_it < max_it)) {
        new_state <- lsoda(current_state,c(0,dt),model,params,rtol=1e-4)[2,2:(2+no_vars-1)]
        current_tol <- sqrt(sum((new_state-current_state)^2)/no_vars)
        current_state <- new_state
        current_it = current_it+1
    }
    rtn = list(current_state,current_it,current_it>max_it)
    names(rtn) <- c("state","iterations","max_it_exceeded")
    rtn
}

# debug(fnVillageLike)
fnVillageLike <- function(villpv,all_ps,sumRes,villindex,villpn,updtsp=FALSE) {

    all_ps[villpn]=villpv

    lnlike <- 0

    minsol = fnSchistAnalyticalSteady(all_ps)
    out <- fnSchistIsAnalytical(minsol$minimum,all_ps)[[2]]

    model_ps_h <- c(  	all_ps["N_H"]-out["I_L"]-out["I_H"],       out["I_L"],   out["I_H"]) / all_ps["N_H"]
    model_ps_c <- c(    all_ps["N_C"]-out["I_C"],       out["I_C"]) / all_ps["N_C"]
    model_ps_d <- c(    all_ps["N_D"]-out["I_D"],       out["I_D"]) / all_ps["N_D"]
    model_ps_p <- c(    all_ps["N_P"]-out["I_P"],       out["I_P"]) / all_ps["N_P"]
    model_ps_w <- c(    all_ps["N_W"]-out["I_W"],       out["I_W"]) / all_ps["N_W"]
    model_ps_r <- c(    all_ps["N_R"]-out["I_R"],       out["I_R"]) / all_ps["N_R"]
    
    model_ys <- list(model_ps_h, model_ps_c,model_ps_d,model_ps_p,model_ps_w,model_ps_r)
    noSpecies <- length(sumRes)
  
    for (i in 1:noSpecies) {
    	dlnlk <- dmultinom(round(sumRes[[i]][villindex,]),prob=model_ys[[i]],log=TRUE)
        lnlike <- lnlike + dlnlk
        if (updtsp) glbTabVecComps[i] <<- glbTabVecComps[i] + dlnlk
    }
    lnlike
}

fnVillageLikeFit <- function(all_ps_in,allDemog,sumRes,villindex,villpn) {
    vill_level_pars <- vill_assign(all_ps_in,allDemog,villindex)
    if (villpn=="") lnlike <- fnVillageLike(0,vill_level_pars,sumRes,villindex,"",updtsp=TRUE)
    else {
        optResult = optimise(fnVillageLike,c(0,1),tol=0.000001, maximum=TRUE,
            all_ps=vill_level_pars,villindex=villindex,sumRes=sumRes,villpn=villpn)
        lnlike = optResult$objective
        global_villpv[villindex] <<- optResult$maximum
        dummy <- fnVillageLike(optResult$maximum,vill_level_pars,sumRes,villindex,villpn,updtsp=TRUE)
    }
    lnlike
}

schistLike <- function(fitted_ps,fitted_names,all_ps,allDemog,sumRes,villpn) {    
    lnlike <- 0
	glbTabVecComps[] <<- 0  
    fwdth=5
    goodvalue=TRUE  
    no_fitted = length(fitted_ps)
    no_data = dim(sumRes[[1]])[1]
    fit_to_list = 1:(no_data)           
    # fit_to_list = c(1,25,49)          
    # fit_to_list = c(1)
    
    for(i in 1:no_fitted) {
        all_ps[fitted_names[i]]=fitted_ps[i]
        if (fitted_ps[i] < 0) goodvalue=FALSE          
    }
    
    if (goodvalue) {
            for (i in fit_to_list) { 
                    lnlike <- lnlike + fnVillageLikeFit(all_ps,allDemog,sumRes,i,villpn)
            }
    } else {
        lnlike <- -1000000
    }
    
    report <- format(lnlike,digits=fwdth)
    report <- paste(report,"c(")
    for (i in 1:no_fitted) {
        report <- paste(report,format(fitted_ps[i],digits=fwdth),sep="")
        if (i != no_fitted) report <- paste(report,",",sep="")
    }
    report <- paste(report,")",sep="")
    for (i in 1:no_fitted) report <- paste(report,fitted_names[i])
    cat(villpn,report,"\n")
    flush.console()
    lnlike
}

schistLikeProfileScore <- function(pval,pindex,baselnlike,thresh,fitted_ps,ps_to_fit_names,parvector,datDemog,dataSum,villlevelp) {
    fitted_ps[pindex]<-pval
    lnlike <- schistLike(fitted_ps,ps_to_fit_names,parvector,datDemog,dataSum,villlevelp)
    rtn <- (baselnlike-thresh-lnlike)*(baselnlike-thresh-lnlike)
    rtn
}

updatePVector <- function(values,names,pvec) {                                                                        
    no_updates <- length(values)
    for (i in 1:no_updates) pvec[names[i]]=values[i]
    pvec
}

generateModelPlottableResults <- function(df,mod,pms,villlevelp) {   
 
    # df <- datDemog
    # mod <- schistAdjMam
    # pms <- fitted_params
    noData <- dim(df)[1]
    rtn <- as.data.frame(array(0,c(noData)))
    names(rtn)[1] <- "Village"
    rtn[] <- df[,"Baranguay"]
    
    for (i in 1:noData) {
        pms <- vill_assign(pms,df,i)
        pms[villlevelp]=global_villpv[i]
        minsol = fnSchistAnalyticalSteady(pms)
        out <- fnSchistIsAnalytical(minsol$minimum,pms)[[2]]
        rtn$y_L[i] = out["I_L"]/pms["N_H"]
        rtn$y_H[i] = out["I_H"]/pms["N_H"]
        rtn$y_C[i] = out["I_C"]/pms["N_C"]
        rtn$y_D[i] = out["I_D"]/pms["N_D"]
        rtn$y_P[i] = out["I_P"]/pms["N_P"]
        rtn$y_W[i] = out["I_W"]/pms["N_W"]
        rtn$y_R[i] = out["I_R"]/pms["N_R"]
        rtn$vill_p[i]=global_villpv[i]
    }
    rtn
}

lookupPattern <- function(pattern,status,patTable) {
    pattMask <- which(patTable$Pattern==pattern)
    statusMask <- which(patTable$True.Status[pattMask]==status)
    if (length(statusMask) != 1) stop("problem in Pattern Lookup Routine")
    rtn = patTable$median[pattMask[statusMask[1]]]
    rtn
}

makeAdjObsTable <- function(dsum,vecVariableNames) {

    rtn <- as.data.frame(rownames(as.data.frame(unclass(dsum[1][[1]]))))
    names(rtn)[1] <- "Village"
    noRows <- dim(rtn)[1]
    tmp <- array(0,c(noRows))
    suffixes <- c("","_lci","_uci")

    for (j in 1:3) {
    
        for (i in 1:noRows) tmp[i]=dsum[[1]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_L",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[1]][i,3+j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_H",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[2]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_C",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[3]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_D",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[4]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_P",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[5]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_W",suffixes[j],sep="") 

        for (i in 1:noRows) tmp[i]=dsum[[6]][i,j]
        rtn <- cbind(rtn,tmp)
        names(rtn)[dim(rtn)[2]] <- paste("x_R",suffixes[j],sep="")
        
    }
    rtn
    
}

makeAdjObsTableOld <- function(dsum,vecVariableNames) {

    rtn <- as.data.frame(rownames(as.data.frame(unclass(dsum[1][[1]]))))
    names(rtn)[1] <- "Village"
    noRows <- dim(rtn)[1]
    tmp <- array(0,c(noRows))

    for (i in 1:noRows) tmp[i]=dsum[[1]][i,2]/sum(dsum[[1]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_L" 

    for (i in 1:noRows) tmp[i]=dsum[[1]][i,3]/sum(dsum[[1]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_H" 

    for (i in 1:noRows) tmp[i]=dsum[[2]][i,2]/sum(dsum[[2]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_C" 

    for (i in 1:noRows) tmp[i]=dsum[[3]][i,2]/sum(dsum[[3]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_D" 

    for (i in 1:noRows) tmp[i]=dsum[[4]][i,2]/sum(dsum[[4]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_P" 

    for (i in 1:noRows) tmp[i]=dsum[[5]][i,2]/sum(dsum[[5]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_W" 

    for (i in 1:noRows) tmp[i]=dsum[[6]][i,2]/sum(dsum[[6]][i,])
    rtn <- cbind(rtn,tmp)
    names(rtn)[dim(rtn)[2]] <- "x_R" 

    rtn
    
}

makeAdjObsTableOldOld <- function(dsum,listDatA,listPatA,vecStatus,vecVariableNames) {

    rtn <- as.data.frame(rownames(as.data.frame(unclass(dsum[1][[1]]))))
    names(rtn)[1] <- "Village"
    noRows <- dim(rtn)[1]
    
    #noVariables <- length(listDatA)
    #if (noVariables != length(listPatA)) stop("Mismatched lists in makeAdjObsTable")
    
    for (i in 1:noVariables) {
        datA <- listDatA[i][[1]]
        patA <- listPatA[i][[1]]
        colNames <- names(as.data.frame(unclass(patA))) 
        noPatterns <- length(colNames)
        tmp <- array(0,c(noRows))
        tmp_size <- array(0,c(noRows))
        for (j in 1:noRows) {
            for (k in 1:noPatterns) {
                tmp[j] <- tmp[j] + patA[j,colNames[k]]*lookupPattern(colNames[k],vecStatus[i],datA)
                tmp_size[j] <- tmp_size[j]+patA[j,colNames[k]]
            }
        }
        rtn <- cbind(rtn,tmp/tmp_size)
        names(rtn)[dim(rtn)[2]] <- vecVariableNames[i] 
    }
    rtn
}

sr_chart_pos <- function(	xindex,yindex,xn,yn,xlm=0,xrm=0,
							xg=0.04,ybm=0,ytm=0,yg=0.01) {

    x_left_mar = xlm
    x_right_mar = xrm
    x_gap = xg
    x_n_charts = xn
    x_width = (1-x_left_mar-x_right_mar-x_gap*(x_n_charts-1))/x_n_charts
    x_charts_l=rep(0,x_n_charts)
    for (i in 1:x_n_charts) x_charts_l[i]=x_left_mar+(i-1)*x_width+(i-1)*x_gap
    x_charts_r=rep(0,x_n_charts)
    for (i in 1:x_n_charts) x_charts_r[i]=x_charts_l[i]+x_width

    y_top_mar = ytm
    y_bottom_mar = ybm
    y_gap = yg
    y_n_charts = yn
    y_height = (1-y_top_mar-y_bottom_mar-y_gap*(y_n_charts-1))/y_n_charts
    y_charts_b = rep(0,y_n_charts)
    for (i in 1:y_n_charts) y_charts_b[i] = y_bottom_mar+(i-1)*y_height+(i-1)*y_gap
    y_charts_t = rep(0,y_n_charts)
    for (i in 1:y_n_charts) y_charts_t[i] = y_charts_b[i]+y_height

    rtn <- c(x_charts_l[xindex],x_charts_r[xindex],y_charts_b[yindex],y_charts_t[yindex])
    
    rtn

}

sr_schist_multiplot <- function(	rks,yl,pos,lyv,lycis,
									argNew=TRUE,argYlab="",ylabline=4,
									xat=c(0,1,10,20,30,40,49,50),xlab=c("","","","","","","",""),
									yat=c(0.00001,0.0001,0.001,0.01,0.1,1),ylab=c("","","","","","")) {
 
    rankSeries <- rks
    ylims <- yl
    pos <- pos
    listYVals <- lyv
    listYCIs <- lycis
	tlen <- +0.03

    noPoints <- length(rankSeries)
    noSeries <- length(listYVals)
    noCIs <- length(listYCIs)
    if (argNew==TRUE) par(new=argNew) 
    par(fig=pos) 
    plot(1:2,xlim=c(0,noPoints+1),ylim=ylims,type="n",axes=FALSE,log="y",xlab="",ylab="")
    axis(1,at=xat,labels=xlab,tck=tlen)
    axis(2,at=yat,labels=ylab,tck=tlen,las=1)
    axis(3,at=xat,labels=rep("",length(xat)),tck=tlen)
    axis(4,at=yat,labels=rep("",length(yat)),tck=tlen)
    mtext(argYlab,side=2,line=ylabline,padj=0.5)
    for (i in 1:noCIs) for (j in 1:noPoints) 	lines(c(j,j),c(listYCIs[[1]][order(rankSeries)][j],listYCIs[[2]][order(rankSeries)][j]),
    											col="grey")
    for (i in 1:noSeries) points(1:noPoints,listYVals[i][[1]][order(rankSeries)],pch=c(20,20)[i],cex=c(1,1)[i],col=c("grey","red")[i])
}

summeriseCorrectedResults <- function(ltStat,ltPat,ltProb) {

    noSpecies <- length(ltStat)
    rtn <- list(noSpecies)
    noVillages <- dim(ltPat[[1]])[[1]]
    
    for (i in 1:noSpecies) {
        tabPat <- ltPat[[i]]
        noPatterns <- dim(tabPat)[2]
        tabProb <- ltProb[[i]]
        vecStat <- ltStat[[i]]
        noStatus <- length(vecStat)
        rtnItem <- array(0,c(noVillages,length(vecStat)))
        vecPatterns <- names(as.data.frame(unclass(tabPat)))
        for (j in 1:noVillages) {
            personCount <- 0
            for (k in 1:noPatterns) {    
                personCount <- personCount + tabPat[j,k]
                for (l in 1:noStatus) {
                    rtnItem[j,l] <- rtnItem[j,l] + tabPat[j,k]*lookupPattern(vecPatterns[k],vecStat[l],tabProb)
                }
            }
            rtnItem[j,] <- rtnItem[j,] / sum(rtnItem[j,]) * personCount
        }
        rtn[[i]] <- rtnItem
    }
    rtn 
}

calcConfidenceIntervals <- function(lstMeans) {
    lstMeans <- dataSum
    minLowerBound <- 1e-10
    lstMeans <- dataSum
    noSpecies <- length(lstMeans)
    rtn <- list(noSpecies)
    noVillages <- dim(lstMeans[[1]])[[1]] 
    for (i in 1:noSpecies) {
        tabMeans <- lstMeans[[i]]
        noMeans <- dim(tabMeans)[2]
        rtnItem <- array(0,c(noVillages,(noMeans-1)*3))
        for (j in 1:noVillages) {
            N <- sum(tabMeans[j,])
            for (k in 2:noMeans) {
                n <- tabMeans[j,k]
                cis <- binconf(n,N)
                rtnItem[j,(k-2)*3+1]=cis[1]
                rtnItem[j,(k-2)*3+2]=cis[2]
                rtnItem[j,(k-2)*3+3]=cis[3]                                 
            }
        } 
        rtn[[i]] <- rtnItem
    }
    rtn
}


fnSrRunif <- function(min=0, max=1,log=FALSE) {
	rv = runif(1)
	if (log) rtn <- min*10^(rv*(log10(max)-log10(min)))
	else rtn <- min + (max - min) * rv
	rtn 
}
