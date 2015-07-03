rm(list=ls(all=TRUE))

SR_version="linux"

source("schisto_funcs.r")
source("schist_load_data.r")
source("schist_params.r")

fnSrRunif <- function(min=0, max=1,log=FALSE) {
	rv = runif(1)
	if (log) rtn <- min*10^(rv*(log10(max)-log10(min)))
	else rtn <- min + (max - min) * rv
	rtn 
}

datDemog["N_S"] <- round(runif(dim(datDemog)[1],min=1000,max=3000))

# What will eventually be function arguments
of_stem <- "/mnt/windesk/projects/schisto/results/070402/H1_rand"
villlevelp = "beta_MS"
ps_to_fit_names <- c(
                "beta_SM",
                "gamma_L",
                "epsilon_H","omega_H","gamma_H",
                "epsilon_C","omega_C",    
                "epsilon_D","omega_D",    
                "epsilon_P","omega_P",    
                "epsilon_W","omega_W",    
                "epsilon_R","omega_R"   ) 

ranseed <- 120987234
no_starting_values <- 30
range_min <- 1e-5
range_max <- 10
max_iterations <- 10000

array_start_vals <- paste("s_",1:no_starting_values,sep="")
array_end_vals <- paste("e_",1:no_starting_values,sep="")
col_headers <- c(array_start_vals,array_end_vals,"minr","maxr","log","lb","ub")
row_names <- c(ps_to_fit_names,"lnlike","fns","conv")
no_tab_columns <- length(col_headers)
no_tab_rows <- length(row_names)
no_params <- length(ps_to_fit_names)
ps_start <- array(no_params)
outTab <- as.data.frame(array(NA,c(no_tab_rows,no_tab_columns)))
names(outTab) <- col_headers
rownames(outTab) <- row_names

for (i in 1:no_params) {
	outTab[i,"log"]=TRUE
	outTab[i,"minr"]=range_min
	outTab[i,"maxr"]=range_max
}

set.seed(ranseed)
for (i in 1:no_params) {
	for (j in 1:no_starting_values) {
		outTab[i,j] = fnSrRunif(min=outTab[i,"minr"],max=outTab[i,"maxr"],log=outTab[i,"log"]) 
	}
}

write.table(outTab,file=paste(of_stem,"_initial.out"),sep="\t",col.names = NA,row.names = TRUE)

for (i in 1:no_starting_values) {
	strStartCol = paste("s_",i,sep="")
	strEndCol = paste("e_",i,sep="")
	for (j in 1:no_params) ps_start[j]=outTab[ps_to_fit_names[j],strStartCol]
	schistEst <- optim(ps_start,schistLike,control=list(trace=0,fnscale=-1,maxit=max_iterations),
        fitted_names=ps_to_fit_names,all_ps = parvector,allDemog=datDemog,villpn=villlevelp,
        sumRes=dataSum)
    outTab[1:no_params,strEndCol] = schistEst$par
    outTab["lnlike",strEndCol] = schistEst$value
    outTab["fns",strEndCol] = schistEst$counts["function"]
    outTab["conv",strEndCol] = schistEst$convergence
    write.table(outTab,file=paste(of_stem,"_intermediate_par_",i,".out"),sep="\t",col.names = NA,row.names = TRUE)
}

maxCol = paste("e_",which.max(outTab["lnlike",(no_starting_values+1):(2*no_starting_values+1)]),sep="")

for (i in 1:no_params) {
	outTab[i,"lb"] <-   (optimise(schistLikeProfileScore,c(outTab[i,"minr"],outTab[i,maxCol]),tol=outTab[i,maxCol]*0.01,
                   		pindex=i,baselnlike=outTab["lnlike",maxCol],thresh=1.96,fitted_ps=outTab[1:no_params,maxCol],
                    	ps_to_fit_names=ps_to_fit_names,parvector=parvector,datDemog=datDemog,dataSum=dataSum,villlevelp=villlevelp))$minimum
	outTab[i,"ub"] <-   (optimise(schistLikeProfileScore,c(outTab[i,"maxr"],outTab[i,maxCol]),tol=outTab[i,maxCol]*0.01,
                   		pindex=i,baselnlike=outTab["lnlike",maxCol],thresh=1.96,fitted_ps=outTab[1:no_params,maxCol],
                    	ps_to_fit_names=ps_to_fit_names,parvector=parvector,datDemog=datDemog,dataSum=dataSum,villlevelp=villlevelp))$minimum
    
    write.table(outTab,file=paste(of_stem,"_intermediate_ci_",i,".out"),sep="\t",col.names = NA,row.names = TRUE)
	
}

write.table(outTab,file=paste(of_stem,"_final.out"),sep="\t",col.names = NA,row.names = TRUE)
