rm(list=ls(all=TRUE))
require("rgdal")
require("RArcInfo")

# fileroot = "D:\\files\\projects"
# setwd(paste(fileroot,"\\..\\..\\eclipse\\R Source",sep=""))

SR_version <- "windows" # Required for schist_load_data

source("schisto//schisto_funcs.r")
source("schisto//schist_load_data.r")
source("schisto//schist_params.r")

tabPlottableData <- makeAdjObsTable(    dataSumCIs, 
                                        c("x_L","x_H","x_C","x_D","x_P","x_W","x_R")  )
                                        
x_Hu <- tabPlottableData$x_L + tabPlottableData$x_H
tabPlottableData <- cbind(tabPlottableData,x_Hu)

# Little hack to get contrinutions to log likelihood

tab3 <- as.data.frame(array(NA,c(noSpecies,3)))
names(tab3) <- c("H0","H1","H2")
row.names(tab3) <- nameSpecies
 
# Set up the plot window for Figure 2
f2w <- 17.75/cm(1)
f2h <- 20/cm(1)
windows(width=f2w,height=f2h)
par(	mai=c(0,0,0,0), 		# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=TRUE,
		mgp=c(0,0,0),
		par(cex=1))
extractSet = 1:49
com_cex <- 1

# Paste model solution here
villlevelp = ""
colindex <- 1
ps_to_fit_names <- c(
                "beta_SM",
                "beta_MS",
                "gamma_L",
                "epsilon_H","omega_H","gamma_H",
                "epsilon_C","omega_C",    
                "epsilon_D","omega_D",    
                "epsilon_P","omega_P",    
                "epsilon_W","omega_W",    
                "epsilon_R","omega_R"	)  

# From first submission version from 070213
# Second submission version also from 070827/H0_final.out
pvals_MLE <- c(  	
	0.000742797	,
	2.817143262	,
	2.641030098	,
	1.341722414	,
	0.105764958	,
	0.150108557	,
	1.13E-05	,
	0.026904871	,
	0.003867743	,
	0.291019568	,
	0.23693769	,
	0.035223242	,
	0.086910648	,
	0.002726709	,
	0.862361332	,
	0.672911413	
)	

fitted_params = updatePVector(pvals_MLE,ps_to_fit_names,parvector)
output <- schistLike(pvals_MLE,ps_to_fit_names,parvector,datDemog,dataSum,villlevelp)
modelResults <- generateModelPlottableResults(datDemog,schistAdjMam,fitted_params,villlevelp)
tab3[,"H0"] <- glbTabVecComps

lgap <- (3.0/cm(1))/f2w
rgap <- (0.25/cm(1))/f2w
tgap <- (0.75/cm(1))/f2h
bgap <- (1.5/cm(1))/f2h
twnvert <- (0.15/cm(1))/f2h
twnhor <- (0.15/cm(1))/f2w

bxlab <- c("","","10","","30","","49","")
lylab <- c("","0.0001 ","","0.01 ","","1.0 ")

sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.01,1),
	sr_chart_pos(colindex,7,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_L[extractSet],modelResults$y_L[extractSet]),
    list(tabPlottableData$x_L_lci[extractSet],tabPlottableData$x_L_uci[extractSet]),
	argNew=FALSE,argYlab="Humans\n(Light)", 
    yat=c(0.01,0.1,1.0),ylab=c("","0.1 ","1.0 "))
    
mtext(expression(bold(H[0])),side=1,line=-6.5)
    
sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,6,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_H[extractSet],modelResults$y_H[extractSet]),
    list(tabPlottableData$x_H_lci[extractSet],tabPlottableData$x_H_uci[extractSet]),
    argYlab="Humans\n(Heavy)",ylab=lylab)
    
sr_schist_multiplot(tabPlottableData$x_C[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,5,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_C[extractSet],modelResults$y_C[extractSet]),
    list(tabPlottableData$x_C_lci[extractSet],tabPlottableData$x_C_uci[extractSet]),
    argYlab="Cats",ylab=lylab)
    
sr_schist_multiplot(tabPlottableData$x_D[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,4,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_D[extractSet],modelResults$y_D[extractSet]),
    list(tabPlottableData$x_D_lci[extractSet],tabPlottableData$x_D_uci[extractSet]),
    argYlab="Dogs",ylab=lylab)

sr_schist_multiplot(tabPlottableData$x_P[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,3,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_P[extractSet],modelResults$y_P[extractSet]),
    list(tabPlottableData$x_P_lci[extractSet],tabPlottableData$x_P_uci[extractSet]),
    argYlab="Pigs",ylab=lylab)
    
sr_schist_multiplot(tabPlottableData$x_W[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,2,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_W[extractSet],modelResults$y_W[extractSet]),
    list(tabPlottableData$x_W_lci[extractSet],tabPlottableData$x_W_uci[extractSet]),
    argYlab="Water\nBuffalo",ylab=lylab)

sr_schist_multiplot(tabPlottableData$x_R[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,1,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_R[extractSet],modelResults$y_R[extractSet]),
    list(tabPlottableData$x_R_lci[extractSet],tabPlottableData$x_R_uci[extractSet]),
    argYlab="Rats",ylab=lylab,xlab=bxlab)

mtext("Village (ranked by overall species prevalence)",side=1,at=76,cex=1.2,font=2,adj=0.5,line=1.5)
mtext("Prevalence",side=2,at=10^13.5,cex=1.2,font=2,padj=0.5,line=5.5)

# Paste model solution here
villlevelp = "beta_MS"
colindex <- 2
ps_to_fit_names <- c(
                "beta_SM",
                "gamma_L",
                "epsilon_H","omega_H","gamma_H",
                "epsilon_C","omega_C",    
                "epsilon_D","omega_D",    
                "epsilon_P","omega_P",    
                "epsilon_W","omega_W",    
                "epsilon_R","omega_R")  

# From first submission version from 070402/H1 _final.out
# Second submission version also from 070829/H1 _final.out
pvals_MLE <- c(
	0.092976451	,
	2.544063098	,
	0.009959121	,
	0.040031566	,
	0.082462588	,
	0.482567222	,
	0.028707275	,
	0.061965658	,
	0.298159777	,
	2.36E-07	,
	0.039589998	,
	0.062959318	,
	0.002526189	,
	0.043942048	,
	0.385012959	
)	

fitted_params = updatePVector(pvals_MLE,ps_to_fit_names,parvector)
output <- schistLike(pvals_MLE,ps_to_fit_names,parvector,datDemog,dataSum,villlevelp)
modelResults <- generateModelPlottableResults(datDemog,schistAdjMam,fitted_params,villlevelp)

tab3[,"H1"] <- glbTabVecComps

sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.01,1),
	sr_chart_pos(colindex,7,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_L[extractSet],modelResults$y_L[extractSet]),
    list(tabPlottableData$x_L_lci[extractSet],tabPlottableData$x_L_uci[extractSet]),
    yat=c(0.01,0.1,1.0),ylab=c("","",""))
    
mtext(expression(bold(H[1])),side=1,line=-6.5)
    
sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,6,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_H[extractSet],modelResults$y_H[extractSet]),
    list(tabPlottableData$x_H_lci[extractSet],tabPlottableData$x_H_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_C[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,5,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_C[extractSet],modelResults$y_C[extractSet]),
    list(tabPlottableData$x_C_lci[extractSet],tabPlottableData$x_C_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_D[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,4,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_D[extractSet],modelResults$y_D[extractSet]),
    list(tabPlottableData$x_D_lci[extractSet],tabPlottableData$x_D_uci[extractSet]))

sr_schist_multiplot(tabPlottableData$x_P[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,3,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_P[extractSet],modelResults$y_P[extractSet]),
    list(tabPlottableData$x_P_lci[extractSet],tabPlottableData$x_P_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_W[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,2,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_W[extractSet],modelResults$y_W[extractSet]),
    list(tabPlottableData$x_W_lci[extractSet],tabPlottableData$x_W_uci[extractSet]))

sr_schist_multiplot(tabPlottableData$x_R[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,1,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_R[extractSet],modelResults$y_R[extractSet]),
    list(tabPlottableData$x_R_lci[extractSet],tabPlottableData$x_R_uci[extractSet]),
	xlab=bxlab)
 
# Third column 
villlevelp = "beta_SM"
colindex <- 3

ps_to_fit_names <- c(
                "beta_MS",
                "gamma_L",
                "epsilon_H","omega_H","gamma_H",
                "epsilon_C","omega_C",    
                "epsilon_D","omega_D",    
                "epsilon_P","omega_P",    
                "epsilon_W","omega_W",    
                "epsilon_R","omega_R"   )    

# From first submission version from 070326\H2 _final.out
# Second submission version e5 from 070828\H1 _final.out 
pvals_MLE <-  c(	
	0.523463557	,
	0.846554487	,
	0.587316926	,
	1.665173892	,
	1.012228082	,
	2.446661407	,
	0.080466731	,
	1.031010002	,
	0.901763767	,
	0.402658839	,
	0.106252732	,
	0.663879099	,
	0.008249639	,
	0.854909872	,
	1.927729519	
) 
                
fitted_params = updatePVector(pvals_MLE,ps_to_fit_names,parvector)
output <- schistLike(pvals_MLE,ps_to_fit_names,parvector,datDemog,dataSum,villlevelp)
modelResults <- generateModelPlottableResults(datDemog,schistAdjMam,fitted_params,villlevelp)
tab3[,"H2"] <- glbTabVecComps

sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.01,1),
	sr_chart_pos(colindex,7,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_L[extractSet],modelResults$y_L[extractSet]),
    list(tabPlottableData$x_L_lci[extractSet],tabPlottableData$x_L_uci[extractSet]),
    yat=c(0.01,0.1,1.0),ylab=c("","",""))
    
mtext(expression(bold(H[2])),side=1,line=-6.5)
    
sr_schist_multiplot(tabPlottableData$x_Hu[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,6,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_H[extractSet],modelResults$y_H[extractSet]),
    list(tabPlottableData$x_H_lci[extractSet],tabPlottableData$x_H_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_C[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,5,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_C[extractSet],modelResults$y_C[extractSet]),
    list(tabPlottableData$x_C_lci[extractSet],tabPlottableData$x_C_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_D[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,4,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_D[extractSet],modelResults$y_D[extractSet]),
    list(tabPlottableData$x_D_lci[extractSet],tabPlottableData$x_D_uci[extractSet]))

sr_schist_multiplot(tabPlottableData$x_P[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,3,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_P[extractSet],modelResults$y_P[extractSet]),
    list(tabPlottableData$x_P_lci[extractSet],tabPlottableData$x_P_uci[extractSet]))
    
sr_schist_multiplot(tabPlottableData$x_W[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,2,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_W[extractSet],modelResults$y_W[extractSet]),
    list(tabPlottableData$x_W_lci[extractSet],tabPlottableData$x_W_uci[extractSet]))

sr_schist_multiplot(tabPlottableData$x_R[extractSet],c(0.00001,1),
	sr_chart_pos(colindex,1,3,7,xg=twnhor,yg=twnvert,xlm=lgap,xrm=rgap,ytm=tgap,ybm=bgap),
    list(tabPlottableData$x_R[extractSet],modelResults$y_R[extractSet]),
    list(tabPlottableData$x_R_lci[extractSet],tabPlottableData$x_R_uci[extractSet]),
	xlab=bxlab)
	
# Figure 3
require("rgdal")
require("RArcInfo")

strLandscanFile <- "C:/offlinefiles/files/data/landscan/asia03/asia03/w001001.adf"
strVillCoordsFile <- "D:\\files\\projects\\schisto\\data\\barangay_coords\\bgy_coords_sp_r_format_less_b_220.txt"

villCoords = read.table(strVillCoordsFile,header=TRUE,sep="\t")
normVillP <- modelResults$vill_p / max(modelResults$vill_p)

# detail easting, detail westing..., wide eastin, wide westing...
de <- 124.3
dw <- 124.3 + 13/12
dn <- 12.3
ds <- 12.3 - 13/12

we <- 119
ww <- 127 
wn <- 19
ws <- 4

spd <- 	120 	# Squares per degree
mn  <- 	85 		# map north
mw  <- 	34 		# map west

info <- GDALinfo(strLandscanFile)
wide <- readGDAL(strLandscanFile,region.dim=c((wn-ws)*spd,(ww-we)*spd),offset=c((mn-wn)*spd,(we-mw)*spd))
detail <- readGDAL(strLandscanFile,region.dim=c((dn-ds)*spd,(dw-de)*spd),offset=c((mn-dn)*spd,(de-mw)*spd))

f3w <- 8.6/cm(1)
widew <- 3/cm(1)	# Detail width
det_border <- 0.5/cm(1)
f3h <- (f3w-det_border)*(dn-ds)/(dw-de)+det_border
posdet <- c(0.01,0.99-det_border/f3w,0.01,0.99-det_border/f3h) 
poswide <- c(0.99-widew/f3w,0.99,0.99-widew*(wn-ws)/(ww-we)/f3h,0.99)

windows(width=f3w,height=f3h)
par(	mai=c(0,0,0,0), 		# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,
		mgp=c(0,0,0),
		par(cex=1))

par(fig=posdet)
# image(detail,col=grey(0.9-0.8*0:100/100))
image(detail,col="grey")
for (i in 1:4) axis(i,at=c(0,130),tck=0,col="blue",lwd=2)
points(villCoords$x,villCoords$y,cex=normVillP*2+0.1,pch=2,col="red")
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
polygon(c(0,1,1,0),c(0,0,1,1),col="white",border=NA)
par(new=TRUE,fig=poswide)
image(wide,col="grey")
#image(wide,col=grey(0.9-0.8*0:100/100))
for (i in 1:4) axis(i,at=c(0,130),tck=0)
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
lines(1-c(dw-ww,de-ww,de-ww,dw-ww,dw-ww)/(we-ww),c(ds-ws,ds-ws,dn-ws,dn-ws,ds-ws)/(wn-ws),col="blue",lwd=2)
