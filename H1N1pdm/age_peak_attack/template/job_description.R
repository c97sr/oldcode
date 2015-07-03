# Objective for this script is to run just MCMC chain fitting the HKU hospitalization data
rm(list=ls(all=TRUE))
options(error=NULL)
require("odesolve")
require("date")

source(paste("~/Dropbox/svneclipse/H1N1pdm/funcs.R",sep=""))
source("baseline_HK_Hosp_Params.R")
load("data_bins_cats_fit.Rdata")

set.seed(875869043)

mcmcH1N1v1(	
		"params_to_fit.csv",
		"trace.csv",
		baseHK4,
		PdmLikeHKHosp,
		no_samples=50,
		samp_freq=1,
		report_block=10,
		casemix=c(-1,-1,-1,-1),
		obs=agdata_fit,
		vecbins=wkBins_fit,
		agnames=agcats
)		