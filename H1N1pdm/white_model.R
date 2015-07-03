# To dos
# - next: put the age specific  

rm(list=ls(all=TRUE))
options(error=NULL)

require("odesolve")
require("date")
source("../nh1n1/swine_flu_funcs.r")
source("../nh1n1/funcs.r")

# Define the different data structures from which to simulate

noDays <- 20
noTypes <- 2

tmp_params <- data.frame(
		R0=2,
		Tg=2.7
)
params <- as.vector(tmp_params)
names(params) <- names(tmp_params)

matNGM <- matrix(c(1,1,1,1),ncol=noTypes,nrow=noTypes)

matSeed <- matrix(NA,nrow=noDays,ncol=noTypes)
matSeed[,1] <- 1
matSeed[,2] <- 0

vecPop <- as.vector(c(1e10,1e10),mode="numeric")

SimModel(matNGM,matSeed,params,vecPop)

