rm(list=ls())
# options(error=recover)
source("id_inference/funcs.r")

simData <- simSEIR(stochastic=TRUE)
detModel <- simSEIR(stochastic=FALSE)
approx_like_1 <- sum(dpois(simData$ave,detModel$ave,TRUE))
