rm(list=ls(all=TRUE))
source("amoy_time_since_functions.r")

# Next thing requried here is for the input to take generic ParameterUpdate
# file and assemble the mcmc object with correct headings, thin interval.
# Do this with the results from fitsw to the large data set!

sample_file = "D:\\files\\projects\\sars\\amoy\\results\\070504b\\debug_param_samples.out"
output = read.table(sample_file,header=TRUE)
nosamples = dim(output)[1]
onetenth = round(nosamples/10)
voi="q_SSE"
plot(output[,voi],type="l",log="y",ylim=c(0.0001,10))
plot(output[(nosamples-9.5*onetenth):nosamples,voi],type="l")
plot(output[,voi],type="l")

plot(output[(nosamples-3*onetenth):nosamples,"q_SSE"],output[(nosamples-3*onetenth):nosamples,"q_ptp"],type="l")
plotMCMCDensity(output,voi,burnin=10000,vmin=0.001,vmax=100000,nk=20)


