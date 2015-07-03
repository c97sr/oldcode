require("coda")

wk_dir = "D:\\share\\files\\projects\\je2005\\results\\051129b\\"
filename = "051129b_samples.out"
full_fn = paste(wk_dir,filename,sep="")

data <- read.table(file=full_fn,header=TRUE,sep="\t",row.names=NULL)
freq_thin = 10
no_samples = dim(data)[1]
prop_to_drop = 0.6
start = no_samples*prop_to_drop
mc_data = mcmc(data=data[start:no_samples,],thin=freq_thin,start=0)
plot(mc_data,ask=dev.interactive())

crosscorr.plot(mc_data)
summary(mc_data)
rejectionRate(mc_data)

no_samples = dim(data)[1]
require("logspline")
x1 <- logspline(data$omega_p)

plot(x1,log="x",xlim=c(1,500))

# test_breaks
# x <- hist(data$number_m,breaks=)
