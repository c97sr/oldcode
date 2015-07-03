require("coda")

x <- read.csv("~/Dropbox/results/100511/trace.csv")
size <- dim(x)[1]
burnin <- size*0.2
x_sub <- x[burnin:size,]
plot(x_sub$t0,type="l")
plot(x_sub$r,x_sub$p_H_base)

# Stuff suing code
xd <- mcmc(x_sub)
summary(xd)
effectiveSize(xd)
plot(xd)

mle_loc <- which.max(x_sub$lnlike)
x_sub[mle_loc,c("phi_1")]
