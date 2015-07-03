x3=1:999/10

list_grt <- function(list,vals) {
y={}
for (i in 1:length(vals)) y=c(y,sum(list>vals[i]))
y/length(vals)
}

hazard_sur <-function(x,mu,sigma) {
1-(1-exp(-plnorm(x,mu,sigma)))/(1-exp(-1))
}

hazard_ran <-function(n,mu,sigma) {
exp( mu+sigma*qnorm(-log( 1 - runif(n)*(1-exp(-1)) )))
}

plot(x3,hazard_sur(x3,3,0.65),type="l",ylim=c(0,1))
points(x3,list_grt(hazard_ran(length(x3),3,0.65),x3),col="red",pch=".")


