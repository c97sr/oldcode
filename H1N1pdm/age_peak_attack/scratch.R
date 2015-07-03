# All these lines will require other bits of other scripts to be run

bps1 <- (1:35)*2+0.5
bps2 <- (1:70)+0.5
freq_hosp <- (hist(hosp.data$adm.wk,breaks=bps1,plot=FALSE))$counts
freq_sev <- (hist(hosp.data$adm.wk[ind_sev],breaks=bps1,plot=FALSE))$counts
sev_rate <- freq_sev / freq_hosp
plot(sev_rate)
hist(hosp.data$adm.wk,breaks=bps2,plot=TRUE,freq=TRUE)

points(solution$time,solution$dIv1,type="l",col="red")
points(solution$time,solution$dIv2,type="l",col="green")
points(solution$time,solution$dIv3,type="l",col="blue")
points(solution$time,solution$dIv4,type="l",col="magenta")
