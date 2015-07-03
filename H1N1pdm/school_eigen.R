# To dos: 
# - Next : tidy up the data code for plotting the confidence bounds then estimate the correct adjustment for the 
# - growth rate based on the daily data
# - next, setup growth rate

rm(list=ls(all=TRUE))
options(error=NULL)
source("SR_misc.r")

noSch 		<- 1100
schSize		<- 700
popsize 	<- 7000000
propSch		<- noSch*schSize/popsize
hvals 		<- c(1*1,0.5*2.5,1*2.5)
noHStates	<- length(hvals)
theta_cc	<- 0.58/2
theta_cs	<- 0.58/2
theta_ca	<- 1 - theta_cc - theta_cs
theta_aa	<- 0.79
theta_ac	<- 1 - theta_aa
phi			<- 2
beta		<- 0.002729466

ngm <- array(1,dim=c((noSch+1)*noHStates,(noSch+1)*noHStates)) 

for (i in 1:(noSch+1)) {
	for (j in 1:(noSch+1)) {
		for (k in 1:noHStates) {
			for (l in 1:noHStates) {
				val <- beta*hvals[k]
				if (j <= noSch) val <- val * phi
				if (i <= noSch && i==j) val <- val * theta_cs / schSize
				if (i <= noSch && i != j && j <= noSch) val <- val * theta_cs / (popsize * propSch)
				if (i <= noSch && i != j && j == noSch) val <- val * theta_ca / (popsize * (1-propSch))
				if (i == (noSch+1) && j <= noSch) val <- val * theta_ac / (popsize * propSch)
				if (i == (noSch+1) && i == (noSch+1)) val <- val * theta_aa / (popsize * (1-propSch))
				ngm[k+(i-1)*noHStates,l+(j-1)*noHStates] <- val
			}
		}
	}	
}

ev <- eigen(ngm)$value
revs <- ev[Im(ev)==0]
maxrev <- max(Re(revs))


