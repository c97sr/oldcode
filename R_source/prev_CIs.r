rm(list=ls(all=TRUE))

pOptim <- function(p,cdfVal,number,pos) {
    powp <- 10^p
    guess <- pbinom(pos,number,powp,lower.tail=TRUE)
    rtnval <- abs(guess-cdfVal)
    rtnval[1]
}

pOptim <- function(p,P,N,n) {
	guess <- pbinom(n,N,P,lower.tail=TRUE)
	rtn <- P-guess
	rtn
}

optimise(pOptim,interval=c(-5,0),maximum=FALSE,cdfVal=0.025,number=446,pos=0)

topdir = "D:\\tmp\\"
input_file = paste(topdir,"gavin_data.csv",sep="")
output_file = paste(topdir,"gavin_output.csv",sep="")

pD = read.table(input_file,header=TRUE,sep=",")   # poultry data

noData <- dim(pD)[1]

for (i in 1:noData) {
    noPos <- pD$pos[i]
    noSamp = pD$number[i]
    pAve = noPos / noSamp
    if (noPos != 0) searchInt <- c(log(pAve,10)-1,log(pAve,10)+1) else searchInt = c(log(5/noSamp,10)-3,log(5/noSamp,10))
    if (noPos==0) pD$lower_ci[i] <- 0 else {
        tmp <- optimise(pOptim,interval=searchInt,maximum=FALSE,cdfVal=0.975,number=noSamp,pos=noPos)$minimum
        pD$lower_ci[i] <- 10^tmp
    }
    tmp <- optimise(pOptim,interval=searchInt,maximum=FALSE,cdfVal=0.025,number=noSamp,pos=noPos)$minimum
    pD$upper_ci[i] <- 10^tmp
}

write.table(pD,file = output_file, sep=",", row.names=FALSE)
