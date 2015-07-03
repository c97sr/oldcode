#This is my first r script

# Checking distributions
x <- rlnorm(100000,0,0.65)
ed <- density(x)
y <- dlnorm(1:1000/100,0,0.65)
plot(1:1000/100,y,type="l")
points(ed)

amoy_sur <- function(x,mu,sigma) {
	## inthaz = plnorm(x,mu,sigma)
	inthaz = 1/mu*x
	rtnval = exp(-inthaz)
	rtnval
}

lm_func <- function(n,mu,sigma) {
	## draws = rnorm(n,0,1)
	## rtnval = exp(sigma*draws)
	draws = runif(n,0,1)
	rtnval = -log(draws)*mu	
	rtnval
}

list_grt <- function(list,vals) {
	nolist = length(list)
	novals = length(vals)
	rtnval = array(0,c(novals))
	for (i in 1:novals) {
		for (j in 1:nolist) if (list[j] > vals[i]) rtnval[i] <- rtnval[i]+1
	}
	rtnval <- rtnval / nolist
	rtnval
}

lm_list = lm_func(1000,3,0.65)
sur_points <- list_grt(lm_list,0:999/100)
plot(1:1000/100,amoy_sur(1:1000/100,3,0.65),type="n",ylim=c(0,1))
points(0:999/100,sur_points,pch=20)
points(1:1000/100,amoy_sur(1:1000/100,3,0.65),type="l",col="red")

# stuff I'm using for amoy
tmp_data <- read.table(file="D:\\share\\files\\projects\\inference2005\\results\\debug_param_samples.out",header=TRUE,sep="\t")
require("logspline")
par(mfrow=c(2,2))
plot(tmp_data[,4])
d1 <- logspline(tmp_data[100:999,1])
plot(d1,log="x")
d2 <- logspline(tmp_data[100:999,2])
plot(d2,log="x")
d3 <- logspline(tmp_data[100:999,3])
plot(d3,log="x")

# stuff for flu

root = "D:\\share\\files\\projects\\flu2005\\results\\050908\\scenB\\"
stem = "set"
leaf = "_Steven.txt"
pfile = "D:\\share\\files\\projects\\flu2005\\results\\050908\\ParameterList.txt"
params <- read.table(pfile)

flu_processor <- function(root,stem,leaf,policy,s,f) {
	## s for start and f for finish. Inclusive of both.
	if (s >= f) stop("f must be greater than s\n")
	colnames <- c("S_365","S_365_min","S_365_max")
	nocols <- length(colnames)
	norows <- f-s+1
	rtnval = array(0,dim=c(norows,nocols))
	rtnval <- as.data.frame(rtnval)
	names(rtnval) <- colnames 
	for (i in s:f) {
		tab <- read.table(file=paste(root,stem,i,leaf,sep=""),header=TRUE,sep="\t",row.names=NULL)
		names(tab) <- names(tab)[2:length(names(tab))]
		## figure out the stats stuff to go here to interrogate the data
		a = array(0,c(0))
		for (j in 1:nrow(tab)) {
			if (tab$Policy[j]==policy) a <- c(a,tab$Susceptibles_Left_On_Day_ThreeSixFive[j])
		}
		rtnval$S_365[i-s+1] <- mean(a)
		# cat(paste(root,stem,i,leaf,sep=""),"\n")
	}
	rtnval
}


check <- flu_processor(root,stem,leaf,"None",9,10)

leaf = "_Steven.txt"
output_data(root,stem,leaf,9,10)

# other useful stuff
alpha <- 0.5
help.start()
a <- split.screen(c(2,2))
screen(a[1])
hist(tmp_data[,3], xlim = c(1,10), probability=TRUE, axes=FALSE, main="", xlab="")
axis(1,pos=0,at=1:10)
axis(2,pos=1,at=(0:10)/10)
screen(a[2])
hist(tmp_data[,3], xlim = c(1,10), probability=TRUE, axes=FALSE, main="", xlab="")
axis(1,pos=0,at=1:10)
axis(2,pos=1,at=(0:10)/10)
hist(tmp_data[,3], xlim = c(1,10), probability=TRUE, axes=FALSE, main="", xlab="")
axis(1,pos=0,at=1:10)
axis(2,pos=1,at=(0:10)/10)
newfun <- function(data){data^2}
newfun(tmp_data[,3])
b <- (tmp_data[,3]>2)
tmp_data[b,3]
a <- logspline(tmp_data[,3])
root = "c:\\"
paste(root,"BANANA",sep="PEAR")
quantile(1:125)
text(1.6,0.1,"Household based quarantine",pos=4)
points(params$R0[1:100],results[1:100,1]/10000,cex=1,col="red")
plot(params$R0[1:100],res_QT[1:100,1]/10000,cex=1,xlim=c(0,10),type="n")
par(pch=20)
windows(width=7,height=4)
dev.list()