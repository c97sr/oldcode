rm(list=ls(all=TRUE))
options(recover=NULL)
require(date)

# Read in the UK flu data
rawData <- read.csv("/Volumes/NO NAME/data/influenza/uk2011/Data_13012011.csv")
rawData$Monday_of_week <- as.date(as.character(rawData$Monday_of_week))
rawData <- rawData[rawData$Monday_of_week < as.date("05-Jan-2011"),]

axdates <- c(	"15-Nov-2010","22-Nov-2010","29-Nov-2010","15-Nov-2010","06-Dec-2010",
				"13-Dec-2010","20-Dec-2010","27-Dec-2010","03-Jan-2011")

plot(	0:1,
		ylab="No DM pos pdm (+1, log scale)",
		xlab="Monday of week",
		xlim=c(min(as.date(axdates)),max(as.date(axdates))),
		type="n",
		ylim=c(1,2000),
		axes=FALSE,
		log="y"
)

axis(1,at=as.date(axdates),labels=axdates)
axis(2)

vecFields <- c(	"DatamartPanflu.pos_Total","DatamartPanflu.pos_0..05", 
				"DatamartPanflu.pos_05.14","DatamartPanflu.pos_15.44", 
				"DatamartPanflu.pos_45.","DatamartPanflu.pos_NK")
vecNames <- c("Total","0-4","5-14","15-44","45+","NK")
vecCols <- c("black","red","pink","orange","green","grey")
nolines <- length(vecFields)

for (i in 1:nolines) {

	points(rawData$Monday_of_week,rawData[,vecFields[i]]+1,type="l",col=vecCols[i],lwd=2)
	
}

legend(min(as.date(axdates)),y=2000,legend=vecNames,col=vecCols,lwd=2)
		
# 17th Jan 2011
agegroups <- as.vector(c(1119*1000, 2105*1000, 2367*1000, 1260*1000))
rates <- c(8.1,6.5,16.7,24.3)
ave_w <- sum(rates * agegroups / sum(agegroups))  
ave_uw <- sum(rates * c(1,1,1,1) / sum(1,1,1,1))  

# Next - axis labels and axis tagss

# Using matrix with named dimensions
# Using convention [x,y], i.e. first listed dimension is x however the thing prints
xdim <- 1:4
ydim <- c("a","b","c")
mat <- matrix(1:(length(xdim)*length(ydim)),nrow=length(xdim),ncol=length(ydim),dimnames=list(xdim,ydim))
mat[1,"a"]

# A very simple regression model
x <- data.frame(dist=c(2,1,2,3,2),time=c(1,1,2,1,2),within=c(0,0,1,1,0))
glm(x$dist ~ x$time + as.factor(x$within), family=poisson)

N <-1000
b <-0.001
I0 <-1
Di <-24

t <- 0
t_end <- 50
S <- N-I0
dt <- 0.1

while (t <= t_end){
	lambda <- b*(N-S)*S/N
	inf <- rbinom(1,S,(1-exp(-dt*lambda))) 
	rec <- rbinom(1,N-S,(1-exp(-dt/Di)))
	S <- S - inf + rec
	cat(t,t+dt,S,N-S,inf,rec,"\n")
	t <- t + dt
}

strLandscanFile <- "../../data/landscan/asia03/asia03/w001001.adf"
ExcerptLStoASCII(strLandscanFile,"c:/tmp/prd_sim_ascii.txt",112.789,114.486,22.150,23.518)

# datadir <- "D:\\files\\data_ns\\landscan\\africa03"
# infodir <- "D:\\files\\data_ns\\landscan\\africa03\\info\\"
# coveragedir <- "D:\\files\\data_ns\\landscan\\africa03\\afric03"

datadir <- "D:\\files\\data_ns\\landscan\\africa03"
infodir <- "D:\\files\\data_ns\\landscan\\africa03\\info\\"
coveragedir <- "D:\\files\\data_ns\\landscan\\africa03\\afric03"

covnames<-get.namesofcoverages(datadir)
tablenames<-get.tablenames(infodir)

arc <- get.arcdata(datadir,"africa03")

#Display the name of the table and its fields
for(i in 1:length(tablenames[[1]])) {
	print(c("Table: ",tablenames$TableName[i]))
	fields<-get.tablefields(infodir,tablenames$TableName[i])
	print("Fields")
	for(j in 1:length(fields))
	print(fields[[j]][1])
	#Get the data
	if(i==1)
	tabledata<-get.tabledata(infodir,tablenames$TableName[i])
	else
	tabledata<-c(tabledata, get.tabledata(infodir,tablenames$TableName[i]) )
}

get.arcdata(datadir2,"africa03")

dateConvert <- function(x) {
	y <- as.POSIXct(strptime(as.character(x),"%m/%d/%Y"))
	y
}

test <- fnSteadyState(ics,parvector,schistSimple,0.01,10,100)
test
schistSimple(0,test[[1]],parvector)[[2]]["y_H_L"]


quickFn <- function(v,a,b) {
	rtn <- a-(v[1]-b)^2
	rtn
}

optTest <- optim(c(-200),quickFn,a=2,b=3,control=list(trace=6,fnscale=-1))

out2 <- schistLike(ps_to_fit_vals,ps_to_fit_names,parvector,ics,infDat)

end_time <- 340
no_intervals <- 10
times = c((0:no_intervals)/no_intervals*end_time)
out <- lsoda(ics,times,schistSimple,parvector,rtol=1e-4)
out


fnFindStable <- function(y) {
	#undebug(schistSimple)
	ans <- schistSimple(0,y,parvector)[[1]]
	no_vars = length(ans)
	rtn = 0
	tot = 0
	for (i in 1:no_vars) {
		rtn = rtn + ans[[i]]*ans[[i]]
		if (y[i] < 0) rtn=1e10
		tot = tot + y[i]
	}
	if (tot<10) rtn = rtn + 1e10
	rtn
}

optim(ics,fnFindStable)

ans <- schistSimple(0,ics,parvector)
ans
ans <- schistSimple(0,ans[[1]]+ics,parvector)
ans

vecN <- c(52,13,13,13,13,3,3,5)
vecn <- c(39,9,7,4,3,0,0,0)

progs <- list.files()

getwd()

setwd("D:\\tmp")
x <- read.csv(file="malcolm_pca.csv",header=TRUE,row.names=NULL)
noResp <- dim(x)[1]
noQs <- dim(x)[2]-1
y <- prcomp(x[,2:(1+noQs)],retx=TRUE)
x <- cbind(x,PC1=array(NA,c(noResp)))
x <- cbind(x,PC2=array(NA,c(noResp)))
rotMat <- y[[2]]
for (i in 1:noResp) for (j in c("PC1","PC2")) x[i,j] <- sum(x[i,2:(1+noQs)]*rotMat[,j])
plot(	1:2,type="n",
		xlim=c(min(x[,"PC1"]),max(x[,"PC1"])),ylim=c(min(x[,"PC2"]),max(x[,"PC2"])),
		xlab="Principal component", ylab="Secondary component",main="Dummy PCA chart")
text(x[,"PC1"],x[,"PC2"],labels=x[,"Responant"])

# set the first line to the directory where your data file is
# you need to use "\\" where you would normaly type "\" for the directory location
# then copy and paste this script into the R console
# a chart should appear in a different window and the results of the PCA
# as the last few lines in the console.

setwd("D:\\tmp")
x <- read.csv(file="malcolm_pca.csv",header=TRUE,row.names=NULL)
noResp <- dim(x)[1]
noQs <- dim(x)[2]-1
y <- prcomp(x[,2:(1+noQs)],retx=TRUE)
x <- cbind(x,PC1=array(NA,c(noResp)))
x <- cbind(x,PC2=array(NA,c(noResp)))
rotMat <- y[[2]]
for (i in 1:noResp) for (j in c("PC1","PC2")) x[i,j] <- sum(x[i,2:(1+noQs)]*rotMat[,j])
plot(	1:2,type="n",
		xlim=c(min(x[,"PC1"]),max(x[,"PC1"])),ylim=c(min(x[,"PC2"]),max(x[,"PC2"])),
		xlab="Principal component", ylab="Secondary component",main="Dummy PCA chart")
text(x[,"PC1"],x[,"PC2"],labels=x[,"Responant"])
y

# Household multinomial likelihood for immune states
data <- list(c(70,63),c(52,66,71),c(15,25,39,29))
nohhsizes <- length(data)
lnlike <- 0
for (i in 1:nohhsizes) {
	chhs <- data[[i]]
	sizehh <- length(chhs-1)
	satps <- chhs / sum(chhs)
	lnlike <- lnlike + dmultinom(chhs,prob=satps,log=TRUE)
}

# Landscan inset / outset figure
rm(list=ls(all=TRUE))
require("rgdal")
require("RArcInfo")

strLandscanFile <- "D:\\files\\data_ns\\landscan\\asia03\\asia03\\w001001.adf"

strShapeDir <- "D:\\files\\data_ns\\worldsf"
strShapeLayer <- "world_adm0"

world <- readOGR(strShapeDir,layer=strShapeLayer)

# detail easting, detail westing..., wide eastin, wide westing...
de <- 114+7/60
dw <- 112+34/60
dn <- 23+35/60
ds <- 22+32/60

we <- 123 +  4/60
ww <- 108 +  7/60
wn <- 25  + 17/60
ws <- 18  +  3/60

spd <- 	120 	# Squares per degree
mn  <- 	85 		# map north
mw  <- 	34 		# map west

info <- GDALinfo(strLandscanFile)
wide <- readGDAL(strLandscanFile,region.dim=c((wn-ws)*spd,(we-ww)*spd),offset=c((mn-wn)*spd,(ww-mw)*spd))
detail <- readGDAL(strLandscanFile,region.dim=c((dn-ds)*spd,(de-dw)*spd),offset=c((mn-dn)*spd,(dw-mw)*spd))

f3w <- 10/cm(1)
widew <- 3/cm(1)	# Detail width
det_border <- 0.5/cm(1)
f3h <- (f3w-det_border)*(dn-ds)/(de-dw)+det_border
posdet <- c(0.01,0.99-det_border/f3w,0.01,0.99-det_border/f3h) 
poswide <- c(0.99-widew/f3w,0.99,0.99-widew*(wn-ws)/(we-ww)/f3h,0.99)

windows(width=f3w,height=f3h)
par(	mai=c(0,0,0,0), 		# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,
		mgp=c(0,0,0),
		par(cex=1))

par(fig=posdet)
image(detail,col=grey(0.9-0.8*0:100/100))
# image(detail,col="grey")
for (i in 1:4) axis(i,at=c(0,130),tck=0,col="blue",lwd=2)
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
polygon(c(0,1,1,0),c(0,0,1,1),col="white",border=NA)
par(new=TRUE,fig=poswide)
#image(wide,col="grey")
image(wide,col=grey(0.9-0.8*0:100/100))
for (i in 1:4) axis(i,at=c(0,130),tck=0)
par(new=TRUE,fig=poswide)
plot(1:2,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
lines(c(dw-ww,de-ww,de-ww,dw-ww,dw-ww)/(we-ww),c(ds-ws,ds-ws,dn-ws,dn-ws,ds-ws)/(wn-ws),col="blue",lwd=2)

simSEIR <- function(R0=1.8,D_E=1,D_I=1,N=1000,I0=1,
		tps=0:6,reals=1,dt=0.1,stochastic=TRUE) {
	
	noTPts <- length(tps) 
	epsilon <- 1e-10
	
	if (noTPts < 3) 
		stop("simSEIR: tps must be of length 2 or greater")
	if (reals > 1 && stochastic==FALSE) 
		stop("simSEIR: why make multiple realisations of a deterministic model?")
	
	rtn <- array(0,c(noTPts-1,3))
	
	for (i in 1:reals) {
		
		t_cur <- tps[1]
		ind_t <- 2
		t_next <- tps[ind_t]
		
		S <- N-I0
		E <- 0
		I <- I0
		R <- 0
		
		while (ind_t <= noTPts) {
			
			pInf <- 1-exp(-(R0*I)/(D_I*N)*dt)
			pBecInf <- 1-exp(-dt/D_E)
			pRec <- 1 - exp(-dt/D_I)
			
			if (stochastic) {
				nInf <- rbinom(1,S,pInf)
				nBecInf <- rbinom(1,E,pBecInf)
				nRec <- rbinom(1,I,pRec)
			} else {
				nInf <- S*pInf
				nBecInf <- E*pBecInf
				nRec <- I*pRec
			}
			
			S <- S - nInf
			E <- E + nInf - nBecInf
			I <- I + nBecInf - nRec
			R <- R + nRec
			
			rtn[ind_t-1,1] <- rtn[ind_t-1,1] + nInf 
			rtn[ind_t-1,2] <- rtn[ind_t-1,2] + nBecInf 
			rtn[ind_t-1,3] <- rtn[ind_t-1,3] + nRec 
			
			t_cur <- t_cur + dt
			
			if (t_cur > (t_next - epsilon)) {
				t_next <- tps[ind_t]
				ind_t <- ind_t+1
			}
			
		}
	}
	
	list(ave=(rtn/reals))
	
}

require("date")

f1 <- function(t) {
	
	t_d1 	<- 3.5
	t_d2 	<- 3.5
	t_r		<- as.date("03Jul2009")
	t_0		<- as.date("10Jun2009")
	r1		<- log(2)/t_d1
	r2		<- log(2)/t_d2
	
	if (t < t_r) rtn <- exp(r1*(t-t_0))
	else rtn <- exp(r1*(t_r-t_0))*exp(r2*(t-t_r))
	rtn

}

f2 <- function(t) {
	
	s1 		<- 1.0
	s2 		<- 0.05
	t_es1 	<- as.date("27Jun2009")
	t_ss2	<- as.date("10Jul2009")
	
	if (t < t_es1) s_t <- s1
	else if (t > t_ss2) s_t	<- s2
	else s_t <- s2 + (t_ss2 - t) / (t_ss2-t_es1) * (s1 - s2)
	
	rtn <- f1(t)*s_t
	rtn
	
}

drange <- c(as.date("01Jun2009"),as.date("01Aug2009"))
xvals <- drange[1]:drange[2]
yvals <- as.vector(rep(NA,length(xvals)))
for (i in 1:length(xvals)) yvals[i] <- f2(xvals[i])
plot(1:2,type="n",xlim = drange,ylim <- c(0,1000))
points(xvals,yvals)

# Example for 2D plots

x <- runif(10)
y <- runif(10)

hist <- sr2DHist(x,y)
nocols <- 11
cscale <- rainbow(nocols)
breakps <- 0:nocols/nocols
cscale[1] <- rgb(1,1,1)
hist$f <- hist$f / sum(hist$f)
image(x=hist$xm,y=hist$ym,z=hist$f,col=cscale,breaks=breakps)

# List of components for graphics

# Set up the defaul margins for a plot

# Some fragments for colin
likeInc <- function(pvec,freqtab,firstbin=18,lastbin=27) {
	starttime 	<- freqtab$lb[firstbin-1]
	rtn <- 0
	for (i in firstbin:lastbin) {
		rtn <- rtn + dpois(freqtab$f[i],pvec["A"]*exp(pvec["r"]*(freqtab$lb[i]-starttime)),log=TRUE)
	}
	rtn
}

pvec <- as.vector(mode="numeric",2^(1:10))
names(pvec) <- c("A","r")
fb <- 1
lb <- 10

mle1 <- optim(pvec,likeInc,control=list(trace=0,fnscale=-1,maxit=1000),freqtab=freq1,firstbin=fb,lastbin=lb)

# Some work with Kendra 9th July 2010
like <- function(X,pR0,pTg,ak=10,aT=30,aI0=10) {
	X[1,1] <- 0
	lnlike <- 0
	for (i in 1:(T+1)) {
		lnlike <- lnlike + dpois(sum(X[i,]),pR0*sum(X[,i+1]),log=TRUE)
	}
	lnlike
}


k <- 10
T <- 30
R0 <- 3
I0 <- 10
Tg <- 5
pvec <- dpois(1:k,Tg)
pvec <- pvec / sum(pvec)

X <- array(0,dim=c(T+1,T+k+1))
X[1,1] <- I0
for (i in 1:(T+1)) {
	X_i_dot <- rpois(1,sum(X[,i]*R0))
	X[i,(i+1):(i+k)] <- rmultinom(1,X_i_dot,pvec)
}

like(X,2,5)

xplot <- 1:50/10
yplot <- rep(-1,50)
for (i in 1:50) yplot[i] <- like(X,xplot[i],5)

plot(xplot,yplot,type="l")


# Next here: code the full model with variable mixing, infectivity, and susceptibility

# Putting in some simple things about the linear attack rate
tmp <-summarize_models(y,file="~/Dropbox/tmp/small.csv",transpose=TRUE)
