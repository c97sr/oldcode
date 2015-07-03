testMetric <- function(d,hi,max=1240) {abs(d-hi)}

smithMetric <- function(d,hi,max=1240) {

	# cat(d,hi,max,"\n")
	# Put the conditional smith function in here
	# Test the whole shebang on the three strain set
	# Then run for the larger example
	
	abs(d-hi)
}

extractVector <- function(coordTable) {
	# coordTable <- workingMapData
	rows = dim(coordTable)[1]
	dims =  dim(coordTable)[2] - 2
	rtn = array(rows*dims,c(0))	
	for (i in 1:rows) for (j in 1:dims) rtn[(i-1)*dims+j] <- coordTable[i,2+j]
	rtn
}

insertVector <- function(coordVec, coordTable) {
	dims = dim(coordTable)[2] - 2
	rows = dim(coordTable)[1]
	for (i in 1:rows) for (j in 1:dims) coordTable[i,2+j] <- coordVec[(i-1)*dims+j]
	coordTable
}

genDistFunc <- function(cv,st,ns,maxhi,dims=2,metfun=testMetric) {
	nsera <- ns
	noLocs <- length(cv)/dims
	nstra <- noLocs - nsera
	error <- 0
	for (i in 1:nsera) {
		for (j in 1:nstra) {
			sqrdist <- 0
			for (k in 1:dims) sqrdist <- sqrdist + (cv[(i-1)*dims+k] - cv[(j-1)*dims+nsera*dims+k])^2 
			eucdist <- sqrt(sqrdist)
			error <- error + metfun(eucdist,st[i,j+1],maxhi[j])	
		}
	}
	error
}

maskedDistFunc <- function(mcv,fcv,st,ns,maxhi,mask,mf=testMetric) {
	for (i in 1:length(mask)) fcv[mask[i]] <- mcv[i]
	rtn <- genDistFunc(fcv,st,ns,maxhi,metfun=mf)
	# cat(mcv,rtn,"\n")
	# flush.console()
	rtn	
}

pntsMap <- function(vec,tags,colour="black",tagson=TRUE) {
	# noanti <- length(vec)/2 - nosera
	colvec <- rainbow(length(vec)/2)
	if (tags==0) tags <- rep("",length(vec)/2)
	# for (i in 1:nosera) points(vec[(i-1)*2+1],vec[(i-1)*2+2],pch=22,col=colvec[i]) 
	# for (i in (nosera+1):(length(vec)/2)) points(vec[(i-1)*2+1],vec[(i-1)*2+2],pch=23,col=colvec[i]) 
	for (i in 1:(length(vec)/2)) {
		points(vec[(i-1)*2+1],vec[(i-1)*2+2],pch=19,col=colvec[i])
		if (tagson) text(vec[(i-1)*2+1],vec[(i-1)*2+2],tags[i],pos=4)
	}
	
}

cart2polar <- function(x,y) {
	rw <- sqrt(x*x + y*y)
	thetaw <- atan2(y,x)
	list(r=rw,theta=thetaw)
}

polar2cart <- function(r,theta) {
	xw <- r*cos(theta)
	yw <- r*sin(theta)
	list(x=xw,y=yw)
}

cartplustheta <- function(x,y,th) {
	z <- cart2polar(x,y)
	rtn <- polar2cart(z$r,z$theta + th)
	rtn
}

hiMap <- function(file,seed=1234,no_optims=10,max_iters=999999,l_bound=-15,u_bound=15,no_dims=2) {

	# file = "chinapoultry\\actual_data.txt"
	# seed=1234
	# no_optims=10
	# max_iters=999999
	# l_bound=-15
	# u_bound=15
	# no_dims=2
	
	set.seed(seed)
	
	hiData = read.table(file,header=TRUE,sep="\t")
	noRef = dim(hiData)[2]-1
	noSera = dim(hiData)[1]
	maxSera = rep(-999,noRef)
	for (i in 2:(noRef+1)) maxSera[i-1] <- max(hiData[,i])
	x <- rep(0,(noSera+noRef)*no_dims)
	msk <- 3:length(x)
	opvals <- array(-1,c(no_optims))
	opcoords <- array(-999,c(no_optims,length(x)+1))
	theta_target <- 0
	
	for (i in 1:no_optims) {
		y <- array(-999,c(length(x)-2))
		y[3:length(y)] <- runif(length(y)-2)*(u_bound-l_bound)+l_bound
		z <- optim(	y,maskedDistFunc,control=list(trace=0,fnscale=1,maxit=max_iters),
				method = "L-BFGS-B",lower=l_bound,upper=u_bound,
				fcv=x,st=hiData,ns=noSera,mask=msk,maxhi=maxSera,mf=smithMetric)
		cat(i,z$value,"\n")
		flush.console()
		opcoords[i,length(x)+1] <- z$value
		opcoords[i,3:length(x)] <- z$par
		theta_result <- cart2polar(opcoords[i,3],opcoords[i,4])$theta
		if (i==0) theta_target <- theta_result
		theta_delta <- theta_target - theta_result
		for (j in 2:(length(x)/2)) {
			tmp <- cartplustheta(opcoords[i,2*(j-1)+1],opcoords[i,2*(j-1)+2],theta_delta)
			opcoords[i,2*(j-1)+1] <- tmp$x
			opcoords[i,2*(j-1)+2] <- tmp$y		
		}
		opcoords[i,1:2] <- x[1:2]
	}
	opcoords
}

getNamesNoSera <- function(file) {
	x <- read.table(file,header=TRUE,sep="\t")
	noRef <- dim(x)[2]-1
	noSera <- dim(x)[1]
	namelist <- c(as.character(x$Sera),tail(names(x),noRef))
	list(no_ref = noRef, no_sera = noSera, names = namelist)
}

