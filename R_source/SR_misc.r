SR_from_unit <- function(x,min=1,max=100,logbase=10,logflag=FALSE) {

	# Takes a values and put it onto t a unit scale. Assumes a log scale if required. 

	if (logflag) {
		rtn <- min*logbase^(x*(log(max,logbase)-log(min,logbase)))
	} else {
		rtn <- min + (max-min)*x
	}
	rtn
}

SR_to_unit <- function(y,min=1,max=100,logbase=10,logflag=FALSE) {
	if (logflag) {
		rtn <- (log(y,logbase)-log(min,logbase))/(log(max,logbase)-log(min,logbase))
	} else {
		rtn <- (y-min)/(max-min) 
	}
	rtn
}

SR_BinOptim <- function(p,P,N,n) {
	guess <- pbinom(n,N,p,lower.tail=TRUE)
	rtn <- P-guess
	rtn
}

binCI <- function(vecN,vecn,P,min=0) {
	noObs <- length(vecN)
	rtn <- array(NA,c(noObs))
	for (i in 1:noObs) {
		if (vecN[i] > 0) {
			sol <- uniroot(SR_BinOptim,c(0,1),P=P,N=vecN[i],n=vecn[i])
			if (sol$root > min) rtn[i] <- sol$root
			else rtn[i] <- min
		} else rtn[i] <- min
	}
	rtn
}

binCINew <- function(N,n,P,min=0.000001) {
	rtn <- NA
	if (N>0) {
		if (N==n) n <- n-1
		if (n==0 && P > 0.5) rtn <- min
		else rtn <- (uniroot(SR_BinOptim,c(0,1),P=P,N=N,n=n))$root
	}
	rtn
}

binCINew(5,0,0.025,min=0.01)
pbinom(100,100,0.9,lower.tail=TRUE)

SR_updatePVector <- function(values,names,pvec) {                                                                        
    no_updates <- length(values)
    for (i in 1:no_updates) pvec[names[i]]=values[i]
    pvec <- jeUpdateAuxParams(pvec)
    pvec
}

sr_chart_pos <- function(	xindex,yindex,xn,yn,xlm=0,xrm=0,
							xg=0.01,ybm=0,ytm=0,yg=0.01) {

    x_left_mar = xlm
    x_right_mar = xrm
    x_gap = xg
    x_n_charts = xn
    x_width = (1-x_left_mar-x_right_mar-x_gap*(x_n_charts-1))/x_n_charts
    x_charts_l=rep(0,x_n_charts)
    for (i in 1:x_n_charts) x_charts_l[i]=x_left_mar+(i-1)*x_width+(i-1)*x_gap
    x_charts_r=rep(0,x_n_charts)
    for (i in 1:x_n_charts) x_charts_r[i]=x_charts_l[i]+x_width

    y_top_mar = ytm
    y_bottom_mar = ybm
    y_gap = yg
    y_n_charts = yn
    y_height = (1-y_top_mar-y_bottom_mar-y_gap*(y_n_charts-1))/y_n_charts
    y_charts_b = rep(0,y_n_charts)
    for (i in 1:y_n_charts) y_charts_b[i] = y_bottom_mar+(i-1)*y_height+(i-1)*y_gap
    y_charts_t = rep(0,y_n_charts)
    for (i in 1:y_n_charts) y_charts_t[i] = y_charts_b[i]+y_height

    rtn <- c(x_charts_l[xindex],x_charts_r[xindex],y_charts_b[yindex],y_charts_t[yindex])
    
    rtn

}

cm2i <- function(d) {d/2.54}

sr1DHist <- function(vec,min=0,max=1,bins=10,login=FALSE) {
	
	nosamp <- length(vec)
	
	tmpvec <- rep(0,bins)
	rtn <- data.frame(lb=tmpvec,ub=tmpvec,f=tmpvec)
	rtn$lb <- SR_from_unit(0:(bins-1)/bins,min=min,max=max,logflag=login)
	rtn$ub <- SR_from_unit(1:bins/bins,min=min,max=max,logflag=login)
	tmp <- SR_to_unit(vec,min=min,max=max,logflag=login)
	
	for (i in 1:nosamp) {
		bin_id <- floor(tmp[i]*bins)+1
		if (bin_id > 0 && bin_id <= bins) rtn$f[bin_id] <- rtn$f[bin_id] + 1 
	}
	
	rtn
	
}

srPlot1DHist <- function(srh,logarg="",add=FALSE,mainin="",xlabin="",ylabin="") {
	nobins <- dim(srh)[1]
	xbounds <- c(srh$lb[1],srh$ub[nobins])
	ybounds <- c(0,max(srh$f))
	if (!add) plot(1:2,type="n",xlim=xbounds,ylim=ybounds,log=logarg,main=mainin,xlab=xlabin,ylab=ylabin)
	for (i in 1:nobins) {
		polygon(	c(srh$lb[i],srh$lb[i],srh$ub[i],srh$ub[i]),
					c(0,srh$f[i],srh$f[i],0))
	}
}

sr2DHist <- function(	vec_x,vec_y,min_x=0,min_y=0,max_x=1,max_y=1,
						bins_x=10,bins_y=10,loginx=FALSE,loginy=FALSE) {
	
	# Function to take two vectors of parameters and create a 2D histogram
	# Returns a matrix of values and two vectors for the axis
	# The vectors are of length one greater than the dimensions of the axis
	
	nosamp <- length(vec_x)
	if (length(vec_y) != nosamp) stop("sr2DHist: Unequal length vectors")
		
	rtn <- matrix(0,nrow=bins_x,ncol=bins_y)
	
	tmp_x <- SR_to_unit(vec_x,min=min_x,max=max_x,logflag=loginx)
	tmp_y <- SR_to_unit(vec_y,min=min_y,max=max_y,logflag=loginy)
	
	for (i in 1:nosamp) {
		bin_idx <- floor(tmp_x[i]*bins_x)+1
		bin_idy <- floor(tmp_y[i]*bins_y)+1
		if (bin_idx > 0 && bin_idx <= bins_x &&
			bin_idy > 0 && bin_idy <= bins_y ) rtn[bin_idx,bin_idy] <- rtn[bin_idx,bin_idy] + 1 
	}
	
	xbounds <- SR_from_unit(0:(bins_x)/bins_x,min=min_x,max=max_x,logflag=loginx)
	ybounds <- SR_from_unit(0:(bins_y)/bins_y,min=min_y,max=max_y,logflag=loginy)
	
	xmids <- SR_from_unit((1:(bins_x)-0.5)/bins_x,min=min_x,max=max_x,logflag=loginx)
	ymids <- SR_from_unit((1:(bins_y)-0.5)/bins_y,min=min_y,max=max_y,logflag=loginy)

	list(f=rtn,xb=xbounds,yb=ybounds,xm=xmids,ym=ymids)
	
}

ExcerptLStoASCII <- function(lsfile,asciifile,de,dw,ds,dn,lsspd=120,lsmn=85,lsmw=34) {
	# Takes landscan file (lsfile) and produces and ascii grid file (asciifile) 
	# Returns nothing
	# coord agurments are degree east (exe), degree west (exw) ...
	browser()
	detail <- readGDAL(	lsfile,
			region.dim=c((dn-ds)*lsspd,(dw-de)*lsspd),
			offset=c((lsmn-dn)*lsspd,(de-lsmw)*lsspd))
	# writeAsciiGrid(detail,asciifile)
	write.asciigrid(detail,asciifile)
}

eventImage <- function(data,popgrid,ev,st,et,sr,er) {
	noEvents = dim(data)[1]
	rtnval = popgrid
	rtnval$z[]=NA
	xmin = popgrid$x[1]
	xgap = popgrid$x[2]-xmin
	nx = length(popgrid$x)-1
	xmax = xmin+nx*xgap
	ymin = popgrid$y[1]
	ygap = popgrid$y[2]-ymin
	ny = length(popgrid$y)-1
	ymax = ymin + ny*ygap
	for (i in 1:noEvents) {
		if (data$Event[i]==ev && data$Day[i] >= st && data$Day[i] <= et && data$Run[i] >= sr && data$Run[i] <= er) {
			xReal = data$X[i]
			yReal = data$Y[i]
			if (xReal >= xmin && xReal < xmax && yReal >= ymin && yReal < ymax) {
				xcoord = (xReal-xmin)/xgap+1
				# ycoord = ny - (yReal-ymin)/ygap+1
				ycoord = (yReal-ymin)/ygap+1
				if (is.na(rtnval$z[xcoord,ycoord])) rtnval$z[xcoord,ycoord]=1
				else rtnval$z[xcoord,ycoord]=rtnval$z[xcoord,ycoord]+1
			}
		}
	}
	noRuns = er-sr+1
	rtnval$z=rtnval$z/noRuns
	rtnval	
}

decLongLatDist <- function(x1,y1,x2,y2,translate=FALSE) {
	if (translate) {
		x1 <- x1*pi/180
		x2 <- x2*pi/180
		y1 <- y1*pi/180
		y2 <- y2*pi/180
	}
	norecs <- length(x1)
	rtn <- rep(-1,norecs)
	for (i in 1:norecs) {
		if (any(sapply(c(x1[i],x2[i],y1[i],y2[i]),is.na))) {rtn[i] <- NA}
		else {
			acosarg <- sin(y1[i])*sin(y2[i])+cos(y1[i])*cos(y2[i])* cos(x2[i]-x1[i])
			if (acosarg > 1) rtn[i] <- 0
			else rtn[i] <- 6378.7*acos(acosarg)
		}
	}
	rtn
}

sr_load_mcmc <- function(	fn,  
							thinfactor=1) {
	thinfactor <- round(thinfactor)
	od <- read.table(fn,header=TRUE)
	if (thinfactor>1) {
		nodata <- dim(od)[1]
		newsize <- nodata %/% thinfactor
		od <- od[(1:newsize)*thinfactor,]
	}
	cat("sr_load_mcmc finished")
	od
}

fnReadFileAsParamVector <- function(fname) {

	# Takes a filename for a parameter file with first row = "V1"
	# Other rows are "param name in quotes"	<white space> param value <cr>
	# Returns a vector with the correct names

	tmp <- as.data.frame(t(read.table(fname,row.names=1,header=FALSE)))
	rtn <- as.vector(tmp,mode="numeric")
	names(rtn) <- names(tmp)
	
	rtn

}

genStochMeanMedianCIs <- function(tab,vecnames,conf=c(0.025,0.975)) {
	
	# Examples of input below as used for developing
	vecnames <- c("a","b")
	tab <- data.frame(	a=c(1,2,3,4,5,6,7,8,9,10),
						b=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.1),
						c=c(10,20,30,40,50,60,70,80,90,100))
				
	nonames <- length(vecnames)
	for (s in vecnames) if (is.na(match(s,names(tab)))) 
			stop("genStochMeanMedianCIs: vecnames containes a string that's not in tab")
	rtn <- data.frame(mean=rep(0,nonames),median=rep(0,nonames),lb=rep(0,nonames),ub=rep(0,nonames))
	row.names(rtn) <- vecnames
	
	for (s in vecnames) {
		rtn[s,"mean"] 	<- mean(tab[,s])
		quants			<- quantile(tab[,s],probs=c(conf[1],0.5,conf[2]))
		rtn[s,"lb"]		<- quants[1]
		rtn[s,"median"]	<- quants[2]
		rtn[s,"ub"]		<- quants[3]
	}
	
	rtn
	
}
