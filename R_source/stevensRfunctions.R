# Copyright Steven Riley (ste.riley@gmail.com)

# This file is part of the library idsource.

# idsource is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

# require("date")
# require("odesolve")
# require("seqinr")

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

LoadCleanCHPData <- function(filename) {
	
	# Designed to take and clean a specific snapshot of data used early in the analysis
	
	#   index  	conf_date 	age 	sex admission  onset    	class 		country infector_index	arrival
	#	1 		01/05/2009  25  	M 	30/04/2009 28/4/2009 	Imported  	Mexico               	30/04/2009
	
	data <- read.csv(filename)
	norecs <- dim(data)[1]
	data[,"conf_date"] <- as.date(as.character(data[,"conf_date"]),order="dmy")
	data[,"admission"] <- as.date(as.character(data[,"admission"]),order="dmy")
	data[,"onset"] <- as.date(as.character(data[,"onset"]),order="dmy")
	data[,"arrival"] <- as.date(as.character(data[,"arrival"]),order="dmy")
	
	data <- cbind(data,imponset=rep(FALSE,norecs))
	data <- cbind(data,impimp=rep(FALSE,norecs))
	
	# Extra things for Max and Eric
	data <- cbind(data,maxclass=rep(FALSE,norecs))
	data <- cbind(data,maxonset=rep(NA,norecs))
	data <- cbind(data,maxonreport=data$conf_date-as.date("1Jan2009")+1)
	
	for (i in 1:norecs) {
		
		# Clean onset and import class
		if (is.na(data$onset[i])) data$imponset[i] <- TRUE
		if (!is.na(data$onset[i])) data$maxonset[i] <- data$onset[i]-as.date("1Jan2009")+1
		
		# Define age classes
		tmp <- data$age[i]
		if (is.na(tmp) || tmp < 0) data$ageclass[i] <- "NK"
		else if (tmp < 12) data$ageclass[i] <- "01yc"
		else if (tmp < 19) data$ageclass[i] <- "02oc"
		else if (tmp < 54) data$ageclass[i] <- "03ya"
		else if (tmp < 110) data$ageclass[i] <- "04oa"
		else stop("Problem with age allocation")
		
		# define age class 2
		tmp <- data$age[i]
		if (is.na(tmp) || tmp < 0) data$ageclass[i] <- "NK"
		else if (tmp < 12) data$ageclass2[i] <- "y"
		else if (tmp < 19) data$ageclass2[i] <- "o"
		else if (tmp < 110) data$ageclass2[i] <- "a"
		else stop("Problem with age allocation")
		
		# define simp import classes
		tmp <- data$class[i]
		if (tmp %in% c("Imported")) data$simpimp[i] <- "Imported"
		else if (tmp %in% c("Linked","NotLinked")) data$simpimp[i] <- "NotImported"
		else if (tmp %in% c("")) {
			data$simpimp[i] <- "NK"
			data$impimp[i] <- TRUE
		}
		else stop(paste("Problem with simplified importation class",tmp))
		data$ageimpclass[i] <- paste(data$simpimp[i],"_",data$ageclass[i],sep="")
		
		# define simp import classes
		tmp <- data$class[i]
		if (tmp %in% c("Imported")) data$maxclass[i] <- 2
		else if (tmp %in% c("","Linked","NotLinked","Undetermined")) data$maxclass[i] <- 1
		else stop("Problem with max's importation class")
		
	}
	
	data
	
}

AddZeros1DTable <- function(tab) {
	days <- as.numeric(names(tab))
	start <- min(days)
	end <- max(days)
	rtn <- vector(mode="numeric",length=end-start+1)
	for (d in start:end) {
		lookup <- match(d,days)
		if (is.na(lookup)) rtn[d-start+1] <- 0
		else rtn[d-start+1] <- tab[lookup]
	}
	list(days=start:end,vals=rtn)
}

MatrixWithZeroRows <- function(tab) {
	
	# Assumes row names are integers and adds zero rows
	
	maxday <- max(as.numeric(row.names(tab)))
	minday <- min(as.numeric(row.names(tab)))
	nodays <- maxday - minday + 1
	
	rtn <- matrix(data=0,nrow=nodays,ncol=dim(tab)[2])
	
	for (day in minday:maxday) 
		if (as.character(day) %in% row.names(tab)) 
			rtn[day-minday+1,] <- tab[as.character(day),]
	
	rtn
	
}

SimModel <- function(matNGM,matSeed,vecPop,vecParams) {
	
	# Take a set of seeds of different types by onset day and simulates an outbreak
	# Define input and output. Next generation matrix defined such that m[i,j] is the 
	# expected number of type j generated by one of type i.
	
	# To be used for the white-style model
	
	nodays <- dim(matSeed)[1]
	notypes <- dim(matSeed)[2]
	
	# Loop to generate infectious contacts on a given future day
	rtn <- matSeed
	rtn[] <- NA
	
	for (i in 1:notypes) {
		for (j in 1:nodays) {
			
		}
	}
	
	rtn
	
}

FigDataOnsetSeedLocal <- function(	filename, 
		matrix,setscolstop=c(c(1),c(2,3)),
		setscolsbot=c(c(1),c(2,3)),
		colors=c("red","blue"),
		widthcm=10,
		heightcm=10) {
	
	# Make figure 2B
	# - Assumes that ventilation.r has been run
	# - who made the 
	
	# browser()
	
	cwidth 	<- widthcm
	cheight <- heightcm
	
	pdf(filename, width=cwidth/cm(1), height=cheight/cm(1))
	
	# Gaps
	g1 <- 2.5 
	g2 <- 2.5
	g3 <- 0.25
	g4 <- 0.25
	
	xdates <- c("27Apr2009","4May2009","11May2009","18May2009","25May2009",
			"1Jun2009","8Jun2009","15Jun2009","22Jun2009")
	
	axticks <- as.date(xdates)
	axstart <- min(axticks)
	axend	<- max(axticks)
	axnames	<- c("1","2","3","4","5","6","7","8")
	ytickv	<- c(0,10,20,30,40,50)
	ytickn	<- ytickv
	ymax 	<- max(ytickv)
	ytickvB	<- c(0,10,20,30,40,50)
	yticknB	<- ytickvB 
	ymaxB 	<- max(ytickvB)
	lineylab <- 3.75
	
	# browser()
	
	# Spaces for the plots
	A <- (cwidth - g1 - g3)  
	B <- (cheight - g2 - 2*g4) / 2 
	
	# Positions of the different parts
	posA <- c(g1/cwidth,(g1+A)/cwidth,(g2+B+g4)/cheight,(g2+2*B+g4)/cheight)
	posB <- c(g1/cwidth,(g1+A)/cwidth,g2/cheight,(g2+B)/cheight)
	
	# Parameter settings for all parts of the figure
	par(cex=0.8)
	par(mgp=c(2.5,1,0))
	# par(tcl=0.25)
	par(mai=(c(0,0,0,0)))
	par(cex.axis=0.9)
	
	# Part A 
	par(fig=posA)
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=c(axstart,axend),
			ylim=c(0,ymax))
	
	axis(2,las=1,at=ytickv,labels=ytickn)
	mtext("Incidence",side=2,line=lineylab)
	axis(1,at=axticks,labels=rep("",length(axticks)))
	text(axstart,max(ytickv),"A",font=2,adj=c(0,1))
	
	legend(as.date("15May2009"),max(ytickv)*3/4,
			c(		"> 54",
					"19 - 54",
					"< 19"),
			lty=c(1,1,1),
			lwd=c(8,8,8),
			col=c("blue","green","red"),
			cex=0.75, 
			bty="n",
	)
	
	# Part B 
	par(fig=posB,new=TRUE)
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=c(axstart,axend),
			ylim=c(0,ymaxB))
	
	axis(2,las=1,at=ytickvB,labels=yticknB)
	mtext("Prevalence",side=2,line=lineylab)
	mtext("(per 100,000)",side=2,line=lineylab-1)
	axis(1,at=axticks,labels=rep("",length(axticks)),las=3)
	
	text(axstart,max(ytickvB),"B",font=2,adj=c(0,1))
	
	dev.off()
	
}

LikePoissonOnsetReport <- function(ps,vecOnsets,vecReports,vecImputes) {
	rtn <- 0
	mu <- ps[1]
	novals <- length(vecOnsets)
	for (i in 1:novals) {
		if (!vecImputes[i]) rtn <- rtn + dpois(vecReports[i]-vecOnsets[i],mu,log=TRUE)
		if (rtn < -1e100 || rtn > 1e100) browser()
	}
	rtn
}

ImputeOnsets <- function(mu,vecImputes,vecReports,vecOnsets) {
	noPoints <- length(vecImputes)
	rtn <- vecOnsets
	for (i in 1:noPoints) if (vecImputes[i]) rtn[i] <- vecReports[i]-rpois(1,mu)
	rtn
}

ImputeImport <- function(p_imp,vecImputes,vecSimpimp) {
	noPoints <- length(vecImputes)
	rtn <- vecSimpimp
	for (i in 1:noPoints) if (vecImputes[i]) {
			if (runif(1) < p_imp) rtn[i] <- "Imported"
			else rtn[i] <- "NotImported"
		}
	rtn
}

HKLikeInc <- function(pvec,inc) {
	nodays <- length(inc)
	rtn <- sum(dpois(inc,pvec["A"]*exp(pvec["r"]*((1:nodays)-1)),log=TRUE))
	rtn
}

lambda_U <- function(t,p,v) {
	if (t < p["tc_U"]) R0 <- p["R0_U1"]
	else R0 <- p["R0_U2"]
	rtn <- R0/p["D_I"]*(v["I_MU"]+v["I_UU"])/(p["N_U"]-p["P_U"]*p["N_T"]*p["D_S"])
	rtn
}

lambda_M <- function(t,p,v) {
	if (t < p["tc_M"]) R0 <- p["R0_M1"]
	else R0 <- p["R0_M2"]
	rtn <- R0/p["D_I"]*(v["I_MM"]+v["I_UM"])/(p["N_M"]-(1-p["P_U"])*p["N_T"]*p["D_S"])
	rtn
}

likeInc <- function(pvec,freqtab,firstbin=18,lastbin=27) {
	starttime 	<- freqtab$lb[firstbin-1]
	rtn <- 0
	for (i in firstbin:lastbin) {
		rtn <- rtn + dpois(freqtab$f[i],pvec["A"]*exp(pvec["r"]*(freqtab$lb[i]-starttime)),log=TRUE)
	}
	rtn
}

postProcLineList <- function(td,order="dmy",ageclass=10) {
	
	setinclude 	<- c(	"Arizona Case","CONFIRMED Colorado Case","Confirmed","PROBABLE",
			"Probable","confirmed","probable","CONFRIMED","P")
	
	td 	<- td[td$CaseStatus %in% setinclude,]
	
	norecords <- dim(td)[1]
	
	rtn <- as.data.frame(cbind(
					cleanOnset=rep(0,norecords),
					cleanInitialReportDate=rep(0,norecords),
					cleanHospPres=rep(0,norecords),
					cleanHosp=rep(0,norecords),
					cleanICU=rep(0,norecords),
					cleanOnsetWk=rep(0,norecords),
					cleanContact=rep(0,norecords),
					cleanAge=rep(0,norecords),
					cleanAG=rep(0,norecords),
					cleanAG2=rep(0,norecords),
					cleanAgeClass=rep(0,norecords),
					cleanMexTrav=rep(0,norecords)
			))
	
	for (i in 1:norecords) {
		
		# Clean up the onset date
		tmp <- as.date(as.character(td$Onset[i]),order=order)
		if (is.na(tmp)) rtn$cleanOnset[i] <- 0
		else rtn$cleanOnset[i] <- tmp
		
		# Some obvious manual corrections should be wihtout "as.date" in the if condition
		if (rtn$cleanOnset[i]==as.date("5Jan2009")) rtn$cleanOnset[i] 	<- as.date("1May2009")
		if (rtn$cleanOnset[i]==as.date("5Feb2009")) rtn$cleanOnset[i] 	<- as.date("2May2009")
		if (rtn$cleanOnset[i]==as.date("5Mar2009")) rtn$cleanOnset[i] 	<- as.date("3May2009")
		if (rtn$cleanOnset[i]==as.date("24Feb2009")) rtn$cleanOnset[i] 	<- 0
		if (rtn$cleanOnset[i]==as.date("27Apr3009")) rtn$cleanOnset[i] 	<- as.date("27Apr2009")
		if (rtn$cleanOnset[i]==as.date("25Apr2029")) rtn$cleanOnset[i] 	<- as.date("25Apr2009")
		if (rtn$cleanOnset[i]==as.date("9Apr2025")) rtn$cleanOnset[i] 	<- as.date("27Apr2009")
		if (rtn$cleanOnset[i]==as.date("9Apr2027")) rtn$cleanOnset[i] 	<- as.date("25Apr2009")
		if (rtn$cleanOnset[i]==as.date("5Apr2009")) rtn$cleanOnset[i] 	<- as.date("4May2009")
		if (rtn$cleanOnset[i]==as.date("5Jan1909")) rtn$cleanOnset[i] 	<- as.date("1May2009")
		if (rtn$cleanOnset[i]==as.date("22Apr1909")) rtn$cleanOnset[i] 	<- as.date("22Apr2009")
		if (rtn$cleanOnset[i]==as.date("24Apr1909")) rtn$cleanOnset[i] 	<- as.date("24Apr2009")
		if (rtn$cleanOnset[i]==as.date("1May1909")) rtn$cleanOnset[i] 	<- as.date("1May2009")
		if (rtn$cleanOnset[i]==-18506) rtn$cleanOnset[i] 	<- as.date("2May2009")
		if (rtn$cleanOnset[i]==15826) rtn$cleanOnset[i] 	<- as.date("2May2009")
		if (rtn$cleanOnset[i]==18053) rtn$cleanOnset[i] <- as.date("6May2009")
		if (rtn$cleanOnset[i]==18083) rtn$cleanOnset[i] <- as.date("7May2009")
		if (rtn$cleanOnset[i]==18114) rtn$cleanOnset[i] <- as.date("8May2009")
		if (rtn$cleanOnset[i]==18145) rtn$cleanOnset[i] <- as.date("9May2009")
		if (rtn$cleanOnset[i]==18174) rtn$cleanOnset[i] <- 0
		if (rtn$cleanOnset[i]==18175) rtn$cleanOnset[i] <- as.date("10May2009")
		
		# Enter the weekstart value
		if (rtn$cleanOnset[i] != 0) {
			rtn$cleanOnsetWk[i] <- rtn$cleanOnset[i] - (rtn$cleanOnset[i] - as.date("5Jan2009")) %% 7 
		} else {
			rtn$cleanOnsetWk[i] <- 0
		}
		
		# Clean up the date of report
		tmp <- as.date(as.character(td$InitialReportDate[i]))
		if (is.na(tmp)) rtn$cleanInitialReportDate[i] <- 0
		else rtn$cleanInitialReportDate[i] <- tmp
		if (rtn$cleanInitialReportDate[i]==as.date("5Jun2009")) rtn$cleanInitialReportDate[i] 	<- as.date("6May2009")		
		
		# Allocate hospitalized status
		tmp <- trimSpace(as.character(td$Hospitalized[i]))
		if (tmp %in% 
				c(	"Y","Yes","y","Y 5/3 (discharged 5/4)","Y (discharged)",
						"Y (disch)","Y 5/2 (discharged 4/3)","Y ","YES",
						"Yes, Pregnant, MILD illness","Y-readmitted 4/30")) 
			rtn$cleanHosp[i] <- "Y"
		else if (tmp %in% c("N","N(ER)","N (ER)","No","n","NO","NN")) rtn$cleanHosp[i] <- "N"
		else if (tmp %in% c("","U")) rtn$cleanHosp[i] <- "NK" 
		else stop("Encountered a value of Hospitalized that is not recognized >>",tmp,"<< at ",i)
		
		# Allocate ICU status
		tmp <- trimSpace(as.character(td$ICU[i]))
		if (tmp %in% c("Y","Yes")) rtn$cleanICU[i] <- "Y"
		else if (tmp %in% c("N","No","n") || rtn$cleanHosp[i] == "N") rtn$cleanICU[i] <- "N"
		else if (	is.na(tmp) || 
				tmp %in% c("","U","NA","Unknown")) rtn$cleanICU[i] <- "NK" 
		else {
			browser()
			stop("Encountered a value of ICU that is not recognized >>",tmp,"<< at ",i)
		}	
		# Clean contact status
		tmp <- trimSpace(as.character(td$ContactofConfCase[i]))
		if (tmp %in% c("Y","Y ","Yes","Yes, Probable","Y (child of 09V01540)")) rtn$cleanContact[i] <- "Y"
		else if (tmp %in% c("N","No","n")) rtn$cleanContact[i] <- "N"
		else if (tmp %in% c(""," ","U","U ","u","Unknown")) rtn$cleanContact[i] <- "NK" 
		else stop("Encountered a value of clean contact that is not recognized >>",tmp,"<< at ",i)
		
		# Clean Age
		tmp_a <- as.character(td$age[i])
		tmp_d <- as.character(td$dob[i])
		if (!is.na(tmp_a) && tmp_a!="" && !as.numeric(tmp_a) < 0) { 
			rtn$cleanAge[i] <- as.numeric(tmp_a)
		} else {
			tmpdob <- as.date("1Jan2200")
			if (!is.na(as.date(tmp_d,order="dmy"))) tmpdob <- as.date((tmp_d),order="dmy")
			else if (!is.na(as.date(tmp_d,order="mdy"))) tmpdob <- as.date((tmp_d),order="mdy")
			tmponset <- as.date("24Apr2009")
			if (rtn$cleanOnset[i]!= 0) tmponset <- rtn$cleanOnset[i]
			if (tmpdob < tmponset) rtn$cleanAge[i] <- (tmponset - tmpdob) / 365.25
			else rtn$cleanAge[i] <- -1
		}
		
		# Allocate age status
		if (rtn$cleanAge[i] <0) rtn$cleanAG[i] <- "NK"
		else if (rtn$cleanAge[i] < 19) rtn$cleanAG[i] <- "c"
		else if (rtn$cleanAge[i] < 55) rtn$cleanAG[i] <- "a"
		else rtn$cleanAG[i] <- "e"
		rtn$cleanAgeClass[i] <- rtn$cleanAge[i] %/% ageclass
		if (rtn$cleanAge[i] <0) rtn$cleanAG2[i] <- "NK"
		else if (rtn$cleanAge[i] < 5) rtn$cleanAG2[i] <- "00t04"
		else if (rtn$cleanAge[i] < 18) rtn$cleanAG2[i] <- "05t17"
		else if (rtn$cleanAge[i] < 64) rtn$cleanAG2[i] <- "18t64"
		else rtn$cleanAG2[i] <- "64p"
		
		# Clean contact status
		tmp <- trimSpace(as.character(td$TravelMexico[i]))
		if (tmp %in% c("Y","yes","Y (3/27 - 4/3) - Cancun","Y (3/27 - 4/5) - Cancun")) rtn$cleanMexTrav[i] <- "Y"
		else if (tmp %in% c(	"N","n","N (Arizona)","N (Costa Rica)","N (travel to DC)",
				"N (San Francisco)","N (travel to DC)","No","M","N (DC)",
				"N (San Diego 4/13-23)","no")) rtn$cleanMexTrav[i] <- "N"
		else if (tmp %in% c("","U","u","Unknown")) 
			rtn$cleanMexTrav[i] <- "NK" 
		else stop("Encountered a value of cleanMexTrav that is not recognized >>",tmp,"<< at ",i)
		
	}
	
	rtn
	
}

venModel <- function(t,y,p) {
	
	rtn=array(0,c(10))
	names(rtn) 	<- c("S","Ia","Is","Ih","Iv","dS","dIa","dIs","dIh","dIv")
	names(y) 	<- c("S","Ia","Is","Ih","Iv","dS","dIa","dIs","dIh","dIv")
	lambda 		<- p["beta"]*(y["Ia"]+y["Is"])/p["N"]
	rtn["S"] 	<- - lambda*y["S"]
	rtn["Ia"] 	<- lambda*y["S"]*(1-p["p_R"]) - p["gamma"]*y["Ia"]
	rtn["Is"] 	<- lambda*y["S"]*p["p_R"] - p["gamma"]*y["Is"]
	rtn["Ih"] 	<- p["gamma"]*y["Is"]*p["p_H_base"]*(1-p["p_I"]) - p["gamma_h"]*y["Ih"]
	rtn["Iv"] 	<- p["gamma"]*y["Is"]*p["p_H_base"]*p["p_I"] - p["gamma_v"]*y["Iv"]
	
	rtn["dS"] 	<- lambda*y["S"]									# infection incidence
	rtn["dIa"] 	<- lambda*y["S"]*(1-p["p_R"])						# asymptomatic infection incidence
	rtn["dIs"] 	<- lambda*y["S"]*p["p_R"]							# symptomatic infection incidence
	rtn["dIh"] 	<- p["gamma"]*y["Is"]*p["p_H_base"]*(1-p["p_I"])	# hospitalized infection incidence
	rtn["dIv"] 	<- p["gamma"]*y["Is"]*p["p_H_base"]*p["p_I"]		# ventilated case incidence
	rtn["dD"] 	<- p["gamma_v"]*y["Iv"]*p["p_D"]					# incidence death
	
	list(rtn)
	
}

usParams <- function(){
	
	c(		N 			= 1,
			seed 		= 1,
			r 			= 0.26,
			Tg 			= 1.8,
			phi_1		= 2.0,
			phi_2		= 1.0,
			p_R 		= 0.86,
			p_H_base 		= 0.031,
			p_H_base_1 		= 0.02489,
			p_H_base_2 		= 0.03575,
			p_H_base_3 		= 0.09338,
			p_I 		= 0.13,
			p_I_1 		= 0.1364,
			p_I_2 		= 0.1111,
			p_I_3 		= 0.2222,
			gamma_v 	= 1/5.2,									# Gowardman
#			p_I_1 		= 0.1364/2,
#			p_I_2 		= 0.1111/2,
#			p_I_3 		= 0.2222/2,
#			gamma_v 	= 1/9.0,									# Gowardman
			gamma_h 	= 1/7.0,									# Feagan		
			t0 			= as.date("12Apr2009"),
			p_1			= (82640086-21514358/5)/304059724,			# From http://www.census.gov/popest/national/asrh/NC-EST2008-sa.html
			p_2			= (170378099 - 4/5*21514358)/304059724,
			p_3			= 72555897/304059724,
			mixmatindex	= 1
	)
}

modelSIR <- function(pname=NULL,pvals=NULL,venParams=usParams()) {
	
	if (length(pname) != length(pvals)) stop("Error: modelSIR: the two vectors must be the same length")
	
	if (length(pname) > 0) {
		for (i in 1:length(pname)) {
			if (pname[i] %in% names(venParams) == FALSE) stop("Unknown parameter passed to a model SIR")
			else venParams[pname[i]] <- pvals[i] 
		}
	}
	
	venParams["R0"] = 1 + venParams["r"]*venParams["Tg"]
	venParams["gamma"] = 1 / venParams["Tg"]
	
	venParams["S0"] = venParams["N"] - venParams["seed"]
	venParams["Is0"] = venParams["seed"]
	
	venParams["beta"] = venParams["gamma"]*venParams["R0"]	
	initial_conditions <- c(venParams["S0"],0,venParams["Is0"],0,0,0,0,0,0,0,0)
	names(initial_conditions) <- c("S","Ia","Is","Ih","Iv","dS","dIa","dIs","dIh","dIv","dD")
	
	# Debugging from here
	
	no_invals 	<- 364
	no_years 	<- 1
	time_points <- ((0:no_invals)*364*no_years)/no_invals
	timestep 	<- time_points[2]-time_points[1]
	last_index 	<- length(time_points+1)
	
	solution <- lsoda(initial_conditions,time_points,venModel,venParams,atol=1e-80)
	
	solution[1,"dS"] <- 0
	for (i in 2:last_index) solution[i,"dS"] <- (solution[i,"dIs"] - solution[i-1,"dIs"]) / timestep
	
	solution[,"time"] <- solution[,"time"] + venParams["t0"] 
	
	solution
	
}

MixingCAE <- function() {
	
	# Mixing matrix based on the UK data from mossang
	# Age groups here are 0-18, 19-54, 55+
	
	rtn <- t(array(c(	0.58,	0.36,	0.06,
							0.23,	0.63,	0.14,
							0.15,	0.52,	0.33	
					), dim=c(3,3)))
	
	rtn
	
}

MixingCCA <- function() {
	
	# Mixing matrix based on the UK data from mossang
	# Age groups here are 0-11, 12-18, 18+
	
	rtn <- t(array(c(		0.43,	0.10,	0.46,
							0.08,	0.53,	0.39,
							0.09,	0.11,	0.80	
					), dim=c(3,3)))
	
	rtn
	
}

# Scale mixing value
ScaleMixingValueOldOld <- function(val,scale) {val[1,1] + scale * (val - val[1,1])}

ScaleMixingValue <- function(val,scale) {
	val^scale
}

ScaleMixingValueOld <- function(val,scale) {
	nrow <- (dim(val))[1]
	ncol <- (dim(val))[2]
	rtn <- matrix(0,nrow=nrow,ncol=ncol)
	for (i in 1:nrow) {
		for (j in 1:ncol) {
			if (i==j && i==1) rtn[i,j] <- val[i,j]
			else if (val[i,j] > val[1,1]) rtn[i,j] <- val[1,1] + scale * (val[i,j] - val[1,1])
			else rtn[i,j] <- val[1,1] - scale * (val[1,1] - val[i,j])
		}
	}
	if (min(rtn) < 1e-10) stop("Scale vaulues giving negative mixing is not allowed")
	rtn
}


MixingCA3 <- function(scale=1) {
	
	# Mixing matrix based on the UK data from mossang
	# Age groups here are 0-20, 20+
	# Still using three age groups just for convenience
	
	# 7.8	2.15
	# 5.63	8.210909091
	# Up to here
	
	rtn <- t(array(c(		7.8,	2.15,	2.15,
							3.57,	4.10,	4.10,
							3.57,	4.10,	4.10	
					), dim=c(3,3)))
	
	rtn <- ScaleMixingValue(rtn,scale)
	
	rtn
	
}

MixingCCAA <- function() {
	
	# Mixing matrix based on the UK data from mossang
	# Age groups here are as per lipsitch severity, 0-4, 5-17, 18-64, 65+
	
	# These data are worked out from mossong {Mossong, 2008 p01379} in the file
	# age_mixing.xls
	
	# Note that the m=t(array(c(...),dim=(.,.))) gives m[i,j] in the same way that R would print m
	# without the transpose t(.), the i and j are switched
	
	rtn <- t(array(c(
							0.216,	0.198,	0.560,	0.026,
							0.029,	0.585,	0.366,	0.020,
							0.040,	0.173,	0.725,	0.062,
							0.023,	0.151,	0.616,	0.209	
					),c(4,4)))
	
	rtn
	
}

MixingNull <- function() {
	
	# Null mixing matrix
	
	rtn <- t(array(c(	1,		1,		1,
							1,		1,		1,
							1,		1,		1	
					), dim=c(3,3)))
	rtn
	
}

MixingNullFour <- function() {
	
	# Null mixing matrix
	
	rtn <- t(array(c(	1,	1,	1,	1,
							1,	1,	1,	1,
							1,	1,	1,	1
					), dim=c(4,4)))
	
	rtn
	
}

ageVenModel <- function(t,y,p) {
	
	# Set the number of state variables and a "skip" auxilliary variable
	noAges 	<- 3
	noStates <- 6
	noAux 	<- 3 
	skp <- noAux + noStates
	
	# Set up the main return functions
	rtn <- array(0,c(noAges*(noAux+noStates)))
	
	# Generate the normalized mixing matrix
	if (p["mixmatindex"]==1) mixm <- MixingNull()
	else if (p["mixmatindex"]==2) mixm <- MixingCAE()
	else if (p["mixmatindex"]==3) mixm <- MixingCCA()
	else if (p["mixmatindex"]==4) mixm <- MixingCA3(p["mixsense"])
	else stop("Code for mixing matrix not found.")
	
	# Scale beta for sesonality
	if ((t %% 364) > p["dur_seas"]) beta <- p["beta"]*(1-p["amp_seas"])
	else beta <- p["beta"]
	
	# Generate the age specific vector for number of people, probability of hospitalization and probability of being severe
	vecN 	<- p["N"]*array(c(p["p_1"],p["p_2"],p["p_3"])/(p["p_1"]+p["p_2"]+p["p_3"]))
	vecpH	<- array(c(p["p_H_base_1"],p["p_H_base_2"],p["p_H_base_3"]))
	vecpI	<- array(c(p["p_I_1"],p["p_I_2"],p["p_I_3"]))
	vecPhi 	<- array(c(p["phi_1"],p["phi_2"],p["phi_3"]))
	
	# browser()
	
	# If aging turned on, generate the aging in and aging out terms for each state
	ageing <- matrix(0,nrow=noAges,ncol=noStates)
	if (p["aging_on"]==1) {
		ageing[1,1] 		<- c(	p["mu_3"]*sum(y[((3-1)*skp+1):((3-1)*skp+6)]) - p["mu_1"]*y[(1-1)*skp+1])
		ageing[1,2:noAges] 	<- - p["mu_1"]*y[(1-1)*skp+2:noAges]
		ageing[2,1:noAges]	<- + p["mu_1"]*y[(1-1)*skp+1:noAges] - p["mu_2"]*y[(2-1)*skp+1:noAges]
		ageing[3,1:noAges]	<- + p["mu_2"]*y[(2-1)*skp+1:noAges] - p["mu_3"]*y[(3-1)*skp+1:noAges]
	}
	
	# Start the main loop
	for (i in 1:noAges) {
		
		# Set the force of infection to zero and start a loop through j infecting classes
		lambda <- 0
		for (j in 1:noAges) {
			
			# Each j class mixes a proportion minm[i,j] of its time with the i class
			# (j-1)*skp+2 and (j-1)*skp+3] pick out the infectious classes
			lambda <- lambda + 	mixm[i,j]*(y[(j-1)*skp+2]+y[(j-1)*skp+3]) / vecN[j] 
			
		}
		
		# The force of infection is adjusted according to the susceptibility of the ith class
		# Different ngm matrices require different parameterizations
		lambda <- lambda *  beta * vecPhi[i]
		
		# Eqn for susceptibles
		rtn[(i-1)*skp+1]	<- -(lambda*y[(i-1)*skp+1]+p["trickle"]/3) + p["gamma_R"]*y[(i-1)*skp+6] + ageing[i,1]
		
		# Eqn for asymp infections
		rtn[(i-1)*skp+2]	<- (lambda*y[(i-1)*skp+1]+p["trickle"]/3)*(1-p["p_R"]) - p["gamma"]*y[(i-1)*skp+2] + ageing[i,2]
		
		# Eqn for symp infections
		rtn[(i-1)*skp+3]	<- (lambda*y[(i-1)*skp+1]+p["trickle"]/3)*p["p_R"] - p["gamma"]*y[(i-1)*skp+3] + ageing[i,3]
		
		# Eqn for hosp normal ward
		rtn[(i-1)*skp+4]	<- p["gamma"]*y[(i-1)*skp+3]*vecpH[i]*(1-vecpI[i]) - p["gamma_h"]*y[(i-1)*skp+4] + ageing[i,4]
		
		# Eqn for hosp ICU ward
		rtn[(i-1)*skp+5]	<- p["gamma"]*y[(i-1)*skp+3]*vecpH[i]*vecpI[i] - p["gamma_v"]*y[(i-1)*skp+5] + ageing[i,5]
		
		# Eqn for recovered
		rtn[(i-1)*skp+6]	<- p["gamma_v"]*y[(i-1)*skp+5] + p["gamma_h"]*y[(i-1)*skp+4] + 
				p["gamma"]*y[(i-1)*skp+3]*(1-vecpH[i]) + p["gamma"]*y[(i-1)*skp+2] - 
				p["gamma_R"]*y[(i-1)*skp+6] + ageing[i,6]
		
		# infection incidence
		rtn[(i-1)*skp+7] 	<- lambda*y[(i-1)*skp+1]							
		
		# asymptomatic infection incidence
		rtn[(i-1)*skp+8] 	<- lambda*y[(i-1)*skp+1]*(1-p["p_R"])
		
		# symptomatic infection incidence
		rtn[(i-1)*skp+9] 	<- lambda*y[(i-1)*skp+1]*p["p_R"]
		
	}
	
	# if (t>500) browser()
	
	list(rtn)
	
}

ageVenModelFour <- function(t,y,p) {
	
	# Set the number of state variables and a "skip" auxilliary variable
	noAges 	<- 4
	noStates <- 5
	noAux 	<- 4 
	skp <- noAux + noStates
	
	# Set up the main return functions
	rtn <- array(0,c(noAges*(noAux+noStates)))
	
	# Generate the normalized mixing matrix
	if (p["mixmatindex"]==1) mixm <- MixingNullFour()
	else if (p["mixmatindex"]==2) mixm <- MixingCCAA()
	else stop("Code for mixing matrix not found for ageVenModelFour.")
	
	# Implement school holidays reduction in transmission
	if ((t + p["t0"]) > p["t_st"]) mixm[2,2] <- mixm[2,2] * p["delta_sm"] 
	
	# Generate the age specific vector for number of people, probability of hospitalization and probability of being severe
	vecN 	<- p["N"]*array(c(p["p_1"],p["p_2"],p["p_3"],p["p_4"])/(p["p_1"]+p["p_2"]+p["p_3"]+p["p_4"]))
	vecpH	<- array(p["p_H_base"]*c(p["p_H_base_1"],p["p_H_base_2"],p["p_H_base_3"],p["p_H_base_4"]))
	vecpI	<- array(c(p["p_I_1"],p["p_I_2"],p["p_I_3"],p["p_I_4"]))
	vecPhi	<- array(c(p["phi_1"],p["phi_2"],p["phi_3"],p["phi_4"]))
	
	# Start the main loop
	for (i in 1:noAges) {
		
		# Set the force of infection to zero and start a loop through j infecting classes
		lambda <- 0
		for (j in 1:noAges) {
			
			# Each j class mixes a proportion minm[i,j] of its time with the i class
			# (j-1)*skp+2 and (j-1)*skp+3] pick out the infectious classes
			lambda <- lambda + 	mixm[i,j]*(y[(j-1)*skp+2]+y[(j-1)*skp+3]) / vecN[j] 
			
		}
		
		# The force of infection is adjusted according to the susceptibility of the ith class
		lambda <- lambda * p["beta"] * vecPhi[i]	
		
		# Eqn for susceptibles
		rtn[(i-1)*skp+1]	<- - lambda*y[(i-1)*skp+1]
		
		# Eqn for asymp infections
		rtn[(i-1)*skp+2]	<- lambda*y[(i-1)*skp+1]*(1-p["p_R"]) - p["gamma"]*y[(i-1)*skp+2]
		
		# Eqn for symp infections
		rtn[(i-1)*skp+3]	<- lambda*y[(i-1)*skp+1]*p["p_R"] - p["gamma"]*y[(i-1)*skp+3]
		
		# Eqn for hosp normal ward
		rtn[(i-1)*skp+4]	<- p["gamma"]*y[(i-1)*skp+3]*vecpH[i]*(1-vecpI[i]) - p["gamma_h"]*y[(i-1)*skp+4]
		
		# Eqn for hosp ICU ward
		rtn[(i-1)*skp+5]	<- p["gamma"]*y[(i-1)*skp+3]*vecpH[i]*vecpI[i] - p["gamma_v"]*y[(i-1)*skp+5]
		
		# infection incidence
		rtn[(i-1)*skp+6] 	<- lambda*y[(i-1)*skp+1]							
		
		# asymptomatic infection incidence
		rtn[(i-1)*skp+7] 	<- lambda*y[(i-1)*skp+1]*(1-p["p_R"])
		
		# symptomatic infection incidence
		rtn[(i-1)*skp+8] 	<- lambda*y[(i-1)*skp+1]*p["p_R"]
		
		# ICU incidence
		rtn[(i-1)*skp+9] 	<- p["gamma"]*y[(i-1)*skp+3]*vecpH[i]*vecpI[i]	
		
	}
	
	list(rtn)
	
}

MakeNgm <- function(mix,params,vecN) {
	
	# Returns a next generation matrix with variable susceptibility by age group 
	# Designed for multiple levels of mixing
	
	rtn <- 1/params["gamma"]*t(array(c(	
							params["phi_1"]*mix[1,1],					params["phi_1"]*(vecN[1]*mix[1,2]/vecN[2]),	params["phi_1"]*(vecN[1]*mix[1,3]/vecN[3]),
							vecN[2]*mix[2,1]/vecN[1],					mix[2,2],									vecN[2]*mix[2,3]/vecN[3],
							params["phi_2"]*(vecN[3]*mix[3,1]/vecN[1]),	params["phi_2"]*(vecN[3]*mix[3,2]/vecN[2]),	params["phi_2"]*(mix[3,3])
					), dim=c(3,3)))
	rtn
}

MakeNgmFour <- function(mix,params,vecN) {
	
	# Returns a next generation matrix with variable susceptibility by age group 
	# Designed for multiple levels of mixing
	
	rtn <- 1/params["gamma"]*t(array(c(	
							params["phi_1"]*mix[1,1],					params["phi_1"]*vecN[1]*mix[1,2]/vecN[2],	params["phi_1"]*vecN[1]*mix[1,3]/vecN[3],				params["phi_1"]*vecN[1]*mix[1,4]/vecN[4],
							params["phi_2"]*vecN[2]*mix[2,1]/vecN[1],	params["phi_2"]*mix[2,2],					params["phi_2"]*vecN[2]*mix[2,3]/vecN[3],			params["phi_2"]*vecN[2]*mix[2,4]/vecN[4],
							params["phi_3"]*vecN[3]*mix[3,1]/vecN[1],	params["phi_3"]*vecN[3]*mix[3,2]/vecN[2],	params["phi_3"]*mix[3,3],														params["phi_3"]*vecN[4]*mix[3,4]/vecN[4],
							params["phi_4"]*vecN[4]*mix[4,1]/vecN[1],	params["phi_4"]*vecN[4]*mix[4,2]/vecN[2],	params["phi_4"]*vecN[4]*mix[4,3]/vecN[3],			params["phi_4"]*mix[4,4]
					), dim=c(4,4)))
	
	rtn
}

likeSusEigen <- function(vecS,params,vecN,mix,cm,ngm) {
	
	# Likelihood function for different susceptibilities given mixing matrix and observed infections
	
	if (vecS[1] < 0 || vecS[2] < 0) rtn <- -1e100
	else {
		params["gamma"] <- 1
		params["phi_1"] <- vecS[1]
		params["phi_2"] <- vecS[2]
		params["phi_3"] <- 1
		mG <- ngm(mix,params,vecN)
		ev  <- Re(eigen(mG)$vector[,1])
		if (ev[1] < 0) ev <- ev * -1
		rtn <- dmultinom(cm,prob=ev,log=TRUE)
	}
	rtn
}

likeSusEigenFour <- function(vecS,params,vecN,mix,cm,ngm) {
	
	# Likelihood function for different susceptibilities given mixing matrix and observed infections
	
	if (vecS[1] < 0 || vecS[2] < 0 || vecS[3] < 0) rtn <- -1e100
	
	else {
		
		params["gamma"] <- 1
		params["phi_1"] <- vecS[1]
		params["phi_2"] <- vecS[2]
		params["phi_3"] <- vecS[3]
		params["phi_4"] <- 1
		mG 	<- ngm(mix,params,vecN)
		ev  <- Re(eigen(mG)$vector[,1])
		if (ev[1] < 0) ev <- ev * -1
		rtn <- dmultinom(cm,prob=ev,log=TRUE)
		
	}
	
	rtn
	
}

SirModel3AgeClasses <- function(pname=NULL,pvals=NULL,casemix=c(1,1,1),vp=usParams(),NGM=MakeNgm) {
	
	# Calling function for the age structured ode model
	# Returns table of all the useful output variables
	# pname is a list of parameter names to change
	# pvals is the list of values to assigne to them
	# case mix is the proportion of individuals in each age class
	# vp is the parameter set to be used
	# NGM is the function to be used to make the next generation matrix
	
	# Check that the length of pname is the same as that of pvec
	if (length(pname) != length(pvals)) stop("Error: modelSIR: the two vectors must be the same length")
	
	# Add some extra parameters which are fully constrained by the other parameters
	# Allow only r or R0 to be present, not both
	if (!is.na(vp["R0"]) && !is.na(vp["r"])) stop("Both R0 and r cannot be specified in this model")
	if (is.na(vp["R0"]) && !is.na(vp["r"])) vp["R0"] = 1 + vp["r"]*vp["Tg"]
	
	# For non-null parameter vectors, change the appropriate parameters
	if (length(pname) > 0) {
		for (i in 1:length(pname)) {
			if (pname[i] %in% names(vp) == FALSE) stop("Unknown parameter passed to a model SIR")
			else vp[pname[i]] <- pvals[i] 
		}
	}
	
	# rate of recovery
	vp["gamma"] = 1 / vp["Tg"]
	
	# Assign the different mixing matrices based on the value of the appropriate parameter
	if (vp["mixmatindex"]==1) mm <- MixingNull()
	else if (vp["mixmatindex"]==2) mm <- MixingCAE()
	else if (vp["mixmatindex"]==3) mm <- MixingCCA()
	else if (vp["mixmatindex"]==4) mm <- MixingCA3(vp["mixsense"])
	else stop("Code for mixing matrix not found.")
	
	# Set up the population vector for the NGM calcuations
	vecP <- c(vp["p_1"],vp["p_2"],vp["p_3"]) / (vp["p_1"]+vp["p_2"]+vp["p_3"]) 
	vN <- c(vp["N"]*vecP[1],vp["N"]*vecP[2],vp["N"]*vecP[3])
	
	# Obtain the values for relative susceptibility, fitted using the NGM data
	# DEBUG_LINE likeSusEigen(c(2,1,1),vp,vN,mm,casemix,ngm=NGM)
	
	if (vp["fitsus"]==1) {
		bfSus <- optim(	c(2,1),likeSusEigen,control=list(trace=0,fnscale=-1,maxit=10000),
				params=vp, vecN=vN, mix=mm, cm=casemix, ngm=NGM)			
		vp["phi_1"] <- bfSus$par[1] # 2.2
		vp["phi_2"] <- bfSus$par[2] # 0.34
	}
	
	# Need to assign phi_3 even if not fitting susceptibility
	vp["phi_3"] <- 1
	
	# Browser here to check susceptibility parameters
	# DEBUG_LINE browser()
	
	# Defined the next generation matrix with, implicitly, beta set = 1
	mG <- NGM(mm,vp,vN)
	esol <- eigen(mG)
	
	if (vp["fitsus"]==1) {	
		# Double check that the eigen vector matches that implied used to contrain susceptibility
		evec <- esol$vectors[,1]
		bench <- casemix
		evec <- evec/sum(evec)*sum(bench)
		if (abs(sum(evec-bench)) > 1) stop("Problem matching the eigen vector to initial age-specific incidence") 
	}
	
	# Calculate the eigen values, the real eigen values and then the largest real eigen value
	ev  <- esol$value
	revs <- ev[Im(ev)==0]
	maxrev <- max(Re(revs))
	
	# Define beta using the basic reproductive number and the largest eigen value
	vp["beta"] <- vp["R0"] / maxrev
	
	# Set the initial conditions assuming that the seed is split evenly between age classed
	initial_conditions <- c(
			vN[1]-vp["seed"]/3,0,vp["seed"]/3,0,0,0,0,0,0,
			vN[2]-vp["seed"]/3,0,vp["seed"]/3,0,0,0,0,0,0,
			vN[3]-vp["seed"]/3,0,vp["seed"]/3,0,0,0,0,0,0	)
	
	# Define the names for the output variables
	vecnames <- c()
	for (i in 1:3) vecnames <- c(vecnames,
				paste(	c("S","Ia","Is","Ih","Iv","R","dS","dIa","dIs"),i,sep="")	)
	names(initial_conditions) <- vecnames
	
	# Set run parameters for the model
	no_invals 	<- vp["nodts"]
	no_years 	<- vp["noyears"]
	time_points <- ((0:no_invals)*364*no_years)/no_invals
	timestep 	<- time_points[2]-time_points[1]
	last_index 	<- length(time_points+1)
	
	# Solve the differential equations
	# browser()
	solution <- lsoda(initial_conditions,time_points,ageVenModel,vp,atol=1e-80)
	
	# Add some auxilliary variables and calculate their values
	solution <- cbind(solution, 
			dS 				= rep(NA,length(time_points)),
			dI1 			= rep(NA,length(time_points)),
			dI2 			= rep(NA,length(time_points)),
			dI3 			= rep(NA,length(time_points)),
			Iv 				= rep(NA,length(time_points)),
			Ih 				= rep(NA,length(time_points)),
			N1 				= rep(NA,length(time_points)),
			N2 				= rep(NA,length(time_points)),
			N3 				= rep(NA,length(time_points)),
			N 				= rep(NA,length(time_points)))
	
	solution[1,"dI1"] <- 0
	solution[1,"dI2"] <- 0
	solution[1,"dI3"] <- 0
	for (i in 2:last_index) {
		solution[i,"dI1"] <- (solution[i,"dIs1"] - solution[i-1,"dIs1"]) / timestep 
		solution[i,"dI2"] <- (solution[i,"dIs2"] - solution[i-1,"dIs2"]) / timestep 
		solution[i,"dI3"] <- (solution[i,"dIs3"] - solution[i-1,"dIs3"]) / timestep 
	} 	
	solution[,"dS"] <- solution[,"dI1"] + solution[,"dI2"] + solution[,"dI3"]
	for (i in 1:last_index) {
		solution[i,"Iv"] <- solution[i,"Iv1"] + solution[i,"Iv2"] + solution[i,"Iv3"]
		solution[i,"Ih"] <- solution[i,"Ih1"] + solution[i,"Ih2"] + solution[i,"Ih3"]
	}
	
	solution[,"N1"] <- rowSums(solution[,paste(c("S","Ia","Is","Ih","Iv","R"),c(1),sep="")])
	solution[,"N2"] <- rowSums(solution[,paste(c("S","Ia","Is","Ih","Iv","R"),c(2),sep="")])
	solution[,"N3"] <- rowSums(solution[,paste(c("S","Ia","Is","Ih","Iv","R"),c(3),sep="")])
	solution[,"N"] <- rowSums(solution[,paste(c("N"),1:3,sep="")])
	
	# Debug lines below for tuning demographic parameters
	# browser()
	# plot(solution[,"N"],type="l",ylim=c(0,max(solution[,"N"])))
	# points(solution[,"N1"],type="l",col="red")
	# points(solution[,"N2"],type="l",col="blue")
	# points(solution[,"N3"],type="l",col="green")	
	
	solution[,"time"] <- solution[,"time"] + vp["t0"] 
	
	# Return the solution and the value of the parameters
	list(sol=solution,par=vp)
	
}

SirModelFourAgeClasses <- function(pname=NULL,pvals=NULL,casemix=c(1,1,1,1),vp=usParams(),NGM=MakeNgmFour) {
	
	# Calling function for the age structured ode model, four age class model
	# This probably needs to be called with number of age classes as an argument?
	# Returns table of all the useful output variables
	# pname is a list of parameter names to change
	# pvals is the list of values to assigne to them
	# case mix is the proportion of individuals in each age class
	# vp is the parameter set to be used
	# NGM is the function to be used to make the next generation matrix
	
	# Check that the length of pname is the same as that of pvec
	if (length(pname) != length(pvals)) stop("Error: modelSIR: the two vectors must be the same length")
	
	# For non-null parameter vectors, change the appropriate parameters
	if (length(pname) > 0) {
		for (i in 1:length(pname)) {
			if (pname[i] %in% names(vp) == FALSE) stop("Unknown parameter passed to a model SIR")
			else vp[pname[i]] <- pvals[i] 
		}
	}
	
	# basic reproductive number
	vp["R0"] = 1 + vp["r"]*vp["Tg"]
	
	# rate of recovery
	vp["gamma"] = 1 / vp["Tg"]
	
	# Assign the different mixing matrices based on the value of the appropriate parameter
	if (vp["mixmatindex"]==1) mm <- MixingNullFour()
	else if (vp["mixmatindex"]==2) mm <- MixingCCAA()
	else stop("Code for mixing matrix not found in four age class model.")
	
	# Set up the population vector for the NGM calcuations
	vecP <- c(vp["p_1"],vp["p_2"],vp["p_3"],vp["p_4"]) / (vp["p_1"]+vp["p_2"]+vp["p_3"]+vp["p_4"]) 
	vN <- c(vp["N"]*vecP[1],vp["N"]*vecP[2],vp["N"]*vecP[3],vp["N"]*vecP[4])
	
	if (vp["fitsus"]!=0) {
		
		# Obtain the values for relative susceptibility, fitted using the NGM data
		# likeSusEigen(c(2,1,1),vp,vN,mm,casemix,ngm=NGM)
		bfSus <- optim(	c(2,1,1),likeSusEigenFour,control=list(trace=0,fnscale=-1,maxit=10000),
				params=vp, vecN=vN, mix=mm, cm=casemix, ngm=NGM)			
		vp["phi_1"] <- bfSus$par[1] # 2.2
		vp["phi_2"] <- bfSus$par[2] # 0.34
		vp["phi_3"] <- bfSus$par[3]
		vp["phi_4"] <- 1
	}
	
	# Browser here to check susceptibility parameters
	# DEBUG_LINE browser()
	
	# Defined the next generation matrix with, implicitly, beta set = 1
	mG <- NGM(mm,vp,vN)
	esol <- eigen(mG)
	
	# Double check that the eigen vector matches that implied used to contrain susceptibility
	evec <- esol$vectors[,1]
	bench <- casemix
	evec <- evec/sum(evec)*sum(bench)
	if (vp["fitsus"]!=0 && sum(abs(evec-bench)) > 0.01) stop(sum(abs(evec-bench)),"0.001","Problem matching the eigen vector to initial age-specific incidence") 
	
	# Calculate the eigen values, the real eigen values and then the largest real eigen value
	ev  <- esol$value
	revs <- ev[Im(ev)==0]
	maxrev <- max(Re(revs))
	
	# Define beta using the basic reproductive number and the largest eigen value
	vp["beta"] <- vp["R0"] / maxrev
	
	# Set the initial conditions assuming that the seed is split evenly between age classed
	initial_conditions <- c(
			vN[1]-vp["seed"]/4,0,vp["seed"]/4,0,0,0,0,0,0,
			vN[2]-vp["seed"]/4,0,vp["seed"]/4,0,0,0,0,0,0,
			vN[3]-vp["seed"]/4,0,vp["seed"]/4,0,0,0,0,0,0,
			vN[4]-vp["seed"]/4,0,vp["seed"]/4,0,0,0,0,0,0	)
	
	# Define the names for the output variables
	vecnames <- c()
	for (i in 1:4) vecnames <- c(vecnames,
				paste(	c("S","Ia","Is","Ih","Iv","dS","dIa","dIs","dIv"),i,sep="")	)
	names(initial_conditions) <- vecnames
	
	# Set run parameters for the model
	no_invals 	<- 364
	no_years 	<- 1
	# Adjust here
	# browser() # XXXX editing here
	# rounded_start_day <- round(vp[])
	time_points	<- 0:(round((vp["tf"]-vp["t0"]))+1) 
	# solution[1,"time"] <- solution[1,"time"] + vp["t0"] %% 1 
	
	# time_points <- ((0:no_invals)*364*no_years)/no_invals
	timestep 	<- time_points[2]-time_points[1]
	last_index 	<- length(time_points+1)
	
	# Solve the differential equations
	solution <- lsoda(initial_conditions,time_points,ageVenModelFour,vp,atol=1e-80)
	
	# Add some auxilliary variables and calculate their values
	solution <- cbind(solution, 
			dS 				= rep(NA,length(time_points)),
			dI1 			= rep(NA,length(time_points)),
			dI2 			= rep(NA,length(time_points)),
			dI3 			= rep(NA,length(time_points)),
			dI4 			= rep(NA,length(time_points)),
			Iv 				= rep(NA,length(time_points)),
			Ih 				= rep(NA,length(time_points)),
			dIv 			= rep(NA,length(time_points))
	)
	
	solution[1,"dI1"] <- 0
	solution[1,"dI2"] <- 0
	solution[1,"dI3"] <- 0
	solution[1,"dI4"] <- 0
	
	for (i in 1:last_index) {
		solution[i,"Iv"] <- solution[i,"Iv1"] + solution[i,"Iv2"] + solution[i,"Iv3"] + solution[i,"Iv4"]
		solution[i,"Ih"] <- solution[i,"Ih1"] + solution[i,"Ih2"] + solution[i,"Ih3"] + solution[i,"Ih4"]
		solution[i,"dIv"] <- solution[i,"dIv1"] + solution[i,"dIv2"] + solution[i,"dIv3"] + solution[i,"dIv4"]
	}
	
	for (i in last_index:2) {
		solution[i,"dI1"] <- (solution[i,"dIs1"] - solution[i-1,"dIs1"]) / timestep 
		solution[i,"dI2"] <- (solution[i,"dIs2"] - solution[i-1,"dIs2"]) / timestep 
		solution[i,"dI3"] <- (solution[i,"dIs3"] - solution[i-1,"dIs3"]) / timestep 
		solution[i,"dI4"] <- (solution[i,"dIs4"] - solution[i-1,"dIs4"]) / timestep 
		solution[i,"dIv"] <- (solution[i,"dIv"] - solution[i-1,"dIv"]) / timestep 
	}
	
	solution[,"dS"] <- solution[,"dI1"] + solution[,"dI2"] + solution[,"dI3"] + solution[,"dI4"]
	solution[,"time"] <- solution[,"time"] + round(vp["t0"]) 
	
	# Return the solution and the value of the parameters
	list(sol=solution,par=vp)
	
}

srDailyInc <- function(vecInc,sd=0,ed=10) {
	nobins 	<- (ed - sd) + 1
	inc 	<- vector(mode="numeric",nobins)
	inc[] 	<- 0
	for (i in 1:length(vecInc)) {
		tmp <- vecInc[i]
		if (tmp >= sd && tmp <= ed) {
			inc[tmp - sd + 1] <- inc[tmp - sd + 1] + 1 
		}	
	}
	inc
}

CIsForPlot <- function(vecN, vecn, ratelb=0.00001) {
	
	noObs <- length(vecn)
	if (noObs != length(vecN)) stop("CIsForPlot: Vectors need to be the same size.")
	rtn <- data.frame(pt=rep(NA,noObs),ub=rep(NA,noObs),lb=rep(NA,noObs))
	
	for (i in 1:noObs) {
		
		tmp_n <- vecn[i]
		tmp_N <- vecN[i]
		if (tmp_N >= 1) {
			rtn$pt[i] <- tmp_n / tmp_N
			if (tmp_n == tmp_N) rtn$ub[i] <- 1
			else rtn$ub[i] <- binCINew(tmp_N,tmp_n,0.025,min=0.000001)
			rtn$lb[i] <- binCINew(tmp_N,tmp_n,0.975,min=lbrates)
		}	
		
	}
	
	rtn	
	
}

plotstack <- function(sol,ac,col,var="Iv",mult=1,thresh=1,const=0) {
	
	ts <- dim(sol)[1]
	x <- c(sol[,"time"],rev(sol[,"time"]))
	
	if (ac==1) y <- c((sol[,paste(var,"1",sep="")]),rep(0,ts))
	else if (ac==2) y <- c((sol[,paste(var,"1",sep="")] + sol[,paste(var,"2",sep="")]),rev((sol[,paste(var,"1",sep="")])))
	else if (ac==3) y <- c((sol[,paste(var,"1",sep="")] + sol[,paste(var,"2",sep="")] + sol[,paste(var,"3",sep="")]),rev((sol[,paste(var,"1",sep="")]+sol[,paste(var,"2",sep="")])))
	else stop("Age class not specified.")
	
	polygon(x,y*mult+const,col=col,border=NA)
	
}

agerate <- function(sol,mult=c(1,1,1),thresh=100) {
	ts 	<- dim(sol)[1] 
	rtn <- 		1/ sol[,"dS"] * (
				sol[,"dI1"]*mult[1] +
				sol[,"dI2"]*mult[2] +
				sol[,"dI3"]*mult[3] )
	
	for (i in 1:ts) if (sol[i,"dS"] < thresh) rtn[i] <- NA
	rtn 
	
	
}


poisLike <- function(vecModTimes, vecModPrevs, vecObsTimes,vecObsPrevs) {
	
	# Takes four vectors.
	# The first two are the same length and are the times and prevalences from the model.
	# The second two are the same size and are the times and prevalances from the observations.
	# The function assumes that the obs arose from poisson distributions with the model mean.
	# NAs in the vector of observed prevalences are ignored.
	# Assumes that every time point in obs is in model, but not vice versa
	# Assumes that times are monotonically increasing
	
	# An epsilpon value for testing equivalence of doubles
	eps <- 1e-10
	
	# Initialize the return value
	rtn <- 0
	
	# Measure the length of the two vectors and initialize the pointer to the current model time
	intNumberMod 	<- length(vecModTimes)
	intNumberObs 	<- length(vecObsTimes)
	intPointerMod 	<- 1
	
	# Start the main loop increment by observation time
	for (i in 1:intNumberObs) {
		dblCurrentObsTime <- vecObsTimes[i]
		while (vecModTimes[intPointerMod] < dblCurrentObsTime) intPointerMod <- intPointerMod + 1
		dblCurrentModTime <- vecModTimes[intPointerMod]
		if (abs(dblCurrentModTime - dblCurrentObsTime) > eps) {
			browser()
			stop("poisLike: Preconditions not satisfied.")
		}
		intCurrentObsPrev <- vecObsPrevs[i]
		dblCurrentModPrev <- vecModPrevs[intPointerMod]
		if (!is.na(intCurrentObsPrev)) {
			if (intCurrentObsPrev < eps && dblCurrentModPrev > eps) stop("poisLike: Zero model prevalence and non-zero observation.")
			rtn <- rtn + dpois(intCurrentObsPrev,dblCurrentModPrev,log=TRUE)
		}
	}
	
	# Return the log "likelihood" 
	rtn	
	
}

PdmSolveLikeV1 <- function(pvals,pname,pbase,casemix,ngm,likefunc,obs) {
	
	# Generate the model solution
	tmp 		<- SirModelFourAgeClasses(pname=pname,pvals=pvals,
			casemix=casemix,vp=pbase,NGM=ngm)		 
	solution 	<- tmp$sol
	endparams	<- tmp$par
	
	ModelTimes <- as.vector(solution[,"time"])
	ModelPrevs <- as.vector(solution[,"Ih"] + solution[,"Iv"])*4380000
	ObsTimes <- obs[,"Date"]
	ObsPrevs <- as.vector(obs[,"Qld_Hosp"])
	
	rtn <- poisLike(ModelTimes,ModelPrevs,ObsTimes,ObsPrevs)
	
	list(like=rtn,sol=solution,par=endparams)
	
}

PdmSolveLikeV2 <- function(pvals,pname,pbase,casemix,ngm,obs) {
	
	# Generate the model solution
	tmp 		<- SirModelFourAgeClasses(pname=pname,pvals=pvals,casemix=casemix,vp=pbase,NGM=ngm)		 
	solution 	<- tmp$sol
	endparams	<- tmp$par
	
	# Define the other variables needed
	incData 	<- obs$inc
	noObs 		<- length(incData)
	incModel 	<- vector(mode="numeric",length=noObs)
	minday 		<- solution[1,"time"] 
	
	for (i in 1:noObs) {
		
		# Extract incidence of ICU from the solution
		incModel[i] <- sum(solution[(obs$date_start[i]-minday+1):(obs$date_finish[i]-minday+1),"dIv"])
		
	}
	
	rtn <- sum(dpois(incData,incModel*5000000,log=TRUE))	
	
	list(like=rtn,sol=solution,par=endparams)
	
}

PdmLikeHKHosp <- function(pvals,pname,pbase,casemix,obs,vecbins,agnames) {
	
	# Generate the model solution
	tmp 		<- SirModelFourAgeClasses(pname=pname,pvals=pvals,casemix=casemix,vp=pbase)		 
	solution 	<- as.data.frame(tmp$sol)
	endparams	<- tmp$par
	
	incModel <- extFourAgeICU(solution,vecbins,agnames)
	
	rtn <- sum(dpois(obs,incModel,log=TRUE))	
	
	list(like=rtn,sol=solution,par=endparams)
	
}


fnProposeParamUpdates <- function(ptab,fmask=1:(dim(ptab)[1])) {
	
	# Takes biological parameters, a mask which is a list of the parameters
	# being fitted, a table of the fitted parameters and their ranges and
	# the number of parameters to be fitted. 
	# Returns a proposed vector of parameters.
	# Comment in and out bouncing and cyclical boundary conditions
	# for the random walk
	
	bps <- ptab[,"val"]
	rtn <- bps
	nofit <- length(rtn)
	
	for (i in 1:nofit) {
		
		# Set up and transform to unit scale
		rv <- runif(1)
		rv <- (rv-0.5)* ptab[i,"step"]
		x <- bps[fmask[i]]
		x <- SR_to_unit(x,min=ptab[i,"min"],max=ptab[i,"max"],logflag=ptab[i,"log"])
		x <- x + rv
		
		# Bouncing boundary conditons
		# if (x < 0) x <- -x	
		# if (x > 1) x <- 2 - x
		
		# Cyclical boundary conditions
		if (x < 0) x <- 1 + x	
		if (x > 1) x <- x - 1
		
		# Test for errors and return to originl scales
		if (x < 0 || x > 1) stop("problem here")		
		rtn[fmask[i]] <- SR_from_unit(x,min=ptab[i,"min"],max=ptab[i,"max"],logflag=ptab[i,"log"])
		
	}
	
	rtn
}

fnProposeParamUpdatesSingle <- function(ptab,index) {
	
	# Takes biological parameters, a mask which is a list of the parameters
	# being fitted, a table of the fitted parameters and their ranges and
	# the number of parameters to be fitted. 
	# Returns a proposed vector of parameters.
	# Comment in and out bouncing and cyclical boundary conditions
	# for the random walk
	
	bps <- ptab[,"val"]
	rtn <- bps
	if (index < 0 || index > dim(ptab)[1]) stop("index must be less than or equal to number of parameters")
	
	# Set up and transform to unit scale
	rv <- runif(1)
	rv <- (rv-0.5)* ptab[index,"step"]
	x <- bps[index]
	x <- SR_to_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])
	x <- x + rv
	
	# Cyclical boundary conditions
	if (x < 0) x <- 1 + x	
	if (x > 1) x <- x - 1
	
	# Test for errors and return to originl scales
	if (x < 0 || x > 1) stop("problem here")		
	rtn[index] <- SR_from_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])
	
	rtn
}

mcmcH1N1v1 <- function(	fnParams,
		fnOutput,
		params,
		funcLike,
		no_samples=10,
		samp_freq=1,
		report_block=1000,
		tune_block=100,
		tune_min=20,
		tune_max=80,
		tune_factor=2,
		propsingle=TRUE,
		dbgaccept=FALSE,
		...) {
	
	# Set up the mcmc chain
	x 				<- read.csv(fnParams,row.names=1)
	fittedparams 	<- rownames(x)
	no_params 		<- length(fittedparams)
	nosamples 		<- no_samples
	sampfreq		<- samp_freq
	nosamples 		<- nosamples - nosamples %% sampfreq
	norecorded		<- report_block
	nonparheadings	<- c("sample","lnlike")
	nnph			<- length(nonparheadings)
	chain 			<- array(0,dim=c(norecorded,2*length(fittedparams)+nnph))
	colnames(chain) <- c(nonparheadings,fittedparams,paste(fittedparams,"_step",sep=""))
	popvec			<- c(params["p_1"],params["p_2"],params["p_3"],params["p_4"])
	popvecn			<- popvec / sum(popvec)
	
	# Set up tuning characteristics of the main chain
	propcount <- vector(no_params,mode="numeric")
	propcount[]=0
	acceptcount <- vector(no_params,mode="numeric")
	acceptcount[]=0
	
	previoussol <- funcLike(pvals=x[,"val"],pname=fittedparams,pbase=params,...)
	
	lnlike <- previoussol$like
	
	recordno <- 1
	first_write	<- TRUE
	for (i in 0:(nosamples-1)) {
		
		if (i %% sampfreq == 0) {
			chain[recordno,"lnlike"] 	<- lnlike
			chain[recordno,"sample"] 	<- i+1
			chain[recordno,(nnph+1):((nnph+1)+no_params-1)] <- x[,"val"]
			chain[recordno,(nnph+no_params+1):((nnph+1)+2*no_params-1)] <- x[,"step"]
			recordno <- recordno + 1
			if (recordno > norecorded) {
				if (first_write) {
					write.table(chain,file=fnOutput,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
					first_write <- FALSE
				} else {
					write.table(chain,file=fnOutput,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
				}
				recordno <- 1
			}
		}
		
		chosenparam <- ceiling(runif(1)*no_params)
		prop_x <- fnProposeParamUpdatesSingle(x,chosenparam)
		# prop_x <- fnProposeParamUpdates(x)
		
		
		currentsol <- funcLike(pvals=prop_x,pname=fittedparams,pbase=params,...)
		
		prop_lnlike <- currentsol$like 	
		
		diff_like <- prop_lnlike - lnlike
		
		if (diff_like > 0) accept <- TRUE
		else if (exp(diff_like) > runif(1)) accept <- TRUE
		else accept <- FALSE
		
		propcount[chosenparam] <- propcount[chosenparam]+1
		
		if (accept || dbgaccept) {
			x[,"val"] <- prop_x
			lnlike <- prop_lnlike
			previoussol <- currentsol
			acceptcount[chosenparam] <- acceptcount[chosenparam]+1
		}	
		
		if (propcount[chosenparam]==tune_block) {
			if (acceptcount[chosenparam] < tune_min) x[chosenparam,"step"] <- x[chosenparam,"step"]/tune_factor
			else if (acceptcount[chosenparam] > tune_max) {
				x[chosenparam,"step"] <- x[chosenparam,"step"]*tune_factor
				if (x[chosenparam,"step"] > 1) x[chosenparam,"step"] <- 1
			}
			propcount[chosenparam] <- 0
			acceptcount[chosenparam] <- 0
		}
		
	}
	
	if (recordno > 1) write.csv(chain[1:(recordno-1),],file=fnOutput,row.names=FALSE)
	
}

genMeanMedianCIs <- function(tab,vecnames,conf=c(0.025,0.975)) {
	
	# Examples of input below as used for developing
	# vecnames <- c("a","b")
	# tab <- data.frame(	a=c(1,2,3,4,5,6,7,8,9,10),
	#		b=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.1),
	#		c=c(10,20,30,40,50,60,70,80,90,100))
	
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

CharToDate <- function(x,order="dmy") {as.date(as.character(x),order=order)}

LoadAndCleanHKSeroInd <- function(filename) {
	
	# reads in the csv file of the participant data from the serosurvey and does
	# any necessary post processing
	
	filename <- "/Volumes/NO NAME/data/influenza/hk_serosurvey/fu_attendence_sr_noident.csv"
	
	rtn <- read.csv(filename)
	
	no_rows <- dim(rtn)[1]
	
	# Fix easy dates
	rtn$base_int_date <- CharToDate(rtn$base_int_date,order="mdy")
	rtn$dob <- CharToDate(rtn$dob,order="mdy")
	
	# Fix recruitment dates for each houehold
	rtn <- cbind(rtn,hh_rec_d=rep(NA,no_rows))
	for (i in 1:no_rows) rtn$hh_rec_d[i] <- FixRecruitmentDate(i,rtn)
	
	# Add age group
	rtn$AG <- rep(1,no_rows)
	for (i in 1:no_rows) {
		tmp <- rtn$dob[i]
		if (!is.na(tmp)) {
			age <- round((as.date("06Jul2009") - tmp) / 365.25)
			if (age > 110 || age < 0) browser() #("problem with age group assignment")
			if (age < 5) rtn$AG[i] <- 2
			else if (age < 19) rtn$AG[i] <- 3
			else if (age < 65) rtn$AG[i] <- 4
			else rtn$AG[i] <- 5
		}
	}
	
	# Let the 
	
	rtn
	
}

LoadAndCleanHKSeroSymp <- function(filename) {
	
	# reads in the csv file of the symptom data from the serosurvey and does
	# any necessary post processing
	# filename <- "/Volumes/NO NAME/data/influenza/hk_serosurvey/fu_symptoms_sr.csv"
	
	rtn <- read.csv(filename)
	rtn$contact_d <- CharToDate(rtn$contact_d)
	rtn$onset_d <- CharToDate(rtn$onset_d)
	rtn
	
}

FixRecruitmentDate <- function(i,tab=c()) {
	
	# i <- 3
	# tab <- ind_data
	
	rtn <- 0
	if (!is.na(tab$base_int_date[i])) rtn <- tab$base_int_date[i]
	else {
		for (j in (i-tab$hh_index[i]):(i-tab$hh_index[i]+tab$hh_size[i]-1))
			if (!is.na(tab$base_int_date[j])) rtn <- tab$base_int_date[j]		
	}
	rtn
	
}

LookUpMain <- function(hh,ind,tab) {
	i <- 1
	max <- dim(tab)[1]
	while (tab$hh_id[i] != hh || tab$hh_index[i] != ind && i <= max) i <- i + 1
	if (i > max) stop("error in LookupMainList")
	i
}

unit_to_log <- function(v,lb,ub) {
	rtnval = lb*10^((log10(ub)-log10(lb))*v)
	rtnval
}

unit_to_linear <- function(v,lb,ub) {
	rtnval = lb + (ub-lb)*v
	rtnval
}

hyper_vector <- function(n,lb,ub,log=FALSE) {
	wk = array(runif(n,0,1),c(n))
	for (i in 1:n) {
		wk[i] = (i-1)/n+wk[i]/n
		if (log) wk[i] = unit_to_log(wk[i],lb,ub) 
		else wk[i] = unit_to_linear(wk[i],lb,ub)
	}
	ranrank = array(runif(n,0,1),c(n))
	index = order(ranrank)
	rtnval = array(0,c(n))
	for (i in 1:n) rtnval[i]=wk[index[i]]
	rtnval
}

procEFlu <- function(filename) {
	
	# Read in the file
	raw <- read.csv(filename)
	
	# Get list of indexes of ever PICU, NICU, ICU
	# index.only.severe <- (raw$hosp.ever.nic == 3) | (raw$hosp.ever.icu == 3)
	index.hosp <- (raw$hosp.ipas.inst.cd != 0) 
	hosp <- raw[index.hosp,] 	# only severe
	noRows <- dim(hosp)[1]
	
	# Make up a matrix of all symptoms and admission date
	# This would be nicer as a list of symptom field strings, but not worth effort and would be slower 
	matSymp <- matrix(nrow=noRows,ncol=12)
	matSymp[,1] <- as.date(as.character(hosp$sym.fever.onset),order="dmy")
	matSymp[,2] <- as.date(as.character(hosp$sym.pneumonia.onset),order="dmy")
	matSymp[,3] <- as.date(as.character(hosp$sym.sore.onset),order="dmy")
	matSymp[,4] <- as.date(as.character(hosp$sym.acute.respiratory.onset),order="dmy")
	matSymp[,5] <- as.date(as.character(hosp$sym.head.onset),order="dmy")
	matSymp[,6] <- as.date(as.character(hosp$sym.myalgia.onset),order="dmy")
	matSymp[,7] <- as.date(as.character(hosp$sym.chills.onset),order="dmy")
	matSymp[,8] <- as.date(as.character(hosp$sym.fatigue.onset),order="dmy")
	matSymp[,9] <- as.date(as.character(hosp$sym.vomit.onset),order="dmy")
	matSymp[,10] <- as.date(as.character(hosp$sym.diarr.onset),order="dmy")
	matSymp[,11] <- as.date(as.character(hosp$sym.other.onset),order="dmy")
	vecAdm <- as.date(as.character(hosp$hosp.adm.dtm),order="dmy")
	vecMinSymp <- apply(cbind(matSymp,vecAdm),1,min,na.rm=TRUE)
	maskTooEarly <- vecMinSymp < 18000
	vecMinSymp[maskTooEarly] <- vecAdm[maskTooEarly]
	
	# Setup the return dataframe
	rtn <- as.data.frame(cbind(
					id  = hosp$patient.key,
					dob = as.date(as.character(hosp$dob),order="dmy"),
					adm = as.date(as.character(hosp$hosp.adm.dtm),order="dmy"),
					dth = (hosp$hosp.dischg.cd == 1),
					sev = (hosp$hosp.ever.nic == 3) | (hosp$hosp.ever.icu == 3),
					first.onset = vecMinSymp))
	
	# Calculate age and assign age groups
	rtn$age <- (rtn$adm - rtn$dob) / 365.25
	for (i in 1:noRows) {
		if (rtn$age[i] > 0 && rtn$age[i] < 5) rtn$ag1[i] <- "1_0t4"
		else if (rtn$age[i] < 18) rtn$ag1[i] <- "2_5t17"
		else if (rtn$age[i] < 65) rtn$ag1[i] <- "3_18t64"
		else if (rtn$age[i] < 110) rtn$ag1[i] <- "4_64plus"
		else {
			browser()
			stop("problem with age group assignemnt 98237492")
		}
		
		if (rtn$age[i] > 2 && rtn$age[i] < 19) rtn$ag2[i] <- "1_3t18"
		else if (rtn$age[i] < 49) rtn$ag2[i] <- "2_19t48"
		else if (rtn$age[i] < 65) rtn$ag2[i] <- "3_49t64"
		else if (rtn$age[i] < 110) rtn$ag2[i] <- "4_65plus"
		else {
			browser()
			stop("problem with age group assignemnt 298374283")
		}
		
		if (rtn$age[i] > 2 && rtn$age[i] < 20) rtn$ag20[i] <- "1_3t19"
		else if (rtn$age[i] < 40) rtn$ag20[i] <- "2_20t39"
		else if (rtn$age[i] < 60) rtn$ag20[i] <- "3_40t59"
		else if (rtn$age[i] < 110) rtn$ag20[i] <- "4_60plus"
		else {
			browser()
			stop("problem with age group assignemnt 928374289")
		}
		
		
	}
	
	# Designate admission weeks
	# st_w_1 <- as.date("29Dec2008") # start week one. Needs checking
	# rtn$adm.wk <- floor((rtn$adm - st_w_1) / 7)
	
	rtn
	
}

plotDiscInc <- function(matInc,vecBins,vecCols) {
	
	# Puts a plot of rectangles on an alreday established plotting surface
	# to be used for incidence
	
	norows <- dim(matInc)[1]
	noweeks <- dim(matInc)[2]
	nobounds <- length(vecBins)
	
	if (length(vecCols)!=norows) stop("matrix and colour vector don't match in 092384209")
	if (noweeks!=nobounds - 1) stop("matrix and bounds vector don't match in 092384209")
	
	for (i in 1:(nobounds-1)) {
		height <- 0
		for (j in 1:norows) {
			tmpinc <- matInc[j,i]
			if (tmpinc > 0) {
				polygon(	c(vecBins[i],vecBins[i+1],vecBins[i+1],vecBins[i]),
						c(height, height, height+tmpinc, height+tmpinc), col=vecCols[j])
				height <- height + tmpinc
			}
		}
	}
	
}

plotHalfHalf <- function(vecMod,vecData,vecBins,vecUp95=NULL,vecLow95=NULL,color="black",linewidth=0.5) {
	
	# Puts a plot of rectangles on an alreday established plotting surface
	# to be used for individual age category incidences to compare 
	
	noweeks <- length(vecMod)
	nobounds <- length(vecBins)
	
	if (length(vecData)!=noweeks) stop("matrix and colour vector don't match in 092384209")
	if (noweeks!=nobounds - 1) stop("matrix and bounds vector don't match in 092384209")
	
	for (i in 1:(nobounds-1)) {
		datheight <- vecData[i]
		polygon(	c(vecBins[i+1]-(vecBins[i+1]-vecBins[i])/2,vecBins[i+1],vecBins[i+1],vecBins[i+1]-(vecBins[i+1]-vecBins[i])/2),
				c(0, 0, datheight, datheight), col=color, border=color,lwd=linewidth)
	}
	
	for (i in 1:(nobounds-1)) {
		modheight <- vecMod[i]
		polygon(	c(vecBins[i],vecBins[i+1]-(vecBins[i+1]-vecBins[i])/2,vecBins[i+1]-(vecBins[i+1]-vecBins[i])/2,vecBins[i]),
				c(0, 0, modheight, modheight), col=NA, border="black",lwd=linewidth)
	}
	
}

extFourAgeICU <- function(sol,vecBins,agnames,varnames=c("dIv1","dIv2","dIv3","dIv4")) {
	
	# Extracts incidence of ICU admission from sol
	# Assumes that sol has a time step of one day and that vecBins are evenly spaced
	# agnames are the names of the age groups and there must be four of them
	# assumes that vecBins is at least of length 2 and is montonically increasing
	
	if(abs(sol$time[2]-sol$time[1]-1.0)>1e-10) stop("arg sol must have dt=1")
	
	noagegroups <- length(agnames)
	noWeeks <- length(vecBins)-1
	rtn <- matrix(0,nrow = noagegroups,ncol=noWeeks)
	row.names(rtn) <- agnames
	nosoltimes <- dim(sol)[1]
	
	binsindex <- 1
	maxbinsindex <- length(vecBins)
	eps <- 1e-10
	currentinc <- as.vector(c(0,0,0,0),mode="numeric")
	
	# Need to scroll through solution time
	for (i in 1:(nosoltimes-1)) {
		
		# Assign current time and bin boundaries
		curtime <- sol$time[i]
		curmin <- vecBins[binsindex]
		curmax <- vecBins[binsindex+1]
		
		# If current time falls within the current bin, add the incidence
		if ((curtime > curmin - eps) && (curtime < curmax + eps)) {
			for (j in 1:length(varnames)) {
				rtn[agnames[j],binsindex] <- rtn[agnames[j],binsindex] + sol[i+1,varnames[j]] - sol[i,varnames[j]] 
			}
		} 
		
		# Check what the next time is and update bins index if needed
		nexttime <- sol$time[i+1]
		if (nexttime > (curmax - eps) && (binsindex < maxbinsindex - 1)) binsindex <- binsindex + 1
		
	}
	
	rtn	
	
}


# Routine to read the used in the CID paper and make equivalent fields
readJoeEFlu <- function(fnFile) {
	
	browser()
	raw <- read.csv(fnFile)
	
}

multimatch <- function(x,vec) {
	rtn <- NULL
	size <- length(vec)
	start <- 1	
	while (start <= size && !is.na(match(x,vec[start:size]))) {
		current <- match(x,vec[start:size]) + start - 1
		rtn <- c(rtn,current)
		start <- current + 1
	}
	rtn
}

post_proc_sero <- function(baseline,laboratory1,laboratory2,symptoms,sympdiary,questionnaire,recruitment,base_date=as.date("04Jul2009")) {
	
	# This function takes the baseline individual data, lab data and symptom reports as input
	# It does all the required post-processing, and adds the data onto the baseline return 
	# - generating a tested / not tested field
	# - generating an age group field
	# - generating a reported symptoms/ not reported symptoms field...	
	
	# Put in the tested field
	norowbaseline <- dim(baseline)[1]
	baseline <- cbind(baseline,pop_recruit=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,recruit_source=rep(-1,norowbaseline))
	baseline <- cbind(baseline,hh_recruit=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,hh_telsource=rep(-1,norowbaseline))
	baseline <- cbind(baseline,tested1=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,tested2=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,final_tested=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,hh_tested=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,test1_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,test2_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,final_fourfold=rep(-1,norowbaseline))
	baseline <- cbind(baseline,b20=rep(NA,norowbaseline))
	baseline <- cbind(baseline,f40=rep(NA,norowbaseline))
	baseline <- cbind(baseline,b40=rep(NA,norowbaseline))
	baseline <- cbind(baseline,f20=rep(NA,norowbaseline))
	baseline <- cbind(baseline,pre_titre=rep(-1,norowbaseline))
	baseline <- cbind(baseline,pre_titre_agg=rep(-1,norowbaseline))
	baseline <- cbind(baseline,post_titre=rep(-1,norowbaseline))
	baseline <- cbind(baseline,post_titre_agg=rep(-1,norowbaseline))
	baseline <- cbind(baseline,age=rep(-1,norowbaseline))
	baseline <- cbind(baseline,ag1=rep(0,norowbaseline))
	baseline <- cbind(baseline,ag10=rep(0,norowbaseline))
	baseline <- cbind(baseline,ag20=rep(0,norowbaseline))
	baseline <- cbind(baseline,sy_ph_crude=rep(-1,norowbaseline))
	baseline <- cbind(baseline,sy_ph_tpu=rep(-1,norowbaseline))
	baseline <- cbind(baseline,child_hh=rep(-1,norowbaseline))
	
	# New SR edits April 30
	# Will do these manipulations out of the setup file now
	
	# baseline <- cbind(baseline,bl_wk=rep(-1,norowbaseline))
	# baseline <- cbind(baseline,fu_wk=rep(-1,norowbaseline))
	
	baseline <- cbind(baseline,blood.pair=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,symp.fever.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.ILI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,symp.ARI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,diary.fever.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.ILI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,diary.ARI.rp=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.reported=rep(FALSE,norowbaseline))
	baseline <- cbind(baseline,quest.fever=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.ILI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,quest.ARI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.fever=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.ILI=rep(-1,norowbaseline))
	baseline <- cbind(baseline,any.ARI=rep(-1,norowbaseline))
	
	# Some rows added by SR to consolidate Danny's edits for the regression model
	baseline  <-cbind(baseline,profession2=rep(-1,norowbaseline))
	baseline  <-cbind(baseline,education2=rep(-1,norowbaseline))
	
	norowsymptoms <- dim(symptoms)[1]
	symptoms <- cbind(symptoms,id=rep(NA,norowsymptoms))
	symptoms <- cbind(symptoms,fever.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,cough.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,sputum.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,sore_t.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,running_n.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,shortness_b.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,chills.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,headache.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,vomit.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,diarrehea.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,myalgia.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ILI=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ILI.rp=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ARI=rep(-1,norowsymptoms))
	symptoms <- cbind(symptoms,ARI.rp=rep(-1,norowsymptoms))
	symptoms$id <- paste(symptoms$houseid,symptoms$sub_no,sep="-")
	
	sympdiary <- sympdiary[,1:13] #remove unused columns
	sympdiary <- sympdiary[1:1962,] #remove unused rows
	norowdiary <- dim(sympdiary)[1]
	sympdiary <- cbind(sympdiary,temperature.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,temperature_up.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,temperature_low.clean=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,fever=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,fever.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,chills=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,chills.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,headache=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,headache.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sore_t=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sore_t.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,cough=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,cough.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sputum=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,sputum.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,stuffy_n=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,stuffy_n.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,running_n=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,running_n.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,myalgia=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,myalgia.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,vomit=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,vomit.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,diarrehea=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,diarrehea.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ILI=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ILI.rp=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ARI=rep(-1,norowdiary))
	sympdiary <- cbind(sympdiary,ARI.rp=rep(-1,norowdiary))
	
	questionnaire <- questionnaire[1:861,] #remove unused rows
	norowquest <- dim(questionnaire)[1]
	questionnaire <- cbind(questionnaire,subid=rep(NA,norowquest))
	questionnaire <- cbind(questionnaire,ILI=rep(-1,norowquest))
	questionnaire <- cbind(questionnaire,ARI=rep(-1,norowquest))
	
	# Correct some formats
	baseline[,"base_ind_date"] <- as.date(as.character(baseline[,"base_ind_date"]),order="dmy")
	baseline[,"fu_ind_date"] <- as.date(as.character(baseline[,"fu_ind_date"]), order="dmy")
	baseline[,"base_dob"] <- as.date(as.character(baseline[,"base_dob"]),order="dmy")
	baseline$age <- (base_date - baseline$base_dob)/ 365.25
	
	baseline$blood.pair <- ifelse (baseline$base_blood == "Yes" & baseline$fu_blood == "Yes", 1, 0)
	
	str1<-substr(questionnaire$id,1,5)
	str2<-substr(questionnaire$id,6,6)
	str0 <- noquote(paste("S0",str1))
	str1.new <- gsub( "[^[:alnum:]]", "", str0)
	questionnaire$subid <- paste(str1.new,str2,sep="-")
	
	temp.tmp <- as.numeric(levels(sympdiary$temp)[sympdiary$temp])
	
	# Determine prevalence of ILI and ARI on phone symptoms
	for (m in 1:norowsymptoms) {
		if (symptoms$fever[m]==1 && (symptoms$cough[m]==1 || symptoms$sore_t[m]==1)) 
			symptoms$ILI[m] <- 1
		else symptoms$ILI[m] <- 0
		
		symp.sum <- symptoms$fever[m]+symptoms$cough[m]+symptoms$sputum[m]+symptoms$sore_t[m]+symptoms$running_n[m]+symptoms$myalgia[m]
		if (symp.sum>=2) symptoms$ARI[m] <- 1
		else symptoms$ARI[m] <- 0
		
	}
	
	
	# Set symptoms on phone reporting during reporting period
	for (m in 1:norowsymptoms) {
		if (!is.na(symptoms$record[m]) && symptoms$record[m] == 1) {
			record_start <- m
			record_finish <- m -1 + symptoms$record.num[m]
			for (k in record_start: record_finish) {
				if (symptoms$fever[k]==1) symptoms$fever.rp[m] <- 1
				else symptoms$fever.rp[m] <- 0
				if (symptoms$cough[k]==1) symptoms$cough.rp[m] <- 1
				else symptoms$cough.rp[m] <- 0
				if (symptoms$sputum[k]==1) symptoms$sputum.rp[m] <- 1
				else symptoms$sputum.rp[m] <- 0
				if (symptoms$sore_t[k]==1) symptoms$sore_t.rp[m] <- 1
				else symptoms$sore_t.rp[m] <- 0
				if (symptoms$running_n[k]==1) symptoms$running_n.rp[m] <- 1
				else symptoms$running_n.rp[m] <- 0
				if (symptoms$shortness_b[k]==1) symptoms$shortness_b.rp[m] <- 1
				else symptoms$shortness_b.rp[m] <- 0
				if (symptoms$chills[k]==1) symptoms$chills.rp[m] <- 1
				else symptoms$chills.rp[m] <- 0
				if (symptoms$headache[k]==1) symptoms$headache.rp[m] <- 1
				else symptoms$headache.rp[m] <- 0
				if (symptoms$vomit[k]==1) symptoms$vomit.rp[m] <- 1
				else symptoms$vomit.rp[m] <- 0
				if (symptoms$diarrehea[k]==1) symptoms$diarrehea.rp[m] <- 1
				else symptoms$diarrehea.rp[m] <- 0
				if (symptoms$myalgia[k]==1) symptoms$myalgia.rp[m] <- 1                   
				else symptoms$myalgia.rp[m] <- 0
				if (symptoms$ILI[k]==1) symptoms$ILI.rp[m] <- 1
				else symptoms$ILI.rp[m] <- 0
				if (symptoms$ARI[k]==1) symptoms$ARI.rp[m] <- 1
				else symptoms$ARI.rp[m]}
		}
	}                  
	
	
	# Perform computation on symptoms diary
	temp.tmp <- as.numeric(levels(sympdiary$temp)[sympdiary$temp])
	#sympdiary$temp_up and sympdiary$temp_low are numerics, instead of factor
	
	for (j in 1:norowdiary) {    
		
		# Convert fahrenheit to celsius
		if (is.na(temp.tmp[j])) sympdiary$temperature.clean[j] <- -1
		else {if (temp.tmp[j] >= 45) 
				sympdiary$temperature.clean[j] <- (temp.tmp[j] - 32) *5/9
			else sympdiary$temperature.clean[j] <- temp.tmp[j]}
		
		if (is.na(sympdiary$temp_up[j])) sympdiary$temperature_up.clean[j] <- -1
		else {if (sympdiary$temp_up[j] >= 45) 
				sympdiary$temperature_up.clean[j] <- (sympdiary$temp_up[j] - 32) *5/9
			else sympdiary$temperature_up.clean[j] <- sympdiary$temp_up[j]}
		
		if (is.na(sympdiary$temp_low[j])) sympdiary$temperature_low.clean[j] <- -1
		else {if (sympdiary$temp_low[j] >= 45) 
				sympdiary$temperature_low.clean[j] <- (sympdiary$temp_low[j] - 32) *5/9
			else sympdiary$temperature_low.clean[j] <- sympdiary$temp_low[j]}
		
		# Assign symptoms on diary
		if (sympdiary$temperature.clean[j] >= 37.5) sympdiary$fever[j] <- 1        
		else {if (sympdiary$temperature.clean[j] == -1)
			{if (sympdiary$temperature_up.clean[j] == -1 && sympdiary$temperature_low.clean[j] == -1) 
					sympdiary$fever[j] <- -1
				else {if (sympdiary$temperature_low.clean[j] >= 37.5) sympdiary$fever[j] <- 1
					else sympdiary$fever[j] <- 0}
			}
			else sympdiary$fever[j] <- 0
		}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$chills[j] <- -1 
		else {tmp <- regexpr("1",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$chills[j] <- 1
			else if (tmp[1] == -1) sympdiary$chills[j] <- 0
			else sympdiary$chills[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$headache[j] <- -1
		else {tmp <- regexpr("2",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$headache[j] <- 1
			else if (tmp[1] == -1) sympdiary$headache[j] <- 0        
			else sympdiary$headache[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$sore_t[j] <- -1
		else {tmp <- regexpr("3", sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$sore_t[j] <- 1
			else if (tmp[1] == -1) sympdiary$sore_t[j] <- 0
			else sympdiary$sore_t[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$cough[j] <- -1
		else {tmp <- regexpr("4",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$cough[j] <- 1
			else if (tmp[1] == -1) sympdiary$cough[j] <- 0
			else sympdiary$cough[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$sputum[j] <- -1
		else {tmp <- regexpr("5",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$sputum[j] <- 1
			else if (tmp[1] == -1) sympdiary$sputum[j] <- 0
			else sympdiary$sputum[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$stuffy_n[j] <- -1
		else {tmp <- regexpr("6",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$stuffy_n[j] <- 1
			else if (tmp[1] == -1) sympdiary$stuffy_n[j] <- 0
			else sympdiary$stuffy_n[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$running_n[j] <- -1
		else {tmp <- regexpr("7",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$running_n[j] <- 1
			else if (tmp[1] == -1) sympdiary$running_n[j] <- 0
			else sympdiary$runnning_n[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$myalgia[j] <- -1
		else {tmp <- regexpr("8",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$myalgia[j] <- 1
			else if (tmp[1] == -1) sympdiary$myalgia[j] <- 0
			else sympdiary$myalgia[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$vomit[j] <- -1
		else {tmp <- regexpr("9",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$vomit[j] <- 1
			else if (tmp[1] == -1) sympdiary$vomit[j] <- 0
			else sympdiary$vomit[j] <- -1}
		
		if (is.na(sympdiary$symptoms[j])) sympdiary$diarrehea[j] <- -1
		else {tmp <- regexpr("10",sympdiary$symptoms[j])
			if (tmp[1] > 0) sympdiary$diarrehea[j] <- 1
			else if (tmp[1] == -1) sympdiary$diarrehea[j] <- 0
			else sympdiary$diarrehea[j] <- -1}
		
		
		# Determine prevalence of ILI and ARI on symptoms diary
		if (sympdiary$fever[j]==1 && (sympdiary$cough[j]==1 || sympdiary$sore_t[j]==1)) 
			sympdiary$ILI[j] <- 1
		else sympdiary$ILI[j] <- 0
		
		diary.sum <- sympdiary$fever[j]+sympdiary$cough[j]+sympdiary$sputum[j]+sympdiary$sore_t[j]+sympdiary$running_n[j]+sympdiary$myalgia[j]
		if (diary.sum>=2) sympdiary$ARI[j] <- 1
		else sympdiary$ARI[j] <- 0
		
	}
	
	# Set symptoms on symptoms diary during reporting period
	for (j in 1:norowdiary) {
		if (!is.na(sympdiary$record.clean[j]) && sympdiary$record.clean[j] == 1) {
			record_start <- j
			record_finish <- j -1 + sympdiary$record.num[j]
			for (k in record_start: record_finish) {
				if (sympdiary$fever[k]==1) sympdiary$fever.rp[j] <- 1
				else sympdiary$fever.rp[j] <- 0
				if (sympdiary$chills[k]==1) sympdiary$chills.rp[j] <- 1
				else sympdiary$chills.rp[j] <- 0
				if (sympdiary$headache[k]==1) sympdiary$headache.rp[j] <- 1
				else sympdiary$headache.rp[j] <- 0
				if (sympdiary$sore_t[k]==1) sympdiary$sore_t.rp[j] <- 1
				else sympdiary$sore_t.rp[j] <- 0
				if (sympdiary$cough[k]==1) sympdiary$cough.rp[j] <- 1
				else sympdiary$cough.rp[j] <- 0
				if (sympdiary$sputum[k]==1) sympdiary$sputum.rp[j] <- 1
				else sympdiary$sputum.rp[j] <- 0
				if (sympdiary$stuffy_n[k]==1) sympdiary$stuffy_n.rp[j] <- 1
				else sympdiary$stuffy_n.rp[j] <- 0
				if (sympdiary$running_n[k]==1) sympdiary$running_n.rp[j] <- 1
				else sympdiary$running_n.rp[j] <- 0
				if (sympdiary$myalgia[k]==1) sympdiary$myalgia.rp[j] <- 1
				else sympdiary$myalgia.rp[j] <- 0
				if (sympdiary$vomit[k]==1) sympdiary$vomit.rp[j] <- 1
				else sympdiary$vomit.rp[j] <- 0
				if (sympdiary$diarrehea[k]==1) sympdiary$diarrehea.rp[j] <- 1
				else sympdiary$diarrehea.rp[j] <- 0
				if (sympdiary$ILI[k]==1) sympdiary$ILI.rp[j] <- 1
				else sympdiary$ILI.rp[j] <- 0
				if (sympdiary$ARI[k]==1) sympdiary$ARI.rp[j] <- 1
				else sympdiary$ARI.rp[j] <- 0}
		}
	}          
	
	for (r in 1:norowquest) {
		
		#Correct some typos from inconsistencies
		if (questionnaire$id[r] == 900902) questionnaire$no_doctor[r] <- 1
		
		# Determine prevalence of ILI and ARI on contact questionnaire 
		if (questionnaire$fever.clean[r]==1 && (questionnaire$cough.clean[r]==1 || questionnaire$sore_t.clean[r]==1))
			questionnaire$ILI[r] <- 1
		else questionnaire$ILI[r] <- 0
		
		quest.sum <- questionnaire$fever.clean[r]+questionnaire$cough.clean[r]+questionnaire$sputum.clean[r]+questionnaire$sore_t.clean[r]+questionnaire$running_n.clean[r]+questionnaire$myalgia.clean[r]
		if (quest.sum>=2) questionnaire$ARI[r] <- 1
		else questionnaire$ARI[r] <- 0
		
	}
	
	## -- Kendra's editing ends -- ##
	for (i in 1:norowbaseline) {
		
		# Some FU interview dates are in year 2000
		# KW added on 4 June 2010
		basedate1 <- as.date("01Jan1900"); basedate1 <- as.numeric(basedate1)
		basedate2 <- as.date("31Dec1900"); basedate2 <- as.numeric(basedate2)
		
		if (!is.na(baseline$fu_ind_date[i])) {
			if (baseline$fu_ind_date[i] <= basedate2 && baseline$fu_ind_date[i] >= basedate1) {
				tmp <- as.numeric(baseline$fu_ind_date[i])
				tmp <- baseline$fu_ind_date[i] + 40177 #convert those in 1900 to 2010
				baseline$fu_ind_date[i] <- as.date(as.character(tmp, order="dmy"))
			}
		}
		
		
		# Tested or not during test1
		if (baseline$ind_id[i] %in% laboratory1$Our_ref) baseline$tested1[i] <- TRUE
		
		# Result of the testing
		if (baseline$tested1[i]) {           
			lab1row <- match(baseline$ind_id[i],laboratory1$Our_ref)
			if (is.na(lab1row)) stop("There shouldn't be an NA here.")
			baseline$b20[i] <- laboratory1$Base20[lab1row]
			baseline$f40[i] <- laboratory1$Post40[lab1row]
			baseline$b40[i] <- laboratory1$Base40[lab1row]
			baseline$f20[i] <- laboratory1$Post20[lab1row]
			if (baseline$b20[i]==0 || baseline$b40[i]==0) baseline$pre_titre[i] <- 0
			if (baseline$b20[i]==1) baseline$pre_titre[i] <- 20
			if (baseline$b40[i]==1) baseline$pre_titre[i] <- 40
			if (baseline$f20[i]==0 || baseline$f40[i]==0) baseline$post_titre[i] <- 0
			if (baseline$f20[i]==1) baseline$post_titre[i] <- 20
			if (baseline$f40[i]==1) baseline$post_titre[i] <- 40
			if (baseline$b20[i]==0 && baseline$f40[i]==1) baseline$test1_fourfold[i] <- 1
			else baseline$test1_fourfold[i] <- 0
			
		}
		
		## -- Kendra edited: -- ##
		
		# Tested or not
		if (baseline$ind_id[i] %in% laboratory2$ind_id) baseline$tested2[i] <- TRUE
		
		# Result of the testing
		if (baseline$tested2[i]) {           
			lab2row <- match(baseline$ind_id[i],laboratory2$ind_id)
			if (is.na(lab2row)) stop("There shouldn't be an NA here.")
			baseline$test2_fourfold[i] <- laboratory2$Kpos[lab2row]
			baseline$pre_titre[i] <- laboratory2$Base_titre[lab2row]
			baseline$post_titre[i] <- laboratory2$Post_titre[lab2row]
		}
		
		baseline$pre_titre_agg[i] <- baseline$pre_titre[i]
		baseline$post_titre_agg[i] <- baseline$post_titre[i]
		
		if (baseline$pre_titre[i]==0 || baseline$pre_titre[i]==5) baseline$pre_titre_agg[i] <- 10
		if (baseline$post_titre[i]==0 || baseline$post_titre[i]==5) baseline$post_titre_agg[i] <-10
		
		if (baseline$tested1[i]==TRUE || baseline$tested2[i]==TRUE) baseline$final_tested[i] <- TRUE
		
		if (baseline$tested1[i]==TRUE && baseline$tested2[i]==FALSE) baseline$final_fourfold[i] <- baseline$test1_fourfold[i]
		if (baseline$tested1[i]==FALSE && baseline$tested2[i]==TRUE) baseline$final_fourfold[i] <- baseline$test2_fourfold[i] 
		if (baseline$tested1[i]==TRUE && baseline$tested2[i]==TRUE) baseline$final_fourfold[i] <- baseline$test2_fourfold[i]
		
		#Correct some formats
		if (baseline$ind_id[i] == "S090433-0") baseline$base_ind_date[i] <- as.date("27Aug2009")
		
		if (baseline$ind_id[i] == "S090002-0") {baseline$base_dob[i] <- as.date("04Oct2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090008-0") {baseline$base_dob[i] <- as.date("11Feb2004")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090033-0") {baseline$base_dob[i] <- as.date("10Apr2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090048-2") {baseline$base_dob[i] <- as.date("10Aug2001")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090056-0") {baseline$base_dob[i] <- as.date("01Apr2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090061-4") {baseline$base_dob[i] <- as.date("10Oct1918")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090066-3") {baseline$base_dob[i] <- as.date("01Jul1914")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090076-0") {baseline$base_dob[i] <- as.date("12Jun2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090093-0") {baseline$base_dob[i] <- as.date("12Feb2009")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090093-1") {baseline$base_dob[i] <- as.date("12Jan2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090094-2") {baseline$base_dob[i] <- as.date("01Jan1928")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090095-0") {baseline$base_dob[i] <- as.date("09Sep2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090097-0") {baseline$base_dob[i] <- as.date("01Apr2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090097-1") {baseline$base_dob[i] <- as.date("01Apr2006")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090102-3") {baseline$base_dob[i] <- as.date("01Jan1924")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-0") {baseline$base_dob[i] <- as.date("01Feb2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-1") {baseline$base_dob[i] <- as.date("04Feb2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090105-2") {baseline$base_dob[i] <- as.date("05Jul2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090116-0") {baseline$base_dob[i] <- as.date("09Jun2000")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090117-4") {baseline$base_dob[i] <- as.date("06Apr1927")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090117-5") {baseline$base_dob[i] <- as.date("26Nov1921")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090125-0") {baseline$base_dob[i] <- as.date("06Jul2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-0") {baseline$base_dob[i] <- as.date("01Nov2008")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-1") {baseline$base_dob[i] <- as.date("01May2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090128-2") {baseline$base_dob[i] <- as.date("01Dec2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090143-0") {baseline$base_dob[i] <- as.date("04Dec2002")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090144-3") {baseline$base_dob[i] <- as.date("01Jan1924")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090152-0") {baseline$base_dob[i] <- as.date("07Jan2003")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090157-0") {baseline$base_dob[i] <- as.date("07Nov2007")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090168-0") {baseline$base_dob[i] <- as.date("06Dec2007")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090173-3") {baseline$base_dob[i] <- as.date("06Feb1921")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090178-0") {baseline$base_dob[i] <- as.date("01Jan2005")   
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090194-1") {baseline$base_dob[i] <- as.date("04Aug1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090202-2") {baseline$base_dob[i] <- as.date("01Jan1911") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090203-1") {baseline$base_dob[i] <- as.date("06Aug2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090204-2") {baseline$base_dob[i] <- as.date("29Jul1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090207-3") {baseline$base_dob[i] <- as.date("01Jan1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090213-0") {baseline$base_dob[i] <- as.date("07Feb2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090222-1") {baseline$base_dob[i] <- as.date("07Dec2000")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090245-0") {baseline$base_dob[i] <- as.date("09Apr2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090249-1") {baseline$base_dob[i] <- as.date("01Sep1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090266-1") {baseline$base_dob[i] <- as.date("03Oct2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090267-4") {baseline$base_dob[i] <- as.date("01Jan1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090275-0") {baseline$base_dob[i] <- as.date("06Apr2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090283-3") {baseline$base_dob[i] <- as.date("01Jan1929")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090285-0") {baseline$base_dob[i] <- as.date("10Feb2009") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090293-0") {baseline$base_dob[i] <- as.date("29Aug2007") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090300-1") {baseline$base_dob[i] <- as.date("11Feb2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090332-1") {baseline$base_dob[i] <- as.date("10Aug1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090332-2") {baseline$base_dob[i] <- as.date("15Apr1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090339-2") {baseline$base_dob[i] <- as.date("01Mar1924") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090334-1") {baseline$base_dob[i] <- as.date("01Jan1926") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090363-0") {baseline$base_dob[i] <- as.date("01Jan1908") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090371-4") {baseline$base_dob[i] <- as.date("01Jan1920")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090372-0") {baseline$base_dob[i] <- as.date("10Jul2000") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090373-2") {baseline$base_dob[i] <- as.date("01Jan1922")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090379-0") {baseline$base_dob[i] <- as.date("01May2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090379-1") {baseline$base_dob[i] <- as.date("01May2005")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090401-1") {baseline$base_dob[i] <- as.date("01Jan1920") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090412-1") {baseline$base_dob[i] <- as.date("04May2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090413-0") {baseline$base_dob[i] <- as.date("01Jan2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090415-0") {baseline$base_dob[i] <- as.date("30Aug2003") #updated on 7 June 2010
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090418-0") {baseline$base_dob[i] <- as.date("01Jan2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090418-1") {baseline$base_dob[i] <- as.date("01Jan2002") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090430-0") {baseline$base_dob[i] <- as.date("06Apr2006") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}    
		if (baseline$ind_id[i] == "S090430-1") {baseline$base_dob[i] <- as.date("09Apr2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090454-0") {baseline$base_dob[i] <- as.date("06Apr2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090464-1") {baseline$base_dob[i] <- as.date("01Jan1921") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090470-4") {baseline$base_dob[i] <- as.date("01Sep1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090476-4") {baseline$base_dob[i] <- as.date("01Jan1925") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090486-0") {baseline$base_dob[i] <- as.date("01Jan2004") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090513-1") {baseline$base_dob[i] <- as.date("05Jan2001") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090518-3") {baseline$base_dob[i] <- as.date("01Jan1913") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090525-0") {baseline$base_dob[i] <- as.date("11Jun2001")
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090525-3") {baseline$base_dob[i] <- as.date("01Jan1928") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090531-4") {baseline$base_dob[i] <- as.date("01Jan1917") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090550-1") {baseline$base_dob[i] <- as.date("08May2000") #updated on 7 June 2010
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090556-0") {baseline$base_dob[i] <- as.date("04Dec2008")  
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090583-1") {baseline$base_dob[i] <- as.date("06Nov2005") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090584-2") {baseline$base_dob[i] <- as.date("29May1926") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		if (baseline$ind_id[i] == "S090587-0") {baseline$base_dob[i] <- as.date("04Aug2003") 
			baseline$age[i] <- (base_date - baseline$base_dob[i])/ 362.5}
		
		if (baseline$ind_id[i] == "S090324-2") {
			baseline$fu_ind_date[i] <- as.date("29Dec2009") #added on 06 July 2010
			baseline$fu_attendance[i] <- "Yes"
			baseline$fu_blood[i] <- "Yes"}
		
		# Add by KMW on 1 June 2010
		if (baseline$ind_id[i] == "S090008-0") baseline$profession[i] <- 12
		if (baseline$ind_id[i] == "S090008-1") baseline$profession[i] <- 12
		if (baseline$ind_id[i] == "S090358-0") baseline$profession[i] <- 3
		
		## -- Kendra's editing ends -- ##
		
		# Set age and age groups
		if (is.na(baseline[i,"age"])) baseline$ag1[i] <- -1
		else if (baseline[i,"age"] < 19) baseline[i,"ag1"] <- 1
		else if (baseline[i,"age"] <= 48) baseline[i,"ag1"] <- 2
		else if (baseline[i,"age"] < 65) baseline[i,"ag1"] <- 3
		else if (baseline[i,"age"] < 110) baseline[i,"ag1"] <- 4
		else stop("Problem with the age group ag1 allocation")
		
		#Add by KMW on 8 June 2010
		if (is.na(baseline[i,"age"])) baseline$ag10[i] <- -1
		else if (baseline[i,"age"] < 10) baseline$ag10[i] <- 1
		else if (baseline[i,"age"] < 20) baseline$ag10[i] <-  2
		else if (baseline[i,"age"] < 30) baseline$ag10[i] <- 3
		else if (baseline[i,"age"] < 40) baseline$ag10[i] <- 4
		else if (baseline[i,"age"] < 50) baseline$ag10[i] <- 5
		else if (baseline[i,"age"] < 60) baseline$ag10[i] <- 6
		else if (baseline[i,"age"] < 70) baseline$ag10[i] <- 7
		else if (baseline[i,"age"] < 80)  baseline$ag10[i] <- 8
		else if (baseline[i,"age"] < 110)  baseline$ag10[i] <- 9
		
		if (is.na(baseline[i,"age"])) baseline$ag20[i] <- -1
		else if (baseline[i,"age"] < 20) baseline[i,"ag20"] <- 1
		else if (baseline[i,"age"] < 40) baseline[i,"ag20"] <- 2
		else if (baseline[i,"age"] < 60) baseline[i,"ag20"] <- 3
		else if (baseline[i,"age"] < 110) baseline[i,"ag20"] <- 4
		else stop("Problem with ag20 allocation")
		
		## -- Kendra edited: -- ##
		# Reported symptoms or not through the phone
		if (baseline$ind_id[i] %in% symptoms$id) baseline$symp.reported[i] <- TRUE
		
		# Result of the reporting
		sympindex <- match(baseline$ind_id[i],symptoms$id)
		if (baseline$symp.reported[i]) {
			if (is.na(sympindex)) stop("There shouldn't be an NA in sympindex among reported")
			baseline$symp.fever.rp[i] <- symptoms$fever.rp[sympindex]
			baseline$symp.ILI.rp[i] <- symptoms$ILI.rp[sympindex]
			baseline$symp.ARI.rp[i] <- symptoms$ARI.rp[sympindex]          
		}
		
		# Reported symptoms or not on symptoms diary 
		if (baseline$ind_id[i] %in% sympdiary$id) baseline$diary.reported[i] <- TRUE
		
		# Result of the reporting
		diaryindex <- match(baseline$ind_id[i], sympdiary$id)
		if (baseline$diary.reported[i]) {
			if (is.na(diaryindex)) stop("There shouldn't be an NA in diaryindex among reported")
			baseline$diary.fever.rp[i] <- sympdiary$fever.rp[diaryindex]
			baseline$diary.ILI.rp[i] <- sympdiary$ILI.rp[diaryindex]
			baseline$diary.ARI.rp[i] <- sympdiary$ARI.rp[diaryindex]
		}
		
		
		# Reported symptoms or not on contact questionniare
		if (baseline$ind_id[i] %in% questionnaire$subid) baseline$quest.reported[i] <- TRUE
		
		# Result of the reporting
		questindex <- match(baseline$ind_id[i], questionnaire$subid)
		if (baseline$quest.reported[i]) {
			if (is.na(questindex)) stop("There shouldn't be an NA in questindex among reported")
			baseline$quest.fever[i] <- questionnaire$fever[questindex]
			baseline$quest.ILI[i] <- questionnaire$ILI[questindex]
			baseline$quest.ARI[i] <- questionnaire$ARI[questindex]
		}
		
		if (baseline$symp.fever.rp[i]==1 || baseline$diary.fever.rp[i]==1 || baseline$quest.fever[i]==1) 
			baseline$any.fever[i] <- 1
		else baseline$any.fever[i] <- 0
		
		if (baseline$symp.ILI.rp[i]==1 || baseline$diary.ILI.rp[i]==1 || baseline$quest.ILI[i]==1)
			baseline$any.ILI[i] <- 1
		else baseline$any.ILI[i] <- 0
		
		if (baseline$symp.ARI.rp[i]==1 || baseline$diary.ARI.rp[i]==1 || baseline$quest.ARI[i]==1) 
			baseline$any.ARI[i] <- 1
		else baseline$any.ARI[i] <- 0
		
		
		#Where did POP recruit this household
		if (baseline$ind_id[i] %in% recruitment$ind_id) baseline$pop_recruit[i] <- recruitment$Match[i]
		recruitindex <- match(baseline$ind_id[i], recruitment$ind_id)
		if (baseline$pop_recruit[i]) baseline$recruit_source[i] <- recruitment$Source[recruitindex]
		
		
		## -- Kendra's editing ends -- ##
		
		
		# Most crude phone symptoms
		if (!is.na(sympindex)) baseline[i,"sy_ph_crude"] <- 1
		else baseline[i,"sy_ph_crude"] <- 0
		
		if (!is.na(sympindex) && symptoms[sympindex,"self_r"] == 1) baseline[i,"sy_ph_tpu"] <- 1
		else baseline[i,"sy_ph_tpu"] <- 0
		
		# Copying in Danny's reclassifications
		if (baseline[i,"profession"] %in% c(1,2)) baseline[i,"profession2"] <- 1
		else if (baseline[i,"profession"] %in% c(3,4)) baseline[i,"profession2"] <- 2
		else if (baseline[i,"profession"] %in% c(5,6)) baseline[i,"profession2"] <- 3
		else if (baseline[i,"profession"] %in% c(7,8,9)) baseline[i,"profession2"] <-4
		else if(baseline[i,"profession"] %in% c(10,11,13,14)) baseline[i,"profession2"] <- 5
		else if(baseline[i,"profession"] %in% c(12)) baseline[i,"profession2"] <- 6
		
#		# date of clinic visit, use the median date as the cutoff date 	
#		baseline  <-cbind(baseline,date_visit=rep(-1,norowbaseline))
#	
#		# cute off for every two weeks
#
#		if (baseline[i,"base_ind_date"] < as.date("4jul9")+14) baseline[i,"date_visit"] <- 1
#		else if (baseline[i,"base_ind_date"] < as.date("18jul9")+14) baseline[i,"date_visit"] <-2
#		else if (baseline[i,"base_ind_date"] < as.date("1Aug9")+14) baseline[i,"date_visit"] <-3
#		else if (baseline[i,"base_ind_date"] < as.date("14Aug9")+14) baseline[i,"date_visit"] <-4
#		else if (baseline[i,"base_ind_date"] < as.date("28Aug9")+28) baseline[i,"date_visit"] <-5
#		#else if(baseline[i,"base_ind_date"] < as.date("11Sep2009")+14) baseline[i,"date_visit"] <-6
		
		#if (is.na(baseline[i,"education"])) baseline$education2[i] <- -1
		if (baseline[i,"education"] ==0) baseline[i,"education2"] <- 0
		else if (baseline[i,"education"] ==1) baseline[i,"education2"] <- 1
		else if (baseline[i,"education"] ==2) baseline[i,"education2"] <- 2
		else if (baseline[i,"education"] ==3) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==4) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==5) baseline[i,"education2"] <- 3
		else if(baseline[i,"education"] ==6) baseline[i,"education2"] <- 4
		else stop("Problem with the grouping education level")
		
	}          
	
	## More editing from Kendra
	# Household stuff
	for (i in 1:norowbaseline) {
		
		#Check whether someone in the household has had their paired samples tested
		if (baseline[i,"hh_index"]==0) {
			hh_index_start <- i
			hh_index_finish <- i -1 + baseline[i,"hh_size"]
		}
		
		# Set presence child_household
		if (1 %in% baseline[hh_index_start:hh_index_finish,"ag1"]) baseline[i,"child_hh"] <- 1
		else if (-1 %in% baseline[hh_index_start:hh_index_finish,"ag1"]) {
			if (1 %in% baseline[hh_index_start:hh_index_finish,"clean_is_child"]) 
				baseline[i,"child_hh"] <- 1
			else 
				baseline[i,"child_hh"] <- 2
		}
		else baseline[i,"child_hh"] <- 2
		
		if (TRUE %in% baseline[hh_index_start: hh_index_finish,"final_tested"]) baseline$hh_tested[i] <- TRUE
		
		if (baseline$ind_id[i] == "S090546-x") baseline$hh_tested[i] <- FALSE
		
		# Where did the household recruit from
		if (TRUE %in% baseline[hh_index_start: hh_index_finish, "pop_recruit"]) baseline$hh_recruit[i] <- TRUE
		if (1 %in% baseline[hh_index_start: hh_index_finish, "recruit_source"]) baseline$hh_telsource[i] <- 1
		else if (2 %in% baseline[hh_index_start: hh_index_finish, "recruit_source"]) baseline$hh_telsource[i] <- 2
		
		
	}
	
	# Set up year week of recruitment SR editing here
	
	#baseline$ag1 <- relevel(factor(baseline$ag1),ref="2")
	
	list(base=baseline,labs1=laboratory1,labs2=laboratory2,symp=symptoms,diary=sympdiary,quest=questionnaire,recruit=recruitment)
	
}


sev_g <- function(k,bi,fi,q1,q2,q3) {
	if (identical(bi-k,2)) 			rtn <- 1-q2
	else if (identical(bi-k,1)) 	rtn <- 1-q1
	else if (k >= bi && k < fi-2) 	rtn <- 1
	else if (identical(fi-k,2)) 	rtn <- q2
	else if (identical(fi-k,1)) 	rtn <- q1
	else 							rtn <- 0
	as.numeric(rtn)
}

sev_r <- function(p,index,ya,b,f,N,q1,q2,minw=1,maxw=length(ya)) {
	rtn <- 0
	bi <- b[index]
	fi <- f[index]
	for (i in minw:maxw) rtn <- rtn + sev_g(i,bi,fi,q1,q2)*ya[i]
	rtn <- rtn / p / N
	as.numeric(rtn)
}

sev_like <- function(p,q1,q2,xa,ba,fa,ya,N,noobs=length(xa),verbose=FALSE) {
	rtn <- 0
	for (i in 1:noobs) {
		prob <- sev_r(p,i,ya,ba,fa,N,q1,q2)
		if (prob > 1) rtn <- rtn - 1e10
		else {
			if (xa[i]==1) rtn <- rtn + log(prob)
			else rtn <- rtn + log(1-prob)
		}
	}
	if (verbose) {
		cat("p ",p," val ",rtn,"\n")
		flush.console()
	}
	rtn
}

sev_like_ci <- function(p,target,q1,q2,xa,ba,fa,ya,N,noobs=length(xa)) {
	raw <- sev_like(p,q1,q2,xa,ba,fa,ya,N)
	abs(raw-(target-1.96))
}

sim_study <- function(fnSimData,vecPs,vecAgs,matWkSev,vecWkBins,q1,q2) {
	
	# Function to simulate a sero study
	# Takes a study protocol, a vector of age specific probabilities, 
	# Not sure that this is working at all XXXX
	
	if (is.character(fnSimData)) {
		
		rtn <- read.csv(fnSimData)
		rtn$base_ind_date <- as.date(as.character(rtn$base_ind_date))
		rtn$fu_ind_date <- as.date(as.character(rtn$fu_ind_date))
		
	} else {
		
		rtn <- fnSimData
		
	}
	
	noSubs <- dim(rtn)[1]
	noWks <- length(vecWkBins)
	if (noWks != dim(matWkSev)[2]+1) stop("problem matching wk bins with severity mtrix")
	noAgs <- length(vecPs)
	if (noAgs != length(vecAgs)) stop("problem matching age group prob vector and ag size vector")
	stFirstWk <- vecWkBins[1]
	if (rtn$base_ind_date[1] < stFirstWk) stop("likely a problem with the date format in a study protocol file")
	
	# Cycle through each invidiual and each week
	for (i in 1:noSubs) {
		base_wk <- floor((rtn$base_ind_date[i]-stFirstWk)/7)+1
		fu_wk <- floor((rtn$fu_ind_date[i]-stFirstWk)/7)+1
		ag <- rtn$ag20[i]
		infs <- 0
		if (base_wk - 2 > 0) infs <- infs + matWkSev[ag,base_wk-2]* (1 - q2) / vecPs[ag]
		if (base_wk - 1 > 0) infs <- infs + matWkSev[ag,base_wk-1]* (1 - q1) / vecPs[ag]
		for (j in base_wk:(fu_wk-3)) infs <- infs + matWkSev[ag,j] / vecPs[ag]
		infs <- infs + matWkSev[ag,fu_wk-2] / vecPs[ag] * q2
		infs <- infs + matWkSev[ag,fu_wk-1] / vecPs[ag] * q1
		ind_p <- infs / vecAgs[ag]
		if (ind_p > 1) stop("problem with simulation")
		if (runif(1) < ind_p) rtn$final_fourfold[i] <- 1
		else rtn$final_fourfold[i] <- 0
	}
	
	rtn
	
}

summarize_models <- function(
		lstModels,
		outfunc=exp,
		writetab=TRUE,
		file="modsum.csv",
		sigdigits=3,
		transpose=FALSE) {
	
	# Figure out the number of models
	nomods <- length(lstModels) 
	
	# Make a vector of all coefficients
	allCoeffs <- c()
	for (i in 1:nomods) {
		
		# Select current model results
		mod <- lstModels[[i]]
		
		# Get a list of variables
		vars <- names(mod$coefficients)
		novars <- length(vars)
		
		# Go through each variabel and add it if its not already in the list
		for (j in 1:novars) {
			
			# Get the variable name
			curname <- vars[j]
			
			# Test for the presence of the variable in the master list
			var_present <- (curname %in% allCoeffs)
			
			# If not in the list add it
			if (!(var_present)) allCoeffs <- c(allCoeffs,curname)
			
			# Close the for loop for j
		}
		
		# Close the for loop for i	
	}
	
	# Define the data structures used to extract the information from the models 	
	noCoeffs <- length(allCoeffs)
	matPointEst <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	matLB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	matUB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
	vecAIC <- vector(mode="numeric",length=nomods)
	
	# Loop back though the models and the coeffciients to populate the data structures
	for (i in 1:nomods) {
		
		# Select current model results
		mod <- lstModels[[i]]
		cis <- confint.default(mod)
		
		# Get a list of variables
		vars <- names(mod$coefficients)
		novars <- length(vars)
		
		# Record the AIC
		vecAIC[i] <- mod$aic
		
		# Go through each variabel and add it if its not already in the list
		for (j in 1:novars) {
			
			# Get the variable name
			curname <- vars[j]
			
			# Extract the point estimate and confidence intervals for the parameters
			matPointEst[i,curname] 	<- mod$coefficients[curname]  
			matLB[i,curname] 		<- cis[curname,1]  
			matUB[i,curname] 		<- cis[curname,2]  
			
			# Close the for loop for j
		}
		
		# Close the for loop for i	
	}
	
	# If selected, write a nicely formatted csv table for the parameters and models
	if (writetab) {
		
		if (transpose) {
			
			# Declare the output string
			strTable <- ""
			
			# Put in the first header row
			strTable <- paste(strTable,"Parameter",sep="")
			for (i in 1:noCoeffs) strTable <- paste(strTable,",",allCoeffs[i],",",allCoeffs[i],sep="")
			strTable <- paste(strTable,",AIC\n",sep="")
			
			# Put in the second header row
			strTable <- paste(strTable,"Model",sep="")
			for (i in 1:noCoeffs) strTable <- paste(strTable,",PE,CI",sep="")
			strTable <- paste(strTable,",AIC\n",sep="")
			
			# Output individual model lines, starting with coefficient loop
			for (i in 1:nomods) {
				
				# Pull the name of the current coefficient
				# curname <- allCoeffs[i]
				
				# Put in the name of the coefficient
				strTable <- paste(strTable,i,sep="")
				
				# Cycle through the tables looking at the different models
				for (j in 1:noCoeffs) {
					
					# Itentify the current coefficient
					curname <- allCoeffs[j]
					
					# Put in the point estimates and confidence intervals for each parameter / model combination
					curPE <- signif(outfunc(matPointEst[i,curname]),digits=sigdigits)
					curLB <- signif(outfunc(matLB[i,curname]),digits=sigdigits)
					curUB <- signif(outfunc(matUB[i,curname]),digits=sigdigits)
					
					# Paste in the parameter values and the confidence intervals
					if (is.na(curPE)) {
						
						# Put in the entry for NA results
						strTable <- paste(strTable,",","-",",","-",sep="")
						
					} else {
						
						# Put in the entry for non NA results
						strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
						
					}
					
					# End j loop for coefficients
				}
				
				# Add the AIC at the end of the line, with a return
				mod <- lstModels[[i]]
				curAIC <- round(mod$aic,digits=1)
				strTable <- paste(strTable,",",curAIC,"\n",sep="")
				
				# End the i for loop for models
			}
			
			# End the if clause for transpose
		} else {
			
			# Declare the output string
			strTable <- ""
			
			# Put in the first header row
			strTable <- paste(strTable,",Model 1",sep="")
			if (nomods>1) for (i in 2:nomods) strTable <- paste(strTable,",,Model ",i,sep="")
			strTable <- paste(strTable,"\n",sep="")
			
			# Put in the second header row
			if (nomods>1) for (i in 1:nomods) strTable <- paste(strTable,",Estimate,(95% CI)",sep="")
			strTable <- paste(strTable,"\n",sep="")
			
			# Output individual coefficient lines, starting with coefficient loop
			for (i in 1:noCoeffs) {
				
				# Pull the name of the current coefficient
				curname <- allCoeffs[i]
				
				# Put in the name of the coefficient
				strTable <- paste(strTable,curname,sep="")
				
				# Cycle through the tables looking at the different models
				for (j in 1:nomods) {
					
					# Put in the point estimates and confidence intervals for each parameter / model combination
					curPE <- signif(outfunc(matPointEst[j,curname]),digits=sigdigits)
					curLB <- signif(outfunc(matLB[j,curname]),digits=sigdigits)
					curUB <- signif(outfunc(matUB[j,curname]),digits=sigdigits)
					
					# Paste in the parameter values and the confidence intervals
					if (is.na(curPE)) {
						
						# Put in the entry for NA results
						strTable <- paste(strTable,",","-",",","-",sep="")
						
					} else {
						
						# Put in the entry for non NA results
						strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
						
					}
					
					# End model for loop
				}
				
				# Return at the end of the line
				strTable <- paste(strTable,"\n",sep="")
				
				# End for for coeffs	
			}
			
			
			# Write the row name for the AICs	
			strTable <- paste(strTable,"AIC",sep="")
			
			# Start the for loop for the AICs
			for (i in 1:nomods) {
				
				# Get the current model
				mod <- lstModels[[i]]
				
				# Format the AIC for the current model
				curAIC <- round(mod$aic,digits=1)
				
				# Write the value and the space
				strTable <- paste(strTable,",",curAIC,",",sep="")
			}
			
			# Return at the end of the AIC line
			strTable <- paste(strTable,"\n",sep="")
			
			# End else statement for transpose of 
		}
		
		# Write the string to the selected file
		cat(strTable,file=file)
		
	}
	
	# Return 
	list(pe=matPointEst,lb=matLB,ub=matUB,aic=vecAIC)
	
}

# A function to make a rolling average of p ff based on age of width width
# Doesn't assume that x is sorte, but does assume it can be ordered
# Uses the median value for the y axis if median is true, otherwise uses mean
windowInc <- function(x,p,width=100,use_median=TRUE) {
	
	# Set up some auxilliary variables
	noMeasures <- length(x)
	noRolling <- noMeasures - width + 1
	if (noRolling < 1) stop("not enough data points for that rolling average")
	
	# Setup return values
	vecX <- vector(mode="numeric",length=noRolling)
	vecY <- vector(mode="numeric",length=noRolling)
	
	# Sort x and p
	x_dash <- x[order(x)]
	p_dash <- p[order(x)]
	
	# Populate the rolling averagesim_study
	
	for (i in width:noMeasures) {
		
		# Setup the loop variabels
		start_window <- i - width + 1
		end_window <- i
		
		# Calculate the rolling average
		vecY[i] <- sum(p_dash[start_window:end_window]) / width
		if (use_median) vecX[i] <- median(x_dash[start_window:end_window])
		else vecX[i] <- mean(x_dash[start_window:end_window])
		
	}
	
	list(x=vecX,y=vecY)
	
}


# An aux routine to compute Binomial 95% Confidence Interval
SR_BinOptim <- function(p,P,N,n) {
	guess <- pbinom(n,N,p,lower.tail=TRUE)
	rtn <- P-guess
	rtn
}

# Computes Binomial 95% Confidence Interval
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

# Auxilliary function for the full set of results
srPieFunc <- function(x,y,vec,rscale=1/10,base=4) {
	agcounts <- vec
	totalcount <- sum(vec)
	cumulativecount <- 0
	for (ag in (1:length(vec))) {
		if (totalcount > 0) circle(x=x,y=y,r= rscale*(log(totalcount,base=base)+1),theta=c(cumulativecount/totalcount*360,(cumulativecount+agcounts[ag])/totalcount*360),col=vecColors[ag])
		cumulativecount <- cumulativecount + agcounts[ag] 
	}	
}

# Make the symptom plot
# Data frame for counts, total number of positives, filename, ylim for chart
plot.symptoms <- function(df, N, filename, ylim_in = c(0,1.0), ms=0.3, ...) {
	
	# Define the number of ...
	fw <- 8.3/cm(1)
	fh <- 8.3/cm(1)
	
	
	# Open the pdf file
	pdf(filename,height=fh,width=fw)
	
	# Set some standard parameter options
	par(mai=(c(0.15*fh,0.175*fw,0,0)))
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	par(cex = 10/12)
	
	# Set up the vectors
	vecLev1 <- names(df)
	noLev1 <- length(vecLev1)
	offvec1 <- (1:noLev1 - 0.5)
	
	vecLev2 <- row.names(df)
	noLev2 <- length(vecLev2)
	offvec2 <- ((0:(noLev2-1))/(noLev2-1)*ms - ms/2)
	
	colvec <- c("red","green","blue","cyan","magenta")
	
	
	# Setup the axes
	plot(1:2,type="n", xlim=c(0,noLev1), ylim=ylim_in, axes=FALSE)
	axis(2,las=1)
	veclabs <- rep(" ",noLev1*2+1)
	for (i in 1:noLev1) veclabs[i*2] <- vecLev1[i]
	axis(1,at=(0:(noLev1*2))/2,labels=veclabs)
	
	# Set up a loop for the main offset
	for (i in 1:noLev1) {
		
		# Set up a loop for the secondary offset
		for (j in 1:noLev2) {
			
			# Plot the points and the confidence intervals
			xoff <- offvec1[i] + offvec2[j]
			n <- df[j,i]
			pest <- n/N
			lb <- binCI(c(N),c(n),0.975)
			ub <- binCI(c(N),c(n),0.025)
			points(c(xoff),c(pest),col=colvec[j],pch=22,bg=colvec[j])
			points(c(xoff,xoff),c(lb,ub),type="l",col=colvec[j])
			
		}
		
	}
	
	legend(0,max(ylim_in),legend=vecLev2,col=colvec,pt.bg=colvec,pch=22,bty="n")
	
	# Close the pdf device
	dev.off()
	
}

# A function to take a character string of a typical serological result and convert it to a raw titre
# Different lookup-tables need to be coded in different ways
# Initially, assumes a full dilution schedule starting at 1:10 
convertRawTitre <- function(charTitre,serotable="full_start_1_10") {
	
	# Define the table
	if (serotable=="full_start_1_10") {
		luTable <- data.frame(	raw=c("<1:20","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560"),
				lu=c(1,2,3,4,5,6,7,8,9)
		)
	}
	
	# Define return vector
	tmp <- luTable[luTable$raw==charTitre,"lu"]
	if (length(tmp) != 1) stop("Wrong value in convertRawTitre")
	
	# Return the correct value
	tmp
	
}

# A function to return the illustrative parameters
abp.params.baseline <- function() {
	c(		N 			= 7000000,
			seed 		= 100,
			trickle		= 1,
			R0 			= 1.8,
			Tg 			= 2.6,
			phi_1		= 1.0,
			phi_2		= 1.0,
			p_R 		= 1.0,
			p_H_base 	= 0.01,
			p_H_base_1 	= 1.0,
			p_H_base_2 	= 1.0,
			p_H_base_3 	= 1.0,
			p_I 		= 1.0,
			p_I_1 		= 1.0,
			p_I_2 		= 1.0,
			p_I_3 		= 1.0,
			gamma_v 	= 1/5.2,									# Gowardman
			gamma_h 	= 1/7.0,									# Feagan		
			t0 			= 0,
			p_1			= 0.2,			
			p_2			= 0.4,
			p_3			= 0.4,
			mixmatindex	= 2,										# CAE from polymod
			mixsense	= 1,
			fitsus		= 0,
			noyears		= 30,
			nodts		= 52*10,
			dur_seas	= 7*4,										# COuld be this low from shamen et al 09
			amp_seas	= 0.47,
			gamma_R		= 1/99999999,
			aging_on	= 0,
			mu_1		= 1/364/19, 
			mu_2		= 1/364/(65-19),
			mu_3		= 1/364/50									# Set to give very slow overall population growth
	)	
}

# A function to take three solutions to a simple model
# and produce a summary chart
abp.fig.illustrate <- function(file="~/Dropbox/tmp/fig1.pdf") {
	
	# Load up the default params
	ill_params <- abp.params.baseline()
	
	# Mass action and uniform susceptibility
	# Model 3, Table 3, Sherman at al 2011
	sol_A <- (SirModel3AgeClasses(
						pname=c("gamma_R","amp_seas","R0","Tg"),
						pvals=c(1/364/3.28,1.09/2.31,2.31,4.18),
						casemix=c(100,40,40),
						vp=ill_params
				))$sol
	
	# Model 3, Table 3, Sherman at al 2011
	sol_B <- (SirModel3AgeClasses(
						pname=c("gamma_R","amp_seas","R0","Tg"),
						pvals=c(1/364/7.77,1.195/2.495,2.495,2.59),
						casemix=c(100,40,40),
						vp=ill_params
				))$sol
	
	
	sol_null <- (SirModel3AgeClasses(
						pname=c("gamma_R","amp_seas","R0"),
						pvals=c(1/364/10,0,0.68),
						casemix=c(100,40,40),
						vp=ill_params
				))$sol
	
	# Open the file to put the chart into
	pdf(file,height=10/cm(1),width=10/cm(1))
	
	# standard figure parameter adjustments
	par(	cex=0.8,
			mgp=c(2.5,1,0),
			mai=(c(0,0,0,0)),
			cex.axis=0.9,
			fig=c(0.2,1.0,0.1,1.0),
			las=1)
	
	# Set up values for the axes
	yaxisvals <- 10^(0:5)
	yaxislabs <- c("1","10","100","1,000","10,000","100,000")
	
	# Plot the x axis
	plot(	1:2,
			#log="y",
			xlab="Time (years)",
			ylab="Incidence (per day)",
			type="n",
			log="y",
			xlim=c(0,ill_params["noyears"]),
			ylim=c(min(yaxisvals),max(yaxisvals)),
			axes=FALSE)
	axis(1)
	axis(2)
	#axis(2,at=yaxisvals,labels=yaxislabs)
	points(sol_A[,"time"]/364,(sol_A[,"dS"]+1),type="l",col="red")
	points(sol_B[,"time"]/364,(sol_B[,"dS"]+1),type="l",col="green")
	points(sol_null[,"time"]/364,(sol_null[,"dS"]+1),type="l",col="black")
	dev.off()
	
}

abp.annualmin <- function(ts,noyears,yl=360,incvar="dS",defaultval=1e100) {
	
	# This function takes output from a standard SR model ts and 
	# calculates the minimum annual incidence for the last noyears
	# of the sample
	# assumes uniform timesteps
	# yl is year length
	# incvar is the variable name for incidence
	
	# browser()
	
	value 	<- vector(mode="numeric",length=noyears)
	value[]	<- defaultval
	
	curind 	<- dim(ts)[1]
	curyear	<- noyears
	curmin	<- 1e100
	dt 		<- ts[curind,"time"] - ts[curind-1,"time"]
	
	nex_yb <- (ts[curind,"time"] %/% yl) * yl
	
	
	while (curyear >= 1 && curind >= 1) {
		if (ts[curind,incvar] < curmin) curmin <- ts[curind,incvar]
		if (ts[curind,"time"] < nex_yb) {
			value[curyear] <- curmin
			curyear <- curyear - 1
			nex_yb <- (ts[curind,"time"] %/% yl) * yl
			curmin <- 1e100
		}
		curind <- curind - 1
	}
	
	value
	
	# End definition of annualmin	
} 


# A funciton for the AB persistence project that calculates the maximum minimum annual value of the last few years of model solutions
abp.hyper <- function(	nohypersamples=1,
								hyperparams = data.frame(params=c("gamma_R","dur_seas"),log=c(0,0),lb=c(1/20,7),ub=c(1/1,7*26)),
								mintakenof = 5,
								noyears = 30,
								lfname = "~/Dropbox/tmp/abp_hyper_log.txt"
								) {	
									
	# Load up default parameters
	ill_params <- abp.params.baseline()								
									
	# Set up hypervector parameter search
	nohyperparams <- dim(hyperparams)[1]
	
	# Set up the hypersquare
	hypersquare <- matrix(0,ncol=nohyperparams,nrow=nohypersamples)
	for (i in 1:nohyperparams) {
		hypersquare[,i] <- hyper_vector(nohypersamples,hyperparams$lb[i],hyperparams$ub[i],log=hyperparams$log[i])
	}
			
	# This is the 
	cat("Parameter sweep started on",date(),"on",(Sys.info())["nodename"],"\n",file=lfname,append=FALSE)
	
	matEachAnnMin <- matrix(-1,nrow=nohypersamples,ncol=mintakenof)
	for (i in 1:nohypersamples) {
		sol_tmp <- (SirModel3AgeClasses(
							pname=as.character(hyperparams$params),
							pvals=hypersquare[i,],
							casemix=c(100,40,40),
							vp=ill_params
					))$sol	
		minvec <- abp.annualmin(sol_tmp,mintakenof)
		matEachAnnMin[i,] <- minvec[]
		cat("Completed",i,"of",nohypersamples,"\n",file=lfname,append=TRUE)
	}
	
	list(matmin=matEachAnnMin,hypersquare=hypersquare)	
	
}

abp.test.run <- function() {
	dfParams = data.frame(params=c("gamma_R","dur_seas"),log=c(0,0),lb=c(1/20,7),ub=c(1/1,7*26))
	system.time(tmp <- abp.hyper(nohypersamples=5,dfParams))
	
	matmin <- tmp$matmin
	hypersquare <- tmp$hypersquare
	nosamples <- dim(matmin)[1]
	browser()
	
}
