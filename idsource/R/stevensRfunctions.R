# Copyright Steven Riley (sr@stevenriley.net)

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

cat( 	"Reading in stevensRfunctions.R","\n",
		"The latest version of this file can be loaded into R with\n",
		"source(\"http://idsource.googlecode.com/svn/trunk/R/stevensRfunctions.R\") or","\n",
		"source(\"http://tinyurl.com/5t7gwnv\")","\n",
		"Revision X of this file can be loaded into R with","\n",
		"source(\"http://idsource.googlecode.com/svn-history/rX/trunk/R/stevensRfunctions.R\")","\n",
		sep="")

# List of three letter project prefixes
# abp - Comparison of the persistance of inflenzas A and B
# apa - Age and peak attack rate

# Translates a value from a unit scale to either a log or a linear scale
# x is the number on the unit scale
# min is the minimum value of the non-unit scale
# max is the maximum value of the non-unit scale
# logbase is the base of the log scale if used
# log is a boolean with TRUE for a log scale and FALSE for a linear scale
# Returns the value on the non-unit scale
SR_from_unit <- function(x,min=1,max=100,logbase=10,logflag=FALSE) {
	if (logflag) {
		rtn <- min*logbase^(x*(log(max,logbase)-log(min,logbase)))
	} else {
		rtn <- min + (max-min)*x
	}
	rtn
}

# Translates a number from a log or linear scale to a unit scale
# x is the number on the non-unit scale
# min is the assumed minimum value of the non-unit scale
# max is the assumed maximum value of the non-unit scale
# logbase is the base of the log scale if used
# log is a boolean for whether the scale is log or linear
# Returns the value on the unit scale
SR_to_unit <- function(y,min=1,max=100,logbase=10,logflag=FALSE) {
	if (logflag) {
		rtn <- (log(y,logbase)-log(min,logbase))/(log(max,logbase)-log(min,logbase))
	} else {
		rtn <- (y-min)/(max-min) 
	}
	rtn
}

binCINew <- function(N,n,P,min=0.000001) {
	rtn <- NA
	if (N>0) {
		if (N==n) n <- n-1
		if (n==0 && P > 0.5) rtn <- min
		else rtn <- (uniroot(srg.ci.binom.optim,c(0,1),P=P,N=N,n=n))$root
	}
	rtn
}

srg.bin.ci <- function(N,n,P,min=0.000001) {
	rtn <- NA
	if (N>0) {
		if (N==n) n <- n-1
		if (n==0 && P > 0.5) rtn <- min
		else rtn <- (uniroot(srg.ci.binom.optim,c(0,1),P=P,N=N,n=n))$root
	}
	rtn
}


# binCINew(5,0,0.025,min=0.01)
# pbinom(100,100,0.9,lower.tail=TRUE)

SR_updatePVector <- function(values,names,pvec) {                                                                        
	no_updates <- length(values)
	for (i in 1:no_updates) pvec[names[i]]=values[i]
	pvec <- jeUpdateAuxParams(pvec)
	pvec
}

# mount -t smbfs //username:userpass@myserver/PUBLIC /smb/public

srg.chart.pos <- function(	xindex,yindex,xn,yn,xlm=0,xrm=0,
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

sr2DHist <- function(vec_x,vec_y,min_x=min(vec_x)-1,min_y=min(vec_y)-1,max_x=max(vec_x)+1,max_y=max(vec_y)+1,
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
	
	list(z=rtn,xb=xbounds,yb=ybounds,xm=xmids,ym=ymids)
	
}

ExcerptLStoASCII <- function(lsfile,asciifile,de,dw,ds,dn,lsspd=120,lsmn=85,lsmw=34) {
	# Takes landscan file (lsfile) and produces and ascii grid file (asciifile) 
	# Returns nothing
	# coord agurments are degree east (exe), degree west (exw) ...
	# require("rgdal")
	require("maptools")
	require("raster")
	tmp <- raster(lsfile)
	ext <- extent(dw,de,ds,dn)
	rtn <- crop(tmp,ext)
	write(rtn,"~/Dropbox/tmp/asciiout.ascii","ascii")
	detail <- readGDAL(	
			lsfile,
			region.dim=c((dn-ds)*lsspd,(dw-de)*lsspd),
			offset=c((lsmn-dn)*lsspd,(de-lsmw)*lsspd))
	# writeAsciiGrid(detail,asciifile)
	write.asciigrid(detail,asciifile)
}

srg.popdens.excerpt <- function(lsfile="~/tmp/lspop2008.flt",x1=-0.4,x2=-0.2,y1=51.8,y2=52) {
	require(rgdal)
	require(raster)
	world <- raster(lsfile)
	ext <- extent(x1,x2,y1,y2)
	rtn <- crop(world,ext)
	rtn
}

# Function to make an R object out of a subset of an acii grid
# Runs from main home
# Need to use this to make a very small subset around some of the villages
# And then test the neighbourhood calculation
# source("stevensRfunctions.R")
# x <- srg.read.subset.asciigrid()
# pdf(file="~/Dropbox/tmp/Rout.pdf")
# image(x)
# dev.off()
# quartz()
srg.read.subset.asciigrid <- function(	
		filename="~/Dropbox/svneclipse/fluscape/data/landscan/prd_sim_ascii_n_44384655.txt",
		south=22.2,
		west=114.1,
		north=22.5,
		east=114.3,
		writeout=TRUE, 
		outputfile="~/Dropbox/tmp/srg.readmask.asciigrid.output.R") {
	
	require(sp)
	
	large <- read.asciigrid(filename,as.image=TRUE)
	nox <- length(large$x)
	noy <- length(large$y)
	dimz <- dim(large$z)
	if (nox < 2) stop("Problem in srg.read.subset.asciigrid")
	if (noy < 2) stop("Problem in srg.read.subset.asciigrid")
	dx <- large$x[2] - large$x[1]
	dy <- large$y[2] - large$y[1]
	minx <- large$x[1]
	maxx <- large$x[nox]
	miny <- large$y[1]
	maxy <- large$y[noy]
	
	# Calculate index bounds
	min_ix <- floor((west - minx) / dx) 
	max_ix <- ceiling((east - minx) / dx) 
	min_iy <- floor((south - miny) / dy) 
	max_iy <- ceiling((north - miny) / dy) 
	
	# Check for index errors
	if (min_ix > max_ix || min_iy > max_iy || min_ix < 1 || max_ix > nox || min_iy < 1 || max_iy > noy) stop("Index problem in srg.read.subset.asciigrid")
	
	# Make output object
	x_rtn <- large$x[min_ix:max_ix]
	y_rtn <- large$y[min_iy:max_iy]
	nox_rtn <- length(x_rtn)
	noy_rtn <- length(y_rtn)
	z_rtn <- matrix(nrow=nox_rtn,ncol=noy_rtn)
	z_rtn[] <- large$z[min_ix:max_ix,min_iy:max_iy]
	
	# Take the subset
	image_saved_from_srg_readmask_asciigrid <- list(x=x_rtn,y=y_rtn,z=z_rtn)
	
	# Write file and return image if required
	if (writeout) save(image_saved_from_srg_readmask_asciigrid,file=outputfile)
	image_saved_from_srg_readmask_asciigrid
	
}

# Function to make a population density image from grump data that is 
# 1 degree larger on all sides than the extreme values in the study population
# fsc.pop.image.as.R.object.from.ascii()
fsc.pop.image.as.R.object.from.ascii <- function(	asciifile="~/Dropbox/tmp/asup00g.asc",
		locationfile="~/Dropbox/svneclipse/fluscape/data/Location_Density_bad.csv",
		outputfile="~/Dropbox/tmp/fsc_op_image_as_R_object_from_grump.R",
		margindecdeg=0.5) {
	
	# Load up the loc density files
	loc <- read.csv(locationfile)
	south_in <- min(loc$LOC_Lat,na.rm = TRUE) - margindecdeg
	north_in <- max(loc$LOC_Lat,na.rm = TRUE) + margindecdeg
	west_in <- min(loc$LOC_Long,na.rm = TRUE) - margindecdeg
	east_in <- max(loc$LOC_Long,na.rm = TRUE) + margindecdeg
	
	rtn <- srg.read.subset.asciigrid(
			filename=asciifile,
			south=south_in,
			west=west_in,
			north=north_in,
			east=east_in,
			writeout=TRUE,
			outputfile=outputfile) 
	
}

# Function to make a population density image from landscan data that is 
# 1 degree larger on all sides than the extreme values in the study population
fsc.pop.image.as.asciigrid.from.landscan <- function(	
		landscanfile="/home/sriley/Dropbox/dataplain/gis/landscan/asia03/asia03/w001001.adf",
		locationfile="/home/sriley/Dropbox/svneclipse/fluscape/data/Location_Density_bad.csv",
		outputfile="/home/sriley/Dropbox/tmp/fsc_op_asciigrid_from_landscan.R",
		margindecdeg=0.5) {
	
	# Load up the loc density files
	loc <- read.csv(locationfile)
	south <- min(loc$LOC_Lat,na.rm = TRUE) - margindecdeg
	north <- max(loc$LOC_Lat,na.rm = TRUE) + margindecdeg
	west <- min(loc$LOC_Long,na.rm = TRUE) - margindecdeg
	east <- max(loc$LOC_Long,na.rm = TRUE) + margindecdeg
	
	ExcerptLStoASCII(landscanfile,outputfile,east,west,south,north,lsspd=120,lsmn=85,lsmw=34)
	
}

# Make 3 different calls to the lanscan data to extract areas A, B and C from the data 
fsc.pop.image.as.asciigrid.from.landscan.gravity <- function(
		margin=0,
		landscanfile="/Users/sriley/Dropbox/dataplain/gis/landscan/asia03/asia03/w001001.adf",
		outputstem=paste("~/Dropbox/tmp/fsc_gravity_area_",margin,".R",sep="")) {
	
	# Load up the loc density files
	south <- 23.03 + margin
	north <- 23.6 + margin
	west <- 113.2 + margin
	east <- 113.85 + margin
	
	ExcerptLStoASCII(landscanfile,outputfile,east,west,south,north,lsspd=120,lsmn=85,lsmw=34)
	
}


# Load landscan
# x <- fsc.load.pop.image()
# image(x,xlim=c(113,113.1),lty="blank")
# quartz()
# image(x)
fsc.load.pop.image <- function(file="~/Dropbox/svneclipse/fluscape/data/fsc_landscan_subset_small.R") {
	load(file)
	image_saved_from_srg_readmask_asciigrid
}

fsc.make.good.cols.breaks <- function(ncols=10,max=20^5) {
	scale <- (0:ncols)/ncols*log10(max)
	bks <- c(-9999,10^scale)
	cls <- c("black",heat.colors(ncols))
	list(breaks=bks,colours=cls)
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
				ycoord = ny - (yReal-ymin)/ygap+1
				# ycoord = (yReal-ymin)/ygap+1
				if (is.na(rtnval$z[xcoord,ycoord])) rtnval$z[xcoord,ycoord]=1
				else rtnval$z[xcoord,ycoord]=rtnval$z[xcoord,ycoord]+1
			}
		}
	}
	noRuns = er-sr+1
	rtnval$z=rtnval$z/noRuns
	rtnval	
}

srg.decLongLatDist <- function(x1,y1,x2,y2,translate=FALSE) {
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
		if (rtn < -1e100 || rtn > 1e100) stop("error akjshdkajhdakj")
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

srg.plotstack <- function(sol,ac,col,var="Iv",mult=1,thresh=1,const=0) {
	
	ts <- dim(sol)[1]
	x <- c(sol[,"time"],rev(sol[,"time"]))
	
	if (ac==1) y <- c((sol[,paste(var,"1",sep="")]),rep(0,ts))
	else if (ac==2) y <- c((sol[,paste(var,"1",sep="")] + sol[,paste(var,"2",sep="")]),rev((sol[,paste(var,"1",sep="")])))
	else if (ac==3) y <- c((sol[,paste(var,"1",sep="")] + sol[,paste(var,"2",sep="")] + sol[,paste(var,"3",sep="")]),rev((sol[,paste(var,"1",sep="")]+sol[,paste(var,"2",sep="")])))
	else stop("Age class not specified.")
	
	polygon(x,y*mult+const,col=col,border=NA)
	
}

srg.plotstack.dash <- function(sol,col,var="Iv",varbelow=NULL,mult=1,thresh=1,const=0,timevar="time") {
	
	ts <- dim(sol)[1]
	x <- c(sol[,timevar],rev(sol[,timevar]))
	
	lowery <- rep(const,ts)
	if (!is.null(varbelow)) for (varb in varbelow) {
			lowery <- lowery + mult*sol[,varb]
		}
	uppery <- lowery + mult*sol[,var]
	
	y <- c(lowery,rev(uppery))
	
	polygon(x,y,col=col,border=NA)
	
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


srg.param.propose.update <- function(ptab,fmask=1:(dim(ptab)[1])) {
	
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
	if (index < 0 || index > dim(ptab)[1]) {
		stop("index must be less than or equal to number of parameters")
	}
	
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
		# prop_x <- srg.param.propose.update(x)
		
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
			if (age > 110 || age < 0) stop("problem with age group assignment")
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

srg.hyper.vector <- function(n,lb,ub,log=FALSE) {
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
			stop("problem with age group assignemnt 98237492")
		}
		
		if (rtn$age[i] > 2 && rtn$age[i] < 19) rtn$ag2[i] <- "1_3t18"
		else if (rtn$age[i] < 49) rtn$ag2[i] <- "2_19t48"
		else if (rtn$age[i] < 65) rtn$ag2[i] <- "3_49t64"
		else if (rtn$age[i] < 110) rtn$ag2[i] <- "4_65plus"
		else {
			stop("problem with age group assignemnt 298374283")
		}
		
		if (rtn$age[i] > 2 && rtn$age[i] < 20) rtn$ag20[i] <- "1_3t19"
		else if (rtn$age[i] < 40) rtn$ag20[i] <- "2_20t39"
		else if (rtn$age[i] < 60) rtn$ag20[i] <- "3_40t59"
		else if (rtn$age[i] < 110) rtn$ag20[i] <- "4_60plus"
		else {
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

annualmin <- function(ts,noyears,yl=360,incvar="dS",defaultval=1e100) {
	
	# This function takes output from a standard SR model ts and 
	# calculates the minimum annual incidence for the last noyears
	# of the sample
	# assumes uniform timesteps
	# yl is year length
	# incvar is the variable name for incidence
	
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
	
}

# Routine to read the used in the CID paper and make equivalent fields
readJoeEFlu <- function(fnFile) {
	
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
	# Not sure that this is working at all 
	
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

srg.summarise.glm <- function(
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
 vecDEX <- vector(mode="numeric",length=nomods)
	
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
		vecDEX[i] <- (1-mod$deviance/mod$null.deviance)*100
  
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
			strTable <- paste(strTable,",AIC,DEX\n",sep="")
			
			# Put in the second header row
			strTable <- paste(strTable,"Model",sep="")
			for (i in 1:noCoeffs) strTable <- paste(strTable,",PE,CI",sep="")
			strTable <- paste(strTable,",AIC,DEX\n",sep="")
			
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
    curDEX <- round((1-mod$deviance/mod$null.deviance)*100,digits=1)
				strTable <- paste(strTable,",",curAIC,",",curDEX,"\n",sep="")
				
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
			if (nomods>1) for (i in 1:nomods) {
     strTable <- paste(strTable,",Estimate,(95% CI)",sep="")
   }
			strTable <- paste(strTable,"\n",sep="")
			
			# Output individual coefficient lines, starting with coefficient loop
			for (i in 1:noCoeffs) {
				
				# Pull the name of the current coefficient
				curname <- allCoeffs[i]
				
				# Put in the name of the coefficient
				strTable <- paste(strTable,curname,sep="")
				
				# Cycle through the tables looking at the different models
				for (j in 1:nomods) {
					
					# Put in the point estimates and confidence intervals for each 
     # parameter / model combination
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

   # Write the row name for the DEXs	
   strTable <- paste(strTable,"DEX",sep="")
   
   # Start the for loop for the DEX
   for (i in 1:nomods) {
     
     # Get the current model
     mod <- lstModels[[i]]
     
     # Format the AIC for the current model
     curDEX <- round((1-mod$deviance/mod$null.deviance)*100,digits=1)
     
     # Write the value and the space
     strTable <- paste(strTable,",",curDEX,",",sep="")
   }
 
   # Return at the end of the DEX line
   strTable <- paste(strTable,"\n",sep="")
   
			# End else statement for transpose of 
		}
		
		# Write the string to the selected file
		cat(strTable,file=file)
		
	}
	
	# Return 
	data.frame(pe=matPointEst,lb=matLB,ub=matUB,aic=vecAIC,dex=vecDEX)
	
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

srg.make.binci.table <- function(numerator,denominator,rounding=3) {
	nosamps <- length(numerator)
	ub <- srg.ci.binom(denominator,numerator,P=0.025)
	lb <- srg.ci.binom(denominator,numerator,P=0.975)
	pe <- numerator/denominator
	round(data.frame(n=numerator,N=denominator,pe=pe,lb=lb,ub=ub),rounding)
}

# Computes Binomial 95% Confidence Interval
srg.ci.binom <- function(vecN,vecn,P,min=0) {
    binom.optim <- function(p,P,N,n) {
    guess <- pbinom(n,N,p,lower.tail=TRUE)
    rtn <- P-guess
    rtn
  }
  
	noObs <- length(vecN)
	rtn <- array(NA,c(noObs))
	for (i in 1:noObs) {
		if (vecN[i] > 0) {
			sol <- uniroot(binom.optim,c(0,1),P=P,N=vecN[i],n=vecn[i])
			if (sol$root > min) rtn[i] <- sol$root
			else rtn[i] <- min
		} else rtn[i] <- min
	}
	rtn
}

# An auxillary function used to make binomial confidence bounds
# returns the difference between a "guess" at the cumulative binomial probability for a number adn the real value

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
			ub <- srg.ci.binom(c(N),c(n),0.025)
			points(c(xoff),c(pest),col=colvec[j],pch=22,bg=colvec[j])
			lb <- srg.ci.binom(c(N),c(n),0.975)
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
srg.convertRawTitre <- function(charTitre,serotable="full_start_1_10") {
	
	# Define the table
	if (serotable=="full_start_1_10") {
		luTable <- data.frame(	raw=c("<1:20","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560"),
				lu=c(1,2,3,4,5,6,7,8,9)
		)
	} else if (serotable=="mahon_ben_1_10") {
		luTable <- data.frame(	raw=c("5","10","40","80","160","320","640","1280","2560","5120","10240"),
				lu=c(1,2,3,4,5,6,7,8,9,10,11))
	}
	
	# Define return vector
	tmp <- luTable[luTable$raw==charTitre,"lu"]
	if (length(tmp) != 1) stop("Wrong value in convertRawTitre")
	
	# Return the correct value
	tmp
	
}

FigPersistenceCompare <- function(solutionA, solutionB, solutionN,
		file="fig1.pdf") {
	
	
	# Plot the illustartive scenarios
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
			xlim=c(0,ill_params["noyears"]),
			ylim=c(1,20),
			axes=FALSE)
	axis(1)
	axis(2)
	#axis(2,at=yaxisvals,labels=yaxislabs)
	points(solutionA[,"time"]/364,(solutionA[,"dS"]+1),type="l",col="red")
	points(solutionB[,"time"]/364,(solutionB[,"dS"]+1),type="l",col="green")
	points(solutionN[,"time"]/364,(solutionN[,"dS"]+1),type="l",col="black")
	dev.off()
	
}

readSpacialInc <- function(fldir,stem,suffix=".csv") {
	
	# Cycle thorugh the fldirectory once to get a list of continuous filenames, index from 0
	no_files <- 0
	while (length(dir(fldir,pattern=paste(stem,no_files,suffix,sep="")))>0) no_files <- no_files+1
	if (no_files == 0) stop("No files found in readSpacialInc")
	
	# Open the first file to get the dimensions
	tmptab <- read.csv(paste(fldir,stem,0,suffix,sep=""),header=FALSE)
	dims <- dim(tmptab)
	nots <- dims[1]
	nopops <- dims[2]-3
	nox <- tmptab[1,2]
	noy <- tmptab[1,3]
	if (nopops != nox*noy) stop("Declared dimensions and nox, noy not consistent in readSpacialInc")
	
	# Define the return array
	rtn <- array(-1,c(no_files,nots,nox,noy))
	
	# Read in the files to the array
	current_file <- 0
	while (current_file < no_files) {
		
		# Read in one file
		tmptab <- read.csv(paste(fldir,stem,current_file,suffix,sep=""),header=FALSE)
		for (i in 1:nots) {
			for (j in 1:nox) {
				for (k in 1:noy)
					rtn[current_file+1,i,j,k] <- tmptab[i,3+(j-1)*noy+k]
			}
		}
		
		# Increment loop
		current_file <- current_file + 1
		
	}
	
	# Give back the main file 4D array
	rtn
	
}

# Function to generate a single set of x,y for incidence from a simulation realization
genSingInc <- function(index,x,dt) {
	
	dimx <- dim(x)
	if (length(dimx) != 4) stop("Incorrect dimensions in genSingleIncidence")
	notimepoints <- dimx[2] - 1
	noreals <- dimx[1]
	if (index < 0 || index > noreals) stop("Index of realization out of bounds in genSingleIncidence")
	
	rtnx <- vector(mode="numeric",length=notimepoints)
	rtny <- vector(mode="numeric",length=notimepoints)
	
	for (i in 1:notimepoints) {
		rtny[i] = sum(x[index,i+1,,]) - sum(x[index,i,,])
		rtnx[i] = dt/2 + (i-1)*dt
	}
	
	list(x=rtnx,y=rtny)
	
}

# Model to generate multiple realizations of a ver
gmm.stochSIRS <- function(	nrs=10,dt=0.25,dtm=1,tmax=100,
		N=1200000,R0=1,D_I=2.6,D_R=360*4,beta_s=0.000001,I0=0) {
	
	# Set some aux variables
	epsilon <- 1e-20
	
	# Define beta in terms of R0
	beta <- R0/D_I
	
	# Set up the return variables using the same conventions as the C++ program
	noSimTs <- round(tmax/dt) + 1
	noMeasTs <- round(tmax/dtm) + 1
	rtn <- array(dim=c(nrs,noMeasTs))
	
	# Start the realization loop
	for (i in 1:nrs) {
		
		# Initialize the loop variables
		sus <- N-I0
		inf <- I0
		rec <- 0
		cuminc <- 0
		nextMeasureIndex <- 1
		
		# Start the time variable
		for (j in 1:noSimTs) {
			
			# Set up the time step variables
			currentTime <- (j-1)*dt
			
			# Record incidence if passed a measure time
			if ((currentTime - (nextMeasureIndex-1)*dtm + epsilon) > 0) {
				rtn[i,nextMeasureIndex] <- cuminc
				nextMeasureIndex <- nextMeasureIndex + 1
			}
			
			# Calculate the probabilities of the different events
			lambda <- beta * inf / N + beta_s
			prob_inf <- 1-exp(-lambda*dt)
			prob_rec <- 1-exp(-dt/D_I)
			prob_bsa <- 1-exp(-dt/D_R)
			
			# Draw the appropriate random variables
			no_inf <- rbinom(1,sus,prob_inf)
			no_rec <- rbinom(1,inf,prob_rec)
			no_bsa <- rbinom(1,rec,prob_bsa)
			
			# Update the state variables
			sus <- sus + no_bsa - no_inf
			inf <- inf + no_inf - no_rec
			rec <- rec + no_rec - no_bsa
			cuminc <- cuminc + no_inf
			
		}
		
	}
	
	rtn
	
}

srg.pair.generation <- function(vecNames = c("a","b","c","d","e")) {
	n <- length(vecNames)
	sizertn <- n*(n-1)/2
	rtn <- data.frame(first=rep(NA,sizertn),second=rep(NA,sizertn))
	currentrow <- 1
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			rtn[currentrow,"first"] <- as.character(vecNames[i])
			rtn[currentrow,"second"] <- as.character(vecNames[j])
			currentrow <- currentrow + 1
		}
	}
	rtn
}

specific.regression.models <- function(dat) {
	mod1 <- glm( formula = dat$com_result ~ dat$age+factor(dat$child_p),family = binomial(logit))
	mod2 <- glm( formula = dat$com_result ~ dat$age+factor(dat$child_p) + dat$n_location,family = binomial(logit))
	mod3 <- glm( formula = dat$com_result ~ dat$age+factor(dat$child_p) + dat$n_location + dat$A_4,family = binomial(logit))
	list(mod1,mod2,mod3)
}

plot.explanitory.variables.hk.contact <- function(df) {
	pos <- df[df$com_result == 1,]
	neg <- df[df$com_result == 0,]
	cpres <- df[df$child_p == 1,]
	cprespos <- cpres[cpres$com_result == 1,]
	cpresneg <- cpres[cpres$com_result == 0,]
	cnpres <- df[df$child_p == 2,]
	cnprespos <- cnpres[cnpres$com_result == 1,]
	cnpresneg <- cnpres[cnpres$com_result == 0,]
	plot(neg$age,jitter(neg$n_location),col="black",xlim=c(0,80),pch=46,cex=5)
	points(pos$age,jitter(pos$n_location),col="red",pch=46,cex=5)
}

load.hk.cohort.pig <- function(file="~/Dropbox/projects/influenza/hk_cohort_swine_herd/age_seralevels_sr.csv") {
	rtn <- read.csv(file)
	rtn
}

transform.hk.cohort.pig.single <- function(x) {
	rtn <- x
	if (rtn < 20) rtn <- 10
	rtn <- log2(rtn/10)
	rtn
}

transform.hk.cohort.pig.all <- function(df) {
	allnames <- names(df)
	namestodo <- setdiff(allnames,c("age"))
	rtn <- df
	norows <- dim(rtn)[1]
	for (var in namestodo) {
		for (i in 1:norows) {
			rtn[i,var] <- transform.hk.cohort.pig.single(rtn[i,var])
		}
	}
	rtn
}

# A function to take a data frame of individual level pre and post titres for different strains
# to look at the possible association with age
run.poisson.regressions <- 	function(	
		df,
		strstem=c("CA42009","X4167","X1110","X1304","Xns29","X1559","X201","g112","X400599")
) {
	
	# Set up the output table
	outputrows <- length(strstem)
	rtn <- 	data.frame(	
			point.pre=rep(-1,outputrows),lb.pre=rep(-1,outputrows),ub.pre=rep(-1,outputrows),
			point.post=rep(-1,outputrows),lb.post=rep(-1,outputrows),ub.post=rep(-1,outputrows)
	)
	
	# Run the different models
	for (i in 1:outputrows) {
		curname <- strstem[i]
		mod.pre <- glm(as.formula(paste(curname,".pre ~ age",sep="")), family=poisson, dat = df)
		mod.post <- glm(as.formula(paste(curname,".post ~ age",sep="")), family=poisson, dat = df)
		cis.pre <- confint.default(mod.pre)
		cis.post <- confint.default(mod.post)
		
		# Extract the point estimate and confidence intervals for the parameters
		rtn[curname,"point.pre"] 	<- mod.pre$coefficients["age"]  
		rtn[curname,"lb.pre"] 		<- cis.pre["age",1]  
		rtn[curname,"ub.pre"] 		<- cis.pre["age",2]  
		rtn[curname,"point.post"] 	<- mod.post$coefficients["age"]  
		rtn[curname,"lb.post"] 		<- cis.post["age",1]  
		rtn[curname,"ub.post"] 		<- cis.post["age",2]  
		
		# Close the for loop	
	}
	
	# Return the table of model results
	rtn
	
# Close the function block 
} 

# A re-code of the simple age peak model
# Set up fast parameters 
apa.setup.params.spq <- function(base=NULL,pnames=NULL,pvals=NULL) {
	
	# Initialize to defaul values if no base option given
	if (identical(base,NULL)) {
		rtn <- vector(mode="numeric",length=0)
		rtn["R0"] 		<- 1.4
		rtn["Tg"] 		<- 2.6
		rtn["De"]		<- 0.001
		rtn["nc"] 		<- 0.2
		rtn["I0c"] 		<- 100/7000000*rtn["nc"]
		rtn["I0a"] 		<- 100/7000000*(1-rtn["nc"])
		rtn["m_aa"]		<- 1
		rtn["m_ac"]		<- 1
		rtn["m_ca"]		<- 1
		rtn["phi"]		<- 1
		rtn["delta"]	<- 1
		rtn["or"]		<- -1
		# Else start from the base values
	} else {
		rtn <- base
	}
	
	# Check for suggested parameter name changes and change if necessary
	if (!identical(pnames,NULL)) {
		noPs <- length(pnames)
		if (noPs!=length(pvals)) stop("Problem with lengths in setup.params.spq")
		rtnnames <- names(rtn)
		for (i in 1:noPs) {
			if (!(pnames[i] %in% rtnnames)) stop("Subbing a parameter that isn't present in setup.params.spq")
			rtn[pnames[i]] <- pvals[i]
		}
	}
	
	# Correct the duration of infectiousness for a non-zero exposure period
	if (rtn["De"] > rtn["Tg"]) stop("Must have De < Tg in apa.setup.params.spq")
	rtn["Di"] <- rtn["Tg"] - rtn["De"]
	
	
	# Define the next generation matrix
	ngm_dev_beta <- rtn["Di"] * t(matrix(c(	1*rtn["nc"],								(rtn["m_ac"])^rtn["delta"]*rtn["nc"],
							(rtn["m_ca"])^rtn["delta"]*rtn["phi"]*(1-rtn["nc"]),		(rtn["m_aa"])^rtn["delta"]*rtn["phi"]*(1-rtn["nc"])), nrow=2,ncol=2))			
	esol 		<- eigen(ngm_dev_beta)		
	R0_dev_beta <- max(esol$values)
	rtn["beta"] <- rtn["R0"] / R0_dev_beta 
	ev  <- Re(esol$vector[,1])
	if (ev[1] < 0) ev <- ev * -1
	ev <- ev/sum(ev)
	rtn["or"] <- ev[1]/rtn["nc"] 
	
	# Return result
	rtn
	
}

apa.setup.state.vector <- function(p) {
	rtn <- vector(mode="numeric",length=0)
	rtn["S"] <- 1 - p["I0c"] - p["I0a"]
	rtn["I"] <- p["I0c"] + p["I0a"]
	rtn["cumInc"] <- 0
	rtn
}

apa.setup.state.vector.2state <- function(p) {
	rtn <- vector(mode="numeric",length=0)
	rtn["Sc"] <- p["nc"] - p["I0c"]
	rtn["Ic"] <- p["I0c"]
	rtn["Sa"] <- (1-p["nc"])-p["I0a"]
	rtn["Ia"] <- p["I0a"]
	rtn["cumIncA"] <- 0
	rtn["cumIncC"] <- 0
	rtn
}

apa.setup.states.2state.SEIR <- function(p) {
	rtn <- vector(mode="numeric",length=0)
	rtn["Sc"] <- p["nc"] - p["I0c"]
	rtn["Ic"] <- p["I0c"]
	rtn["Ec"] <- 0
	rtn["Ea"] <- 0
	rtn["Sa"] <- (1-p["nc"])-p["I0a"]
	rtn["Ia"] <- p["I0a"]
	rtn["cumIncA"] <- 0
	rtn["cumIncC"] <- 0
	rtn
}


# Make the equations
apa.simple.peak.attack.2state <- function(t,y,p){
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	infc 			<- 	p["beta"] * y["Sc"] * ( y["Ic"] + y["Ia"])
	infa 			<- 	p["beta"] * y["Sa"] * ( y["Ic"] + y["Ia"])
	rtn["Sc"] 		<- - infc
	rtn["Ic"] 		<- + infc - 1 / p["Tg"] * y["Ic"]
	rtn["Sa"] 		<- - infa
	rtn["Ia"] 		<- + infa - 1 / p["Tg"] * y["Ia"]
	rtn["cumIncA"] 	<- infa
	rtn["cumIncC"] 	<- infc
	list(rtn)
}

# Make the equations
apa.spa.2state.mix <- function(t,y,p){
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	infc 			<- 	p["beta"] * y["Sc"] * (     1     * y["Ic"] + p["m_ac"] * y["Ia"])
	infa 			<- 	p["beta"] * y["Sa"] * ( p["m_ca"] * y["Ic"] + p["m_aa"] * y["Ia"])
	rtn["Sc"] 		<- - infc
	rtn["Ic"] 		<- + infc - 1 / p["Tg"] * y["Ic"]
	rtn["Sa"] 		<- - infa
	rtn["Ia"] 		<- + infa - 1 / p["Tg"] * y["Ia"]
	rtn["cumIncA"] 	<- infa
	rtn["cumIncC"] 	<- infc
	list(rtn)
}

# Make the equations
apa.spa.2state.sus  <- function(t,y,p){
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	infc 			<- 	p["beta"] * y["Sc"] * (     1     * y["Ic"] + (p["m_ac"])^p["delta"] * y["Ia"])
	infa 			<- 	p["beta"] * p["phi"] * y["Sa"] * ( (p["m_ca"])^p["delta"] * y["Ic"] + (p["m_aa"])^p["delta"] * y["Ia"])
	rtn["Sc"] 		<- - infc
	rtn["Ic"] 		<- + infc - 1 / p["Tg"] * y["Ic"]
	rtn["Sa"] 		<- - infa
	rtn["Ia"] 		<- + infa - 1 / p["Tg"] * y["Ia"]
	rtn["cumIncA"] 	<- infa
	rtn["cumIncC"] 	<- infc
	list(rtn)
}

# Make the equations
apa.spa.2state.SEIR  <- function(t,y,p){
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	infc 			<- 	p["beta"] * y["Sc"] * (     1     * y["Ic"] + (p["m_ac"])^p["delta"] * y["Ia"])
	infa 			<- 	p["beta"] * p["phi"] * y["Sa"] * ( (p["m_ca"])^p["delta"] * y["Ic"] + (p["m_aa"])^p["delta"] * y["Ia"])
	rtn["Sc"] 		<- - infc
	rtn["Ec"]		<- + infc - y["Ec"] / p["De"]
	rtn["Ic"] 		<- + y["Ec"] / p["De"] - 1 / p["Di"] * y["Ic"]
	rtn["Sa"] 		<- - infa
	rtn["Ea"]		<- + infa - y["Ea"] / p["De"]
	rtn["Ia"] 		<- + y["Ea"] / p["De"] - 1 / p["Di"] * y["Ia"]
	rtn["cumIncA"] 	<- infa
	rtn["cumIncC"] 	<- infc
	list(rtn)
}

# Make the equations
apa.spa.2state.sus  <- function(t,y,p){
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	infc 			<- 	p["beta"] * y["Sc"] * (     1     * y["Ic"] + (p["m_ac"])^p["delta"] * y["Ia"])
	infa 			<- 	p["beta"] * p["phi"] * y["Sa"] * ( (p["m_ca"])^p["delta"] * y["Ic"] + (p["m_aa"])^p["delta"] * y["Ia"])
	rtn["Sc"] 		<- - infc	
	rtn["Ic"] 		<- + infc - 1 / p["Tg"] * y["Ic"]
	rtn["Sa"] 		<- - infa
	rtn["Ia"] 		<- + infa - 1 / p["Tg"] * y["Ia"]
	rtn["cumIncA"] 	<- infa
	rtn["cumIncC"] 	<- infc
	list(rtn)
}

# Make the equations
apa.simple.peak.attack 	<- function(t,y,p) {
	rtn 			<- vector(mode="numeric",length=length(y))
	names(rtn) 		<- names(y)
	inf 			<- 	 p["beta"] * y["S"] * y["I"]
	rtn["S"] 		<- - inf
	rtn["I"] 		<- + inf - 1 / p["Tg"] * y["I"]
	rtn["cumInc"] 	<-   inf
	list(rtn)
}

# Extract the attack rate
apa.extract.CAR 	<- function(odeout,p) {
	lts <- dim(odeout)[1]
	c(p["nc"]-odeout[lts,"Sc"],1-p["nc"]-odeout[lts,"Sa"])
}

# Extract incidence
apa.extract.inc <- function(odeout,p) {
	nots <- dim(odeout)[1]
	rtn <- matrix(nrow = nots,ncol=4,dimnames=list(1:nots,c("time","children","adults","total")))
	rtn[1,] <- 0
	rtn[,"time"] <- odeout[,"time"] 
	rtn[2:nots,"children"] <- c(odeout[2:nots,"cumIncC"]-odeout[1:(nots-1),"cumIncC"])
	rtn[2:nots,"adults"] <- c(odeout[2:nots,"cumIncA"]-odeout[1:(nots-1),"cumIncA"])
	rtn[,"total"] <- rtn[,"children"]+rtn[,"adults"]
	rtn
}

# Run the code and make a simple plot
apa.main.simple.peak.attack <- function() {
	require("odesolve")
	p 	 <- apa.setup.params.spq()
	s3 	 <- apa.setup.state.vector.2state(p)
	sol3 <- lsoda(s3,0:180,apa.spa.2state.mix,p)
	inc <- apa.extract.inc(sol3,p)
	plot(sol3[,"time"],inc[,3],type="l")
}

# Generate latin hypercube samples of CAR, Peak attackrate and over representation
apa.lhs <- function(nosamps=1,vecpars=c("R0"),veclb=c(2),vecub=c(3),veclog=c(FALSE),constpnames=c(),constpvals=c()) {
	
	# Include required libraries
	require("odesolve")
	
	# Define the output matrix
	novars = length(vecpars)
	rtn <- matrix(ncol= novars + 3, nrow = nosamps,dimnames=list(1:nosamps,c(vecpars,"or","peak","CAR")))
	
	# Generate the parameter values
	for (i in 1:novars) {
		rtn[,i] <- srg.hyper.vector(nosamps,veclb[i],vecub[i],veclog[i])
	}
	
	# Solve each of the cases
	for (i in 1:nosamps) {
		
		# Run each of the parameter sets
		p 	 			<- apa.setup.params.spq(pnames=c(vecpars,constpnames),pvals=c(rtn[i,1:novars],constpvals))
		initc 			<- apa.setup.state.vector.2state(p)
		sol 			<- lsoda(initc,0:360,apa.spa.2state.sus,p)
		inc 			<- apa.extract.inc(sol,p)
		rtn[i,novars+1] <- p["or"]
		rtn[i,novars+2] <- max(inc[,"total"])
		rtn[i,novars+3] <- sum(apa.extract.CAR(sol,p))
		
	}
	
	# Return the complete matrix
	rtn
	
	# A test line
	# x <- apa.lhs(p,10,c("m_aa","m_ca","m_ac"),c(0.01,0.01,0.01),c(1,1,1),veclog=c(TRUE,TRUE,TRUE))
}

# Function to generate a simple legend strip to be included in a multipart figure
srg.simple.image.legend <- function(min=0,max=1,colpalette=heat.colors(10),logscale=FALSE,axlabs=NULL) {
	number=length(colpalette)
	imagedat <- matrix(nrow=number,ncol=2)
	xvals <- c(0,1,2)
	if (logscale) {
		min <- log10(min)
		max <- log10(max)
	}
	yvals <- (0:(number))/(number)*(max-min)+min
	for (i in 1:number) imagedat[i,] <- (yvals[i+1]+yvals[i])/2
	plot(1:2,type="n",axes=FALSE,xlim=c(0,2),ylim=c(min,max),xlab="",ylab="")
	image(xvals,yvals,t(imagedat),add=TRUE,col=colpalette)
	if (is.null(axlabs)) {
		axis(2,las=1)
	} else if (logscale==TRUE) {
		axis(2,las=1,at=log10(axlabs),las=1,labels=axlabs)
	} else {
		axis(2,las=1,at=axlabs,las=1,labels=axlabs)
	}
	
}

# Function to generate a simple legend strip to be included in a multipart figure
srg.simple.image.legend.v2 <- function(breaks=0:10/10,colpalette=heat.colors(length(breaks)-1),logscale=FALSE,axlabs=NULL) {
	number=length(colpalette)
	imagedat <- matrix(nrow=number,ncol=2)
	xvals <- c(0,1,2)
	yvals <- breaks
	for (i in 1:number) imagedat[i,] <- (yvals[i+1]+yvals[i])/2
	plot(1:2,type="n",axes=FALSE,xlim=c(0,2),frame.plot=TRUE,ylim=c(min(breaks),max(breaks)),xlab="",ylab="")
	image(xvals,yvals,t(imagedat),add=TRUE,col=colpalette)
	axis(2,las=1)	
}

apa.gen.lhs.results <- function(
		nsamps=10,
		params = c("m_aa","m_ca","m_ac"),
		zmin=0.01,
		zmax=1,
		graphdir="~/Dropbox/projects/influenza/age_peak_attack/Current/figs/",
		outputfile="apa_lhs.csv") { 
	
	# This thing needs to make three separate figures for the three mixing values and then 
	noparams <- length(params)
	logzvals <- TRUE
	x <- apa.lhs(nsamps,params,rep(zmin,noparams),rep(zmax,noparams),veclog=c(logzvals,logzvals,logzvals))
	write.csv(x,paste(graphdir,outputfile,sep=""),row.names=FALSE)
	
}

apa.plot.lhs.results <- function(
		ncols=100,
		params = c("m_aa","m_ca","m_ac"),
		zmin=0.01,
		zmax=1,
		graphdir="~/Dropbox/projects/influenza/age_peak_attack/Current/figs/",
		lhsfile="apa_lhs.csv",
		chtJpeg=TRUE,
		cex=1.0) { 
	
	# This thing needs to make three separate figures for the three mixing values and then 
	logzvals <- TRUE
	x <- read.csv(paste(graphdir,lhsfile,sep=""))
	noparams <- length(params)
	colindex = heat.colors(ncols)
	nsamps = dim(x)[1]
	
	height <- 15
	width <- 10.09
	
	if (chtJpeg) {
		jpeg(filename=paste(graphdir,"scatter_lhs.jpg",sep=""),width=width,height=height,units="cm",res=600)
	} else {
		pdf(file=paste(graphdir,"scatter_lhs.pdf",sep=""),useDingbats=FALSE,width=width/cm(1),height=height/cm(1))
	}
	
	# Settings for all parts of the figure
	par(mai=(c(0,0,0,0)))
	nonaxis_cex <- 11/12
	axis_cex <- 10/12
	par(cex.axis=axis_cex)
	xlm=2.5/width
	xrm=2.5/width
	xg=0.5/width
	ybm=1.5/height
	yg=0.25/height
	ytm=0.25/height
	listpos <- list(
			srg.chart.pos(1,3,1,3,xlm=xlm,xrm=xrm,xg=xg,ybm=ybm,ytm=ytm,yg=yg),
			srg.chart.pos(1,2,1,3,xlm=xlm,xrm=xrm,xg=xg,ybm=ybm,ytm=ytm,yg=yg),
			srg.chart.pos(1,1,1,3,xlm=xlm,xrm=xrm,xg=xg,ybm=ybm,ytm=ytm,yg=yg))
	vecnew <- c(FALSE,TRUE,TRUE)
	xlabs <- c("","","Over-representation")
	ylabs <- c("","Peak attack rate","")
	letterlabs <- c(expression(bold(A) ~~ mu[aa]),expression(bold(B) ~~ mu[ca]),expression(bold(C) ~~ mu[ac]))
	
	# Plot the different charts
	for (i in 1:length(params)) {
		if (logzvals) {
			vals <- log10(x[,i])
			lb <- log10(zmin)
			ub <- log10(zmax)
			zvals <- (vals-lb) / (ub-lb)*(ncols-1)+1
		}
		else zvals = (x[,i] - zmin)/(zmax-zmin)*(ncols-1)+1
		yvals <- x[,5]
		par(fig=listpos[[i]],new=vecnew[i])
		plot(x[,4],yvals,pch=19,col=colindex[zvals],ylim=c(0,0.025),xlim=c(0,5),axes=FALSE,cex=cex)
		axis(2,at=c(0,0.005,0.01,0.015,0.02,0.025),las=1,mgp=c(0,1,0))
		mtext(ylabs[i],side=2,line=3)
		text(5,y=0.025,labels=letterlabs[i],adj=c(1,1))
		if (i==3) {
			axis(1,at=0:5,las=1,labels=0:5,mgp=c(0,0.5,0))
			mtext(xlabs[i],side=1,line=1.5)
		} else {
			axis(1,at=0:5,las=1,labels=rep("",6),mgp=c(0,0.5,0))
		}
	}
	
	par(fig=c(0.9,0.975,0.50,(srg.chart.pos(1,3,1,3,xlm=xlm,xrm=xrm,xg=xg,ybm=ybm,ytm=ytm,yg=yg)[4])),new=TRUE)
	srg.simple.image.legend(min=zmin,max=zmax,colpalette=colindex,logscale=logzvals,axlabs=c(0.01,0.02,0.05,0.1,0.2,0.5,1.0))
	
	# Turn off the device
	dev.off() 
	
}


apa.read.grump.china <- function(strFile = "~/Dropbox/tmp/chn_gpwfe_pcount_ascii_25/chnp10ag.asc") {
	require("sp")
	rtn <- read.asciigrid(strFile,proj4string="WGS84")
	rtn
}

eov.severity.setup.states <- function(p) {
	rtn <- vector(mode="numeric",length=0)
	rtn["S"] 	<- 1
	rtn["IA"] 	<- 0
	rtn["IB"] 	<- 0
	rtn["R"] 	<- 0
	rtn["cumA"] <- 0
	rtn["cumB"] <- 0
	rtn
}

# Look at baseline A
# Show all infections and severe case separately
# Show sesitivity to timing of seeding

eov.severity.model <- function(t,y,p){
	rtn 					<- vector(mode="numeric",length=length(y))
	names(rtn) 				<- names(y)
	betat 					<- 	p["R0"] * ( 1 - p["A"] * cos(2 * pi * t / 364) ) / p["TG"]
	infA					<- p["alpha_A"] + betat * y["S"] * y["IA"]
	infB					<- p["alpha_B"] + betat * p["f"] * y["S"] * y["IB"]
	rtn["S"] 				<- - infA - infB + y["R"] / p["TR"]
	rtn["IA"] 				<- + infA - y["IA"] / p["TG"]
	rtn["IB"] 				<- + infB - y["IB"] / p["TG"]
	rtn["R"] 				<- (y["IA"] + y["IB"]) / p["TG"] - y["R"] / p["TR"] 
	rtn["cumA"] 			<- + infA
	rtn["cumB"] 			<- + infB
	list(rtn)
}

# Function to make a 
eov.severity.make.inc <- function(sol,p) {
	
	# Gather aux data
	dims <- dim(sol)
	notps <- dims[1]
	
	# Set up the return object
	rtn <- data.frame(t=sol[,"time"],incA=rep(0,notps),incB=rep(0,notps),sincA=rep(0,notps),sincB=rep(0,notps))
	
	# Calculate the incidences
	rtn[2:notps,c("incA","incB")] <- sol[2:notps,c("cumA","cumB")] - sol[1:(notps-1),c("cumA","cumB")] 
	rtn[,"sincA"] <- rtn[,"incA"] * p["p"] 
	rtn[,"sincB"] <- rtn[,"incB"] * p["p"] * p["g"] 
	
	# Return values
	rtn
	
}

eov.severity.mod.plot <- 	function(	
		psToChange=c("f","g"),
		vecvals = c(20:10/10,rep(10,11)),
		xvals = (0:(52))*7,
		xplotrange = c(min(xvals),max(xvals)),
		makeplot = TRUE
) {
	
	# Add required packages
	# require(odesolve)
	require(deSolve)
	
	# Change the code below to run a loop and add the points to the chart
	
	nops <- length(psToChange)
	nosets <- length(vecvals)/nops
	if (nosets * nops != length(vecvals)) stop("vecvals and psToChange to consistent")
	vals <- matrix(data=vecvals,ncol=nops,nrow=nosets)
	cols <- rev(heat.colors(nosets))
	cols[6] <- "black"
	notps <- length(xvals)
	
	# Set up the variables for the loop 
	yvals <- matrix(ncol=nosets,nrow=notps)	
	
	# Loop through to make model runs
	for (i in 1:nosets) {
		
		# Update the parameters
		p_pb <- eov.severity.setup.params(pnames=psToChange,pvals=vals[i,])
		
		# Generate a single set of yval data
		# Needs to output the individual incidences and the number of severe cases
		inc <- eov.severity.runonce(p_pb,xvals)
		
		# Take results and add in the correct values to yvals
		yvals[,i] <- inc[,"sincA"]+inc[,"sincB"]
		
	}
	
	# Make the plots
	if (makeplot) {
		plot(1:2, type="n",xlim=xplotrange,ylim=c(0,max(yvals)))
		for (i in 1:nosets) points(xvals,yvals[,i],col=cols[i],type="l")
	}
	
	# Return the values
	list(yvals=yvals)
	
}

# Make up a new 
eov.severity.plot.single <- function(
		xvals = (0:(52))*7,
		pnames=NULL,
		pvals=NULL,
		plotdir="~/Dropbox/projects/virulence/observe_flu_severity/") {
	
	# Make up the parameter set
	ps <- eov.severity.setup.params(pnames=pnames,pvals=pvals)
	inc <- eov.severity.runonce(ps,xvals)
	
	# Start the plot file
	pdf(file=paste(plotdir,"single_run_log.pdf",sep=""),height=15/cm(1),width=15/cm(1))	
	
	# Set the chart margins
	xmin=0.1
	xmax=0.9
	ymin=0.1
	ymax=0.95
	
	incmin=1
	incmax=1000000
	
	# Plot resulting chart
	par(
			fig=c(xmin,xmax,ymin,ymax),
			mai=c(0,0,0,0),
			mgp=c(2,1,0),
			new=FALSE
	)
	plot(1:2,type="n",log="y",axes=FALSE,xlim=c(min(xvals),max(xvals)),ylim=c(incmin,incmax))
	axis(1,at=c((0:12)*30,364))
	axis(2,at=10^(0:6))
	points(xvals,10^6*(inc[,"incA"])+1,type="l",xpd=TRUE)
	points(xvals,10^6*(inc[,"incB"])+1,type="l",xpd=TRUE)
	
	srg.plotstack.dash(inc,"red",var="incA",const=1,timevar="t",mult=10^6)
	srg.plotstack.dash(inc,"blue",var="incB",varbelow=c("incA"),const=1,timevar="t",mult=10^6)
	
	mtext("Time (days)",side=1,line=2)
	mtext("Infection incidence (colours, log scale)",side=2,line=2)
	mtext("Severe case incidence (black line, log scale)",side=4,line=2)
	
	
	# Plot resulting chart
	par(
			fig=c(xmin,xmax,ymin,ymax),
			mai=c(0,0,0,0),
			mgp=c(2,1,0),
			new=TRUE
	)
	plot(1:2,type="n",log="y",axes=FALSE,xlim=c(min(xvals),max(xvals)),ylim=c(1,500))
	axis(4)
	points(xvals,9*10^6*(inc[,"sincA"]+inc[,"sincB"])+1,type="l",xpd=TRUE,col="black",lwd=2)
	
	dev.off()
	
}


# Run a single solution of the emergence model for eov
eov.severity.runonce <- function(ps,xvs) {
	
	require("odesolve")
	
	# This needs to handle the time junction a little better	
	s_pb 	<- eov.severity.setup.states(ps)
	maxts <- max(xvs)
	ps["alpha_A"] <- ps["alpha"]
	ps["alpha_B"] <- 0
	
	if (ps["tB"] < maxts) {
		
		ts_partA <- c(xvs[xvs <= ps["tB"]],ps["tB"])
		sol_pb_1	<- lsoda(s_pb,ts_partA,eov.severity.model,ps,atol=1e-20)
		nvars <- dim(sol_pb_1)[2]-1
		laststep <- dim(sol_pb_1)[1]
		restart <- sol_pb_1[laststep,2:(1+nvars)]
		ps["alpha_A"] <- 0
		ps["alpha_B"] <- ps["alpha"]
		ts_partB <- c(ps["tB"],xvs[xvs > ps["tB"]])
		sol_pb_2	<- lsoda(restart,c(ps["tB"],xvs[xvs > ps["tB"]]),eov.severity.model,ps,atol=1e-20)
		sol_pb		<- rbind(sol_pb_1[1:(dim(sol_pb_1)[1]-1),],sol_pb_2[2:(dim(sol_pb_2)[1]),]) 
		
	} else {
		
		sol_pb <- lsoda(s_pb,xvs,eov.severity.model,ps,atol=1e-20)
		
	}
	
	inc_pb 		<- eov.severity.make.inc(sol_pb,ps)
	
	inc_pb
	
}

eov.severity.setup.params <- function(base=NULL,pnames=NULL,pvals=NULL) {
	
	# Initialize to defaul values if no base option given
	if (identical(base,NULL)) {
		rtn <- vector(mode="numeric",length=0)
		rtn["R0"] 		<- 2.32
		rtn["TG"] 		<- 2.6
		rtn["TR"]		<- 5.35*364
		rtn["alpha"] 	<- 1e-6/14
		rtn["alpha_A"] 	<- 0
		rtn["alpha_B"] 	<- 0
		rtn["A"]		<- 0.51
		rtn["tB"]		<- 80
		rtn["p"]		<- 1e-4
		rtn["f"]		<- 1.0
		rtn["g"]		<- 1.0
		
		# Else start from the base values
	} else {
		rtn <- base
	}
	
	# Check for suggested parameter name changes
	# Make this bit a function
	if (!identical(pnames,NULL)) {
		noPs <- length(pnames)
		if (noPs!=length(pvals)) stop("Problem with lengths in read.severity.setup.params")
		rtnnames <- names(rtn)
		for (i in 1:noPs) {
			if (!(pnames[i] %in% rtnnames)) stop("Subbing a parameter that isn't present in read.severity.setup.params")
			rtn[pnames[i]] <- pvals[i]
		}
	}
	
	# Do any other required posthoc comparisons
	
	# Return result
	rtn
	
}

# Make all the plots and 
eov.severity.main <- function(plotdir="~/Dropbox/projects/virulence/observe_flu_severity/") {
	
	eov.severity.plot.single(pnames=c("f","g","tB"),pvals=c(1.5,10,80))
	
	pdf(file=paste(plotdir,"pandemic_raw.pdf",sep=""),height=12/cm(1),width=12/cm(1))
	vals <- eov.severity.mod.plot()
	dev.off()
	
	pnameslong <- c("f","g","tB","alpha")
	pvalslong <- c(1.5,2.0,20*364,1e-9)
	
	pdf(file=paste(plotdir,"seasonal_wide_raw.pdf",sep=""),height=10/cm(1),width=14/cm(1))
	vals <- eov.severity.mod.plot(	
			psToChange=pnameslong,
			vecvals=pvalslong,
			xvals = (0:(40*12))*(364/12))
	dev.off()
	
	pdf(file=paste(plotdir,"seasonal_narrow_raw.pdf",sep=""),height=8/cm(1),width=12/cm(1))
	vals <- eov.severity.mod.plot(	
			psToChange=pnameslong,
			vecvals=pvalslong,
			xvals = (0:(40*12))*(364/12),
			xplotrange = c(17*12,22*12)*(364/12))
	dev.off()
	
}

# Function to return a list of the ids of the subset of individuals taken for malik's eid study
eid.pigstrains.subset <- function() {
	
	rtn <- c(	"S090191-0"	,
			"S090113-1"	,
			"S090191-1"	,
			"S090454-4"	,
			"S090418-2"	,
			"S090241-0"	,
			"S090051-0"	,
			"S090118-0"	,
			"S090222-1"	,
			"S090451-0"	,
			"S090353-1"	,
			"S090477-1"	,
			"S090060-1"	,
			"S090024-4"	,
			"S090191-2"	,
			"S090242-1"	,
			"S090184-4"	,
			"S090150-0"	,
			"S090218-3"	,
			"S090481-2"	,
			"S090413-4"	,
			"S090222-3"	,
			"S090313-0"	,
			"S090313-1"	,
			"S090484-1"	,
			"S090513-2"	,
			"S090538-0"	,
			"S090556-3")
	rtn
	
}

# Function to take the pp_data object from the severity script and add on a variable for whether that sample was in 
# the eid pig sera substudy
eid.pigstrains.addvar <- function(base,studyids) {
	
	rtn <- base
	rtn$eidsub <- rtn$ind_id %in% eid.pigstrains.subset()
	rtn
	
}

# Example figure for on.
plot.example.for.on <- function() {
	
	dummyy = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
	dummyuci = dummyy + 0.09
	dummylci = dummyy - 0.09
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=c(0,5),
			ylim=c(0.01,1),
			log="y")
	axis(1,at=0:5,las=3,labels=c("","A","B","C","D",""))
	axis(2)
	mtext("No",side=1,line=2)
	mtext("onset",side=1,line=3)
	
	offset <- 0.15
	xvals <- c(1-offset,1+offset,2-offset,2+offset,3-offset,3+offset,4-offset,4+offset)
	points(xvals,dummyy)
	for (i in 1:length(dummyy)) lines(c(xvals[i],xvals[i]),c(dummylci[i],dummyuci[i]))
	
	# points(3.5,ICU_no_onset_pt,col="Black")
	# lines(c(3.5,3.5),c(ICU_no_onset_ub,ICU_no_onset_lb),col="Black")
	
}

# source("~/Dropbox/svneclipse/idsource/R/stevensRFunctions.R")
# apa.baseline()
apa.baseline <- function(	strFile="~/Dropbox/projects/influenza/age_peak_attack/Current/figs/illustrate_hetero_base.pdf",
		lbinc=0.0067,
		ubinc=0.0098) {
	
	# Define functions that are largely private
	
	# A little routine for adding a grey polygon
	addpolygon <- function() {
		polygon(c(axstart,axend,axend,axstart),c(lbinc,lbinc,ubinc,ubinc),col="grey",border=NA)
		lines(c(axstart,axend),c(lbinc,lbinc),lty=2)
		lines(c(axstart,axend),c(ubinc,ubinc),lty=2)
	}
	
	# A function to print out sum summaries of the incidence results
	repvals <- function(incX,p,lab="YY") {
		cat(lab," Peak attack rate",max(incX[,"total"]),"\n")
		cat("   Over-representation",p["or"],"\n")
		cat("   Cumulative attack rate: children ",sum(incX[,"children"]),", adults ",sum(incX[,"adults"]),", and total ",sum(incX[,"total"]),"\n\n")
	}
	
	require("date")
	require("odesolve")
	
	# Sol A is the mass action baseline
	pA 	 <- apa.setup.params.spq()
	#sA 	 <- apa.setup.state.vector.2state(pA)
	#solA <- lsoda(sA,0:180,apa.spa.2state.sus,pA)
	sA 	 <- apa.setup.states.2state.SEIR(pA)
	solA <- lsoda(sA,0:180,apa.spa.2state.SEIR,pA)
	incA <- apa.extract.inc(solA,pA)
	incA[,"time"] <- incA[,"time"] + as.date("1Apr2009")
	repvals(incA,pA,lab="A")
	
	# Sol B is adjusted for age-specific susceptibility
	pB 	 <- apa.setup.params.spq(pnames=c("phi"),pvals=c(0.351))
	#sB 	 <- apa.setup.state.vector.2state(pB)
	#solB <- lsoda(sB,0:180,apa.spa.2state.sus,pB)
	sB 	 <- apa.setup.states.2state.SEIR(pB)
	solB <- lsoda(sB,0:180,apa.spa.2state.SEIR,pB)
	incB <- apa.extract.inc(solB,pB)
	incB[,"time"] <- incB[,"time"] + as.date("1Apr2009") 
	repvals(incB,pB,"B")
	
	# Adjust for polymod-reported mixing frequencies
	pC 	 <- apa.setup.params.spq(pnames=c("phi","m_ca","m_ac","m_aa"),pvals=c(0.351,0.28,0.26,0.38))
	#sC 	 <- apa.setup.state.vector.2state(pC)
	#solC <- lsoda(sC,0:180,apa.spa.2state.sus,pC)
	sC 	 <- apa.setup.states.2state.SEIR(pC)
	solC <- lsoda(sC,0:180,apa.spa.2state.SEIR,pC)
	incC <- apa.extract.inc(solC,pC)
	incC[,"time"] <- incC[,"time"] + as.date("1Apr2009")
	repvals(incC,pC,"C")
	
	# Setup plot auxes
	axticks <- as.date(c(
					"1April2009","1May2009","1Jun2009","1Jul2009","1Aug2009"))
	
	axstart <- min(axticks)
	axend	<- max(axticks)
	axnames	<- c("0","1","2","3","4")	
	ytickvA	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02)
	yticknA	<- c("0.0%","","0.5%","","1.0%","","1.5%","","2.0%")
	ymaxA 	<- max(ytickvA)
	
	# Total widths
	fwidth <- 15
	fheight <- 7
	
	# Ratios of plottable areas
	xr 	<- c(1,1,1)
	
	# Gaps in cms divded by total widths
	yg1 <- 1.5 / fheight
	yg2 <- 0.5 / fheight
	yg2b <- 1.5 / fwidth
	xg1 <- 2.0 / fwidth
	xg2 <- 0.5 / fwidth
	pwidth <- (1-xg1-3*xg2)/3
	pheight <- (1-yg1-yg2)
	
	# Define plottable areas
	posA <- c(xg1,xg1+pwidth,yg1,yg1+pheight)
	posB <- c(xg1+pwidth+xg2,xg1+2*pwidth+xg2,yg1,yg1+pheight)
	posC <- c(xg1+2*pwidth+2*xg2,xg1+3*pwidth+2*xg2,yg1,yg1+pheight)	
	
	# Open the pdf file
	pdf(strFile,height=fheight/cm(1),width=fwidth/cm(1),useDingbats=FALSE)
	
	# Settings for all parts of the figure
	par(mai=(c(0,0,0,0)))
	par(mgp=c(2.5,1,0))
	nonaxis_cex <- 11/12
	axis_cex <- 10/12
	par(cex.axis=axis_cex)
	
	# First chart setup
	par(fig=posA)
	plot(1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxA),ylab="")
	axis(2,las=1,at=ytickvA,labels=yticknA,hadj=0.8)
	text(axstart,max(ytickvA),"A",font=2,adj=c(0,1))
	mtext("Daily incidence",side=2,line=2.7,cex=nonaxis_cex)
	axis(1,at=axticks,labels=axnames,line=0,padj=-0.8)
	
	# Plot the lines
	addpolygon()
	srg.plotstack.dash(incA,"red",var="children",varbelow=NULL)
	srg.plotstack.dash(incA,"blue",var="adults",varbelow=c("children"))
	points(incA[,"time"],incA[,"total"],type="l",col="black",lwd=2,lty=1)	
	
	# Second chart setup
	par(fig=posB,new=TRUE)
	plot(1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxA),ylab="")
	text(axstart,max(ytickvA),"B",font=2,adj=c(0,1))
	axis(1,at=axticks,labels=axnames,line=0,padj=-0.8)
	mtext("Time (months)",side=1,line=1.5,cex=nonaxis_cex)
	
	# Plot the lines
	addpolygon()
	srg.plotstack.dash(incB,"red",var="children",varbelow=NULL)
	srg.plotstack.dash(incB,"blue",var="adults",varbelow=c("children"))
	points(incB[,"time"],incB[,"total"],type="l",col="black",lwd=2,lty=1)	
	
	# Third chart setup
	par(fig=posC,new=TRUE)
	plot(	1:2,type="n",axes=FALSE,xlim=c(axstart,axend),ylim=c(0,ymaxA),ylab="")
	axis(1,at=axticks,labels=axnames,line=0,padj=-0.8)
	text(axstart,max(ytickvA),"C",font=2,adj=c(0,1))
	
	# Plot the lines
	addpolygon()
	srg.plotstack.dash(incC,"red",var="children",varbelow=NULL)
	srg.plotstack.dash(incC,"blue",var="adults",varbelow=c("children"))
	points(incC[,"time"],incC[,"total"],type="l",col="black",lwd=2,lty=1)	
	
	# Add a legend
	legend(as.date("07April2009"),max(ytickvA)*3.5/4,
			c(	"Children",
					"Adults",
					"Range 2009 ANZ"),
			lty=c(1,1,1),
			lwd=c(8,8,8),
			col=c("red","blue","grey"),
			cex=axis_cex*0.8, 
			bty="n",
	)
	
	dev.off()
	
}

apa.sense.one <- function() {
	
	require("date")
	
	# Define the illustrative parameters
	ill_params <- 
			c(		N 			= 7000000,
					seed 		= 100,
					trickle		= 1.0,
					r 			= (1.4-1.0)/2.6,
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
					t0 			= as.date("1Apr2009"),
					p_1			= 0.2,			
					p_2			= 0.4,
					p_3			= 0.4,
					mixmatindex	= 4,										# CAE from polymod
					mixsense	= 0,
					fitsus		= 0,
					noyears		= 1,
					nodts		= 52*7,
					dur_seas	= 7*13,
					amp_seas	= 0.10,
					gamma_R		= 1/99999999,
					aging_on	= 0,
					mu_1		= 1/364/19, 
					mu_2		= 1/364/(65-19),
					mu_3		= 1/364/50									# Set to give very slow overall population growth
			)
	
	illRateICU <- 10/100000
	
	
	popsize <- ill_params["N"]
	
	# Runs for the sensitivity analysis figure
	x_lims 		<- c(0,1.2)		# r
	y_lims 		<- c(0.5,5)		# p_R * p_I * p_H_base
	no_bins_x	<- 4
	no_bins_y	<- 4			# often either 2 or 10
	x_size 		<- (x_lims[2]-x_lims[1])/no_bins_x
	y_size		<- (y_lims[2]-y_lims[1])/no_bins_y
	x_bounds	<- x_lims[1]+(0:no_bins_x) * x_size
	y_bounds	<- y_lims[1]+(0:no_bins_y) * y_size
	x_mids		<- x_lims[1] + (1:no_bins_x) * x_size - x_size / 2
	y_mids		<- y_lims[1] + (1:no_bins_y) * y_size - y_size / 2
	
	# Make the array to be plotted
	rvals <- c((1.2-1.0)/2.6,(1.4-1.0)/2.6,(1.8-1.0)/2.6)
	# rvals <- c((1.4-1.0)/2.6)
	noRs <- length(rvals)
	chtData <- array(dim=c(no_bins_x,no_bins_y,noRs))
	chtData_pi <- array(dim=c(no_bins_x,no_bins_y,noRs))
	for (i in 1:no_bins_x) {
		for (j in 1:no_bins_y) {
			# generate solution
			for (k in 1:noRs) {
				sol_tmp <- (SirModel3AgeClasses(
									pname=c("r","N","phi_1","mixsense","fitsus"),
									pvals=c(rvals[k],popsize,y_mids[j],x_mids[i],0),
									casemix=c(20,40,40),
									vp=ill_params))$sol
				# q[i,j,k] <- sum(sol_tmp[,"dS"])/popsize
				chtData[i,j,k] <- illRateICU*max(sol_tmp[,"dS"])/popsize*100000
			}
			cat(paste("Completed",j + (i-1)*no_bins_y," of ",no_bins_x*no_bins_y,"\n"))
			flush.console()
		}
	}
	
	# Generate the latin hypercube samples for the scatter plot
	rep_bnds <- c(1,5)
	mix_bnds <- c(0,5)
	no_samples <- 10
	set.seed(1234)
	hv_rep <- srg.hyper.vector(no_samples,rep_bnds[1],rep_bnds[2])
	hv_mix <- srg.hyper.vector(no_samples,mix_bnds[1],mix_bnds[2])
	out_ar <- matrix(nrow=no_samples,ncol=noRs)
	out_pi <- matrix(nrow=no_samples,ncol=noRs)
	out_phi <- matrix(nrow=no_samples,ncol=noRs)
	
	# Set up the output variables
	for (i in 1:no_samples) {
		for (j in 1:noRs) {
			sol_tmp <- SirModel3AgeClasses(
					pname=c("r","N","phi_1","mixsense","fitsus"),
					pvals=c(rvals[j],popsize,1,hv_mix[i],1),
					casemix=c(20*hv_rep[i],40,40),
					vp=ill_params)
			out_ar[i,j] <- sum((sol_tmp$sol)[,"dS"])/popsize
			out_pi[i,j] <- illRateICU*max((sol_tmp$sol)[,"dS"])/popsize
			out_phi[i,j] <- sol_tmp$par["phi_1"]
		}
		cat(paste("Completed",i," of ",no_samples,"\n"))
		flush.console()
	}
	
	axticks <- as.date(c(
					"1April2009","1May2009","1Jun2009","1Jul2009","1Aug2009"))
	
	axstart <- min(axticks)
	axend	<- max(axticks)
	axnames	<- c("0","1","2","3","4")
	# ytickvA	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175)
	# yticknA	<- c("0.0%","","0.5%","","1.0%","","1.5%","")
	# ymaxA 	<- max(ytickvA)
	# ytickvB	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)
	# yticknB	<- c("0.0%","","0.5%","","1.0%","","1.5%")
	# ymaxB 	<- max(ytickvB)
	# ytickvC	<- c(0,0.0025,0.005,0.0075,0.01,0.0125)
	# yticknC	<- c("0.0%","","0.5%","","1.0%","")
	# ymaxC 	<- max(ytickvC)
	
	ytickvA	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175)
	yticknA	<- c("0.0%","","0.5%","","1.0%","","1.5%","")
	ymaxA 	<- max(ytickvA)
	ytickvB	<- c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)
	yticknB	<- c("0.0%","","0.5%","","1.0%","","1.5%")
	ymaxB 	<- max(ytickvB)
	ytickvC	<- c(0,0.0025,0.005,0.0075,0.01,0.0125)
	yticknC	<- c("0.0%","","0.5%","","1.0%","")
	ymaxC 	<- max(ytickvC)
	
	# Total widths
	xw <- 15
	yw <- 13
	
	# Ratios of plottable areas
	yra <- c(ymaxC,ymaxB,ymaxA)
	yrb <- c(1,1)
	xr 	<- c(5,5,0.5)
	
	# Gaps in cms divded by total widths
	yg1 <- 1.5 / yw
	yg2 <- 0.5 / yw
	yg2b <- 1.5 / yw
	xg1 <- 2.0 / xw
	xg2 <- 1.6 / xw
	xg3 <- 1.0 / xw
	xg4 <- 0.25 / xw
	
	# Sizes of plottable areas
	ysa <- yra / sum(yra) * (1-length(yra)*yg2-yg1) 
	ysb <- yrb / sum(yrb) * (1-yg2b-yg1-yg2) 
	xs	<- xr  / sum(xr) * (1-xg2-xg3-xg1-xg4)
	
	# Define plottable areas
	posA <- c(xg1,xg1+xs[1],yg1+ysa[1]+yg2+ysa[2]+yg2,yg1+ysa[1]+yg2+ysa[2]+yg2+ysa[3])
	posB <- c(xg1,xg1+xs[1],yg1+ysa[1]+yg2,yg1+ysa[1]+yg2+ysa[2])
	posC <- c(xg1,xg1+xs[1],yg1,yg1+ysa[1])
	posD <- c(xg1+xs[1]+xg2,xg1+xs[1]+xg2+xs[2],yg1+ysb[1]+yg2b,yg1+ysb[1]+yg2b+ysb[2])
	posE <- c(xg1+xs[1]+xg2,xg1+xs[1]+xg2+xs[2],yg1,yg1+ysb[1])
	posF <- c(xg1+xs[1]+xg2+xs[2]+xg3,xg1+xs[1]+xg2+xs[2]+xg3+xs[3],yg1+ysb[1]+yg2b,yg1+ysb[1]+yg2b+ysb[2])
	posG <- c(xg1+xs[1]+xg2+xs[2]+xg3,xg1+xs[1]+xg2+xs[2]+xg3+xs[3],yg1,yg1+ysb[1])
	
	pdf("~/Dropbox/projects/influenza/age_peak_attack/Current/figs/illustrate_hetero_sens_one.pdf",height=yw/cm(1),width=xw/cm(1),
			useDingbats=FALSE)
	
	# Settings for all parts of the figure
	# Settings for the default spaces around the plot areas
	par(mai=(c(0,0,0,0)))
	# Not sure
	par(mgp=c(2.5,1,0))
	# Scaling of axis fonts
	nonaxis_cex <- 11/12
	axis_cex <- 10/12
	par(cex.axis=axis_cex)	
	
	# Plot the heat chart
	par(fig=posD)
	
	zlimits <- c(0.1,0.6)
	zlimits2 <- c(0.0,0.5)
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=x_lims,
			ylim=y_lims)
	
	axis(2,las=1,hadj=0.3)
	axis(1,las=1,padj=-0.8)
	
	cscale1 <- rev(heat.colors(130))
	
	image(x_mids,y_mids,chtData[,,1],add=TRUE,zlim=zlimits,col=cscale1)
	
	cntlevs <- c(0.2,0.25,0.3,0.4,0.5)
	contour(x_mids,y_mids,chtData[,,1],
			add=TRUE,
			levels=cntlevs,
			labels=as.character(cntlevs)
	)
	
	
	mtext("Susceptibility",side=2,line=1.5,cex=nonaxis_cex)
	mtext("Mixing",side=1,line=1.5,cex=nonaxis_cex)
	text(x_lims[1],y_lims[2],"D",font=2,adj=c(0,1))
	
	# Plot the legend the heat chart
	legend_resolution <- 50
	start <- zlimits[1]
	end <- zlimits[2]
	width <- end - start
	binsize <- width / legend_resolution
	zmids <- start + binsize / 2 + binsize * (0:(legend_resolution-1))
	legend_data <- matrix(zmids,nrow=legend_resolution,ncol=2)
	
	par(fig=posF,new=TRUE)
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=c(0,1),
			ylim=zlimits)
	
	image(c(0.25,0.75),zmids,t(legend_data),add=TRUE,zlim=zlimits,col=cscale1)
	axis(2,las=1,hadj=0.7,pos=-0.3)
	
	# Make the scatter plot of the dependency on over-representation 
	par(fig=posE,new=TRUE)
	plot(1:2,type="n",axes=FALSE,xlim=rep_bnds,ylim=zlimits2)
	axis(1,las=1,padj=-1.0)
	axis(2,las=1,hadj=0.7)
	mtext("Cumulative attack rate",side=2,line=2.0,cex=nonaxis_cex)
	mtext("Over-representation",side=1,line=1.5,cex=nonaxis_cex)
	text(rep_bnds[1],zlimits2[2],"E",font=2,adj=c(0,1))
	
	# Plot the legend to the scatterplot
	colres <- 100
	cscale2 <- cm.colors(colres)
	legend_resolution <- 50
	start <- 0
	end <- 2
	width <- end - start
	binsize <- width / legend_resolution
	zmids <- start + binsize / 2 + binsize * (0:(legend_resolution-1))
	legend_data <- matrix(zmids,nrow=legend_resolution,ncol=2)
	
	vec_col <- as.vector(1:no_samples)
	for (i in 1:no_samples) vec_col[i] <- cscale2[round((out_phi[i]-start)/(end-start)*colres)]
	for (i in 1:noRs) {
		points(hv_rep,out_pi[,i]*100000,col=vec_col,pch=19)
		points(hv_rep,out_pi[,i]*100000,pch=1)
	}
	
	# Plot the legend for the scatter plot
	par(fig=posG,new=TRUE)
	
	plot(	1:2,
			type="n",
			axes=FALSE,
			xlim=c(0,1),
			ylim=c(start,end))
	
	image(c(0.25,0.75),zmids,t(legend_data),add=TRUE,zlim=c(start,end),col=cscale2)
	axis(2,las=1,hadj=0.7,pos=-0.3)
	
	dev.off()
	
}

# Read in the main study data from fluscape
fsc.read.data <- function(fsc.root="~/Dropbox/svneclipse/fluscape") {
	loc <- read.csv(paste(fsc.root,"/data/Location_Info_V1.csv",sep=""),na.strings=c("NA","\\N"))
	hh <- read.csv(paste(fsc.root,"/data/HouseHolds_V1.csv",sep=""),na.strings=c("NA","\\N"))
	part <- read.csv(paste(fsc.root,"/data/Participants_V1.csv",sep=""),na.strings=c("NA","\\N"))
	cont <- read.csv(paste(fsc.root,"/data/Contacts_V1.csv",sep=""),na.strings=c("NA","\\N"))
	his <- read.csv(paste(fsc.root,"/data/HI_Results_V1.csv",sep=""),na.strings=c("NA","\\N"))
	list(loc=loc,hh=hh,part=part,cont=cont,his=his)
}

# Some R commands to play with the data
fsc.explore <- function(fluscapedir="~/Dropbox/svneclipse/fluscape/",plotscatter=FALSE) {
	
	# load.and.merge is defined in GeneralUtility.r in fluscape source/R directory
	source(paste(fluscapedir,"source/R/GeneralUtility.r",sep=""))
	x <- load.and.merge.part.V1(topdir=fluscapedir)
	if (plotscatter) plot(x$PART_AGE,jitter(fsc.scale.titre(x$HI.H1N1.2009.PDM)))
	x
	
}

# A function to put neutralization titres onto a non-negative integer scale
fsc.scale.titre <- function(titre) {
	if (titre < 1) rtn <- 0
	else rtn <- log2(titre/5)
	rtn
}

# Takes vectors of location, household and participant number and returns a vector of participant table numbers
# Change this to a routine that uses merge to make up the dataset as per Justin's function
fsc.gen.part.index <- function(tabPart,location,household,ind) {
	noparts <- dim(tabPart)[1]
	nolookups <- length(location)
	rtn <- vector(mode=numeric,length=nolookups)
	if (length(household)!=nolookups || length(ind!=nolookups)) stop("Problem in call to fsc.gen.part.index")
	currentlookup <- 1
	while (currentlookup < nolookups) {
		currentpart <- 1
		while (	!(	identical(tabPart[currentpart,"LOC_ID"],location[currentlookup]) &&
					identical(tabPart[currentpart,"HH_ID"],location[currentlookup]) &&
					identical(tabPart[currentpart,"PARTICIPANT_ID"],location[currentlookup]))) currentpart <- currentpart + 1
	}
}


# To get the grid coords of a (long,lat) (x,y) value from a asciigrid as an image
# Needs some default arguments so that it doesn't recalcuate min max etc every time
srg.popimage.dblGetCoords <- function(x,y,pg) {
	minx <- pg$x[1]
	miny <- pg$y[1]
	nx <- length(pg$x)
	ny <- length(pg$y)
	dx <- pg$x[2] - minx
	dy <- pg$y[2] - miny
	
	if ( 	x < minx || x > minx + (nx-1)*dx 	||
			y < miny || y > miny + (ny-1)*dy	||		
			any(is.na(x),is.na(y))) rtn <- c(-1,-1)
	
	else {
		
		rx <- floor(1 + (x - minx) / dx)
		ry <- floor(1 + (y - miny) / dy)
		rtn <- c(rx,ry)
	}
	rtn
}

# To get the mid point of a cell, same caveats as function above re need for optimaization
srg.popimage.dblMidCell <- function(xi,yi,pg) {
	
	minx <- pg$x[1]
	miny <- pg$y[1]
	nx <- length(pg$x)
	ny <- length(pg$y)
	dx <- pg$x[2] - minx
	dy <- pg$y[2] - miny
	
	if ( 	(xi < 1 || xi > nx) ||
			(yi < 1 || yi > ny)		) rtn <- c(9999,9999)
	else {
		rx <- minx + dx * (xi-1) + dx/2.0
		ry <- miny + dy * (yi-1) + dy/2.0
		rtn <- c(rx,ry)
	}
	
	rtn
	
}

# Calculate distance on a grid between two squares 
srg.popimage.gridDist <- function(x1,y1,x2,y2,pg) {
	coord1 <- srg.popimage.dblMidCell(x1,y1,pg)
	coord2 <- srg.popimage.dblMidCell(x2,y2,pg)
	dist <- srg.decLongLatDist(coord1[1],coord1[2],coord2[1],coord2[2],translate=TRUE)
	dist 
}

## To calculate the integrated connectivity for a set of longitude and latitudes, based on the
## - latitude of locations (of length n)
## - longitude of locations (must be of length n)
## - population numbers as an image
## - maximum radius to consider r_max
## - vector of offsets length m for offset power k(x) = 1 / (1+(offset/x^nnpower))
## - vector of non-negative power (must also be length m) for k(x) = 1 / (1+(offset/x^nnpower))
#fsc.calc.connect <- function(loclat,loclong,popimage,r_max,vecoffsets,vecpower,dbuniform=FALSE,dbflat=FALSE) {
#	
#	
#	# Validate the inputs
#	nolocs <- length(loclat)
#	if (!identical(nolocs,length(loclong))) stop("fsc.calc.connect: lats and longs must be of the same length")
#	nooffsets <- length(vecoffsets)
#	nopowers <- length(vecpower)
#	
#	# Set up the variables to be returned
#	connectivity <- array(data=-1,dim=c(nolocs,nooffsets,nopowers))
#	
#	# Define some variables that will be needed for a few things
#	nox <- length(popimage$x)
#	noy <- length(popimage$y)
#	if (dbuniform) popimage$z[] <- 1
#	
#	# Check that none of the boundary squares are too close to the sample squares
#	for (i in 1:nox) {
#		# Not yet implemented
#	}
#	
#	# Cycle through twice to build a list of the x and y for all cells that have a non-zero interaction
#	nononezero <- 0
#	allpossdistances <- array(data=-1,dim=c(nox,noy,nolocs))
#	cat("Calculating all possible distances\n")
#	pb <- txtProgressBar(min = 0, max = nox, style = 3)
#	for (i in 1:nox) {
#		for (j in 1:noy) {
#			if (!is.na(popimage$z[i,j]) && popimage$z[i,j] > 0) {
#				for (k in 1:nolocs) {
#					# sqcentre <- srg.popimage.dblMidCell(i,j,popimage)
#					sqcentre <- c(popimage$x[i],popimage$y[j])
#					x1 <- sqcentre[1]
#					y1 <- sqcentre[2]
#					x2 <- loclong[k]
#					y2 <- loclat[k]
#					tmpdist <- srg.decLongLatDist(x1,y1,x2,y2,translate=TRUE)
#					if (tmpdist < r_max) {
#						allpossdistances[i,j,k] <- tmpdist
#						nononezero <- nononezero + 1
#					}
#				}
#			} 
#		}
#		setTxtProgressBar(pb, i)
#	}
#	# Probably works pretty well above here, but needs attention below
#	
#	
#	# Cycle back through the same set assigning xcoords, ycoods, locations and distances to the 
#	reducedset <- matrix(nrow=nononezero,ncol=3,dimnames=list(1:nononezero,c("location","number","distance"))) 
#	currentrow <- 1
#	cat("Populating reduced set\n")
#	pb <- txtProgressBar(min = 0, max = nox, style = 3)
#	for (i in 1:nox) {
#		for (j in 1:noy) {
#			for (k in 1:nolocs) {
#				if (allpossdistances[i,j,k] > 0) {
#					reducedset[currentrow,"location"] <- k
#					reducedset[currentrow,"distance"] <- allpossdistances[i,j,k]
#					reducedset[currentrow,"number"] <- popimage$z[i,j]					
#					currentrow <- currentrow + 1
#				}
#			}
#		}
#		setTxtProgressBar(pb, i)
#	}
#	
#	# table(reducedset[,"location"],reducedset[,"number"])
#
#	# Loop to calculate the cumulative neighbourhood size for the specified parameters
#	for (i in 1:length(vecoffsets)) {
#		off <- vecoffsets[i]
#		for (j in 1:length(vecpower)) {
#			pow <- vecpower[j]
#			veclocs <- vector(length=nolocs)
#			rowno <- 1
#			while (rowno <= nononezero) {
#				dist 			<- reducedset[rowno,"distance"]
#				nopeople 		<- reducedset[rowno,"number"]
#				loc 			<- reducedset[rowno,"location"]
#				tmp 			<- 1/(1+(dist/off)^pow)
#				veclocs[loc] 	<- veclocs[loc] + tmp
#				rowno 			<- rowno + 1
#			}
#			connectivity[,i,j] <- veclocs[]
#		}
#	}
#	
#	# Return matrix of connectivity
#	list(con=connectivity,offs=vecoffsets,pows=vecpower)
#	
#}

fsc.gen.connect.dists.new <- function(
		loclat,										# Vector of latitudes of origin points							
		loclong,									# Vector of longitudes of origin points
		locid,										# Vector of unique ID integers for origins 
		popimage,									# R image of population density length(popimage$x) == dim(popimage$z)[1]
		rast,										# Raster image of population density
		r_max=30,									# The maximum distance to consider
		dbuniform=FALSE,							# Debug option to ignore popimage and assume uniform density == 1
		centreorigins=TRUE,							# Option to shift origin locations to their nearest cell centre
		csvfile="~/Dropbox/tmp/popdistandsize.csv"	# Name of the output file to write (may be large)
) {
	
	# Validate the inputs
	nolocs <- length(loclat)
	if (!identical(nolocs,length(loclong))) stop("fsc.calc.connect: lats and longs must be of the same length")
	
	# Define some variables that will be needed for a few things
	nox <- length(popimage$x)
	noy <- length(popimage$y)
	if (dbuniform) popimage$z[] <- 1
	
	# Check that none of the boundary squares are too close to the sample squares
	for (i in 1:nox) {
		# Not yet implemented
	}
	
	# Check that origins are within popimage and shift to cell centre if required
	if (centreorigins) {
		dx <- popimage$x[2]-popimage$x[1]
		dy <- popimage$y[2]-popimage$y[1]
		minx <- popimage$x[1] - dx/2
		miny <- popimage$y[1] - dy/2
		for (i in 1:nolocs) {
			# XXX up to here
			indx <- ceiling((loclong[i] - minx) / dx)
			indy <- ceiling((loclat[i] - miny) / dy)
			loclong[i] <- minx + indx * dx - dx / 2
			loclat[i] <- miny + indy * dy - dy / 2
		}
	}
	
	# Cycle through twice to build a list of the x and y for all cells that have a non-zero interaction
	nononezero <- 0
	allpossdistances <- array(data=-1,dim=c(nox,noy,nolocs))
	cat("Calculating all possible distances\n")
	pb <- txtProgressBar(min = 0, max = nox, style = 3)
	for (i in 1:nox) {
		for (j in 1:noy) {
			if (!is.na(popimage$z[i,j]) && popimage$z[i,j] > 0) {
				for (k in 1:nolocs) {
					# sqcentre <- srg.popimage.dblMidCell(i,j,popimage)
					sqcentre <- c(popimage$x[i],popimage$y[j])
					x1 <- sqcentre[1]
					y1 <- sqcentre[2]
					x2 <- loclong[k]
					y2 <- loclat[k]
					tmpdist <- srg.decLongLatDist(x1,y1,x2,y2,translate=TRUE)
					if (tmpdist < r_max) {
						allpossdistances[i,j,k] <- tmpdist
						nononezero <- nononezero + 1
					}
				}
			} 
		}
		setTxtProgressBar(pb, i)
	}
	cat("\n")
	
	# Cycle back through the same set assigning xcoords, ycoods, locations and distances to the 
	reducedset <- matrix(nrow=nononezero,ncol=3,dimnames=list(1:nononezero,c("location","number","distance"))) 
	currentrow <- 1
	cat("Populating reduced set\n")
	pb <- txtProgressBar(min = 0, max = nox, style = 3)
	for (i in 1:nox) {
		for (j in 1:noy) {
			for (k in 1:nolocs) {
				if (allpossdistances[i,j,k] > 0) {
					reducedset[currentrow,"location"] <- locid[k]
					reducedset[currentrow,"distance"] <- allpossdistances[i,j,k]
					reducedset[currentrow,"number"] <- popimage$z[i,j]					
					currentrow <- currentrow + 1
				}
			}
		}
		setTxtProgressBar(pb, i)
	}
	
	# Return matrix of connectivity
	if (!is.null(csvfile)) {
		write.csv(reducedset,file=csvfile,row.names=FALSE)
	} 
	
	reducedset
	
}

# Function to take a dataframe of the distance from sites of interest to other populations and to
# return an effective neighbourhood size for each location
fsc.calc.neighbourhood <- function(vecLocs,dfconn,offset=10,power=5) {
	rtn <- vector(length=length(vecLocs),mode="numeric")
	maxloc <- max(vecLocs)
	lookup <- vector(length=maxloc,mode="numeric")
	for (i in 1:length(vecLocs)) lookup[vecLocs[i]] <- i
	norows <- dim(dfconn)[1]
	currentrow <- 1
	while (currentrow <= norows) {
		location <- dfconn[currentrow,"location"]
		rtnrow <- lookup[location]
		number <- dfconn[currentrow,"number"]
		distance <- dfconn[currentrow,"distance"]
		kernel <- 1 / ( 1 + ( distance / offset ) ^ power )
		rtn[rtnrow] <- rtn[rtnrow] + number * kernel
		currentrow <- currentrow + 1
	}
	rtn
}


# Function to generate a set of parameters for the serolevel ode model
srg.index.ode.params <- function() {
	rtn <- c(imax=2,alpha1=2,alpha2=3)
	rtn
}

# Function to test efficient indexing for ODE models
# The current version of this should make up a 2-state decay model 
srg.index.ode.mod <- function(t,y,p,vecAgeSusArg,blValidateArgs=FALSE) {
	indi <- 1
	imax <- p["imax"]
	rtn <- vector(length=length(y))
	while (indi <= imax) {
		rtn[indi] <- -y[indi]*vecAgeSusArg[indi]
		indi <- indi+1
	}
	list(rtn)
}

# Wrapper for the lsoda function 
# Note the need for atol=1e-40
srg.index.ode.solution <- function() {
	require(deSolve)
	p <- srg.index.ode.params()
	times <- 0:20
	ics <- c(100,200)
	vecAgeSus <- c(p["alpha1"],p["alpha2"])
	sol <- lsoda(ics,times,srg.index.ode.mod,p,vecAgeSusArg=vecAgeSus,atol=1e-40)
}

# Calculate time 
srg.ave.speed <- function(dist,disthr) {
	rtn <- 60* (dist / disthr)
	cat(floor(rtn),":",(rtn-floor(rtn))*60,"\n",sep="")
	rtn
}

# Parameter function for sero-level model
# Requires input functions to generate the matrix functions for the model
# For quantities such as relative infectivity, relative susceptibility and the like
# All of those are just arrays that are built by the vector of scalar parameters
# And the input functions taken here
# In the call to the ODE system, we will send sca as the parameters and the matrices as additional argumants
# This function also needs to make and solve the NGM to give the correct R0 - beta relationship
isl.params.funcs <- function() {
	
	# Defone a set of scalars. 
	# Needs a generic update parameter function here
	# Implement based from one of the other parameter functions above
	sca <- c(
			R0 = 10,
			Tg = 2.6,
			alpha = 1/365,
			maxi = 2,
			maxj = 1,
			maxk = 1,
			maxa = 5,
			maxX = 3,
			N = 1000000,
			seed = 100,
			trickle = 0
	)
	
	# Implement a full R0 calculation when possible
	sca["beta"] <- sca["R0"]/sca["Tg"]
	
	# Create 'lookup' arrays for each of the aux functions
	matM <- isl.make.m.null(sca)
	arrf <- isl.make.f.null(sca)
	arrg <- isl.make.g.simple(sca)
	arrh <- isl.make.h.simple(sca)
	
	# Return lookup arrays for state variables
	tmp <- isl.make.S.I.lookups(sca)
	arrSlu <- tmp$S
	arrIlu <- tmp$I
	arrCIlu <- tmp$CI
	
	sca["novars"] <- max(arrIlu)
	
	list(ps=sca,mat.M=matM,arr.f=arrf,arr.g=arrg,arr.h=arrh,arr.S.lu=arrSlu,arr.I.lu=arrIlu,arr.CI.lu=arrCIlu)
	
}

# Simple function to return a null mixing matrix
isl.make.m.null <- function(ps) {
	
	# Check to see that the parameter argumants are consistent with this version of the function
	# Not needed in this version. Left in because the identical syntax is right
	# if (!identical(as.numeric(ps["maxa"]),2)) stop("Function isl.make.m.2state designed for maxa=2")
	
	# Assign values. Need to update for correct polymopd data later
	rtn <- array(dim=c(ps["maxa"],ps["maxa"]))
	rtn[,] <- 1
	
	# Return result
	rtn
}

# Simple function to return a null infectivity array of the correct dimensons
isl.make.f.null <- function(ps) {
	
	# Define an array of the correct size
	rtn <- array(dim=c(ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"]))
	
	# Assign values to the array
	rtn[,,,] <- 1
	
	# Return the result
	rtn
	
}

# Function to return the most simple immune state boost:
# Current titre to infecting virus is raised by a single 2-fold doubling unless already at the max value
# All other titres left unaffected
# The indeces of this matrix in the methods are X,a,l,m,n,i,j,k
# l,m,n are the initial states and i,j,k are the final states
# Therefore, f(X,a,l,m,n,i,j,k) is the proportion of individuals infected by strain X in age group a
# Who start from state l,m,n and end up in state i,j,k
# Therefore, summing over all i,j,k for any l,m,n should give 1!
isl.make.g.simple <- function(ps) {
	
	# Declare an array to be returned and pray you don't run out of memory
	rtn <- array(0,dim=c(	ps["maxX"],ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"],
					ps["maxi"],ps["maxj"],ps["maxk"]))
	
	# Open all the variable loops
	for (X in 1:ps["maxX"]) {
		for (a in 1:ps["maxa"]) {
			for (l in 1:ps["maxi"]) {
				for (m in 1:ps["maxj"]) {
					for (n in 1:ps["maxk"]) {
						
						# Assign values
						if (identical(as.numeric(X),1)) {
							if (l < ps["maxi"]) rtn[X,a,l,m,n,l+1,m,n] <- 1
							else rtn[X,a,l,m,n,l,m,n] <- 1
						} else if (identical(as.numeric(X),2)) {
							if (m < ps["maxj"]) rtn[X,a,l,m,n,l,m+1,n] <- 1
							else rtn[X,a,l,m,n,l,m,n] <- 1
						} else if (identical(as.numeric(X),3)) {
							if (n < ps["maxk"]) rtn[X,a,l,m,n,l,m,n+1] <- 1
							else rtn[X,a,l,m,n,l,m,n] <- 1
						} else stop("problem with values of X in isl.make.g.simple")
						
						# Close n loop 
					}
					
					# Close m loop
				}
				
				# Close l loop
			}
			
			# Close a loop
		}
		
		# Close X loop
	}
	
	# Return the array
	rtn
	
}

# Takes a scaler set of parameters and returns an array with values for h
# h(X,a,i,j,k) is the susceptibility of individuals to strain X dependent on their
# age and their antibody levels. Values are referenced to h(X,1,1,1,1)=1, i.e.
# those in the youngest age group with no detectable titres always have a relative
# suceptibility of 1
isl.make.h.simple <- function(ps) {
	
	# Define the array to be returned, the default value is to be suscptible
	rtn <- array(1,dim=c(ps["maxX"],ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"]))
	
	# Cycle through the indices adding in protected individuals
	for (X in 1:ps["maxX"]) {
		for (a in 1:ps["maxa"]) {
			for (i in 1:ps["maxi"]) {
				for (j in 1:ps["maxj"]) {
					for (k in 1:ps["maxk"]) {
						
						# Assign values
						if (identical(as.numeric(X),1)) {
							if (i > 1) rtn[X,a,i,j,k] <- 0
						} else if (identical(as.numeric(X),2)) {
							if (j > 1) rtn[X,a,i,j,k] <- 0
						} else if (identical(as.numeric(X),3)) {
							if (k > 1) rtn[X,a,i,j,k] <- 0
						} else stop("Problem with values of X in isl.make.h.simple")
						
						# Close k loop 
					}
					
					# Close j loop
				}
				
				# Close i loop
			}
			
			# Close a loop
		}
		
		# Close X loop
	}
	
	# Return the array
	rtn
	
}

# Function to make a lookup array for the linearized index of the isl model
# Can also be used to make list of column names for the solution table
isl.make.S.I.lookups <- function(ps,givenames=FALSE) {
	
	# Declare the array and the counter
	rtnS <- array(dim=c(ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"]))
	rtnI <- array(dim=c(ps["maxX"],ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"]))
	rtnCumInf <- array(dim=c(ps["maxX"],ps["maxa"],ps["maxi"],ps["maxj"],ps["maxk"]))
	
	counter <- 1
	colnames <- c("time")
	
	# Cycle through the array and assign numbers
	for (a in 1:ps["maxa"]) {
		for (i in 1:ps["maxi"]) {
			for (j in 1:ps["maxj"]) {
				for (k in 1:ps["maxk"]) {
					
					rtnS[a,i,j,k] <- counter
					counter <- counter + 1
					colnames <- c(colnames,paste("S_a_",a,"_i_",i,"_j_",j,"_k_",k,sep=""))
					
					# Close k loop 
				}
				
				# Close j loop
			}
			
			# Close i loop
		}
		
		# Close a loop
	}
	
	# Cycle through the indices adding in protected individuals
	for (X in 1:ps["maxX"]) {
		for (a in 1:ps["maxa"]) {
			for (i in 1:ps["maxi"]) {
				for (j in 1:ps["maxj"]) {
					for (k in 1:ps["maxk"]) {
						
						rtnI[X,a,i,j,k] <- counter
						counter <- counter + 1
						colnames <- c(colnames,paste("I_X_",X,"_a_",a,"_i_",i,"_j_",j,"_k_",k,sep=""))
						
						# Close k loop 
					}
					
					# Close j loop
				}
				
				# Close i loop
			}
			
			# Close a loop
		}
		
		# Close X loop
	}
	
	# Cycle through the indices adding in protected individuals
	for (X in 1:ps["maxX"]) {
		for (a in 1:ps["maxa"]) {
			for (i in 1:ps["maxi"]) {
				for (j in 1:ps["maxj"]) {
					for (k in 1:ps["maxk"]) {
						
						rtnCumInf[X,a,i,j,k] <- counter
						counter <- counter + 1
						colnames <- c(colnames,paste("CInf_X_",X,"_a_",a,"_i_",i,"_j_",j,"_k_",k,sep=""))
						
						# Close k loop 
					}
					
					# Close j loop
				}
				
				# Close i loop
			}
			
			# Close a loop
		}
		
		# Close X loop
	}
	
	if (givenames) rtn <- colnames
	else rtn <- list(S=rtnS,I=rtnI,CI=rtnCumInf)
	
	rtn
	
}

# Takes parameters and a lookup array and returns the correct index for a model
isl.S.lu <- function(ps,lu,a,i,j,k) {
	if (min(c(a,i,j,k)) < 1) rtn <- 0 
	else if (a > ps["maxa"] || i > ps["maxi"] || j > ps["maxj"] || k > ps["maxk"]) rtn <- 0
	else rtn <- lu[a,i,j,k]
	rtn
}

# Takes parameters and a lookup array and returns the correct index for a model
# This one works for cuminf as well as for infectious states
isl.I.lu <- function(ps,lu,X,a,i,j,k) {
	if (min(c(X,a,i,j,k)) < 1) rtn <- 0 
	else if (X > ps["maxX"] || a > ps["maxa"] || i > ps["maxi"] || j > ps["maxj"] || k > ps["maxk"]) rtn <- 0
	else rtn <- lu[X,a,i,j,k]
	rtn
}

# Function to give population weights
isl.pop.weights.null <- function(ps) {
	rtn <- as.vector(rep(1,length=ps["maxa"]))
	rtn <- rtn / sum(rtn) 
	rtn 
}

isl.popweights.5.uk.2003 <- function(ps) {
	if (ps["maxa"]!=5) stop("max number of age classes must be 5")
	rtn <- as.vector(c(3009303,	8731641, 19889635,	12707659,	8453916))
	rtn <- rtn / sum(rtn) 
	rtn 
}

# Naive initial prevalences of different states
isl.initial.prev.null <- function(ps,a,i,j,k) {
	if (i==1 && j==1 && k==1) rtn <- 1
	else rtn <- 0
	rtn
} 

# Initial prevalences for a 5 state model as per Johnson et. al.
# Values estimated by eye or using paper text - not yet with data grab or equivalent
isl.initial.prev.johnson.2009 <- function(ps,a,i,j,k) {
	if (ps["maxa"]!=5) stop("this initialization function is for a 5 age state model")
	if (ps["maxi"]!=2) stop("this initialization function is for a 2 i state model")
	if (ps["maxj"]!=1) stop("this initialization function is for a 1 j state model")
	if (ps["maxk"]!=1) stop("this initialization function is for a 1 k state model")
	if (a==1) {
		prev <- 25/58
		if (i==2) rtn <- prev
		else rtn <- 1-prev
	} else if (a==2) {
		prev <- 0.25
		if (i==2) rtn <- prev
		else rtn <- 1-prev
	} else if (a==3) {
		prev <- 0.18
		if (i==2) rtn <- prev
		else rtn <- 1-prev
	} else if (a==4) {
		prev <- 0.20
		if (i==2) rtn <- prev
		else rtn <- 1-prev
	} else if (a==5) {
		prev <- 0.30
		if (i==2) rtn <- prev
		else rtn <- 1-prev
	} else stop("Problem in ils.initial.prev.johnson.2009")
}

# Setup the initial conditions with everyone susceptible
# Eventually, will want to develop this a little more with some individuals 
isl.make.ics <- function(ps,arrSlu,arrIlu,arrCIlu,fnage=isl.pop.weights.null,fnInitialPrevs=ils.initial.prev.null) {
	
	# Note that numeric vectors are initiated to 0
	rtn	 	<- vector(mode="numeric",length=max(arrCIlu))
	agevec 	<- fnage(ps)
	seedvec <- agevec * ps["seed"] / ps["N"]
	for (a in 1:ps["maxa"]) {
		
		# Put in some intitial sero states 
		for (i in 1:ps["maxi"]) for (j in 1:ps["maxj"]) for (k in 1:ps["maxk"]) rtn[isl.S.lu(ps,arrSlu,a,i,j,k)] <- (agevec[a] - seedvec[a])*fnInitialPrevs(ps,a,i,j,k)
		for (X in 1:ps["maxX"]) rtn[isl.I.lu(ps,arrIlu,X,a,1,1,1)] <- seedvec[a]
	}
	rtn
}

# Function to be called by lsoda to calculate the derivatives of the state
# variables using all the aux functions and arrays
# Initially, just loop thorugh all the variables and have them decay exponentially!
isl.model <- function(t,y,ps,matM,arrf,arrg,arrh,arrSlu,arrIlu,arrCIlu) {
	
	# Set up vector to be returned
	rtn <- as.vector(y)
	
	# Setup any required auxilliary variables
	# First calculate the force of ifnection per strain epr age group
	matFOI <- array(0,dim=c(ps["maxX"],ps["maxa"]))
	X <- 1
	while (X <= ps["maxX"]) {
		a <- 1
		while (a <= ps["maxa"]) {
			tmpval <- 0
			b <- 1
			while (b <= ps["maxa"]) {
				i <- 1
				while (i <= ps["maxi"]) {
					j <- 1
					while (j <= ps["maxj"]) {
						k <- 1
						while (k <= ps["maxk"]) {	
							tmpval <- tmpval + arrf[b,i,j,k] * y[isl.I.lu(ps,arrIlu,X,b,i,j,k)]
							k <- k+1
							# close k loop
						}
						j <- j + 1	
						# close j loop
					}
					i <- i + 1
					# close i loop
				}
				tmpval <- tmpval * matM[a,b]
				b <- b+1
				# Close the b loop	
			}
			matFOI[X,a] <- tmpval * ps["beta"]
			a <- a + 1
			# Close a loop
		}
		X <- X + 1
		# Close X loop
	}
	
	# Set up the variable for boostin terms
	vecInBoost <- vector(mode="numeric",length=ps["maxX"])
	vecOutBoost <- vector(mode="numeric",length=ps["maxX"])
	vecInf <- vector(mode="numeric",length=ps["maxX"])
	
	# Loop to cycle through the susceptible states
	a <- 1
	while (a <= ps["maxa"]) {
		i <- 1
		while (i <= ps["maxi"]) {
			j <- 1
			while (j <= ps["maxj"]) {
				k <- 1
				while (k <= ps["maxk"]) {
					
					# First X loop to calcuate the boosting values
					vecInBoost[] <- 0
					vecOutBoost[] <- 0
					X <- 1
					while (X <= ps["maxX"]) {	
						l <- 1
						while (l <= ps["maxi"]) {
							m <- 1
							while (m <= ps["maxj"]) {
								n <- 1
								while (n <= ps["maxk"]) {
									
									# Add the immune state changes
									vecInBoost[X] <- vecInBoost[X] + y[isl.I.lu(ps,arrIlu,X,a,l,m,n)] * arrg[X,a,l,m,n,i,j,k]
									vecOutBoost[X] <- vecOutBoost[X] + y[isl.I.lu(ps,arrIlu,X,a,i,j,k)] * arrg[X,a,i,j,k,l,m,n]
									
									# Close n loop
									n <- n+1
								}
								m <- m+1
								# Close m loop
							}
							l <- l +1
							# Close l loop
						}
						# Close the X loop
						X <- X + 1
					}
					
					# Start second X loop to calculate the state variables
					delta <- sum(min(i,2),min(j,2),min(k,2),-3)
					X <- 1
					while (X <= ps["maxX"]) {
						
						# Need to put the decay term for infection terms
						decayterm <- -delta * y[isl.I.lu(ps,arrIlu,X,a,i,j,k)]
						if (i < ps["maxi"]) decayterm <- decayterm + y[isl.I.lu(ps,arrIlu,X,a,i+1,j,k)]
						if (j < ps["maxj"]) decayterm <- decayterm + y[isl.I.lu(ps,arrIlu,X,a,i,j+1,k)]
						if (k < ps["maxk"]) decayterm <- decayterm + y[isl.I.lu(ps,arrIlu,X,a,i,j,k+1)]
						
						# Equation for the I states
						vecInf[X] <- y[isl.S.lu(ps,arrSlu,a,i,j,k)]*arrh[X,a,i,j,k]*(matFOI[X,a]+isl.gamma.basic(ps,X,t))
						rtn[isl.I.lu(ps,arrIlu,X,a,i,j,k)] <-  + vecInf[X] + ps["alpha"]*decayterm - 1/(ps["Tg"])*vecOutBoost[X]
						rtn[isl.I.lu(ps,arrCIlu,X,a,i,j,k)] <- + vecInf[X]
						
						# Close the second X loop
						X <- X + 1
						
					}					
					
					# Decay term for susceptibility term
					decayterm <- - delta * y[isl.S.lu(ps,arrSlu,a,i,j,k)]
					if (i < ps["maxi"]) decayterm <- decayterm +  y[isl.S.lu(ps,arrSlu,a,i+1,j,k)]
					if (j < ps["maxj"]) decayterm <- decayterm +  y[isl.S.lu(ps,arrSlu,a,i,j+1,k)]
					if (k < ps["maxk"]) decayterm <- decayterm +  y[isl.S.lu(ps,arrSlu,a,i,j,k+1)]
					
					# S equations will go here because they are not indexed with X
					rtn[isl.S.lu(ps,arrSlu,a,i,j,k)] <- - sum(vecInf) + ps["alpha"]*decayterm + 1/(ps["Tg"])*sum(vecInBoost)
					
					# Close the k loop
					k <- k+1
					
				}
				j <- j+1
			}
			i <- i+1	
		}
		a <- a+1
	}
	
	# Return a list with the derivatives as the first element
	list(rtn)
	
}

# Reset R session by removing all objects and loading in this file
# Note that because of the remove all command, the lcation of the 
# function file cannot be an argument of this function
srg.reset.session <- function() {
	rm(list=ls(all=TRUE))
	source(paste("~/Dropbox/svneclipse/idsource/","R/stevensRfunctions.R",sep=""))
}

# Seeding function for the isl project
isl.gamma.basic <- function(ps,X,t) {
	if (identical(as.numeric(X),1)) rtn <- ps["trickle"]/ps["N"]
	else rtn <- 0
	rtn
}

# Function to label output for isl project
isl.check.boost <- function(ps,arrg) {
	
	# Open all the variable loops
	for (X in 1:ps["maxX"]) {
		for (a in 1:ps["maxa"]) {
			for (i in 1:ps["maxi"]) {
				for (j in 1:ps["maxj"]) {
					for (k in 1:ps["maxk"]) {
						tmp <- 0
						for (idash in 1:ps["maxi"]) {
							for (jdash in 1:ps["maxj"]) {
								for (kdash in 1:ps["maxk"]) {
									tmp <- tmp + arrg[X,a,i,j,k,idash,jdash,kdash]
									# Close n loop
								}
								# Close m loop
							}
							# Close l loop
						}
						if (abs(tmp - 1) > 1e-10) stop("Problem is isl.check.boost")
						# Close k loop 
					}
					# Close j loop
				}
				# Close i loop
			}
			# Close a loop
		}
		# Close X loop
	}
}

# Example function for Rprof example
srg.rprofex.f <- function(n) {
	vecrand <- runif(n)
	rtn <- srg.rprofex.g(vecrand)
	rtn
}

# Example function for Rprof example
srg.rprofex.g <- function(vecvals) {
	rtn <- mean(vecvals)
	rtn
}

# Example to illustrate the ue of Rprof
srg.example.profile <- function(outdir = "~/Dropbox/tmp/",n=10000000,fnProf="rprofexample.out") {
	
	# Turn on profiling
	Rprof(filename=paste(outdir,fnProf,sep=""))
	
	# Call the example function
	x <- srg.rprofex.f(n)
	
	# Turn off profiling
	Rprof(NULL)
	
	# Print the summary
	summaryRprof(filename=paste(outdir,fnProf,sep=""))
	
}

# Function to plot the distribution of titres once the cumulative data have been calculated
# This is essentially a graded color stack chart
isl.plot.strain.tires <- function(sol,agevec,ps,arrSlu,arrIlu,plotdir="./",plotfile=NULL,h=20,w=10) {
	
	# Set up some standard auxillliary functions 
	times 		<- sol[,"time"]
	tmin 		<- min(times)
	tmax 		<- max(times)
	
	if (!is.null(plotfile)) pdf(file=paste(plotdir,plotfile,sep=""),height=h/cm(1),width=w/cm(1)) 
	
	plot(1:2,type="n",xlim=c(tmin,tmax),ylim=c(0,1))
	
	# Make the normalized age vector cumulative
	noagegroups <- length(agevec)
	cumagevec <- vector(mode="numeric",length=noagegroups)
	for (i in 2:noagegroups) cumagevec[i] <- cumagevec[i-1]+agevec[i-1]
	
	# start loop through age categories
	for (j in 1:noagegroups) {	
		pltmat <- isl.gen.strain.titre.data(sol,ps,arrSlu,arrIlu,agestates=j:j)
		dimsmat 	<- dim(pltmat)
		nostates 	<- dimsmat[1]
		notps 		<- dimsmat[2]
		if (length(times) != notps) stop ("Problem with times and pltmat in isl.plot.strain.tires")
		colscale <- heat.colors(nostates)
		polygon(c(times,rev(times)),cumagevec[j]+c(rep(0,notps),rev(pltmat[1,])),col=colscale[1],border=NA,new=TRUE)
		if (nostates > 1) for (i in 2:nostates) polygon(c(times,rev(times)),cumagevec[j]+c(pltmat[i-1,],rev(pltmat[i,])),col=colscale[i],border=NA)
	}
	
	if (!is.null(plotfile)) dev.off() 
	
}


# Function to plot titres from a single strain
# Really can't alter the titres now
isl.gen.strain.titre.data <- function(sol,ps,arrS,arrI,X=1,agestates=NULL,rowcumulative=TRUE) {
	nots <- dim(sol)[1]
	if (is.null(agestates)) agestates <- 1:ps["maxa"]
	if (X==1) {nostates <- ps["maxi"]; iflag <- 1; jflag <- 0; kflag <- 0}
	else if (X==2) {nostates <- ps["maxj"]; iflag <- 0; jflag <- 1; kflag <- 0}
	else if (X==3) {nostates <- ps["maxk"]; iflag <- 0; jflag <- 0; kflag <- 1}
	else stop("Error in isl.gen.strain.titre.data. X not in {1,2,3}")
	rtn <- matrix(ncol=nots,nrow=nostates)
	for (i in 1:nots) {
		for (j in 1:nostates) {
			val <- 0
			for (a in agestates) {
				if (iflag) {
					for (l in 1:ps["maxj"]) {
						for (m in 1:ps["maxk"]) {
							
							# Need the 1+ below because time is the first column of the results
							val <- val + sol[i,1+isl.S.lu(ps,arrS,a,j,l,m)]
							for (Xi in 1:ps["maxX"]) {
								val <- val + sol[i,1+isl.I.lu(ps,arrI,Xi,a,j,l,m)]
							}
						}
					}
				} else if (jflag) {
					stop("Jflag branch not yet implemented")
				} else if (kflag) {
					stop("Jflag branch not yet implemented")
				} else stop("Viable option selected.")
			}
			rtn[j,i] <- val
		}
	}
	if (rowcumulative && nostates > 1) {
		for (j in 2:nostates) rtn[j,] <- rtn[j-1,] + rtn[j,]
	}
	rtn
}

# Place holder function to put the 
isl.gen.first.results <- function(datadir="./",plotname=NULL) {
	require(deSolve)
	su 	<- isl.params.funcs()
	isl.check.boost(su$ps,su$arr.g)
	ics <- isl.make.ics(	su$ps,su$arr.S.lu,su$arr.I.lu,su$arr.CI.lu,
			fnage=isl.popweights.5.uk.2003,fnInitialPrevs=isl.initial.prev.johnson.2009)
	sol <- lsoda(ics,(0:26)*7,isl.model,su$ps,
			rtol=1e-20,
			matM=su$mat.M,arrf=su$arr.f,arrg=su$arr.g,arrh=su$arr.h,arrSlu=su$arr.S.lu,
			arrIlu=su$arr.I.lu,arrCIlu=su$arr.CI.lu)
	colnames(sol) <- isl.make.S.I.lookups(su$ps,givenames=TRUE)
	write.csv(sol,file=paste(datadir,"isl_debug_dump.csv",sep=""))
	isl.plot.strain.tires(sol,isl.popweights.5.uk.2003(su$ps),su$ps,su$arr.S.lu,su$arr.I.lu,datadir,plotname)
}

# Function to read in expected prevalences and then make output file of sample sizes
hpa.mrc.power <- function(infilestem="~/Dropbox/grants/mrc_flu_sero_2011/expected.prevalences") {
	x <- read.csv(paste(infilestem,".csv",sep=""))
	write.csv(x,file=paste(infilestem,"output_from_hpa_mrc_power",".csv",sep=""))
	0
}

# Function to take a baseline sample, expected baseline prevalence, expected new prevalence and desired accuracy.
# It will then find the smallest sample for the second year that gives at least the accuracy requested
# It also takes a maximum sample size and will return -1 if the required sample size is greater than the maximum
# or it will return -1 if there the base sample is so small for the requested difference 
hpa.mrc.acc2 <- function(na=300,pa=0.07,pb=0.55,nbmax=5000,accuracy=0.1) {
	require(Epi)
	tmp_res <- ci.pd(round(na*pa),round(nbmax*pb),na - round(na*pa), nbmax - round(nbmax*pb),print=FALSE)
	init_accuracy <- tmp_res[7] - tmp_res[6]
	if (init_accuracy > accuracy) {
		rtn <- nbmax
		accuracy <- init_accuracy
	} else {
		rtn <- nbmax
		continue <- TRUE
		while (continue) {
			rtn <- rtn - 1
			tmp_res <- ci.pd(round(na*pa),round(rtn*pb),na - round(na*pa), rtn - round(rtn*pb),print=FALSE)
			current_accuracy <- tmp_res[7] - tmp_res[6]
			if (current_accuracy > accuracy) continue <- FALSE
			if (rtn == 0) stop("Error in hpa.mrc.second.sample")
		}
	}
	list(nb=rtn,acc=accuracy)
}

hpa.mrc.acc1.vec <- function(vecsamples,vecprev,accuracy=0.1) {
	nomeasures <- length(vecsamples)
	if (nomeasures < 2) stop("This routine is designed for more than one year hpa.mrc.acc1.vec")
	vecrtn <- vector(mode="numeric",length=nomeasures-1)
	vecacc <- vector(mode="numeric",length=nomeasures-1)
	for (i in 2:nomeasures) {
		if (i==2) tmp <- hpa.mrc.acc2(vecsamples[1],vecprev[i-1],vecprev[i],vecsamples[i],accuracy)
		else tmp <- hpa.mrc.acc2(vecrtn[i-2],vecprev[i-1],vecprev[i],vecsamples[i],accuracy)
		vecrtn[i-1] <- tmp$nb
		vecacc[i-1] <- tmp$acc
	}
	list(nb=vecrtn,acc=vecacc)
}

# Calculates the geodesic distance between two points specified by
# radian latitude/longitude using the Haversine formula
# Taken from http://r.789695.n4.nabble.com/Geographic-distance-between-lat-long-points-in-R-td3442338.html
srg.geod.distance <- function(long1, lat1, long2, lat2) {
	R <- 6371 # Earth mean radius [km]
	delta.long <- (long2 - long1)
	delta.lat <- (lat2 - lat1)
	a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
	c <- 2 * asin(min(1,sqrt(a)))
	d = R * c
	return(d) # Distance in km
}

# Function to generate a one way hash of a vector of non-unique values
srg.one.way.hash <- function(vecDat,seed) {
	require("hash")
	uniquevals <- names(table(vecDat))
	set.seed(seed)
	hashtab <- hash(uniquevals,runif(length(uniquevals),min=1,max=1e10))
	rtn <- hashtab[vecDat]
	rtn
}

# Function to calculate distances between two postcodes for michelle
srg.calc.dists.for.michelle <- function(	
		datafile="~/Dropbox/tmp/unique_ID_and_postcode_site.csv",
		postcodefilename="~/Dropbox/tmp/postcodes.csv",
		outputfile="~/Dropbox/tmp/unique_ID_and_postcode_site_dist_added.csv",
		sitevar="site",
		pcvar="pc11") {
	
	# Define constants
	deg2rad <- 2 * pi / 360
	hlat1 <- 51.54 * deg2rad
	hlong1 <- -0.246 * deg2rad
	hlat2 <- 51.58 * deg2rad
	hlong2 <- -0.336 * deg2rad
	
	# Define the internal functions needed for the loop
	matchpostcode <- function(pc,i) {
		rowno <- match(pc,tabPc[,"postcode"])
		if (is.na(rowno)) stop("Can't find postcode in srg.calc.dists.for.michelle: ",pc,i)
		rowno
	}
	
	# Calc distance
	calcdist <- function(hosp,pc,i) {
		
		# Check for a null postcode
		if (pc=="NULL" || hosp=="NULL") {
			rtn <- -1
			
			# Otherwise assume that the postciode should work
		} else {
			
			# Assign lat and long of correct hospital
			if (hosp=="CMH") { 	
				hlat <- hlat1
				hlong <- hlong1
			} else if (hosp == "NPH") { 	
				hlat <- hlat2
				hlong <- hlong2 
			} else stop("Problem with hospital code in srg.calc.dists.for.michelle ",hosp,i)
			
			# Find the lat and long of the case
			pcrow <- matchpostcode(pc,i)
			clat <- tabPc[pcrow,"latitude"]
			clong <- tabPc[pcrow,"longitude"]
			
			rtn <- srg.geod.distance(clong*deg2rad,clat*deg2rad,hlong,hlat)
		}
		
		# Return either distance or error flag
		rtn
		
	}
	
	# Read in the postcodes and the data
	tabPc <- read.csv(postcodefilename)
	tabDat <- read.csv(datafile)
	noentries <- dim(tabDat)[1]
	
	# Cycle through the data
	tabDat <- cbind(tabDat,dist=rep(-1,noentries))
	i <- 1
	while (i <= noentries) {
		tmp <- calcdist(tabDat[i,sitevar],tabDat[i,pcvar],i)
		tabDat[i,"dist"] <- calcdist(tabDat[i,sitevar],tabDat[i,pcvar])
		i <- i + 1
	}
	
	# Write the new file
	write.csv(tabDat,file=outputfile)
	
}

# Routine to plot the sensitivity of 
apa.R0phidelta.sens <- function(
		nsamps=50,
		graphdir=	"~/Dropbox/projects/influenza/age_peak_attack/Current/figs/",
		outputfile=	"apa_lhs_R0_phi_delta.R",
		params= 	c("phi","delta"),
		R0vals = 	c(1.1,1.4,1.8,2.6),
		vecmins= 	c(0.1,0.1),
		vecmaxes= 	c(1.0,10),
		veclog= 	c(FALSE,TRUE)
) {
	
	# Define array to be saved
	rtn <- array(dim=c(nsamps,length(params)+3,length(R0vals)),dimnames=list(1:nsamps,c(params,"or","peak","CIAR"),1:length(R0vals)))
	
	# Define the parameters, ranges and scales
	for (i in 1:length(R0vals)) {
		rtn[,,i] <- apa.lhs(
				nsamps,
				params,
				vecmins,
				vecmaxes,
				veclog,
				constpnames=c("R0"),
				constpvals=c(R0vals[i]))
	}
	
	# Write the output file
	save(rtn,file=paste(graphdir,outputfile,sep=""))
	
}

# Make a plot of the R0 phi delta sensitivity
apa.plot.R0phidelta <- function(
		graphdir="~/Dropbox/projects/influenza/age_peak_attack/Current/figs/",
		lhsfile="apa_lhs_R0_phi_delta.R",
		chartfile="scatter_r0phidelta.pdf",
		width=8,
		height=10,
		cex=1.0) { 
	
	# Load up the data from the file
	load(paste(graphdir,lhsfile,sep=""))
	
	# Establish the jpeh chart file 
	pdf(file=paste(graphdir,chartfile,sep=""),width=width/cm(1),height=height/cm(1))
	
	# Make the plot
	par(	mai=(c(0,0,0,0)),
			fig=c(0.2,0.95,0.2,0.95),
			new=FALSE
	)
	plot(0:1,xlim=c(1,4),ylim=c(0,0.05),axes=FALSE,type="n")
	axis(1)
	axis(2)
	cols <- c("red","blue","green")
	for (i in 1:(dim(rtn)[3])) {
		points(rtn[,"or",i],rtn[,"peak",i],col=cols[i],pch=19,cex=cex)
	}
	
	# Turn off the device
	dev.off() 
	
}

# Routine to simulate a very quick epidemic in a small group, with infection, symptoms and deaths
# This routine assumes that the q
vqe.main <- function(
		rng.seed=1234,
		seed=10,
		dmax=180,
		mu_inc=12,
		alpha_inc=1000,
		mu_dinf=2,
		alpha_dinf=5,
		p_inf=0.01,
		p_fever=0.9,
		offset_fever=4,
		mu_fever=4,
		mu_fever_dur=8,
		p_severe=0.075,
		mu_severe_start=5,
		mu_severe_dur=20,
		p_death=0.15,
		Npop=30,
		charfile="./characteristics.csv",
		outfile="./sim_epi.csv",
		hidetinf=TRUE) {
	
	# Reset the random number seed
	set.seed(rng.seed)
	
	# Read in the characteristics and add the extra columns
	# rtn <- read.csv(charfile)
	rtn <- data.frame(id=1:Npop)
	nopeople <- dim(rtn)[1]
	
	# Test some parameter values
	if (seed >= nopeople) stop("See size must be less that size of the populations")
	if (seed < 1) stop("See size must be at least 1 in vqe.main")
	if (mu_inc < 1 || mu_dinf < 1 || mu_fever_dur < 1) stop("Means of waiting time distributions must be greater than 1")
	
	# Set up the return data frame
	rtn <- cbind(rtn,
			t.inf=rep(-1,nopeople),
			t.shed.start=rep(-1,nopeople),
			t.shed.finish=rep(-1,nopeople),
			t.fever.start=rep(-1,nopeople),
			t.fever.stop=rep(-1,nopeople),
			t.severe.start=rep(-1,nopeople),
			t.severe.stop=rep(-1,nopeople),
			b.death=rep(0,nopeople),
			index.infector=rep(0,nopeople),
			generation=-1
	)
	
	vecInfectees <- c()
	while (length(vecInfectees) < seed) {
		choice <- round(runif(1,min=0.5,max=nopeople+0.5))
		if (!(choice %in% vecInfectees)) vecInfectees <- c(vecInfectees,choice)
	}
	vecInfectors <- rep(-1,length(vecInfectees))
	
	# Set preconditions for the loop
	day <- 0
	
	while (day <= dmax) {
		
		# Process pending infections
		if (length(vecInfectees)>0) {
			for (i in 1:length(vecInfectees)) {
				
				# Establish the consequences of infection
				# Time of infection
				rtn[vecInfectees[i],"t.inf"] <- day
				
				# Infector
				rtn[vecInfectees[i],"index.infector"] <- vecInfectors[i]
				if (vecInfectors[i]<0) {rtn[vecInfectees[i],"generation"] <- 1}
				else {rtn[vecInfectees[i],"generation"] <- rtn[vecInfectors[i],"generation"] + 1}
				
				# Generate the start of shedding and its finish
				start <- day+1+floor(srg.rgamma(1,mu_inc,alpha_inc))
				finish <- start+1+floor(srg.rgamma(1,mu_dinf,alpha_dinf))
				rtn[vecInfectees[i],"t.shed.start"] <- start
				rtn[vecInfectees[i],"t.shed.finish"] <- finish
				
				# Decide who had a fever and when it occurred
				# For interest, make the timing of the fever be relative to the start if infectiousness
				if (runif(1) < p_fever) {
					rtn[vecInfectees[i],"t.fever.start"] <- rtn[vecInfectees[i],"t.shed.start"] - offset_fever + rpois(1,mu_fever)  
					rtn[vecInfectees[i],"t.fever.stop"] <- rtn[vecInfectees[i],"t.fever.start"] + rpois(1,mu_fever_dur-1)
					
					# Decide on severe cases
					if (runif(1) < p_severe) {
						
						rtn[vecInfectees[i],"t.severe.start"] <- rtn[vecInfectees[i],"t.fever.start"]+ 1 + rpois(1,mu_severe_start-1)  
						rtn[vecInfectees[i],"t.severe.stop"] <- rtn[vecInfectees[i],"t.severe.start"] + 1 + rpois(1,mu_severe_dur-1)
						
						# Decide on death of severe cases
						if (runif(1) < p_death) {
							rtn[vecInfectees[i],"b.death"] <- 1
						}
						
					}
					
				}
				
			}
			
			vecInfectees <- c()
			vecInfectors <- c()
			
		}
		# Generate the next set of infections
		# Scroll through looking for infectious individuals
		for (i in 1:nopeople) {
			if (day >= rtn[i,"t.shed.start"] && day < rtn[i,"t.shed.finish"] ) {
				for (j in 1:nopeople) {
					if (runif(1) < p_inf && rtn[j,"t.inf"] < 0) {
						if (!(j %in% vecInfectees)) {
							vecInfectees <- c(vecInfectees,j)
							vecInfectors <- c(vecInfectors,i)
						}
					}
				}
			}
		}
		
		day <- day + 1
		
		# Close the time loop 
	}
	
	# Supress infection times in output
	if (hidetinf) rtn[,"t.inf"] <- -1
	
	# Write the output file and return the results
	write.csv(rtn,file=outfile,row.names=FALSE)
	rtn
	
}

# SR's preferred parameterization of the gamma distribution
srg.rgamma <- function(n,mu,alpha) {
	shape <- alpha
	scale <- mu / alpha
	rtn <- rgamma(n,alpha,scale=scale)
	rtn
}

# Function to load up anonymized hospital attendence data
pha.load.data <- function(filename="/Volumes/NO\ NAME/data/hospital/2011_10_28\ anonymised\ data.csv") {
	require(date)
	raw <- read.csv(filename)	
	# reduced <- data.frame(hospnoanon=raw$hospnoanon,age=raw$age,date=as.date(as.character(raw$adate),order="dmy")+(as.date("01/01/2009")-as.date("01/01/1909")),record=raw$record)
	reduced <- data.frame(hospnoanon=raw$hospnoanon,age=raw$age,record=raw$record)
	tmp <- reduced[1:1000,]
	newdftmp <- reshape(tmp,idvar="hospnoanon",v.names=c("age"),timevar="record",direction="wide")
	testtab <- table(tmp$hospnoanon)
	newdf
}

# Routine to put a distance lookup file in a temp directory
# Needs a population image as an argument
# Think the problem is here, but can't find it on my local copy yet
# rm(list=ls(all=TRUE))
# source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# fsc.make.dist.lookup()
fsc.make.dist.lookup <- function(
		popimage=fsc.load.pop.image(),
		fluscapedir="~/Dropbox/svneclipse/fluscape/",
		outdir="~/Dropbox/tmp/",
		locfile="Location_Density_bad.csv",
		outfile = "popdistandsize.csv",
		rmax=30,
		locSub = 1:42) {
	
	# Load up required packages and source files
	source(paste(fluscapedir,"source/R/GeneralUtility.r",sep=""))
	connfile 	<- paste(outdir,outfile,sep="")
	dens 		<- read.csv(paste(fluscapedir,"data/",locfile,sep=""))
	
	# Load up the landscan or the grump data and calculate connecivity for a simple gravity-style model
	# Could we set up an interpolation of the connectivity of each site based on the power law
	rtn <- fsc.gen.connect.dists(dens$LOC_Lat[locSub],dens$LOC_Long[locSub],dens$LOC_ID[locSub],popimage=popimage,r_max=rmax,csvfile=connfile)
	
	rtn
	
}

# Calculate neighbourhood size and correlation between them
fsc.neighbour.corr <- function(vecTheta,dfcon,tab) {
	
	vecNeigh <- fsc.calc.neighbourhood(tab$LOC_ID,dfcon,offset=vecTheta[1],power=vecTheta[2])
	plot(log(vecNeigh),tab$coeff,main=paste("Off ",vecTheta[1],", pow ",vecTheta[2]))
	list(cor(log(vecNeigh),tab$coeff))
	
}


fsc.sandbox <- function(
		fluscapedir="~/Dropbox/svneclipse/fluscape/",
		outdir="~/Dropbox/tmp/",
		connfile = "~/Dropbox/tmp/popdistandsize.csv") {
	
	# Load up required packages and source files
	require(sp)
	source(paste(fluscapedir,"source/R/GeneralUtility.r",sep=""))
	
	# Load the data, density info and merge them
	fsd 	<- load.and.merge.part.V1(topdir=fluscapedir)
	dens 	<- read.csv(paste(fluscapedir,"data/Location_Density_bad.csv",sep=""))
	
	
}

fsc.neighbour.sense <- function(
		vecPowers = c(2,4,6,8,10),
		vecOffsets = c(1,5,10,15),
		fluscapedir="~/Dropbox/svneclipse/fluscape/",
		outdir="~/Dropbox/tmp/",
		connfile = "~/Dropbox/tmp/popdistandsize.csv",
		datafile = "manuscripts/main.outcome/loc.coeff.subset.csv") {
	
	# Load up the data file of lat longs and model coefficients and community connectivities
	coeffs 	<- read.csv(paste(fluscapedir,datafile,sep=""))	
	dfConnect <- read.csv(connfile)
	
	# Need a function to calculate the correlation at each point in a grid 
	cors 		<- matrix(nrow=length(vecOffsets),ncol=length(vecPowers),dimnames=list(vecOffsets,vecPowers))
	
	for (i in 1:length(vecOffsets)) {
		for (j in 1:length(vecPowers)) {
			tmp <- fsc.neighbour.corr(c(vecOffsets[i],vecPowers[j]),dfConnect,coeffs)
			cors[i,j] <- tmp[[1]]
			cat("Completed",j+(i-1)*length(vecPowers),"of",length(vecPowers)*length(vecOffsets),"\r")
			flush.console()
		}
	}
	
	# Write output
	write.csv(cors,file=paste(outdir,"fsc.ns",".out1.",srg.file.tag(),".csv",sep=""))
	
}

# Returns the current date and time in a format to be added to a file output
srg.file.tag <- function() {
	tmp <- date()
	tmp <- gsub(" ","-",tmp)
	tmp <- gsub(":","-",tmp)
	tmp <- paste((Sys.info())["nodename"],tmp,sep="_")
	tmp
}

fsc.compare.95.08 <- function() {
	
	# Set up some directories for the comparison
	fluscapedir <- "~/Dropbox/svneclipse/fluscape/"
	outdir <- "~/Dropbox/tmp/"
	connfile <- "~/Dropbox/tmp/popdistandsize.csv"
	
	require(sp)
	source(paste(fluscapedir,"source/R/GeneralUtility.r",sep=""))
	
	# Load the data, density info and merge them
	fsd_tmp <- load.and.merge.part.V1(topdir=fluscapedir)
	fsd_tmp <- fsd_tmp[fsd_tmp$PART_PROV_SAMPLE==1,]
	fsd 	<- fsd_tmp[is.na(fsd_tmp$PART_AGE)==FALSE,]
	fsd 	<- cbind(fsd,log2_95=sapply(fsd$HI.H3N2.1995,function(x){ifelse(x<9,0,log2(x/5))}))
	fsd		<- cbind(fsd,log2_08=sapply(fsd$HI.H3N2.2008,function(x){ifelse(x<9,0,log2(x/5))}))
	fsd		<- cbind(fsd,log2_79=sapply(fsd$HI.H3N2.1979,function(x){ifelse(x<9,0,log2(x/5))}))
	fsd 	<- cbind(fsd,age_cat=sapply(fsd$PART_AGE,function(x){findInterval(x,(0:10)*10)}))
	
	# xy plot by age group
	pdf(file="~/Dropbox/tmp/Rout1.pdf",width=25/cm(1),height=15/cm(1))
	layout(matrix(1:10,ncol=5,byrow=TRUE))
	for (i in 1:10) {
		ageindices <- fsd$age_cat == i
		plot(jitter(fsd$log2_08[ageindices]),jitter(fsd$log2_95[ageindices]),main=paste("Age group ",i),xlim=c(-1,8),ylim=c(-1,8))
	}
	dev.off()
	
	# Jitterplot by age
	pdf(file="~/Dropbox/tmp/Rout2.pdf",width=25/cm(1),height=15/cm(1))
	layout(matrix(1:3,ncol=3,byrow=TRUE))
	plot(fsd$PART_AGE,jitter(fsd$log2_79),main="1979-like")
	plot(fsd$PART_AGE,jitter(fsd$log2_95),main="1995-like")
	plot(fsd$PART_AGE,jitter(fsd$log2_08),main="2008-like")
	dev.off()
	
	# Box plot for the different strains
	pdf(file="~/Dropbox/tmp/Rout3.pdf",width=25/cm(1),height=15/cm(1))
	layout(matrix(1:3,ncol=3,byrow=TRUE))
	boxplot(log2_79 ~ age_cat,data=fsd,main="1979-like")
	boxplot(log2_95 ~ age_cat,data=fsd,main="1995-like")
	boxplot(log2_08 ~ age_cat,data=fsd,main="2008-like")
	dev.off()
	
	# Histogram of the different ages
	pdf(file="~/Dropbox/tmp/Rout4.pdf",width=25/cm(1),height=15/cm(1))
	hist(fsd$PART_AGE,breaks=(0:10)*10)
	dev.off()
	
}

# Function to load up the Hong Kong Contact survey basic data
hkc.load.data.and.add <- function(filename="~/Dropbox/shares/on_contact/data_and_code/task_in_ic1/data762.csv") {
	rtn <- read.csv(filename)
	rtn <- cbind(rtn,ag1=hkc.add.agegroup(rtn))
	rtn
}

hkc.add.agegroup <- function(dat,thresh=c(6,20,64)) {
	thresh <- c(0,thresh,999)
	rtn <- findInterval(dat$age,thresh,rightmost.closed=TRUE)
	rtn
}

hkc.gen.rev.cum <- function(vecX) {
	maxX <- max(vecX)
	sortX <- sort(vecX)
	rtn <- vector(type="numeric",length=maxX)
	currentval <- maxX
	for (i in 1:maxX) {
		if (sortX[i] < i) currentval <- currentval
		rtn[i] <- currentval
		
	} 
}

# Generic function to make extra charts and functions of interest for the HK contacts study
hkc.sandbox <- function(fn="/Volumes/NO NAME/data/influenza/hk_serosurvey/contacts_final_762.csv") {
	
	# Next: make a scatter plot of the data
	dat <- read.csv(fn)
	plot(dat$age,dat$contact_min,log="y")
	table(dat$contact_min)
	
}

srg.chart.pos <- function(	xindex,yindex,xn,yn,xlm=0,xrm=0,
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

sps.event.image <- function(data,popgrid,ev,st,et,sr,er) {
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

sps.plot.snapshots <- function( 
		fnPopdata = "~/Dropbox/svneclipse/spatialsim/spatialsim/data/prd_sim_ascii_n_44384655.txt",
		fnEpidata = "~/Dropbox/svneclipse/spatialsim/spatialsim/runtemplates/rfcid_rep/output/flusars_2_pset_1_Events.out",
		fnOutfile = "~/Dropbox/tmp/sev_flu_cont",
		chttimes = c(2,4,6,8,10,12,14,16,18,20,22)) {
	
	# Declare required libraries
	require("sp")
	require("fields")
	
	# Read in the ascii file
	popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
	
	# Read in the epi data as per my event file format
	data <- read.table(file=fnEpidata,header=TRUE)
	data$X <- data$X * 180 / pi
	data$Y <- data$Y * 180 / pi
	
	winwidth <- 20/cm(1)
	winheight <- 10/cm(1)
	resolution=300
	popaxis <- grey(0.9-0.6*0:100/100)
	infaxis <- tim.colors(256)
	ximax=4
	yimax=2
	chtlet <- c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)","p","q")
	
	
	# Start the pdf file
	pdf(file=paste(fnOutfile,".pdf",sep=""),width=winwidth,height=winheight)
	
	par(cex=0.8)
	par(mai=c(0/cm(1),0/cm(1),0/cm(1),0/cm(1)))
	
	for (yi in yimax:1) {
		for (xi in 1:ximax) {	
			if (xi==1 && yi==yimax) blnew=FALSE 
			else blnew=TRUE
			cur_pos <- srg.chart.pos(xi,yi,ximax,yimax,xg=0.05/winwidth,yg=0.05/winheight) 
			par(fig=cur_pos,new=blnew)
			cur_index <- xi+(yimax-yi)*ximax
			epiImage <- sps.event.image(data,popgrid,0,0,chttimes[cur_index],0,0)
			image(popgrid,col=popaxis,axes=FALSE)
			image(epiImage,col=infaxis,add=TRUE)
			mtext(paste(chtlet[cur_index],chttimes[cur_index],"days"),side=3,line=-1)
		}	
	}
	
	dev.off()
	
}

# Function to decipher a specific format of sequence coding and to look at the distances between the pairs
hhe.dist.and.char <- function(dnaseq,writedist=FALSE,outfile="~/Dropbox/tmp/hhe_dist_and_char.csv") {
	
	require("ape")
	
	# Define some housekeeping variables
	intSeqLength <- length(dnaseq[[1]])
	intNoSequences <- length(dnaseq)
	allDist <- dist.dna(dnaseq,model="N",pairwise.deletion=FALSE,as.matrix=TRUE)
	# hist(allDist)
	
	# Set up a dataframe for all the sequences
	# For each we should know: Name | Current Study (0|1) | Household | Member | Visit
	dfSeqChars <- as.matrix(cbind(	
					current_study=rep(0,intNoSequences),
					managua_study=rep(0,intNoSequences),
					household=rep(0,intNoSequences),
					member=rep(0,intNoSequences),
					visit=rep(0,intNoSequences))
	)
	
	for (i in 1:intNoSequences) {
		
		# Break the string
		strBits <- (strsplit(names(dnaseq[i]),"/"))[[1]]
		
		# Check for membership of studies
		if (length(strBits) > 1 && strBits[2] == "Hong_Kong") dfSeqChars[i,"current_study"] <- 1
		if (length(strBits) > 1 && strBits[2] == "Managua") dfSeqChars[i,"managua_study"] <- 1
		
		# Assign additional details for current study
		if (dfSeqChars[i,"current_study"] > 0) {
			
			# Split the
			if (length(strBits) < 3) stop("Inconsistancy in setUpDetails: 92837429")
			details <- (strsplit(strBits[3],"m|v"))[[1]]
			if (length(details) != 3) stop("Inconsistancy in setUpDetails: 83737429")
			dfSeqChars[i,"household"] <- as.numeric(details[1])
			dfSeqChars[i,"member"] <- as.numeric(details[2])
			dfSeqChars[i,"visit"] <- as.numeric(details[3])
			
		}
		
	}
	
	if (writedist) write.table(allDist,file=outfile,row.names=names(dnaseq),col.names=names(dnaseq),sep=",")
	
	list(lengthdna=intSeqLength,noseq=intNoSequences,distances=allDist,chars=dfSeqChars)
	
}

hhe.make.hists <- function(binarySD,fnOutput="~/Dropbox/tmp/hhe_make_hists_1.csv") {
	
	# Read off the length of the sequences and calculate the distances between each sequence
	x <- hhe.dist.and.char(binarySD)
	
	# Construct the frequency tables
	max_dist <- max(x$distances)
	tabMain <- 	as.matrix(cbind(	
					host=rep(0,max_dist+1),
					house=rep(0,max_dist+1),
					study=rep(0,max_dist+1),
					study_to_out=rep(0,max_dist+1),
					outgroup=rep(0,max_dist+1)
			))
	
	debugcount <- 0
	
	# Cycle though all distinct pair combinations and assign distances
	for (i in 1:(x$noseq-1)) {
		
		# Second member of the pair
		for (j in (i+1):x$noseq) {
			
			debugcount <- debugcount + 1
			
			# Distance between pairs
			intCurrentDistance <- x$distances[i,j] + 1
			
			# Test that i and j are members of the current study
			if (	x$chars[i,"current_study"] == 1 && x$chars[j,"current_study"] == 1) {
				
				# Test for the same household
				if (x$chars[i,"household"] == x$chars[j,"household"]) {
					
					# Test for the same individual
					if (x$chars[i,"member"] == x$chars[j,"member"]) {
						
						# Increment distance measure for within host
						tabMain[intCurrentDistance,"host"] <- tabMain[intCurrentDistance,"host"] + 1 
						
						# Close if for same individual open else	
					} else {
						
						# Increment distance measure for within household
						tabMain[intCurrentDistance,"house"] <- tabMain[intCurrentDistance,"house"] + 1 
						
						# Close else for same individual	
					}
					
					
					# Close the if for same household	
				} else {
					
					# Increment distance measure for within study
					tabMain[intCurrentDistance,"study"] <- tabMain[intCurrentDistance,"study"] + 1 
					
					# Close else for same household	
				}
				
				# Close the if for both current study and open else	
			} else {
				
				# Test if either were in the current study
				if (x$chars[i,"current_study"] == 1 || x$chars[j,"current_study"] == 1) {
					
					# Incremenet the between study count
					tabMain[intCurrentDistance,"study_to_out"] <- tabMain[intCurrentDistance,"study_to_out"] + 1
					
					# Endif of either being in the current study
				} else {
					
					# This case must be that neither are in the current study
					tabMain[intCurrentDistance,"outgroup"] <- tabMain[intCurrentDistance,"outgroup"] + 1
					
					# End else that either are in the current study
				}
				
			}
			
			# Close the j loop
		}
		
		# Close the for i loop
	}
	
	if ((x$noseq * x$noseq - x$noseq) / 2 != sum(tabMain)) stop("Problem with the number of pairs")
	write.csv(tabMain,file=fnOutput)
	
}

hhe.pairwise.compare <- function(seq1,seq2) {
	nosites <- length(seq1)
	if (length(seq2) != nosites) stop("hhe.pairwise.compare: problem")
}

hhe.sandbox <- function() {
	
	require("ape")
	fn2007AllFullTimed <- "~/Dropbox/projects/influenza/dna_household/anon_data/H3N2_all_concat_154.nex"
	rawSeqData <- read.nexus.data(fn2007AllFullTimed)
	binSeqData <- as.DNAbin(rawSeqData)
	x <- hhe.dist.and.char(binSeqData,writedist=TRUE)
	hhe.make.hists(binSeqData)
	dist.dna(binSeqData[1:2],model="N")
	
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
	
	# Load up libraries
	require("deSolve")
	
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
# Next edit to allow the recording of all timepoints for all runs
abp.hyper <- function(	nohypersamples=1,
		hyperparams = data.frame(params=c("gamma_R","dur_seas"),log=c(0,0),lb=c(1/20,7),ub=c(1/1,7*26)),
		mintakenof = 5,
		noyears = 30,
		lfname = "./abp_hyper_log.txt",
		logfile=TRUE
) {	
	
	# Include required libraries
	require("odesolve")
	
	# Load up default parameters
	ill_params <- abp.params.baseline()								
	
	# Set up hypervector parameter search
	nohyperparams <- dim(hyperparams)[1]
	
	# Set up the hypersquare
	dimnamesnext <- list(1:nohypersamples,hyperparams[,1])
	hypersquare <- matrix(0,ncol=nohyperparams,nrow=nohypersamples,dimnames=dimnamesnext)
	for (i in 1:nohyperparams) {
		hypersquare[,i] <- srg.hyper.vector(nohypersamples,hyperparams$lb[i],hyperparams$ub[i],log=hyperparams$log[i])
	}
	
	# This is the 
	if (logfile) cat("Parameter sweep started on",date(),"on",(Sys.info())["nodename"],"\n",file=lfname,append=FALSE)
	
	dimnamesnext <- list(1:nohypersamples,paste("min_yr_",(noyears-(mintakenof-1)):noyears,sep=""))
	
	matEachAnnMin <- matrix(-1,nrow=nohypersamples,ncol=mintakenof,dimnames=dimnamesnext)
	
	for (i in 1:nohypersamples) {
		sol_tmp <- (SirModel3AgeClasses(
							pname=as.character(hyperparams$params),
							pvals=hypersquare[i,],
							casemix=c(100,40,40),
							vp=ill_params
					))$sol
		minvec <- abp.annualmin(sol_tmp,mintakenof)
		matEachAnnMin[i,] <- minvec[]
		if (logfile) cat("Completed",i,"of",nohypersamples,"\n",file=lfname,append=TRUE)
	}
	
	list(matmin=matEachAnnMin,hypersquare=hypersquare)	
	
}

# Function under development to generate a list of value for the 
abp.dev.main <- function(outdir = "~/Dropbox/tmp") {
	dfParams = data.frame(params=c("gamma_R","dur_seas"),log=c(0,0),lb=c(1/20,7),ub=c(1/1,7*26))
	metric <- system.time(tmp <- abp.hyper(nohypersamples=5,dfParams))
	matmin 		<- tmp$matmin
	hypersquare <- tmp$hypersquare
	nosamples 	<- dim(matmin)[1]
}

# Function to simulate multistrain system for da-wei
twf.simulation <- function() {
	strFLU <- function(
			R0_1=1.8,
			R0_2=1.8,
			p_1=0,
			p_2=0,
			alpha_12=0.5,
			t_1=0,
			t_2=50,
			D_I=2.6,
			N=10000,
			D_S=0.25,
			tps=(0:400)/4,
			reals=1,
			dt=0.25,
			stochastic=TRUE,
			plot=TRUE) {
		
		noTPts <- length(tps) 
		epsilon <- 1e-10
		
		if (noTPts < 3) 
			stop("simSEIR: tps must be of length 2 or greater")
		if (reals > 1 && stochastic==FALSE) 
			stop("simSEIR: why make multiple realisations of a deterministic model?")
		
		rtn <- array(0,c(noTPts-1,2))
		
		for (i in 1:reals) {
			
			t_cur <- tps[1]
			ind_t <- 2
			t_next <- tps[ind_t]
			
			S_0 <- round(N*(1-p_1-p_2))
			S_1 <- round(N*p_1)
			S_2 <- N-S_0-S_1
			I_1 <- 0
			I_2 <- 0
			
			while (ind_t <= noTPts) {
				
				if (t_cur > t_1 - epsilon) hSed_1 <- 1/D_S/N 
				else hSed_1 <- 0
				if (t_cur > t_2 - epsilon) hSed_2 <- 1/D_S/N 
				else hSed_2 <- 0
				if (hSed_1 > 0 && hSed_2 > 0) {
					hSed_1 <- hSed_1 / 2
					hSed_2 <- hSed_2 / 2
				}
				
				lambda_1 <- R0_1*I_1/D_I/N + hSed_1
				lambda_2 <- R0_2*I_2/D_I/N + hSed_2
				
				pInf_0 	<- 1 - exp(-(lambda_1+lambda_2)*dt)
				pInf_1 	<- 1 - exp(-lambda_1*dt)
				pInf_2 	<- 1 - exp(-lambda_2*dt)
				pRec 	<- 1 - exp(-dt/D_I)
				
				if (stochastic) {
					if (lambda_1 + lambda_2 > 0) {					
						nInf_0 	<- rbinom(1,S_0,pInf_0)
						nInf_01 <- rbinom(1,nInf_0,lambda_1/(lambda_1+lambda_2))
						nInf_02 <- rbinom(1,nInf_0,lambda_2/(lambda_1+lambda_2))
					} else {
						nInf_01 <- 0
						nInf_02 <- 0
					}
					nInf_21 <- rbinom(1,S_2,pInf_1*(1-alpha_12))
					nInf_12 <- rbinom(1,S_1,pInf_2*(1-alpha_12))
					nRec_1 	<- rbinom(1,I_1,pRec)
					nRec_2 	<- rbinom(1,I_2,pRec)
				} else {
					if (lambda_1 + lambda_2 > 0) {					
						nInf_0 	<- S_0*pInf_0
						nInf_01 <- nInf_0*lambda_1/(lambda_1+lambda_2)
						nInf_02 <- nInf_0*lambda_2/(lambda_1+lambda_2)
					} else {
						nInf_01 <- 0
						nInf_02 <- 0
					}
					nInf_21 <- S_2*pInf_1*(1-alpha_12)
					nInf_12 <- S_1*pInf_2*(1-alpha_12)
					nRec_1 	<- I_1*pRec
					nRec_2 	<- I_2*pRec
				}
				
				S_0 <- S_0 - nInf_01 - nInf_02 
				S_1 <- S_1 + nRec_1  - nInf_12
				S_2 <- S_2 + nRec_2  - nInf_21
				I_1 <- I_1 + nInf_01 + nInf_21 - nRec_1
				I_2 <- I_2 + nInf_02 + nInf_12 - nRec_2
				
				rtn[ind_t-1,1] <- rtn[ind_t-1,1] + nInf_01 + nInf_21 
				rtn[ind_t-1,2] <- rtn[ind_t-1,2] + nInf_02 + nInf_12 
				
				t_cur <- t_cur + dt
				
				if (t_cur > (t_next - epsilon)) {
					t_next <- tps[ind_t]
					ind_t <- ind_t+1
				}
				
			}
		}
		
		rtn <- rtn / reals
		
		if (plot) {
			plot(rtn[,1],type="l",col="red",ylim=c(0,max(rtn[,1],rtn[,2])))
			points(rtn[,2],type="l",col="green")
		}
		
		list(ave=rtn)
		
	}
	
	approxPoissLike <- function(strData,strModel) {
		nstr <- dim(strData)[2]
		lnlike <- 0
		for (i in 1:nstr) {
			lnlike <- lnlike + sum(dpois(strData[,i],strModel[,i],log=TRUE))
		}
		lnlike
	}
	
	x <- strFLU(stochastic=TRUE,R0_2=1.8,plot=FALSE)
	param_range <- 1:150/100+1
	like_range <- rep(0,length(param_range))
	for (i in 1:length(param_range)) {
		y <- strFLU(stochastic=FALSE,R0_2=param_range[i],plot=FALSE)
		like_range[i] <- approxPoissLike(x$ave,y$ave)
	}
	plot(param_range,like_range,type="l")
}

psi.sim.modelA <- function(
		R0 = 1.4,
		Tg=2.6,
		t_0=20,
		seed=1,
		N=10000,
		tps=(0:52)*7,
		reals=10,
		dt=0.25,
		stochastic=TRUE,
		plot=TRUE) {
	
	noTPts <- length(tps) 
	epsilon <- 1e-10
	if (seed >= N) stop("psi.sim.modelA seed greater than N")
	
	if (noTPts < 3) 
		stop("simSEIR: tps must be of length 2 or greater")
	if (reals > 1 && stochastic==FALSE) 
		stop("simSEIR: why make multiple realisations of a deterministic model?")
	
	rtn <- array(0,c(noTPts-1,reals))
	
	for (i in 1:reals) {
		
		t_cur <- tps[1]
		ind_t <- 2
		t_next <- tps[ind_t]
		
		sS <- N
		sI <- 0
		sR <- 0
		
		blNotYetSeed <- TRUE
		
		while (ind_t <= noTPts) {
			
			if (t_cur >= t_0 && blNotYetSeed) {
				seedInf <- seed
				blNotYetSeed <- FALSE
			} else {
				seedInf <- 0
			}
			
			lambda <- R0 * sI / Tg / N
			pInf <- 1 - exp(-lambda*dt)
			pRec 	<- 1 - exp(-dt/Tg)
			
			
			if (stochastic) {
				if (lambda > 0 || seedInf > 0 ) {
					nInf <- rbinom(1,sS,pInf) + seedInf
				} else {
					nInf <- 0
				}
				
				nRec <- rbinom(1,sI,pRec)
			} else {
				if (lambda > 0 ) {
					nInf <- sS*pInf + seedInf
				} else {
					nInf <- 0
				}
				
				nRec <- sI*pRec
				
			}
			
			sS <- sS - nInf
			sI <- sI +nInf - nRec
			sR <- sR - nRec
			
			
			rtn[ind_t-1,i] <- rtn[ind_t-1,i] + nInf 
			
			t_cur <- t_cur + dt
			
			if (t_cur > (t_next - epsilon)) {
				t_next <- tps[ind_t]
				ind_t <- ind_t+1
			}
			
		}
	}
	
	if (plot) {
		plot(rtn[,1],type="l",col="grey",ylim=c(0,max(rtn)))
		if (reals > 1) {
			for (j in 2:reals) points(rtn[,j],type="l",col="grey")
		}
	}
	
	list(allruns=rtn)
	
}

irp.proc.year <- function(strYear) {
	lnYear <- nchar(strYear)
	if (lnYear == 4) rtn <- as.numeric(strYear)
	else {
		if (lnYear != 2) stop(paste("irp.proc.year: Year string not the right length:",strYear))
		if (substr(strYear,1,1)=="0") rtn <- as.numeric(strYear)+2000
		else rtn <- as.numeric(strYear)+1900
	}
	rtn
}

irp.process.name <- function(
		strStrain="Type/Host/Place/Year",
		vecKnownHosts=c("Chicken","Swine"),
		vecKnownPlaces=c("Hebei","England")) {
	
	lsParts <- strsplit(strStrain,split="/")
	vecParts <- unlist(lsParts)
	numParts <- length(vecParts)
	
	# Not quite finished
	rtnyear <- irp.proc.year(vecParts[numParts])
	rtntype <- "A"
	rtnhost <- "Human"
	rtnplace <- "China"
	
	list(type=rtntype,host=rtnhost,place=rtnplace,year=rtnyear)
	
}

irp.parse.strainlist <- function(
		strStrList="~/Dropbox/tmp/strain_names.csv",
		strPlaceList,
		strHostList) {
	
	vecStrains <- (read.csv("~/Dropbox/tmp/strain_names.csv",colClasses=c("character"),header=TRUE))[[1]]
	numStrains <- length(vecStrains)
	rtn <- data.frame(year=rep(-1,numStrains))
	for (i in 1:numStrains) {
		curStrain <- vecStrains[i]
		strainData <- irp.process.name(curStrain)
		rtn[i,"year"] <- strainData$year
	}
	
	rtn
}

fnModel <- function(			
		beta=0.1,
		dur=2,
		seed=1,
		N=10,
		dt=1,
		timesteps=30) {
	
	vector_state_infs <- vector(mode="numeric",length=timesteps+1)
	vector_event_infs <- vector(mode="numeric",length=timesteps+1)
	vector_event_recs <- vector(mode="numeric",length=timesteps+1)
	
	vector_state_infs[1] <- seed
	vector_event_infs[1] <- 0
	
	
	for (i in 1:timesteps) {
		
		# Current number of susceptible individuals
		no_cur_sus <- N - vector_state_infs[i]
		
		# Force of infection
		foi <- beta * vector_state_infs[i]
		
		# Probability of events
		pinf <- 1 - exp(-dt*foi)
		prec <- 1 - exp(-dt/dur)
		
		# Generating events and updating the event vectors
		vector_event_infs[i+1] <- 
				rbinom(1,no_cur_sus,pinf)
		vector_event_recs[i+1] <- 
				rbinom(1,vector_state_infs[i],prec)
		
		# Updating the state variables
		vector_state_infs[i+1] <- 
				vector_state_infs[i] + 
				vector_event_infs[i+1] - 
				vector_event_recs[i+1] 
		
	} 
	
	# Return the three vectors from the function with names
	list(	state_inf=vector_state_infs, 
			event_inf=vector_event_infs, 
			event_rec=vector_event_recs)
	
}

fnSets <- function(N=10,timesteps=30,reals=10) {
	
}

# Teaching example of fixed timestep seir model
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# mod1 <- tex.seir.basic(deterministic=TRUE) 
# mod1 <- tex.seir.basic(deterministic=TRUE,N=10000,R0=5,Tg=4,noTimeSteps=35)
# pdf(file="~/Dropbox/tmp/rout1.pdf")
# plot(mod1$inf_inc,type="l")
# dev.off()
tex.seir.basic <- function(	
	  De=1.48,			# Duration latent
		Tg=2.6,				# Generation time
		R0=1.8,				# Basic reproductive number
		N=6800000,		# Population size
		I0=10,				# Initial number infective
		dt=1,			    # Timestep
		R1=1.0,       # Second R value
		t1=999,       # time point of change
		R2=1.0,       # Third R value
		t2=9999,      # time point of change
		noTimeSteps=10,
		deterministic=TRUE) {
	
	# Define the variable to be returned					
	rtn_inf_inc <- vector(mode="numeric",length=noTimeSteps+1)
	
	# Define derived parameters
	Di <- Tg - De
	
	# Set up state variables
	state_S <- N - I0
	state_E <- 0
	state_I <- I0
	state_R <- 0
	
	# Assign the seed infection to time 0
	rtn_inf_inc[1] <- I0
	
	# Start the main loop
	for (i in 1:noTimeSteps) {
		
	  # Set Beta
	  if (i < t1*dt) {
	    beta = R0 / Di
	  } else if (i < t2*dt) {
	    beta = R1 / Di
	  } else {
	    beta = R2 / Di
	  }
	  
		# Calculate hazard rates
		hazInf 		<- beta * state_I / N
		hazBecInf 	<- 1 / De
		hazRec 		<- 1 / Di 
		
		# Calculate probabilities
		pInf 	<- 1 - exp(-dt * hazInf)
		pBecInf <- 1 - exp(-dt * hazBecInf)
		pRec  <- 1 - exp(-dt * hazRec)
		
		# Calculate the expected or simulated number of events
		if (deterministic) {
			nInf 	<- state_S * pInf
			nBecInf <- state_E * pBecInf
			nRec 	<- state_I * pRec
		} else {
			nInf <- rbinom(1,state_S,pInf)
			nBecInf <- rbinom(1,state_E,pBecInf)
			nRec <- rbinom(1,state_I,pRec)
		}
		
		# Update the state variables
		state_S <- state_S - nInf
		state_E <- state_E + nInf - nBecInf
		state_I <- state_I + nBecInf - nRec
		state_R <- state_R - nRec
		
		# Update the output variables
		rtn_inf_inc[i+1] <- nInf
		
	} # End main loop
	
	# Return incidence
	list(inf_inc = rtn_inf_inc, time=(0:noTimeSteps)*dt)
	
}

# A simple from-scratch simulator using the wikipedia example
# http://en.wikipedia.org/wiki/Logistic_regression
srg.sim.logistic <- function(
		
		b_0=-5,
		b_1 = 2, 
		b_2 = -1, 
		b_3 = 1.2, 
		vecX_1 = c(0,4,6,8),
		vecX_2 = c(0,0,1,0),
		vecX_3 = c(2,2,0,1)) {
	
	vecZ <- b_0 + b_1 * vecX_1 + b_2 * vecX_2 + b_3 * vecX_3
	logit_x <- 1 / (1+exp(-vecZ))
	
}

# Simple household model to generate outbreaks with
# only a parameter p
adm.houseshold <- function(p=0.5,N=4) {
	rtn <- c()
	inf <- 1
	sus <- N-1
	while (inf > 0 && sus > 0) {
		newinf <- 0
		for (i in 1:inf) {
			for (j in 1:sus) {
				if (runif(1) < p) {
					newinf <- newinf + 1
					sus <- sus - 1
				}
			}
		}
		rtn <- c(rtn,newinf)
		inf <- newinf
	}
	rtn
}

# setwd("/Users/sriley/Dropbox/talks/201203")
adm.sheep.likelyscript <- function() {
	source("farm_functions_v2.R")
	farmdf <- adm.load.sheep.prep.sim()
	plot(farmdf$easting,farmdf$northing)
	
	nofarms <- dim(farmdf)[1]
	s <- adm.seed.sheep(farmdf)
	p <- adm.sheep.params()
	errcode <- adm.apply.seed(s,nofarms)
	gencode <- 0
	while (gencode < 3) gencode <- adm.apply.gen.model(p,nofarms)
	adm.plot.gen(farmdf)
	table(farmdf$g_i)
}

# rm(list=ls(all=TRUE))
# source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# psi.sim.modelA()
psi.sim.modelA <- function(
		R0 = 1.6,
		Tg=2.6,
		t_0=0,
		pC=1,
		seed=10,
		N=10000,
		tps=0:(7*10),
		reals=1,
		dt=0.25,
		stochastic=FALSE,
		plot=TRUE) {
	
	noTPts <- length(tps) 
	epsilon <- 1e-10
	if (seed >= N) stop("psi.sim.modelA seed greater than N")
	
	if (noTPts < 3) 
		stop("psi.sim.modelA: tps must be of length 2 or greater")
	if (reals > 1 && stochastic==FALSE) 
		stop("psi.sim.modelA: why make multiple realisations of a deterministic model?")
	
	rtn <- array(0,c(noTPts-1,reals))
	for (i in 1:reals) {	
		t_cur <- tps[1]
		ind_t <- 2
		t_next <- tps[ind_t]
		
		sS <- N
		sI <- 0
		sR <- 0
		
		blNotYetSeed <- TRUE
		
		while (ind_t <= noTPts) {
			
			if (t_cur >= t_0 && blNotYetSeed) {
				seedInf <- seed
				blNotYetSeed <- FALSE
			} else {
				seedInf <- 0
			}
			
			lambda <- R0 * sI / Tg / N
			pInf <- 1 - exp(-lambda*dt)
			pRec 	<- 1 - exp(-dt/Tg)
			
			if (stochastic) {
				if (lambda > 0 || seedInf > 0 ) {
					nInf <- rbinom(1,sS,pInf) + seedInf
					nObs <- rbinom(1,nInf,pC)
				} else {
					nInf <- 0
					nObs <- 0
				}
				
				nRec <- rbinom(1,sI,pRec)
			} else {
				if (lambda > 0 || seedInf > 0) {
					nInf <- sS*pInf + seedInf
					nObs <- nInf*pC
				} else {
					nInf <- 0
					nObs <- 0
				}
				
				nRec <- sI*pRec
				
			}
			
			sS <- sS - nInf
			sI <- sI +nInf - nRec
			sR <- sR - nRec
			
			
			rtn[ind_t-1,i] <- rtn[ind_t-1,i] + nObs 
			
			t_cur <- t_cur + dt
			
			if (t_cur > (t_next - epsilon)) {
				t_next <- tps[ind_t]
				ind_t <- ind_t+1
			}
			
		}
	}
	
	if (plot) {
		avge_curve = c(0,dim=length(rtn[,1]))
		for (i in 1:length(rtn[,1])) avge_curve[i] = mean(rtn[i,])
		plot(rtn[,1],type="l",col="grey",ylim=c(0,max(rtn)),xlim=c(0,20),ylab="Number of Cases", 
				xlab="Time (Weeks)")
		if (reals > 1) {
			for (j in 2:reals) points(rtn[,j],type="l",col="grey")
		}
		points(avge_curve,type="l",col="blue")
	}
	
	rtn
	
}

psi.ode.modelA <- function(t,y,p) {
	
	rtn <- vector(mode="numeric",length=length(y))
	names(rtn) <- c("S","I","R","intDS")
	
	beta <- p["R0"]/p["Tg"]
	
	rtn["S"] = -beta*y["S"]*y["I"]/p["N"]
	rtn["I"] =  beta*y["S"]*y["I"]/p["N"] - y["I"]/p["Tg"]
	rtn["R"] =  y["I"]/p["Tg"]
	rtn["intDS"] = beta*y["S"]*y["I"]/p["N"]
	
	list(rtn)
	
}

psi.ode.modelB <- function(t,y,p) {
	
	rtn <- vector(mode="numeric",length=length(y))
	names(rtn) <- c("S","Ia","Ic","Ip","R","intDS")
	
	beta <- p["R0"]/p["Tg"]
	
	rtn["S"]  = -beta*y["S"]*(y["Ia"]+y["Ic"]+y["Ip"])/p["N"]
	rtn["Ia"] =  (1-p["pc"])*beta*y["S"]*(y["Ia"]+y["Ic"]+y["Ip"])/p["N"] - y["Ia"]/p["Tg"]
	rtn["Ic"] =  p["pc"]*(1-p["pp"])*beta*y["S"]*(y["Ia"]+y["Ic"]+y["Ip"])/p["N"] - y["Ic"]/p["Tg"]
	rtn["Ip"] =  p["pc"]*p["pp"]*beta*y["S"]*(y["Ia"]+y["Ic"]+y["Ip"])/p["N"] - y["Ip"]/p["Tg"]
	rtn["R"]  =  y["Ia"]/p["Tg"]+ y["Ic"]/p["Tg"] + y["Ip"]/p["Tg"]
	rtn["intDS"] = beta*y["S"]*(y["Ic"]+y["Ip"])/p["N"]
	
	list(rtn)
	
}

psi.mbm.lnlike <- function(mod,dat) {
	tmp 	<- dpois(dat,mod,log=TRUE)
	sum(tmp)	
} 

psi.mbm.lnlike.chisq <- function(mod,dat,nfittedp=1) {
	fit 	<- dpois(dat,mod,log=TRUE)
	sat 	<- dpois(dat,dat,log=TRUE)
	rtn_chisq 	<- 2 * sum((sat - fit))
	rtn_pval	<-  pchisq(rtn_chisq,length(dat)-1-nfittedp,log.p=TRUE)
	list(chisq=rtn_chisq,pval=rtn_pval)	
} 

psi.fn.opt <- function(par,basepar,names,vecdata,model="ModelA") {
	
	tmp <- psi.opt.inc(par,basepar,names,vecdata,model=model)
	incidence <- tmp$inc
	ps <- tmp$ps
	rtn <- psi.mbm.lnlike(incidence,vecdata)
	rtn
	
}

psi.change.params <- function(par,basepar,names) {
	
	rtn <- basepar
	nochanges <- length(names)
	lengthbase <- length(basepar)
	if (length(par) != nochanges) stop("fnOpt") 
	for (i in 1:nochanges) rtn[names[i]] <- par[i]
	if (length(rtn) != lengthbase) stop("fnOpt")
	
	rtn
}

psi.opt.inc <- function(par,basepar,names,vecdata,model="ModelA") {
	
	tmp_ps <- psi.change.params(par,basepar,names)
	
	noweeks <- length(vecdata) # max(fnData[,"X"])
	endday <- noweeks*7
	blankweeks <- floor(tmp_ps["t0"]/7)
	fullweeks <- noweeks - 1 - blankweeks
	endfirstweek <- (blankweeks+1)*7 - tmp_ps["t0"] 
	timepts <- c(0,endfirstweek,endfirstweek+(1:(fullweeks-1))*7)
	notps <- length(timepts)
	
	if (model=="ModelA") {
		ics <- c(tmp_ps["N"]-1,1,0,0)
		names(ics) <- c("S","I","R","intDS")
		solution <- lsoda(ics,timepts,psi.ode.modelA,tmp_ps,atol=1e-80)
	} else if (model=="ModelB") {
		ics <- c(tmp_ps["N"]-1,1,0,0,0,0)
		names(ics) <- c("S","Ia","Ic","Ip","R","intDS")
		solution <- lsoda(ics,timepts,psi.ode.modelB,tmp_ps,atol=1e-80)
	} else stop("fnOpt_make_inc incorrect model selected")
	
	inc <- vector(mode="numeric",length=noweeks)
	inc[(blankweeks+2):noweeks] <- solution[2:notps,"intDS"] - solution[1:(notps-1),"intDS"]
	inc[] <- inc[] + tmp_ps["e_nonflu"]
	list(inc=inc,ps=tmp_ps)
	
}

psi.trim.data.in <- function(longvec,peakrel=1e-10) {
	peakheight <- max(longvec)
	thresh <- peakheight*peakrel
	noobs <- length(longvec)
	first <- 1
	last <- noobs
	while (longvec[first] < thresh && first <= noobs) first <- first + 1
	while (longvec[last] < thresh && last >= 1) last <- last - 1
	longvec[first:last]
}

psi.trim.data.out <- function(longvec,peakrel=1e-10) {
	peakheight <- max(longvec)
	peakloc <- match(peakheight,longvec)
	thresh <- peakheight*peakrel
	noobs <- length(longvec)
	first <- peakloc
	while (longvec[first] > thresh && first >= 1) first <- first - 1
	if (longvec[first] < thresh) first <- first + 1
	last <- peakloc
	while (longvec[last] > thresh && last <= noobs) last <- last + 1
	if (longvec[last] < thresh) last <- last - 1
	if (first >= last) stop("fnTrimDataOutwards")
	longvec[first:last]
}


# setwd("/Users/sriley/Dropbox/shares/pete_psi_epi/ste/")
# rm(list=ls(all=TRUE))
# source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# psi.ill.chat.feb.10.2012()
psi.ill.chat.feb.10.2012 <- function(sim=TRUE,colnumber=3,seed=1234) {
	
	set.seed(seed)
	ili_small_pandemic <- read.csv("../data/ILI-small-pandemic-curves-by-zip5-v1.2.csv",row.names=1,colClasses=c("integer","character",rep("integer",51)))
	
	base_ps <- c(R0 = 1.4, Tg = 2.6, t0 = 0, N = 3000, e_nonflu=0.001, pc=0.25, pp=0)
	base_min <- c(R0 = 0.9,Tg = 0.5,t0 = 0,N = 10,e_nonflu=0.001, pc=0, pp=0)
	
	if (sim) {
		x <- (psi.sim.modelA())$allruns
		trimdata <- c(0,0)
		try <- 1
		maxtry <- 1000
		while (sum(trimdata) < 50 && try < maxtry) {
			trimdata <- psi.trim.data.in(x[,ceiling(runif(1)*dim(x)[2])],peakrel=1e-10)
			try <- try + 1
		}
		if (try == maxtry) stop("problem lakjsdlajdal")
	} else {
		trimdata <- fnTrimDataInwards(ili_small_pandemic[,colnumber],peakrel=1e-10)
	}
	
	base_max <- c(R0 = 20,Tg = 50,t0 = (length(trimdata)-2.001)*7,N = 100000,e_nonflu=100,pc=1,pp=1)
	vecnames <- c("R0","N","t0")
	vecstarts <- c(1.4,200,10)
	vecps <- base_ps[vecnames]
	vecmin <- base_min[vecnames]
	vecmax <- base_max[vecnames]
	
	sol <-  psi.opt.inc(vecps,base_ps,vecnames,trimdata,model="ModelA")
	
	psi.mbm.lnlike(sol$inc,trimdata)
	psi.mbm.lnlike.chisq(sol$inc,trimdata,length(vecps))
	
	plot(sol$inc,type="l",col="green",main=paste("colno",colnumber),ylim=c(0,max(sol$inc,trimdata)))
	points(trimdata,type="l",col="blue")
	
	psi.fn.opt(c(1.4,2.6,21),base_ps,vecnames,trimdata,model="ModelA")
	
	res <- optim(
			vecstarts,
			psi.fn.opt,
			method="L-BFGS-B",
			lower=vecmin,
			upper=vecmax,
			basepar=base_ps,
			names=vecnames,
			vecdata=trimdata,
			model="ModelA",
			control=list(fnscale=-1))
	
	cat(paste(vecnames,res$par,"\n"))
	
	sol2 <- psi.opt.inc(res$par,base_ps,vecnames,trimdata)
	solstart <- psi.opt.inc(vecstarts,base_ps,vecnames,trimdata)
	plot(sol2$inc,type="l",col="green",main=paste("colnumber",colnumber,"lnlike =",res$value),ylim=c(0,max(sol2$inc,trimdata)))
	points(trimdata,type="l",col="blue")
	points(solstart$inc,type="l",col="red")
	
	psi.mbm.lnlike(sol2$inc,trimdata)
	psi.mbm.lnlike.chisq(sol2$inc,trimdata,length(vecps))
	
}

# setwd("/Users/sriley/Dropbox/shares/pete_psi_epi/ste/")
# rm(list=ls(all=TRUE))
# source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# psi.make.R0.sense.plots.A(scenB=c(1.7,10000),scenC=c(1.7,1000),weeksfit=3,fnOut="3wks")
# psi.make.R0.sense.plots.A(scenB=c(1.5,4000),scenC=c(1.7,2000),weeksfit=5,fnOut="5wks")
psi.make.R0.sense.plots.A <- function(
		sim=TRUE,
		colnumber=3,
		seed=1234,
		weeksfit=3,
		scenB=c(1.7,50000),
		scenC=c(1.7,1000),
		fnOut="psi.make.R0.sense.plots.A") {
	
	require(odesolve)
	
	set.seed(seed)
	ili_small_pandemic <- read.csv("../data/ILI-small-pandemic-curves-by-zip5-v1.2.csv",row.names=1,colClasses=c("integer","character",rep("integer",51)))
	
	base_ps <- c(R0 = 1.4, Tg = 2.6, t0 = 0, N =30000, e_nonflu=0.001, pc=0.25, pp=0)
	base_min <- c(R0 = 0.9,Tg = 0.5,t0 = 0,N = 10,e_nonflu=0.001, pc=0, pp=0)
	
	# Need to change the base ps here
	
	if (sim) {
		x <- (psi.sim.modelA())$allruns
		trimdata <- c(0,0)
		try <- 1
		maxtry <- 1000
		while (sum(trimdata) < 50 && try < maxtry) {
			trimdata <- psi.trim.data.in(x[,ceiling(runif(1)*dim(x)[2])],peakrel=1e-10)
			try <- try + 1
		}
		if (try == maxtry) stop("problem lakjsdlajdal")
	} else {
		trimdata <- fnTrimDataInwards(ili_small_pandemic[,colnumber],peakrel=1e-10)
	}
	
	base_max <- c(R0 = 20,Tg = 50,t0 = (length(trimdata)-2.001)*7,N = 100000,e_nonflu=100,pc=1,pp=1)
	vecnames <- c("R0","N")
	vecstarts <- c(1.4,200)
	vecps <- base_ps[vecnames]
	vecmin <- base_min[vecnames]
	vecmax <- base_max[vecnames]
	
	sol <-  psi.opt.inc(vecps,base_ps,vecnames,trimdata,model="ModelA")
	
	psi.mbm.lnlike(sol$inc,trimdata)
	psi.mbm.lnlike.chisq(sol$inc,trimdata,length(vecps))
	
	plot(sol$inc,type="l",col="green",main=paste("colno",colnumber),ylim=c(0,max(sol$inc,trimdata)))
	points(trimdata,type="l",col="blue")
	
	res <- optim(
			vecstarts,
			psi.fn.opt,
			method="L-BFGS-B",
			lower=vecmin,
			upper=vecmax,
			basepar=base_ps,
			names=vecnames,
			vecdata=trimdata[1:weeksfit],
			model="ModelA",
			control=list(fnscale=-1))
	
	cat(paste(vecnames,res$par,"\n"))
	psi.fn.opt(res$par,base_ps,vecnames,trimdata[1:weeksfit],model="ModelA")
	
	psi.fn.opt(scenB,base_ps,vecnames,trimdata[1:weeksfit],model="ModelA")
	
	psi.fn.opt(scenC,base_ps,vecnames,trimdata[1:weeksfit],model="ModelA")
	
	solA <- psi.opt.inc(res$par,base_ps,vecnames,trimdata)
	solB <- psi.opt.inc(scenB,base_ps,vecnames,trimdata)
	solC <- psi.opt.inc(scenC,base_ps,vecnames,trimdata)
	
	# solstart <- psi.opt.inc(vecstarts,base_ps,vecnames,trimdata)
	pdf(file=paste("~/Dropbox/tmp/",fnOut,".pdf",sep=""))	
	
	plot(solA$inc,type="l",col="green",main="",ylim=c(0,1000))
	points(trimdata)
	points(trimdata[1:weeksfit],pch=19,col="blue")
	points(solB$inc,type="l",col="red")
	points(solC$inc,type="l",col="cyan")
	
	dev.off()
	
}

psi.sevcase.gen.ex <- function() {
	set.seed(7348829)
	simdat1 <- psi.trim.data.in(psi.sim.modelA(R0=1.6,pC=0.6,seed=10,stochastic=TRUE))
	simdat2 <- psi.trim.data.in(psi.sim.modelA(R0=1.85,pC=0.12,seed=10,stochastic=TRUE))
	plot(simdat1,ylim=c(0,max(c(simdat1,simdat2))),type="l",col="red")
	points(simdat2,ylim=c(0,max(c(simdat1,simdat2))),type="l",col="blue")
	list(severe=simdat1,mild=simdat2)
}

# var = "R0"
# plot(tmp1[,var],col="red",pch=19,cex=0.3)
# points(tmp2[,var],col="blue",pch=19,cex=0.3)
# hist(tmp2[,var])

# Parameter ordering is: R0, pc, I0
psi.A2.sim.fit.params <- function(
		simdat,
		nsamps=5) {
	
	# Define bounds for parameters
	ptab <- data.frame(
			max 	= c(	5.0,	1.0,	10000	),
			min 	= c(	0.5,	0.01,	1		),
			log		= c(	FALSE,	TRUE,	TRUE	),
			step 	= c(	0.1,	0.1,	0.1		),
			val 	= c(	1.4,	0.01,	10		)
	)
	
	# Simulate one single
	maxtime <- length(simdat)
	
	nps <- dim(ptab)[1]
	
	# Set initial conditions
	curpars 	<- ptab$val
	initsol 	<- psi.sim.modelA(R0=ptab$val[1],pC=ptab$val[2],seed=ptab$val[3],stochastic=FALSE,tps=0:maxtime,plot=FALSE)
	lnlike 		<- psi.mbm.lnlike(initsol,simdat)
	
	# Set up measuring table
	nMcmc 	<- nps * nsamps 
	tabMcmc <- data.frame(
			index=1:nMcmc,
			lnlike=rep(NA,nMcmc),
			R0=rep(NA,nMcmc),
			pc=rep(NA,nMcmc),
			I0=rep(NA,nMcmc))
	
	# Start the sampling loop
	for (i in 1:nsamps) {
		for (j in 1:nps) {
			
			sampno <- (i-1)*nps + j
			
			# Propose solution
			newparams <- fnProposeParamUpdatesSingle(ptab,j)
			
			# Calculate new loglike
			newsol <- psi.sim.modelA(R0=newparams[1],pC=newparams[2],seed=newparams[3],stochastic=FALSE,tps=0:maxtime,plot=FALSE)
			
			# Calculate loglike
			newlnlike <- psi.mbm.lnlike(newsol,simdat)
			
			diff_like <- newlnlike - lnlike
			
			if (diff_like > 0) {
				accept <- TRUE
			} else if (exp(diff_like) > runif(1)) {
				accept <- TRUE 
			} else accept <- FALSE
			
			if (accept) {
				ptab[,"val"] <- newparams
				lnlike <- newlnlike
			}	
			
			# Update the record of the chain
			tabMcmc[sampno,"lnlike"] <- lnlike
			tabMcmc[sampno,"R0"] <- ptab[1,"val"]
			tabMcmc[sampno,"pc"] <- ptab[2,"val"]
			tabMcmc[sampno,"I0"] <- ptab[3,"val"]
			
		}
	}
	
	tabMcmc
	
}

psi.gen.case.sols <- function(tmc,indsols,simdat,stoch=FALSE) {
	ndays <- length(simdat)
	nsols <- length(indsols)
	rtn <- matrix(nrow=ndays,ncol=nsols)
	for (i in 1:nsols) {
		tmp <-psi.sim.modelA(
				R0=tmc[indsols[i],"R0"],
				pC=tmc[indsols[i],"pc"],
				seed=tmc[indsols[i],"I0"],
				stochastic=stoch,
				tps=0:ndays,
				plot=FALSE)
		rtn[,i] <- tmp
	}
	rtn
}

# rm(list=ls(all=TRUE))
# source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# z1 <- psi.gen.case.sols(tmp1,1800-(0:19)*10,simdat1)
# plot(1:20,type="n",ylim=c(0,1000),xlim=c(0,70))
# nreals <- dim(z1)[2]
# for (i in 1:nreals) points(z1[,i],type="l",col="grey")
# points(simdat1,pch=19,col="red")
# w <- psi.run.sev.cases(tps=c(10,20,30))
# save(w,file="/Users/sriley/Dropbox/shares/pete_psi_epi/ste/steres_w.RData")
# w1 <- psi.run.sev.cases(tps=3:43)
# save(w1,file="/Users/sriley/Dropbox/shares/pete_psi_epi/ste/steres_w1.RData")

psi.run.sev.cases <- function(
		nmcsamps = 600,
		tps=c(10,20,30),
		nsols=20) {
	
	# Generate simulated epidemics
	# See function defaults for parameter values
	simcases <- psi.sevcase.gen.ex()
	simdat1 <- simcases$severe
	simdat2 <- simcases$mild
	
	# Set up return variables
	ntps <- length(tps)
	halfsample <- floor(nmcsamps*3/2) 
	ndsim1 <- length(simdat1)
	ndsim2 <- length(simdat2)
	arrsev <- array(dim=c(ndsim1,nsols,ntps))
	arrmild <- array(dim=c(ndsim2,nsols,ntps))
	arrmcsev <- array(dim=c(halfsample,2,ntps))
	arrmcmild <- array(dim=c(halfsample,2,ntps))
	mcsubset <- (nmcsamps*3-halfsample+1):(nmcsamps*3)
	sampsubset <- nmcsamps*3-(0:(nsols-1))*60
	
	
	for (i in 1:ntps) {
		tmpsev <- 	psi.A2.sim.fit.params(simdat1[1:tps[i]],nsamps=nmcsamps)
		tmpmild <- 	psi.A2.sim.fit.params(simdat2[1:tps[i]],nsamps=nmcsamps)
		zsev <- psi.gen.case.sols(tmpsev,sampsubset,simdat1)
		zmild <- psi.gen.case.sols(tmpmild,sampsubset,simdat2)
		arrsev[,,i] <- zsev[]
		arrmild[,,i] <- zmild[]
		
		arrmcsev[,1,i] <- tmpsev[mcsubset,"R0"]  
		arrmcsev[,2,i] <- tmpsev[mcsubset,"pc"]  
		arrmcmild[,1,i] <- tmpmild[mcsubset,"R0"]  
		arrmcmild[,2,i] <- tmpmild[mcsubset,"pc"]  
		
		# debugplot <- function(dbdat=simdat1,dbres=zsev) {
		#	plot(dbdat,ylim=c(0,200),xlim=c(0,70))
		#	for (j in 1:nsols) points(dbres[,j],type="l",col="grey")
		# }
		
		# debugplot(simdat1,arrsev[,,i])
		# debugplot(simdat2,arrmild[,,i])
		# plot(tmpsev$R0)
		# plot(tmpmild$R0)
		
	}
	
	list(	timepoints=tps,
			datsev=simdat1,datmild=simdat2,
			arrsev=arrsev,arrmild=arrmild,
			paramsev=arrmcsev,parammild=arrmcmild)
	
}

# Function to make sample size calculations for generic paired flu protocol
gfp.simp.samp.size <- function(
		vecNs = c(50,100,500,1000,2000,5000),
		vecPs = c(0.001,0.005,0.01,0.05,0.1,0.5)
) {
	
	colnames <- vecNs
	rownames <- vecPs
	nocols <- length(colnames)
	norows <- length(rownames)
	tablemat <- matrix(ncol=nocols,nrow=norows,dimnames=list(rownames,colnames))
	
	for (i in 1:norows) {
		for (j in 1:nocols) {
			N <- vecNs[j]
			p <- vecPs[i]
			no <- round(p*N)
			lb <- binCINew(N,no,0.975)
			ub <- binCINew(N,no,0.025)
			tablemat[i,j] <- max(ub-p,p-lb)
		}
	}
	
	tablemat
	
}


# To generate a list of origin destination pairs from a list of origins and 
# an image of population denstiy

fsc.gen.connect.dists <- function(
		loclat,										# Vector of latitudes of origin points							
		loclong,									# Vector of longitudes of origin points
		locid,										# Vector of unique ID integers for origins 
		popimage,									# R image of population density length(popimage$x) == dim(popimage$z)[1]
		r_max=30,									# The maximum distance to consider
		dbuniform=FALSE,							# Debug option to ignore popimage and assume uniform density == 1
		centreorigins=TRUE,							# Option to shift origin locations to their nearest cell centre
		csvfile="~/Dropbox/tmp/popdistandsize.csv"	# Name of the output file to write (may be large)
) {
	
	# Validate the inputs
	nolocs <- length(loclat)
	if (!identical(nolocs,length(loclong))) stop("fsc.calc.connect: lats and longs must be of the same length")
	
	# Define some variables that will be needed for a few things
	nox <- length(popimage$x)
	noy <- length(popimage$y)
	if (dbuniform) popimage$z[] <- 1
	
	# Check that none of the boundary squares are too close to the sample squares
	for (i in 1:nox) {
		# Not yet implemented
	}
	
	# Check that origins are within popimage and shift to cell centre if required
	if (centreorigins) {
		dx <- popimage$x[2]-popimage$x[1]
		dy <- popimage$y[2]-popimage$y[1]
		minx <- popimage$x[1] - dx/2
		miny <- popimage$y[1] - dy/2
		for (i in 1:nolocs) {
			# XXX up to here
			indx <- ceiling((loclong[i] - minx) / dx)
			indy <- ceiling((loclat[i] - miny) / dy)
			loclong[i] <- minx + indx * dx - dx / 2
			loclat[i] <- miny + indy * dy - dy / 2
		}
	}
	
	# Cycle through twice to build a list of the x and y for all cells that have a non-zero interaction
	nononezero <- 0
	allpossdistances <- array(data=-1,dim=c(nox,noy,nolocs))
	cat("Calculating all possible distances\n")
	pb <- txtProgressBar(min = 0, max = nox, style = 3)
	for (i in 1:nox) {
		for (j in 1:noy) {
			if (!is.na(popimage$z[i,j]) && popimage$z[i,j] > 0) {
				for (k in 1:nolocs) {
					# sqcentre <- srg.popimage.dblMidCell(i,j,popimage)
					sqcentre <- c(popimage$x[i],popimage$y[j])
					x1 <- sqcentre[1]
					y1 <- sqcentre[2]
					x2 <- loclong[k]
					y2 <- loclat[k]
					tmpdist <- srg.decLongLatDist(x1,y1,x2,y2,translate=TRUE)
					if (tmpdist < r_max) {
						allpossdistances[i,j,k] <- tmpdist
						nononezero <- nononezero + 1
					}
				}
			} 
		}
		setTxtProgressBar(pb, i)
	}
	cat("\n")
	
	# Cycle back through the same set assigning xcoords, ycoods, locations and distances to the 
	reducedset <- matrix(nrow=nononezero,ncol=3,dimnames=list(1:nononezero,c("location","number","distance"))) 
	currentrow <- 1
	cat("Populating reduced set\n")
	pb <- txtProgressBar(min = 0, max = nox, style = 3)
	for (i in 1:nox) {
		for (j in 1:noy) {
			for (k in 1:nolocs) {
				if (allpossdistances[i,j,k] > 0) {
					reducedset[currentrow,"location"] <- locid[k]
					reducedset[currentrow,"distance"] <- allpossdistances[i,j,k]
					reducedset[currentrow,"number"] <- popimage$z[i,j]					
					currentrow <- currentrow + 1
				}
			}
		}
		setTxtProgressBar(pb, i)
	}
	
	# Return matrix of connectivity
	if (!is.null(csvfile)) {
		write.csv(reducedset,file=csvfile,row.names=FALSE)
	} 
	
	reducedset
	
}

srg.raster.to.image <- function(rpop) {
	require("raster")
	
	xmin 	<- rpop@extent@xmin
	ymin 	<- rpop@extent@ymin
	xres 	<- res(rpop)[1]
	yres 	<- res(rpop)[2]
	nx 		<- rpop@ncols
	ny 		<- rpop@nrows
	
	tmpx <- xres*(0:(nx-1)) + xmin + xres/2
	tmpy <- yres*(0:(ny-1)) + ymin + yres/2
	tmpz <- as.matrix(flip(t(rpop),'x'))
	
	list(x=tmpx,y=tmpy,z=tmpz)
	
}


# Returns landscan asia data as a raster object
# need a working rgdal
srg.load.landscan.as.raster <- function(
		landscanfile="/Users/sriley/Dropbox/dataplain/gis/landscan/asia03/asia03/w001001.adf") {
	require("raster")
	require("rgdal")
	r <- raster(landscanfile)
	r
}

fsc.load.wide.raster <- function(
		fluscapetopdir = "/Users/sriley/Dropbox/svneclipse/fluscape/") {
	require("raster")
	rtn = raster(paste(fluscapetopdir,"data/landscan/fluscape_wide.gri",sep=""))
	rtn
} 

row.russell.2012.fig1 <- function() {
	a <- 10^-5
	plot(1:2,type="n",log="y",ylim=c(10^-25,10^0),xlim=c(0,5),xlab="Time",ylab="Probability")
	xvals <- (0:100)/100*5
	x <- xvals * 6 / 24
	colset <- c("purple","red","orange","green","blue")
	for (s in 1:5) {
		m <- s
		y <- x^m * a^m*(1-a)^(x*s-m)
		points(xvals,y,type="l",col=colset[s],lwd=2)
	}
	
}

abp.for.cam.machine <- function(
		nsamps=2,
		filestem="ab_persist",
		dfPBounds = data.frame(
				params=c("gamma_R","dur_seas","R0","seed","amp_seas","Tg"),
				log=c(1,0,0,1,0,0),
				lb=c(1/20/364,7,1.1,1,0,1.8),
				ub=c(1/1/364,7*26,3,100,0.9,2.6))) {
	tmp <- abp.hyper(nohypersamples=nsamps,dfPBounds)
	fn1 <- paste(filestem,"_parameters.csv",sep="")
	fn2 <- paste(filestem,"_last_few_troughs.csv",sep="")
	write.csv(tmp$hypersquare,file=fn1,row.names=FALSE)
	write.csv(tmp$matmin,file=fn2,row.names=FALSE)
}

syp.seirsModel <- function(t,y,p) {
	rtn <- array(0,c(5))
	names(rtn) <- c("S","E","I","R","dS")
	names(y) <- c("S","E","I","R","dS")
	lambda <- p["beta"]*y["I"]/(y["S"]+y["E"]+y["I"]+y["R"])
	rtn["S"] <- -lambda*y["S"] + p["omega"]*y["R"]+(1-p["p_i"])*p["gamma"]*y["I"]
	rtn["E"] <- lambda*y["S"]-p["theta"]*y["E"]
	rtn["I"] <-  p["theta"]*y["E"]- p["gamma"]*y["I"]
	rtn["R"] <- p["p_i"]*p["gamma"]*y["I"] - p["omega"]*y["R"]
	rtn["dS"] <- lambda*y["S"]
	list(rtn)
}

# This could have been the last thing I was working on
syp.seiei_CG_Model <- function(t,y,p) {
	rtn <- array(0,c(10))
	names(rtn) <- c("S_C","E_1C","I_1C","E_2C","I_2C","E_H","L_C","R_C","S_H","E_1H","I_1H","E_2H","I_2H","E_H","L_H","R_H")
	names(y) <- c("S_C","E_1C","I_1C","E_2C","I_2C","E_H","L_C","R_C","S_H","E_1H","I_1H","E_2H","I_2H","E_H","L_H","R_H") 
	lambda_C <- p["beta"]*(p["m_LL"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HL"]*y["I_H"]/p["N_H"])
	lambda_G <- p["beta"]*(p["m_LH"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HH"]*y["I_H"]/p["N_H"])
	rtn["S_L"] <- -lambda_L*y["S_L"] + p["omega"]*y["R_L"]+(1-p["p_i"])*p["gamma"]*y["I_L"]
	rtn["E_L"] <- lambda_L*y["S_L"]-p["theta"]*y["E_L"]
	rtn["I_L"] <-  p["theta"]*y["E_L"]- p["gamma"]*y["I_L"]
	rtn["R_L"] <- p["p_i"]*p["gamma"]*y["I_L"] - p["omega"]*y["R_L"]
	rtn["S_H"] <- -lambda_H*y["S_H"] + p["omega"]*y["R_H"]+(1-p["p_i"])*p["gamma"]*y["I_H"]
	rtn["E_H"] <- lambda_H*y["S_H"]-p["theta"]*y["E_H"]
	rtn["I_H"] <-  p["theta"]*y["E_H"]- p["gamma"]*y["I_H"]
	rtn["R_H"] <- p["p_i"]*p["gamma"]*y["I_H"] - p["omega"]*y["R_H"]
	rtn["dH"] <- lambda_H*y["S_H"]
	rtn["dL"] <- lambda_L*y["S_L"]
	list(rtn)
}

syp.fnPoissonLike <- function(x,y) {
	# y are model incidences
	# x are data
	# output is log likelihood of x if y mean of independent poisson variables
	noSamples <- length(y)
	lnlike <- 0
	mintol <- 1e-20
	for (i in 1:noSamples) {
		if (y[i] < mintol && x[i] > mintol) {
			lnlike <- lnlike - 1000000
		} else {
			debugtmp <- dpois(x[i],y[i],log=TRUE)
			lnlike <- lnlike + debugtmp
		}
	}
	lnlike
}

syp.RepPs <- function(allp,pfn,pfv) {
	nofitted <- length(pfv)
	for (i in 1:nofitted) {
		allp[pfn[i]]<-pfv[i]
	}
	allp
}

# Define the actual model
syp.seirs_LH_Model <- function(t,y,p) {
	rtn <- array(0,c(10))
	names(rtn) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
	names(y) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH") 
	lambda_L <- p["beta"]*(p["m_LL"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HL"]*y["I_H"]/p["N_H"])
	lambda_H <- p["beta"]*(p["m_LH"]*y["I_L"]/(p["N"]-p["N_H"])+p["m_HH"]*y["I_H"]/p["N_H"])
	rtn["S_L"] <- -lambda_L*y["S_L"] + p["omega"]*y["R_L"]+(1-p["p_i"])*p["gamma"]*y["I_L"]
	rtn["E_L"] <- lambda_L*y["S_L"]-p["theta"]*y["E_L"]
	rtn["I_L"] <-  p["theta"]*y["E_L"]- p["gamma"]*y["I_L"]
	rtn["R_L"] <- p["p_i"]*p["gamma"]*y["I_L"] - p["omega"]*y["R_L"]
	rtn["S_H"] <- -lambda_H*y["S_H"] + p["omega"]*y["R_H"]+(1-p["p_i"])*p["gamma"]*y["I_H"]
	rtn["E_H"] <- lambda_H*y["S_H"]-p["theta"]*y["E_H"]
	rtn["I_H"] <-  p["theta"]*y["E_H"]- p["gamma"]*y["I_H"]
	rtn["R_H"] <- p["p_i"]*p["gamma"]*y["I_H"] - p["omega"]*y["R_H"]
	rtn["dH"] <- lambda_H*y["S_H"]
	rtn["dL"] <- lambda_L*y["S_L"]
	list(rtn)
}

syp.fnModelLike <- function(v,pnamesfitted,allps,alldata_h,alldata_l,report=FALSE) {
	
	# Set up the required parameters and state variables
	allps <- syp.RepPs(allps,pnamesfitted,v)
	ics <- c(allps["N"]-allps["N_H"],0,0,0,allps["N_H"]-allps["seed"],0,allps["seed"],0,0,0)
	names(ics) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
	
	no_years <- length(alldata_h)
	time_points <- 0:no_years
	timestep <- 1
	solution <- lsoda(ics,time_points,syp.seirs_LH_Model,allps)
	yL <- array(-100,no_years)
	yH <- array(-100,no_years)
	for (i in 1:no_years) {
		yL[i] <- solution[i+1,"dL"]-solution[i,"dL"]
		yH[i] <- solution[i+1,"dH"]-solution[i,"dH"]
	}
	like <- syp.fnPoissonLike(alldata_l,yL)
	like <- like + syp.fnPoissonLike(alldata_h,yH)
	
	if (report) {
		cat(like,allps,names(allps),"\n")
		flush.console()
	}
	
	# Return the likelihood
	like
	
}

# At the moment, just a copy paste of the above function
syp.seiei.mod.like <- function(v,pnamesfitted,allps,alldata_h,alldata_l,report=TRUE,plot=TRUE) {
	
	# v <- c(1000,1)
	# pnamesfitted <- c("N","seed")
	# allps <- syp.seirParams
	# alldata_h <- incData_H
	# alldata_l <- incData_L
	
	allps <- syp.RepPs(allps,pnamesfitted,v)
	
	if (	allps["seed"] > allps["N"] 	||
			allps["seed"] < 1 			||
			FALSE) {
		like <- -10000000
	} else {
		ics <- c(allps["N"]-allps["N_H"],0,0,0,allps["N_H"]-allps["seed"],0,allps["seed"],0,0,0)
		names(ics) <- c("S_L","E_L","I_L","R_L","S_H","E_H","I_H","R_H","dL","dH")
		no_years <- 18
		time_points <- 0:no_years
		timestep <- 1
		solution <- lsoda(ics,time_points,syp.seirs_LH_Model,allps)
		yL <- array(-100,no_years)
		yH <- array(-100,no_years)
		for (i in 1:no_years) {
			yL[i] <- solution[i+1,"dL"]-solution[i,"dL"]
			yH[i] <- solution[i+1,"dH"]-solution[i,"dH"]
		}
		like <- syp.fnPoissonLike(alldata_l,yL)
		like <- like + syp.fnPoissonLike(alldata_h,yH)
	}
	
	if (report) {
		cat(like,allps,names(allps),"\n")
	}
	if (plot) { 
		# To be completed
	}
	flush.console()
	like
}

syp.vis.stage.data <- function(
		dat,
		filename,
		cols = c("red","blue","green","magenta","cyan","orange"),
		ymax = 10000,
		xmin = 1994,
		logarg = "y",
		tofile=TRUE) {
	if (tofile) pdf(filename)
	nocols <- length(names(dat))
	vars <- names(dat)[2:nocols]
	novars <- nocols - 1
	width <- 2
	plot(	1:2,
			ylim=c(1,ymax),
			xlim=c(xmin,2011),
			ylab = "Reported incidence (plus 1)",
			xlab="Year",
			log=logarg,
			axes=FALSE)
	axis(1)
	axis(2)
	for (i in 1:novars) points(dat$year,dat[,vars[i]]+1,type="l",col=cols[i],lwd=width)
	legend(xmin,ymax,vars,col=cols,lwd=width,bty="n")
	if (tofile) dev.off()
}

# At the moment, just a copy paste
syp.seiei.Params <- function() {
	rtn <- c(
			N = 7000000,				# Size of the population
			N_C = 50000,				# Number in high risk group
			seed = 1,					# Number initially infectious
			omega = 1/(26/52),			# Duration immune stage
			theta = 1/(3/52),			# Duration of latent period
			gamma = 1/2,				# Duration of infectiousness
			p_i = 0.2,					# Proportion gaining immunity
			beta = 1,					# Infectiousness
			m_LH = 1,					# Relative infectivity of low risk to high risk 
			m_HH = 1,					# Relative infectivity of hight risk to high risk
			m_LL = 1,					# Relatvie infectivity of low risk to low risk
			m_HL = 1					# Relative infectivity of high risk to low risk
	)
	rtn
}

srg.mod.gamma.rng <- function(n,mu,shp=10) {
	# Simulate from the modified poisson
	rtn <- vector(length=n)
	for (i in 1:n) {
		tmp <- 0
		while (tmp < 1) tmp <- rgamma(1,shape=shp,scale=(mu+1)/shp)
		tmp <- tmp - 1
		rtn[i] <- floor(tmp)
		# rtn[i] <- tmp
	}
	rtn
}

srg.mod.gamma.dens <- function(xvals,mu,shp=10) {
	rtn 	<- vector(length=length(xvals))
	for (i in 1:length(xvals)) {
		# rtn[i] 	<- dgamma(xvals[i]+1,shape=shp,scale=(mu+1)/shp,log=TRUE) - log(1-dgamma(0,shape=shp,scale=(mu+1)/shp,log=FALSE)) 
		rtn[i] <- log(pgamma(xvals[i]+2,shape=shp,scale=(mu+1)/shp,log=FALSE)-pgamma(xvals[i]+1,shape=shp,scale=(mu+1)/shp,log=FALSE)) - log(1-pgamma(1,shape=shp,scale=(mu+1)/shp,log=FALSE))
		if (rtn[i] < -1e100) stop("alksdjalksd problem here")
	}
	exp(rtn)
}

# Generic multinomial distribution for titre results
# hist(srg.titre.rng(10000,0,terraces=5),breaks=-0.5+0:12)
# Next: figure out why this doesn't work
srg.titre.pvec <- function(exp,factor,terraces,tmax) {
	
	# Setup the multinomial vector 
	pvec <- vector(length=tmax+1)
	pvec[] <- 1/factor^terraces
	pvec[exp+1] <- 1
	for (i in 1:terraces) {
		if (exp+1-i > 0) pvec[exp+1-i] <- 1/factor^i
		if (exp+1+i <= tmax) pvec[exp+1+i] <- 1/factor^i
	}
	pvec
	
}


srg.titre.rng <- function(n,exp,factor=10,terraces=5,tmax=10) {
	
	# Setup vector
	pvec <- srg.titre.pvec(exp,factor,terraces,tmax)
	
	# Simulate from the modified poisson
	rtn <- vector(length=n)
	for (i in 1:n) {
		ranvec <- rmultinom(1,1,pvec)
		rtn[i] <- match(1,ranvec)-1
	}
	rtn
}

srg.titre.dens <- function(y,exp,factor=10,terraces=6,tmax=12) {
	
	# Setup vector
	pvec <- srg.titre.pvec(exp,factor,terraces,tmax)
	pvec <- pvec/sum(pvec)
	
	rtn 	<- vector(length=length(y))
	for (i in 1:length(y)) rtn <- pvec[y+1]
	rtn
	
}

# hist(srg.titre.rng(10000,3),breaks=-0.5+0:12)
# plot(0:10,srg.titre.dens(0:10,3),type="l")

agr.times.tables <- function(noqs=10,tts=1:10) {
	startt <- proc.time()[3]
	print("Hello Abbie")
	correct <- 0
	noposs <- length(tts)
	for (i in 1:noqs) {
		a <- ceiling(runif(1)*10)
		b <- ceiling(runif(1)*noposs)
		b <- tts[b]
		c <- a*b
		x <- readline(paste("Question",i,"of",noqs,"What is ",a," times ",b," ? "))
		if (x==c) {
			print("Correct")
			correct <- correct + 1
		} else {
			print("Better luck next time")
		}
	}
	endt <- proc.time()[3]
	print(paste("You got",correct,"out of",noqs,"correct. In",round(endt-startt,1),"seconds."))
}

agr.times.levels <- function(level,noqs=10) {
	if (level=="chimpanzee") agr.times.tables(noqs=noqs,tts=rep(2,10)) 
	else if (level=="parrot") agr.times.tables(noqs=noqs,tts=rep(10,10))
	else if (level=="peakcock") agr.times.tables(noqs=noqs,tts=5)
	else stop("Don't recognize that level")
	
}

# Little function to use to easily make sure that Dropbox has synced your changes
srg.test.sync <- function() {print("today's message banana")}

srg.quick.jpeg.line <- function(vecDat,suffix=1) {
	jpeg(paste("tmp_",suffix,".jpg",sep=""))
	plot(vecDat,type="l")
	dev.off()
}

vtv.michael.plot <- function(datFile="serological_michael.csv") {
	
	# Read in the data
	dat <- read.csv(datFile)
	pheight <- 5
	pwidth <- 8
	pdf("human_strains.pdf",useDingbats=FALSE,height=pheight,width=pwidth)
	# Parameter settings for all parts of the figure
	# par(cex=0.8)
	# par(mgp=c(2.5,1,1))
	# par(tcl=0.25)
	par(mai=(c(0,0,0,0)))
	par(fig=c(0.20,0.99,0.15,0.99))
	datHumanH1 <- dat[dat$type2=="human" & dat$subtype=="H1N1",]
	datHumanH3 <- dat[dat$type2=="human" & dat$subtype=="H3N2",]
	datHumanH2 <- dat[dat$type2=="human" & dat$subtype=="H2N2",]
	plot(
			datHumanH1$year,
			log(datHumanH1$vn.titre/5,2),
			pch=19,col="red",
			xlab="Year",
			ylab="Antibody strength",
			ylim=c(0,8),
			xlim=c(1950,2010),
			axes=FALSE)
	axis(1)
	axis(2,at=0:8,labels=paste("1:",5*2^(0:8)),las=1)
	points(datHumanH1$year,log(datHumanH1$vn.titre/5,2),type="l",col="red")
	points(datHumanH3$year,log(datHumanH3$vn.titre/5,2),pch=19,col="green")
	points(datHumanH3$year,log(datHumanH3$vn.titre/5,2),type="l",col="green")
	points(datHumanH2$year,log(datHumanH2$vn.titre/5,2),pch=19,col="blue")
	abline(h=3)
	grid()
	dev.off()
	
	pdf("animal_strains.pdf",useDingbats=FALSE,height=pheight,width=pwidth)
	# Parameter settings for all parts of the figure
	# par(cex=0.8)
	# par(mgp=c(2.5,1,1))
	# par(tcl=0.25)
	par(mai=(c(0,0,0,0)))
	par(fig=c(0.20,0.99,0.15,0.99))
	datAnimalH2 <- dat[dat$type2=="animal",]
	plot(
			datAnimalH2$year,
			log(datAnimalH2$hi.titer/5,2),
			pch=19,col="red",
			xlab="Year",
			ylab="Antibody strength",
			ylim=c(0,8),
			xlim=c(1975,2012),
			axes=FALSE)
	axis(1)
	axis(2,at=0:8,labels=paste("1:",5*2^(0:8)),las=1)
	# points(datAnimalH2$year,log(datAnimalH2$hi.titer/5,2),type="l",col="red")
	grid()
	abline(h=3)
	dev.off()
	
}

# setwd("~/Dropbox/shares/pete_small_share")
# psi.scatter.jan2013()
psi.scatter.jan2013 <- function() {
	ISR <- 0.02
	par(cex=1.4)
	x <- read.csv("sirMLE_omega1.6.csv") 
	x <- cbind(Rbest=apply(cbind(x$R0,x$R.),1,max),x)
	plot(1:2,type="n",log="x",ylim=c(1,3),xlim=c(0.0001,1),axes=FALSE,
			ylab="Basic Reproductive Number, R0",
			xlab="Severe Infection Rate")
	points(ISR*x[x$Ntotal>10000,"pC"],x[x$Ntotal>10000,"Rbest"])
	axis(1,at=10^(-4:0),labels=c("0.01%","0.1%","1%","10%","100%"))
	axis(2)	
}

# Make a single run from a branching process model
# need to start here and run out the change of starting assumptions
# debug with the multiple runs version
ncv.single.run <- function(
  theta,				# Parameter vector
  detThresh = 1E10, 	# Deterministiuc threshold for expected number of infections 
  maxcuminc = 1E10	# Threshold beyond which incidence is no longer calculated
) {
  
  # Call out the required libraries
  require("date")
  
  # Declare local functions
  sumnb <- function(m,n,k,mu) {
    rtn <- vector(mode="numeric",length=m)
    for (i in 1:m) rtn[i] <- sum(rnbinom(n,size=k,mu=mu))
    rtn
  }
  
  # Extract parameters from names vector
  R0 		<- as.numeric(theta["R0"])
  k 		<- as.numeric(theta["k"])
  Tg 		<- as.numeric(theta["Tg"])
  alpha	<- as.numeric(theta["alpha"])
  pS 		<- as.numeric(theta["pS"])
  pDS 	<- as.numeric(theta["pDS"])
  rIntro 	<- as.numeric(theta["rIntro"])
  if (rIntro < 1) rIntro <- 0
  t0		<- as.integer(theta["t0"])
  tend  	<- as.integer(theta["tend"])
  tstart  <- as.integer(theta["tstart"])
  Tgmax <- as.integer(theta["Tgmax"])
  tmax 	<- tend - tstart
  
  # Check time periods
  if (t0 > as.date("13-Jun-2012")) {
    stop("t0 can't be on or after the first reported symptom onset")
  }
  if (t0 <= tstart) stop("to cannot be less than tmin")
  if (tend < as.date("13-Jun-2012")) stop("tmax must be greater than the first case")
  
  # Define the return structure
  rtn <- matrix(nrow=tmax,ncol=3,dimnames=list(1:tmax,c("Inf","Sev","Rep")))
  rtn[] <- 0
  
  # Generate animal introductions
  indexlivedays <- (t0-tstart):tmax 
  rtn[indexlivedays,"Inf"] <- rpois(runif(length(indexlivedays)),rIntro/365)
  rtn[t0-tstart,"Inf"] <- theta["I0"]
  
  # Main loop to generate additional infections
  skipthrough <- FALSE
  for (i in 2:tmax) {
    # if (skipthrough | sum(rtn[,"Inf"]) >= maxcuminc) {
    if (skipthrough | sum(rtn[,"Inf"]) >= maxcuminc) {
      rtn[i,"Inf"] <- 0 
      skipthrough <- TRUE
    } else {
      ntor 	<- rtn[i-1,"Inf"]
      expinf	<- ntor * R0
      repar_p <- 1.0/(1.0+k/R0)
      repar_r <- R0*(1.0-repar_p)/repar_p
      if (expinf < detThresh) {
        
        # put the nbinom code here
        if (ntor > 0) {
          # nonew <- sumnb(1,ntor,k=k,mu=R0)
          nonew <- rnbinom(1,ntor*repar_r,1-repar_p)
        } else {nonew <- 0}
        
        # tmp <- table(1+rpois(nonew,Tg-1))
        tmp <- tabulate(ceiling(rgamma(nonew,alpha,scale=Tg/alpha)),nbins=Tgmax)
        
      } else {
        nonew <- expinf
        tmp 	<- round(dgamma(0:(Tgmax-1),alpha,scale=Tg/alpha)*expinf)
      }
      if (nonew > 0) {
        for (j in 1:Tgmax) {
          if (i+j <= tmax) rtn[i+j,"Inf"] <- rtn[i+j,"Inf"] + tmp[j]
        }
      }
    }
  }
  
  # Now generate severe and observed cases
  rtn[2:tmax,"Sev"] <- rpois(length(rtn[2:tmax,"Inf"]),pS*rtn[2:tmax,"Inf"])
  rtn[2:tmax,"Rep"] <- rpois(length(rtn[2:tmax,"Rep"]),pDS*rtn[2:tmax,"Sev"])
  
  # Return the time series from the model
  rtn
  
}

# ncv.theta.base(c("R0"),c(1))
ncv.theta.base <- function(pnames=NULL,pvals=NULL) {
  rtn <- array(
    c(1.2,0.16,7,4.88,0.05,1,0,1,
      as.date("01-Jan-2012"),as.date("31-Dec-2011"),as.date("01-Jun-2013"),30),
    dim=c(12),
    dimnames=list(
      c("R0","k","Tg","alpha","pS","pDS","rIntro","I0","t0","tstart","tend","Tgmax"))
  )
  if (length(pnames > 0 & length(pnames)==length(pvals))) {
    for (i in 1:length(pnames)) {
      if (pnames[i] %in% names(rtn)) rtn[pnames[i]]<-pvals[i]
    }
  }
  rtn
}

# Next - implement max incidence
# Repeat runs with plotting if required
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# x <- ncv.multiple.reals(10,makeplot=TRUE)
# plot(apply(x[,"Inf",],c(1),mean),type="l")
ncv.multiple.reals <- function(
  nreals = 10,			
  theta = ncv.theta.base(),
  makeplot=FALSE,		# Should there be a summary plot
  logy = "",			
  plotvar = "Inf",
  maxpropfailure = 100,
  propNoExtinctThresh = 10,
  ...
) {			
  
  tmax <- theta["tend"] - theta["tstart"]
  Tg <- theta["Tg"]
  if (tmax < 2*Tg) stop("tmax must be 2*Tg or greater")
  
  rtn <- array(dim=c(tmax,3,nreals),dimnames=list(1:tmax,c("Inf","Sev","Rep"),1:nreals))
  nosuccess <- 1
  maxfailure <- nreals*maxpropfailure
  nofailure <- 0
  
  while (nosuccess <= nreals) {
    tmp <- ncv.single.run(theta,...)
    if (sum(tmp[(tmax-2*Tg):tmax,]) > 0) {
      rtn[,,nosuccess] 	<- tmp
      nosuccess <- nosuccess + 1
    } else {
      nofailure <- nofailure + 1
    }
  }
  if (nofailure > maxfailure) warning("Too many failures")
  
  # Plot all reals
  if (makeplot) {
    ymax = max(rtn[,plotvar,]+1)
    plot(1:2,type="n",ylim=c(1,ymax),xlim=c(0,tmax),log=logy)
    for (i in 1:nreals) points(rtn[,"Inf",i],type="l")
  }
  
  rtn
  
}

# Latin hypercube sampling and filtering
# system.time(x <- ncv.latin())
ncv.latin <- function(
  nops = 1,
  nreals= 10,	
  theta = ncv.theta.base(),	
  pnames = c("R0"),
  lbs = c(1.19),
  ubs = c(1.21),
  logscale = c(FALSE),
  objectname = "ncv_output",
  filename = "~/Dropbox/tmp/ncv_latin_output_",
  ...
) {
  
  # Set up the timescale
  tmax <- theta["tend"]-theta["tstart"]
  
  # Set up the hypercube samples
  # problem here with the hyper vector of dates
  rtnps <- matrix(nrow=nops,ncol=length(pnames),dimnames=list(1:nops,pnames))
  for (i in 1:length(pnames)) {
    rtnps[,pnames[i]] <- srg.hyper.vector(nops,lbs[i],ubs[i],logscale[i]) 
  }
  
  # Setup the main results array
  vals <- array(dim=c(tmax,3,nreals,nops),dimnames=list(1:tmax,c("Inf","Sev","Rep"),1:nreals,1:nops))
  # loopres <- foreach(i=1:nops, .export=ls(envir=globalenv())) %do% {
  loopres <- foreach(i=1:nops, .export=ls(envir=globalenv())) %do% {
    theta_tmp <- theta
    for (j in 1:length(pnames)) theta_tmp[pnames[j]] <- rtnps[i,j]
    ncv.multiple.reals(nreals,theta_tmp,...)
  }
  for (i in 1:nops) vals[,,,i] <- loopres[[i]]
  
  # Define the return object
  rtn <- list(vals=vals,ps=rtnps)
  
  # Save as an R object
  fname <- paste(filename,objectname,".Rout",sep="")
  assign(objectname,rtn)
  save(list=c(objectname),file=fname)
  
  # Return val
  rtn
  
}

# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# setwd("~/Dropbox/tmp")
# load("~/Dropbox/tmp/ncv_latin_output_ncv_output.Rout")
# load("~/Dropbox/tmp/ncv_latin_output_ncv20130605.Rout")
# ncv.hyper.plots(ncv_output)
ncv.hyper.plots.v2 <- function(ho,noreals=10) {
  
  runs_good_v2 <- function(vec,minval=1000,maxval=10000) { 
    tmp <- ifelse((vec > minval) & (vec < maxval),1,0)
    rtn <- sum(tmp)
    rtn
  }
  
  mean_non_zero <- function(vec) {
    nonzmask <- ifelse(vec > 0,TRUE,FALSE)
    nonz <- vec[nonzmask]
    median(nonz)
  }
  
  # Need to work below here to make sure is OK
  ps <- ho$ps
  if ("pS" %in% colnames(ps)) ps[,"pS"] <- log10(ps[,"pS"])
  # if ("k" %in% colnames(ps)) ps[,"k"] <- log10(ps[,"k"])
  if ("rIntro" %in% colnames(ps)) ps[,"rIntro"] <- log10(ps[,"rIntro"])
  val <- ho$vals
  x1 <- val[,"Inf",,]
  x2 <- apply(x1,c(2,3),sum)
  nomeetcrit <- apply(x2,c(2),runs_good_v2)
  mask <- nomeetcrit > 1
  pmask <- ps[mask,]
  zmask <- nomeetcrit[mask] / noreals
  
  tiff("null.jpg")
  plot(as.data.frame(ps),col="red",pch=22)
  dev.off()
  
  ncol <- 100
  axisprop <- 0.9
  colscale <- gray(0:ncol / ncol * axisprop)
  zval <- log10(nomeetcrit[mask])
  zmin <- 0
  zmax <- log10(noreals/2)
  zval <- ifelse(zval > zmax,zmax,zval)
  zval <- ifelse(zval < zmin,zmin,zval)
  znorm <- as.integer((zval - zmin) / (zmax-zmin) * ncol)
  cexmain <- 4
  pchmain <- 20
  
  
  pdf("legend.pdf",width=30/cm(1),height=30/cm(1),useDingbats=FALSE)
  plot(c(1,1,1),c(1,2,3),pch=pchmain,cex=cexmain,col=colscale[c(1,ncol/2,ncol)])
  dev.off()
  
  # Make the correct colour index and then order according to size
  # Points must be added from top to bottom?
  pos1 <- c(0.1,0.4,0.1,0.9)
  
  
  
  # tiff("criteria_1.tif",width=1024,height=1024,compression="lzw")
  pdf("criteria_1.pdf",width=30/cm(1),height=30/cm(1),useDingbats=FALSE)
    pmaskTmp <- pmask
    if ("Tg" %in% colnames(ps)) pmaskTmp[,"Tg"] <- pmask[,"Tg"]+0.5
    par(fig=pos1)
    plot(as.data.frame(pmaskTmp),pch=pchmain,cex=cexmain,col=colscale[znorm])
  dev.off()
  
}


# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# load("~/Dropbox/tmp/ncv_latin_output_ncv_output.Rout")
# load("~/Dropbox/tmp/ncv_latin_output_ncv20130605.Rout")
# ncv.hyper.plots(ncv_output)
ncv.hyper.plots <- function(ho) {
  
  # Define some local functions to manipulate the runs
  case_last_window <- function(vec,window=30) {
    lvec <- length(vec)
    if (lvec < window) stop("XXX")
    sum(vec[(lvec-30+1):lvec])
  }
  
  runs_good <- function(vec,minval=1,maxval=10,prop=0.05) {
    tmp <- ifelse(vec > minval & vec < maxval,1,0)
    rtn <- FALSE
    if ((sum(tmp)/length(tmp)) > (prop - 1e-5)) {rtn <- TRUE}
    rtn
  }
  
  runs_good_v2 <- function(vec,minval=1,maxval=10) { 
    tmp <- ifelse(vec > minval & vec < maxval,1,0)
    rtn <- sum(tmp)
    rtn
  }
  
  mean_non_zero <- function(vec) {
    nonzmask <- ifelse(vec > 0,TRUE,FALSE)
    nonz <- vec[nonzmask]
    median(nonz)
  }
  
  # Need to work below here to make sure is OK
  ps <- ho$ps
  if ("pS" %in% colnames(ps)) ps[,"pS"] <- log10(ps[,"pS"])
  if ("k" %in% colnames(ps)) ps[,"k"] <- log10(ps[,"k"])
  if ("rIntro" %in% colnames(ps)) ps[,"rIntro"] <- log10(ps[,"rIntro"])
  val <- ho$vals
  x1 <-val[,"Sev",,]
  x2 <- apply(x1,c(2,3),case_last_window)
  x3 <- apply(x2,c(2),mean_non_zero)
  
  mask1 <- apply(x2,c(2),runs_good)
  vecNoConReals <- apply(x2,c(2),runs_good_v2,minval=337,maxval=2735)
  
  y1 <- ps[mask1,]
  y2 <- x3[mask1]
  y3 <- vecNoConReals[mask1]
  
  tiff("null.jpg")
  plot(as.data.frame(ps),col="red",pch=22)
  dev.off()
  
  ncol <- 100
  axisprop <- 0.9
  colscale <- gray(0:ncol / ncol * axisprop)
  zval <- log10(y2)
  zmin <- 0
  zmax <- 1
  zval <- ifelse(zval > zmax,zmax,zval)
  zval <- ifelse(zval < zmin,zmin,zval)
  znorm <- as.integer((zval - zmin) / (zmax-zmin) * ncol)
  cexmain <- 2
  pchmain <- 20
  
  pdf("legend.pdf",width=30/cm(1),height=30/cm(1),useDingbats=FALSE)
  plot(c(1,1,1),c(1,2,3),pch=pchmain,cex=cexmain,col=colscale[c(1,ncol/2,ncol)])
  dev.off()
  
  # Make the correct colour index and then order according to size
  # Points must be added from top to bottom?
  
  tiff("criteria_1.tif",width=1024,height=1024,compression="lzw")
  # pdf("criteria_1.pdf",width=30/cm(1),height=30/cm(1),useDingbats=FALSE)
  plot(as.data.frame(y1),pch=pchmain,cex=cexmain,col=colscale[znorm])
  dev.off()
  
  ncol2 <- 10
  colscale2 	<- rev(gray((1:ncol2)/ncol2))
  zval2 		<- ifelse(y3 < ncol2,y3,ncol2)
  y4 			<- y1[order(zval2),]
  zvalsort 	<- zval2[order(zval2)]
  
  tiff("no_passing_hurdle.tif",width=1024,height=1024,compression="lzw")
  # pdf("criteria_1.pdf",width=30/cm(1),height=30/cm(1),useDingbats=FALSE)
  plot(as.data.frame(y4),pch=pchmain,cex=cexmain,col=colscale2[zvalsort])
  dev.off()
  
}


#bdg(10000,10,1,6)
#bdg <- function(m,n,size,mu) {hist(sumnb(m,n,size,mu)-rpois(m,n*mu))}
#bdg2 <- function(m,size,mu) {hist(rnbinom(m,size=size,mu=mu)-sumnb(m,1,size,mu))}
#hist(rnbinom(10000,size=10,mu=3))
#hist(sumnb(10000,100,10,3))
#hist(rpois(1000,300))


# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(job.20130307.a())
job.20130307.a <- function() {
  
  # This call generates an object of size 83.5 Gb
  # and should take about an hour to run on the MBA
  ncv.latin(
    nops = 2500,
    nreals=20,
    theta = array(
      c(1.2,8,0.05,1,3),
      dim=c(5),
      dimnames=list(c("R0","Tg","pS","pDS","rIntro"))),
    pnames = 	c("R0",		"pS",	"rIntro",	"Tg"),
    lbs = 		c(0.5,		0.01,	0,			2),
    ubs = 		c(1.5,		1,		10,			20),
    logscale = 	c(FALSE,	TRUE,	FALSE,		FALSE),
    objectname = "ncv20130307aXXXXX",
    filename = "~/tmp/ncv_latin_output_"
  )
  
}

# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(job.20130312())
job.20130312 <- function() {
  
  # This call generates an object of size 83.5 Gb
  # and should take about an hour to run on the MBA
  ncv.latin(
    nops = 2279,
    nreals=20,
    theta = array(
      c(1.2,10,8,0.05,1,3),
      dim=c(6),
      dimnames=list(c("R0","k","Tg","pS","pDS","rIntro"))),
    pnames = 	c("R0",		"k",	"pS",	"rIntro",	"Tg"),
    lbs = 		c(0.5,		0.1,	0.01,	0,			2),
    ubs = 		c(1.5,		100,	1,		10,			20),
    logscale = 	c(FALSE,	TRUE,	TRUE,	FALSE,		FALSE),
    objectname = "ncv20130312a",
    filename = "~/tmp/ncv_latin_output_"
  )
  
}

# Up to here
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(x <- job.20130313())
job.20130313 <- function() {
  
  # Load up required packages
  require("date")
  
  # This call generates an object of size 83.5 Gb
  # and should take about an hour to run on the MBA
  ncv.latin(
    nops = 7000,
    nreals=20,
    theta = array(
      c(1.2,10,8,0.05,1,3,as.date("13-Jun-2012"),as.date("13-Mar-2013"),1),
      dim=c(9),
      dimnames=list(c("R0","k","Tg","pS","pDS","rIntro","t0","tend","I0"))),
    pnames = 	c("R0",		"k",	"pS",	"rIntro",	"Tg", 	"t0"),
    lbs = 		c(0.5,		0.1,	0.01,	0.1,		2,		as.date("01-Jan-2012")),
    ubs = 		c(1.5,		100,	1,		100,		20,		as.date("12-Jun-2012")),
    logscale = 	c(FALSE,	TRUE,	TRUE,	TRUE,		FALSE,	FALSE),
    objectname = "ncv20130313xxxx",
    filename = "~/tmp/ncv_latin_output_"
  )
  
}

# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(x <- job.20130314())
job.20130314 <- function() {
  
  # Load up required packages
  require("date")
  
  # This call generates an object of size 83.5 Gb
  # and should take about an hour to run on the MBA
  ncv.latin(
    nops = 1000,
    nreals=20,
    theta = array(
      c(1.2,10,8,0.1,1,3,as.date("13-Apr-2012")),
      dim=c(7),
      dimnames=list(c("R0","k","Tg","pS","pDS","rIntro","t0"))),
    pnames = 	c("R0",		"k"),
    lbs = 		c(0.5,		0.1),
    ubs = 		c(1.5,		100),
    logscale = 	c(FALSE,	TRUE),
    objectname = "ncv20130314",
    filename = "~/tmp/ncv_latin_output_"
  )
  
}

# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(x <- job.20130605())
job.20130605 <- function() {
  
  # Load up required packages
  require("date")
  
  # This call generates an object of size 83.5 Gb
  # and should take about an hour to run on the MBA
  ncv.latin(
    nops = 30,
    nreals=20,
    theta = array(
      c(1.2,10,10,0.1,1,3,as.date("13-Apr-2012")),
      dim=c(7),
      dimnames=list(c("R0","k","Tg","pS","pDS","rIntro","t0"))),
    pnames = 	c("R0",		"k", 	"pS",	"rIntro"),
    lbs = 		c(0.5,		0.1,	0.005,	0),
    ubs = 		c(2.0,		100,	0.5,	10),
    logscale = 	c(FALSE,	TRUE,	TRUE,	FALSE),
    objectname = "ncv20130605",
    filename = "~/Dropbox/tmp/ncv_latin_output_"
  )
  
}

# Generate and plot the mean and individual realizations for 
# infections and severe cases for two different parameter sets
# one red for very high severity and one green for low severity

# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# system.time(x <- ncv.ill.scenarios())
ncv.ill.scenarios <- function() {
  
  pHigh <- c(1.2,10,8,0.05,1,3,as.date("13-Jun-2012"),as.date("13-Mar-2013"),1)
  pLow <- c(1.2,10,8,0.05,1,3,as.date("13-Jun-2012"),as.date("13-Mar-2013"),1)
  
  rHigh <- ncv.multiple.reals(
    nreals = 100,			
    theta = array(pHigh,dim=c(7),
      dimnames=list(c("R0","k","Tg","pS","pDS","rIntro","t0")))) 
  
  plot(1:2,xlim=c(0,dim(rHigh)[1]),ylim=c(0,max(rHigh[,"Sev",])),type="n")
  for (i in 1:dim(rHigh)[3]) points(rHigh[,"Sev",i],type="l")
  
  # Upto here
  
  
}

# Not run below, but this is how to make a skeleton function 
srg.make.rcppgsl.package <- function() {
	require(inline)
	inctxt='
		#include <gsl/gsl_matrix.h>
		#include <gsl/gsl_blas.h>
			'
	bodytxt='
			RcppGSL::matrix<double> M = sM; // create gsl data structures from SEXP
			int k = M.ncol();
			Rcpp::NumericVector n(k); // to store results
			for (int j = 0; j < k; j++) {
			RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
			n[j] = gsl_blas_dnrm2(colview);
			}
			M.free() ;
			return n; // return vector
			'
	foo <- cxxfunction(signature(sM="numeric"), body=bodytxt, inc=inctxt, plugin="RcppGSL")
	## see Section 8.4.13 of the GSL manual: create M as a sum of two outer products
	M <- outer(sin(0:9), rep(1,10), "*") + outer(rep(1, 10), cos(0:9), "*")
	print(foo(M))
	package.skeleton( "mypackage",foo)
}

srg.print.format <- function(x) {
  format(signif(x,3),scientific=FALSE,big.mark=",",drop0trailing=TRUE)
}

