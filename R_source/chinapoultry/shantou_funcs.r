require(date)

srSCont <- function(s1,s2) {
	if (s1=="" || s2=="") rtn <- FALSE
	else {
		head <- strsplit(s1,s2)[[1]][1]
		if (head != s1) rtn <- TRUE
		else rtn <- FALSE
	}
	rtn
}

srSContOr <- function(s1,vecs2) {
	return <- FALSE
	for (s2 in vecs2) return <- return || srSCont(s1,s2)
	return
}

genSummary <- function(	db,
		vecBinBounds=as.date("1Jan2009"):as.date("31Jan2009"),
		hterms=c("H9"),
		tterms=c("CK")	) {
	
	# Takes a database db, a sart date, an approximate end date, a regular time interval dt,
	# a list of definative subsype strings and a list of definative host type stings
	
	# Current to dos in this file
	# - Should deal better with other and nonull categories
	# - Make a few simple charts and then come back to this to get a more complete summary of these data
		
	# Set up the column headings
	nameSummaryVars <- c()
	nohterms <- length(hterms)
	notterms <- length(tterms)
	for (a in c("N",hterms,"nonnull")) for (b in c(tterms,"_other")) nameSummaryVars <- c(nameSummaryVars,paste(b,a,sep="."))
	noSummaryVars <- length(nameSummaryVars)

	# Set up the date ranges
	# approxInterval <- approxEndDate - startDate
	# noTimeBins <- floor((approxEndDate-startDate) / dt)
	# noTimePoints <- noTimeBins + 1
	# endDate <- startDate + dt*noTimeBins
	# sumDat <- as.data.frame(array(0,c(noTimeBins,noSummaryVars)))
	# names(sumDat) <- nameSummaryVars
	# vecBinStart <- (0:(noTimeBins-1))*dt+startDate
	# vecBinEnd <- vecBinStart + dt 
	
	# Set up the date ranges
	noTimeBins 		<- length(vecBinBounds)-1
	vecBinStart 	<- vecBinBounds[1:(noTimeBins)]
	vecBinEnd 		<- vecBinBounds[2:(noTimeBins+1)]
	sumDat 			<- as.data.frame(array(0,c(noTimeBins,noSummaryVars+2)))
	names(sumDat) 	<- c("StartBin","EndBin",nameSummaryVars)
	sumDat[,"StartBin"] <- vecBinStart
	sumDat[,"EndBin"]	<- vecBinEnd
		
	i <- 1
	norows <- dim(db)[1]
	while (i <= norows) {
		if (!is.na(db$Date[i])) {
	
			# Find the correct bin for the data
			# dbin <- floor((db$Date[i]-startDate) / dt) + 1
			if (db$Date[i] < vecBinBounds[1] || db$Date[i] >= vecBinBounds[noTimeBins + 1] ) dbin <- 0
			else {
				dbin <- 1
				while (db$Date[i] >= vecBinEnd[dbin]) dbin <- dbin + 1
			}
			
			# If the time bin falls within the range being considered
			if (dbin != 0) {
				
				# Assign temporary values for the subtype and the host type
				charType1 <- as.character(db$Subtype[i])
				charHost1 <- as.character(db$CleanType[i])
				
				# If I've already classified this as a first of a pair or an only
				# note that the pair checking gets rid of (F) samples
				if (db$Pairs[i]=="First" || db$Pairs[i]=="Only") {
					isother <- TRUE
					for (y in tterms) { 
						if (srSCont(charHost1,y)) {
							isother <- FALSE
							sumDat[dbin,paste(y,"N",sep=".")] <- sumDat[dbin,paste(y,"N",sep=".")]+1
							istestedh <- FALSE
							for (x in hterms) {
								if (db$Pairs[i]=="First") {
									charType2 <- as.character(db$Subtype[i+1])
									if ((srSCont(charType1,x) || srSCont(charType2,x))) {
										sumDat[dbin,paste(y,x,sep=".")] <- sumDat[dbin,paste(y,x,sep=".")] + 1	
										istestedh <- TRUE
									}  
								} else {
									if (srSCont(charType1,x)) {
										sumDat[dbin,paste(y,x,sep=".")] <- sumDat[dbin,paste(y,x,sep=".")] + 1
										istestedh <- TRUE
									} 
									
								}
							}
							if (!istestedh) {
								sumDat[dbin,paste(y,"nonnull",sep=".")] <- sumDat[dbin,paste(y,"nonnull",sep=".")] + 1
							}
						}
					}
					if (isother) sumDat[dbin,paste("_other","N",sep=".")] <- sumDat[dbin,paste("_other","N",sep=".")]+1
					if (db$Pairs[i]=="First") i <- i+1
				}				
			}
		}
		i <- i+1
	}
	
	# Return the summary table 
	sumDat
	
}

extractShantouData <- function(yearsinclude,matchedtf,matchedsub,blocksize=10000,readrows=-1) {
	
	arraysize <- blocksize
	mdtnames <- c("Date","Number","Region","Type","Subtype") # Main data names
	mdt <- as.data.frame(array(NA,c(blocksize,length(mdtnames))))
	names(mdt) <- mdtnames
	

	currentrow <- 1 	
	for (j in 1:length(yearsinclude)) {
		y <- yearsinclude[j]
		cy <- read.csv(file=fnFName(y),header=TRUE,nrows=readrows)
		norows <- dim(cy)[1]
		for (i in 1:norows) {
			if (!is.na(cy$Sample.Date[i])) {
				mdt$Date[currentrow] <- as.date(as.character(cy$Sample.Date[i]),order="mdy")
				mdt$Number[currentrow] <- cy$Field.Number[i]
				mdt$Region[currentrow] <- as.character(cy$Sample.Region[i])
				mdt$Type[currentrow] <- as.character(cy[i,matchedtf[j]])
				mdt$Subtype[currentrow] <- as.character(cy[i,matchedsub[j]])
				currentrow <- currentrow + 1
				if (currentrow %% 1000 == 0) {
					cat("Rows done: ",currentrow," of ",norows,"              \r");
					flush.console()
				}
				if (currentrow > arraysize) {
					tmpframe <- data.frame(array(NA,c(blocksize,length(mdtnames)))) 
					names(tmpframe) <- mdtnames
					rownames(tmpframe) <- currentrow:(currentrow+blocksize-1)
					mdt <- rbind(mdt,tmpframe)
					arraysize <- arraysize + blocksize
				}
			}
		}
	}
	
	cat("\n")
	
	mdtout <- as.data.frame(array(NA,c(currentrow-1,length(mdtnames))))
	mdtout[1:(currentrow-1),] <- mdt[1:(currentrow-1),]
	names(mdtout) <- mdtnames
	mdtout
	
}

checkPairs <- function(x) {
	norows <- dim(x)[1]
	pairvec <- array(0,c(norows))
	i = 1
	while (i < norows) {
		type1 <- as.character(x$Type[i])
		type2 <- as.character(x$Type[i+1])
		head1 <- strsplit(type1,"\\(T\\)")
		head2 <- strsplit(type2,"\\(C\\)")
		if (type1 != ""					&&
			type2 != ""					&&
			head1[[1]][1] != type1 		&& 
			head1[[1]] == head2[[1]] 	&&
			is.na(x$Date[i]) != TRUE	&&
			is.na(x$Date[i+1]) != TRUE	&&
			x$Date[i] == x$Date[i+1]) {
			pairvec[i] <- "First"
			pairvec[i+1] <- "Second"
			i <- i + 2
		} else {
			if ( srSCont(type1,"\\(T\\)") ||
				 srSCont(type1,"\\(C\\)") ) pairvec[i] <- "Only"
			else pairvec[i] <- "None"
			i <- i+1
		}
	}
	if (i == norows) pairvec[norows] <- "Only"
	pairvec
}

cleanTypes <- function(x) {
	norows <- dim(x)[1]
	cleantypes <- array(0,c(norows))
	i = 1
	while (i <= norows) {
		type <- as.character(x$Type[i])
		if (srSContOr(type,c("Ck","CK","Chicken"))) cleantypes[i] <- "Chicken"
		else if (srSContOr(type,c("DK","Dk","Duck","Gs","GS","Goose"))) cleantypes[i] <- "Waterfowl"
		else if (srSContOr(type,c("Quail","Qa"))) cleantypes[i] <- "Quail"
		else cleantypes[i] <- "Other"
		i <- i+1
	}
	cleantypes
}

curPos <- function(index) {
	rtn <- sr_chart_pos(1,index,1,notypes,
			xlm=2.0/cm(1)/wwidth,xrm=0.5/cm(1)/wwidth,xg=0.5/cm(1)/wwidth,
			ybm=2.0/cm(1)/wheight,ytm=0.5/cm(1)/wheight,yg=0.5/cm(1)/wheight)
	rtn
}

# Start of the idfferential equation model is below here
shantouModel <- function(t,y,p) {
	dy <- vector(mode="numeric",length(y))
	dy[1] <- 	-	p["beta"]*y[1]*y[2]
	dy[2] <- 	+	p["beta"]*y[1]*y[2] - p["gamma"]*y[2]
	dy[3] <- 		p["gamma"]*y[2]
	list(dy)
}

			


