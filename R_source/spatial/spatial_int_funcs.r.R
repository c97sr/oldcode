# TODO: Add comment
# 
# Author: sriley
###############################################################################

dblGetCoords <- function(x,y,pg) {
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

dblMidCell <- function(xi,yi,pg) {

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

gridDist <- function(x1,y1,x2,y2,pg) {
	coord1 <- dblMidCell(x1,y1,pg)
	coord2 <- dblMidCell(x2,y2,pg)
	dist <- decLongLatDist(coord1[1],coord1[2],coord2[1],coord2[2],translate=TRUE)
	dist 
}

genDistMatrix <- function(popgrid,sampx,sampy,loud=FALSE) {

	nx 	<- length(popgrid$x)
	ny 	<- length(popgrid$y)
	tps <- data.frame(x=sampx,y=sampy)
	no_origin_pts <- dim(tps)[1]
	
	# Generate aux_lookup for distance values
	aux_lookup <- matrix(0,nx,ny)
	current_index <- 0
	aob <- 0
	
	for (i in 1:no_origin_pts) {
		if (!is.na(tps$x[i]) && !is.na(tps$x[i])) {
			tmp <- dblGetCoords(tps$x[i],tps$y[i],popgrid)
			if (tmp[1] >= 0) {
				if (aux_lookup[tmp[1],tmp[2]]==0) {
					current_index <- current_index+1
					aux_lookup[tmp[1],tmp[2]] <- current_index
				}
			} else {
				aob <- aob+1
			}
		}
	}
	
	mss <- current_index
	aux_ori_lu <- matrix(-1,mss,2)
	for (i in 1:nx) {
		for (j in 1:ny) {
			if (aux_lookup[i,j]>0) aux_ori_lu[aux_lookup[i,j],] <- c(i,j)  
		}
	}
	
	aux_dist <- matrix(-1,mss,nx*ny)
	aux_order_x <- matrix(-1,mss,nx*ny)
	aux_order_y <- matrix(-1,mss,nx*ny)
	
	count <- 0
	
	for (i in 1:nx) {
		for (j in 1:ny) {
			if (aux_lookup[i,j]>0) {
				lookindex <- aux_lookup[i,j] 
				for (k in 1:nx) {
					for (l in 1:ny) {
						aux_order_x[lookindex,(k-1)*ny+l] <- k
						aux_order_y[lookindex,(k-1)*ny+l] <- l
						aux_dist[lookindex,(k-1)*ny+l] <- gridDist(i,j,k,l,popgrid)
					}
				}
				dist_oder <- order(aux_dist[lookindex,])
				aux_order_x[lookindex,] <- aux_order_x[lookindex,dist_oder]
				aux_order_y[lookindex,] <- aux_order_y[lookindex,dist_oder]
				aux_dist[lookindex,] <- aux_dist[lookindex,dist_oder]
				if (loud) {
					count <- count + 1
					cat(paste("Completed ",count," of ", mss,"     \r"))
					flush.console()
				}
			}
		}
	}
	
	aux_vecnn <- array(0,c(mss))
	for (i in 1:mss) {
		minsqdist <- 99999
		for (j in 1:mss) {
			if (i!=j) {
				tmpdist <- gridDist(aux_ori_lu[i,1],aux_ori_lu[i,2],aux_ori_lu[j,1],aux_ori_lu[j,2],popgrid)
				if (tmpdist < minsqdist) minsqdist <- tmpdist
			}
		}
		aux_vecnn[i] <- minsqdist
	}
	
	list(	tablu=aux_lookup,		# Table of ny by nx with zero for most and >0 (index of an origin) where there is an origin
			veclu=aux_ori_lu,		# Vector of (x,y)s in which the position in the vector corresponds to the index of the origin square
			vecnn=aux_vecnn,		# A list of distances between centres of nearest origin squares
			dist=aux_dist,			# Table with nx*ny in one dimension and index of origin in the other. Cell is distance from origin square 
			ox=aux_order_x,			# As above, but cell is x coord of square in question
			oy=aux_order_y,			# similar to above.
			aob=aob					# Number of origin points outside the bounds of the popgrid
	)
}

fnKernHaz <- function(r,d,a) {
	1/(1+(d/r)^a)
}

lnlike_kern <- function(theta,dr,hx,hy,wx,wy,popgrid,popdist,loud=TRUE,rmax=5000) { 
		
	# theta <- c(offset,power,nullavedist)
	nx 	<- length(popgrid$x)
	ny 	<- length(popgrid$y)
	mss <- dim(popdist$dist)[1]
	nod <- length(hx)
	max_r_ind <- round(rmax / dr) + 1
	epsilon <- 1e-3

	# calculate sum nx kx dx for all origin squares
	matsumhaz <- array(0,c(mss,max_r_ind))
	vecsumhaz <- array(0,c(mss))
	for (i in 1:mss) {
		xi <- popdist$veclu[i,1]
		yi <- popdist$veclu[i,2]
		neicount <- 1
		max_nc <- nx*ny
		next_r <- dr
		cur_r <- -1
		dr_ind <- 1
		while (neicount <= nx*ny && cur_r < rmax) {
			aggpeop <- 0
			while (cur_r < next_r && neicount <= nx*ny && cur_r < rmax) {
				cur_r <- popdist$dist[i,neicount]
				cur_x <- popdist$ox[i,neicount]
				cur_y <- popdist$oy[i,neicount]
				sqpop <- popgrid$z[cur_x,cur_y]
				if (!is.na(sqpop)) aggpeop <- aggpeop + sqpop	
				neicount <- neicount + 1
			}
			haz <- fnKernHaz(cur_r+dr/2,theta[1],theta[2])
			matsumhaz[i,dr_ind] <- haz * aggpeop
			dr_ind <- dr_ind + 1
			next_r <- next_r + dr
		}
		vecsumhaz[i] <- sum(matsumhaz[i,])
	}
	
	# Calc the probabilities of each actual choice of workplace destination to get a likelihood
	lnlike <- 0
	i <- 1
	while (i < nod) {
		# Check that origin and destination are valid
		ocoords <- dblGetCoords(hx[i],hy[i],popgrid)
		dcoords <- dblGetCoords(wx[i],wy[i],popgrid)
		if (min(ocoords,dcoords) > 1) {
			# c for current below
			hxc <- ocoords[1]
			hyc <- ocoords[2]
			wxc <- dcoords[1]
			wyc <- dcoords[2]
			or_ind <- popdist$tablu[hxc,hyc]
			distc <- gridDist(hxc,hyc,wxc,wyc,popgrid)
			if (distc < epsilon) distc <- popdist$vecnn[or_ind] / 2
			if (distc > rmax) distc <- rmax - epsilon
			r_ind <- round((distc - distc %% dr)/dr + 1 + epsilon)
			lnlike <- lnlike + log(matsumhaz[or_ind,r_ind])-log(vecsumhaz[or_ind])
		}
		i <- i + 1
	}
	
	if (loud) {
		cat(paste(lnlike,theta[1],theta[2],"\n"))
		flush.console()
	}
	
	lnlike

}

downSampleAscii <- function(ag,samp) {
	nx <- length(ag$x) %/% samp
	ny <- length(ag$y) %/% samp
	rtnx <- array(-1,c(nx))
	rtny <- array(-1,c(ny))
	rtnz <- array(NA,c(nx,ny))
	for (i in 1:nx) rtnx[i] <- mean(ag$x[((i-1)*samp+1):((i-1)*samp+samp)])
	for (i in 1:ny) rtny[i] <- mean(ag$y[((i-1)*samp+1):((i-1)*samp+samp)])
	for (i in 1:((nx-1)*samp)) {
		for (j in 1:((ny-1)*samp)) {
			val <- ag$z[i,j]
			tmpx <- (i-1) %/% samp + 1
			tmpy <- (j-1) %/% samp + 1
			if (! is.na(val)) {
				if (is.na(rtnz[tmpx,tmpy])) rtnz[tmpx,tmpy] <- val
				else rtnz[tmpx,tmpy] <- rtnz[tmpx,tmpy] + val
			}
		}
	}
	list(x=rtnx,y=rtny,z=rtnz)
}

genDistMasks <- function(hx,hy,wx,wy,pg) {
	close <- c()
	far <- c()
	nopairs <- length(hx)
	for (i in 1:nopairs) {
		if (dblGetCoords(wx[i],wy[i],pg)[1] >= 0) close <- c(close,i)
		else far <- c(far,i)
	}
	list(close=close,far=far)
}

lnlike_kern_new <- function(theta,hx,hy,wx,wy,
							smallpopgrid,smallpopdist,bigpopgrid,bigpopdist,
							loud=TRUE,drclose=3,drwide=30) { 
						
	dmask=genDistMasks(hx,hy,wx,wy,smallpopgrid)					
						
	lnlike <- 0
	
	lnlike <- lnlike + lnlike_kern(	theta,drclose,
									hx[dmask$close],hy[dmask$close],
									wx[dmask$close],wy[dmask$close],
									smallpopgrid,smallpopdist,
									loud=FALSE)

	lnlike <- lnlike + lnlike_kern(	theta,drwide,
									hx[dmask$far],hy[dmask$far],
									wx[dmask$far],wy[dmask$far],
									bigpopgrid,bigpopdist,
									loud=FALSE)

	if (loud) {
		cat(paste(lnlike,theta[1],theta[2],"\n"))
		flush.console()
	}
	
	lnlike

}

downSampleAscii <- function(ag,samp) {
	nx <- length(ag$x) %/% samp
	ny <- length(ag$y) %/% samp
	rtnx <- array(-1,c(nx))
	rtny <- array(-1,c(ny))
	rtnz <- array(NA,c(nx,ny))
	for (i in 1:nx) rtnx[i] <- mean(ag$x[((i-1)*samp+1):((i-1)*samp+samp)])
	for (i in 1:ny) rtny[i] <- mean(ag$y[((i-1)*samp+1):((i-1)*samp+samp)])
	for (i in 1:((nx-1)*samp)) {
		for (j in 1:((ny-1)*samp)) {
			val <- ag$z[i,j]
			tmpx <- (i-1) %/% samp + 1
			tmpy <- (j-1) %/% samp + 1
			if (! is.na(val)) {
				if (is.na(rtnz[tmpx,tmpy])) rtnz[tmpx,tmpy] <- val
				else rtnz[tmpx,tmpy] <- rtnz[tmpx,tmpy] + val
			}
		}
	}
	list(x=rtnx,y=rtny,z=rtnz)
}
