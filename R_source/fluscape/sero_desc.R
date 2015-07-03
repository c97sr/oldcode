rm(list=ls(all=TRUE))
options(error=NULL)
require(lattice)


ChartLocSeroData <- function(tab) {
	rtn <- tab
	norows <- dim(tab)[1]
	nocols <- dim(tab)[2]
	for (i in 1:norows) {
		for (j in 3:nocols) {
			if (tab[i,j] < 0) rtn[i,j] <- 0
			else rtn[i,j] <- tab[i,j] / max(tab[1:norows,3:nocols])
			# browser()
		}
	}
	rtn
}

require("date")

tmp <- ChartLocSeroData(dat[dat$Test=="VN",])

# Person Test H3N2.ST.90.2003 H3N2.ST.806.2005 H3N2.ST.904.2008
# 1 P1H01P01   HI         0.03125           0.0625          0.03125
# H1N1.ST.104.2005 H1N1.ST.92.2009

cloud(log(tmp$H3N2.ST.90.2003+1)~log(tmp$H3N2.ST.806.2005+1)*log(tmp$H3N2.ST.904.2008+1))
dat <- read.csv(("~/files/projects/influenza/Fluscape/data/SeraFirst16.csv"))

# Plot for figure 2
windows(width=12/cm(1),height=16/cm(1))
par(	mai=c(0,0.05,0,0.25), 	# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,	
		mgp=c(0,0,0),
		par(cex=1)	)

par(fig=c(0.1,0.9,0.5,1.0))

tvm <- persp(log(tmp$H3N2.ST.806.2005+1),log(tmp$H3N2.ST.904.2008+1),log(tmp$H3N2.ST.90.2003+1),
		theta = 135,phi = 20, r=sqrt(3), expand = 0.5, col = "lightblue",
		ltheta = 120, shade = 0.75, ticktype = "detailed",
		xlab = "X", ylab = "Y", zlab = "Sinc( r )",
		zlim = c(0,0.75),axes=FALSE)

tvm <- persp(log(dFig2_sp_ch,10),log(dFig2_sy_rh,10),dFig2b,theta = 135, 
		phi = 20, r=sqrt(3), expand = 0.5, col = "lightblue",
		ltheta = 120, shade = 0.75, ticktype = "detailed",
		xlab = "X", ylab = "Y", zlab = "Sinc( r )",
		zlim = c(0,0.75),axes=FALSE)

xtp <- c(-5,-4,-3,-2,-1,0) 									# x tick places
xtl <- expression("",10^-4,"",10^-2,"",10^0) 				# x tick labels
for (i in 1:length(xtp)) lines(trans3d(c(xtp[i],xtp[i]),c(0,+0.15),c(0,0),tvm))
for (i in 1:length(xtl)) text(trans3d(xtp[i],+0.15,0,tvm),xtl[i],adj=c(0,0.8),cex=0.8)
text(trans3d(-1.75,+1.2,0,tvm),"Stockpile \n size",cex=1.0)

