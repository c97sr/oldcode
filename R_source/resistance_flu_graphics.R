rm(list=ls(all=TRUE))
fileroot = "D:\\files\\projects"
setwd(paste(fileroot,"\\..\\eclipse\\R Source",sep=""))
popsize <- 6800000

read_joe_freq <- function(fn) {
	dFig3a_tmp  <- read.csv(file=fn,header=FALSE,row.names=NULL)
	rtnz <- array(NA,c(dim(dFig3a_tmp)[1]-1,dim(dFig3a_tmp)[2]-1))
	for (i in 1:dim(rtnz)[1]) for (j in 1:dim(rtnz)[2]) {
		tmp <- dFig3a_tmp[i+1,j+1]
		if (tmp > 1e-10) rtnz[i,j] <- dFig3a_tmp[i+1,j+1]
	}
	rtnx <- array(NA,dim(dFig3a_tmp)[2]-1)
	rtny <- array(NA,dim(dFig3a_tmp)[1]-1)
	f3a_xbw <- dFig3a_tmp[1,3] - dFig3a_tmp[1,2]
	f3a_ybw <- dFig3a_tmp[3,1] - dFig3a_tmp[2,1]
	rtnx[1] <- dFig3a_tmp[1,2]- f3a_xbw / 2
	rtny[1] <- dFig3a_tmp[2,1]- f3a_ybw / 2
	for (i in 2:length(rtnx)) rtnx[i] = rtnx[i-1] + f3a_xbw
	for (i in 2:length(rtny)) rtny[i] = rtny[i-1] + f3a_ybw
	list(rtnx,rtny,rtnz)
}

fFig1a = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig1a.csv",sep="")
fFig1b = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig1b.csv",sep="")
fFig1c_sr = paste(fileroot,"\\influenza\\results\\resistance\\070726\\1e6_ResistantAR.csv",sep="")
fFig1c_lr = paste(fileroot,"\\influenza\\results\\resistance\\070726\\HK_ResistantAR.csv",sep="")
fFig1c_sw = paste(fileroot,"\\influenza\\results\\resistance\\070726\\1e6_AR.csv",sep="")
fFig1c_lw = paste(fileroot,"\\influenza\\results\\resistance\\070726\\HK_AR.csv",sep="")
fFig2_sp_ch = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig2_stockpile_column_headers.csv",sep="")
fFig2_sy_rh = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig2_synergy_row_headers.csv",sep="")
fFig2a = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig2AR.csv",sep="")
fFig2b = paste(fileroot,"\\influenza\\results\\resistance\\070726\\fig2ResistantAR.csv",sep="")
fFig3_sa_p1 = paste(fileroot,"\\influenza\\results\\resistance\\070726\\PolicyA_AR_ResAR_hist_P1.csv",sep="")
fFig3_sa_p2 = paste(fileroot,"\\influenza\\results\\resistance\\070726\\PolicyA_AR_ResAR_hist_P2.csv",sep="")
fFig3_sb_p1 = paste(fileroot,"\\influenza\\results\\resistance\\070726\\PolicyB_AR_ResAR_hist_P1.csv",sep="")
fFig3_sb_p2 = paste(fileroot,"\\influenza\\results\\resistance\\070726\\PolicyB_AR_ResAR_hist_P2.csv",sep="")

dFig1a <- read.table(file=fFig1a,header=FALSE,sep=",",row.names=NULL)
names(dFig1a)<-c("time","wild","resistant")

dFig1b <- read.table(file=fFig1b,header=FALSE,sep=",",row.names=NULL)
dFig1b[,1] <- dFig1b[,1]/6800000
names(dFig1b)<-c("cumAtt","allAtt","resAtt")

dFig1c_sr <- read.table(file=fFig1c_sr,header=FALSE,sep=",",row.names=NULL)
names(dFig1c_sr)<-c("rate","five","fifty","ninetyfive")
dFig1c_lr <- read.table(file=fFig1c_lr,header=FALSE,sep=",",row.names=NULL)
names(dFig1c_lr)<-c("rate","five","fifty","ninetyfive")
dFig1c_sw <- read.table(file=fFig1c_sw,header=FALSE,sep=",",row.names=NULL)
names(dFig1c_sw)<-c("rate","five","fifty","ninetyfive")
dFig1c_lw <- read.table(file=fFig1c_lw,header=FALSE,sep=",",row.names=NULL)
names(dFig1c_lw)<-c("rate","five","fifty","ninetyfive")

dFig2_sp_ch  <-read.table(file=fFig2_sp_ch,header=FALSE,sep=",",row.names=NULL)
fig2_nocols <- dim(dFig2_sp_ch)[2]-1
dFig2_sp_ch <- as.vector(t(dFig2_sp_ch[1,1:fig2_nocols]))
dFig2_sy_rh <- as.vector(read.table(file=fFig2_sy_rh,header=FALSE,sep=",",row.names=NULL))
dFig2_sy_rh <- rev(dFig2_sy_rh[,1])
fig2_norows <- length(dFig2_sy_rh)
dFig2a_tmp <- read.csv(file=fFig2a,header=FALSE,row.names=NULL,flush=TRUE)
dFig2b_tmp <- read.table(file=fFig2b,header=FALSE,sep=",",row.names=NULL)
dFig2a = array(0,c(fig2_norows,fig2_nocols))
dFig2b = array(0,c(fig2_norows,fig2_nocols))
for (i in 1:fig2_nocols) {
	for (j in fig2_norows:1) {
		dFig2a[j,i] <- dFig2a_tmp[j,i]
		dFig2b[j,i] <- dFig2b_tmp[j,i]
	}		
}

df3_sa_p1 <- read_joe_freq(fFig3_sa_p1)
df3_sa_p2 <- read_joe_freq(fFig3_sa_p2)
df3_sb_p1 <- read_joe_freq(fFig3_sb_p1)
df3_sb_p2 <- read_joe_freq(fFig3_sb_p2)

# Plot for figure 1
windows(width=11.5/cm(1),height=10/cm(1))
par(mai=c(0,0,0,0),xaxs="i",yaxs="i",xpd=FALSE)
par(mgp=c(0,0,0))
par(cex=1)

par(fig=c(0.125,0.95,0.10,0.95))
plot(1:2,type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,200),ylim=c(0,0.05))
axis(	1,pos=0,tck=-0.01,at=c(0,50,100,150,200),
		labels=c(0,50,100,150,200),outer=TRUE,padj=0.1)
mtext("Time (days)",side=1,line=1)
axis(	2,pos=0,tck=-0.01,las=1,at=(0:5)/100,hadj=1.2,
		labels=c("0%","1%","2%","3%","4%","5%"),outer=TRUE,line=2)
mtext("Daily Incidence (% of population)",side=2,line=2)
mtext("A",side=2,at=0.05,cex=1.2,font=2,las=1,line=2)

polygon(c(dFig1a[,1],rev(dFig1a[,1])),c(rep(0,dim(dFig1a)[1]),rev(dFig1a[,2]/popsize)),
			col="blue",border=NA)
polygon(c(dFig1a[,1],rev(dFig1a[,1])),c(rep(0,dim(dFig1a)[1]),rev(dFig1a[,3]/popsize)),
			col="red",border=NA)

par(fig=c(0.7,0.95,0.675,0.925),new=TRUE)
plot(1:2,type="n",axes=FALSE,xlab="",ylab="",log="x",xlim=c(1e-5,1),ylim=c(0,0.75),new=FALSE)
axis(	1,pos=0,tck=-0.05,outer=TRUE,at=10^(-5:0),
		label=expression("",10^-4,"",10^-2,"",10^0),line=2,cex.axis=0.8)
mtext("Cumulative incidence\nantivirals started",side=1,line=1.7,cex=0.8)
axis(	2,pos=1e-5,tck=-0.05,outer=TRUE,at=c(0,0.25,0.5,0.75),
		label=expression("0%","25%","50%","75%"),las=1,cex.axis=0.8,hadj=1.2)
mtext("Attack rate",side=2,line=1.9,cex=0.8)
mtext("B",side=2,at=0.75,cex=1,font=2,las=1,line=1.9)
lines(dFig1b[,"cumAtt"],dFig1b[,"allAtt"])
lines(dFig1b[,"cumAtt"],dFig1b[,"resAtt"],lty=2)

par(fig=c(0.7,0.95,0.25,0.5),new=TRUE)
plot(1:2,type="n",axes=FALSE,xlab="",ylab="",log="x",xlim=c(1e-5,1),ylim=c(0,0.75),new=FALSE)
axis(	1,pos=0,tck=-0.05,outer=TRUE,at=10^(-5:0),
		label=expression("",10^-4,"",10^-2,"",10^0),line=2,cex.axis=0.8)
mtext("Emergence rate",side=1,line=1,cex=0.8)
axis(	2,pos=1e-5,tck=-0.05,outer=TRUE,at=c(0,0.25,0.5,0.75),
		label=expression("0%","25%","50%","75%"),las=1,cex.axis=0.8,hadj=1.2)
mtext("Resistant\nattack rate",side=2,line=2.0,cex=0.8)
mtext("C",side=2,at=0.75,cex=1,font=2,las=1,line=1.9)
lines(dFig1c_lr[,"rate"],dFig1c_lr[,"fifty"],col="cyan")
lines(dFig1c_lr[,"rate"],dFig1c_lr[,"five"],col="cyan",lty=2)
lines(dFig1c_lr[,"rate"],dFig1c_lr[,"ninetyfive"],col="cyan",lty=2)
lines(dFig1c_sr[,"rate"],dFig1c_sr[,"fifty"],col="magenta")
lines(dFig1c_sr[,"rate"],dFig1c_sr[,"five"],col="magenta",lty=2)
lines(dFig1c_sr[,"rate"],dFig1c_sr[,"ninetyfive"],col="magenta",lty=2)

# Plot for figure 2
windows(width=12/cm(1),height=16/cm(1))
par(	mai=c(0,0.05,0,0.25), 	# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,	
		mgp=c(0,0,0),
		par(cex=1)	)

par(fig=c(0.1,0.9,0.5,1.0))
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

ytp <- c(-5,-4,-3,-2,-1,0) # x tick places
ytl <- expression("",10^-4,"",10^-2,"",10^0) # x tick labels
for (i in 1:length(ytp)) lines(trans3d(c(0,+0.15),c(ytp[i],ytp[i]),c(0,0),tvm)) 
for (i in 1:length(ytl)) text(trans3d(+0.15,ytp[i],0,tvm),ytl[i],adj=c(1.1,0.8),cex=0.8)
text(trans3d(+1.2,-1.75,0,tvm),"Synergy",cex=1.0)

ztp <- c(0,0.25,0.5,0.75) 					# x tick places
ztl <- expression("0%","25%","50%","75%") 	# x tick labels
for (i in 1:length(ztp)) lines(trans3d(c(0,+0.15),c(-5,-5),c(ztp[i],ztp[i]),tvm)) 
for (i in 1:length(ztp)) lines(trans3d(c(-5,-5),c(0,+0.15),c(ztp[i],ztp[i]),tvm)) 
for (i in 1:length(ztl)) text(trans3d(-5,+0.15,ztp[i],tvm),ztl[i],adj=c(-0.1,0.5),cex=0.8)
text(trans3d(-5,0,0.9,tvm),"Resistant attack\nrate",cex=1.0)
text(trans3d(0,-5,0.9,tvm),"A",cex=1.2,font=2)

par(fig=c(0.1,0.9,0,0.5),new=TRUE)
tvm <- persp(log(dFig2_sp_ch,10),log(dFig2_sy_rh,10),dFig2a,theta = 135, 
		phi = 20, r=sqrt(3), expand = 0.5, col = "lightblue",
      	ltheta = 120, shade = 0.75, ticktype = "detailed",
      	xlab = "X", ylab = "Y", zlab = "Sinc( r )",
      	zlim = c(0,0.75),axes=FALSE)
      	
xtp <- c(-5,-4,-3,-2,-1,0) 									# x tick places
xtl <- expression("",10^-4,"",10^-2,"",10^0) 				# x tick labels
for (i in 1:length(xtp)) lines(trans3d(c(xtp[i],xtp[i]),c(0,+0.15),c(0,0),tvm))
for (i in 1:length(xtl)) text(trans3d(xtp[i],+0.15,0,tvm),xtl[i],adj=c(0,0.8),cex=0.8)
text(trans3d(-1.75,+1.2,0,tvm),"Stockpile \n size",cex=1.0)

ytp <- c(-5,-4,-3,-2,-1,0) # x tick places
ytl <- expression("",10^-4,"",10^-2,"",10^0) # x tick labels
for (i in 1:length(ytp)) lines(trans3d(c(0,+0.15),c(ytp[i],ytp[i]),c(0,0),tvm)) 
for (i in 1:length(ytl)) text(trans3d(+0.15,ytp[i],0,tvm),ytl[i],adj=c(1.1,0.8),cex=0.8)
text(trans3d(+1.2,-1.75,0,tvm),"Synergy",cex=1.0)

ztp <- c(0,0.25,0.5,0.75) 					# x tick places
ztl <- expression("0%","25%","50%","75%") 	# x tick labels
for (i in 1:length(ztp)) lines(trans3d(c(0,+0.15),c(-5,-5),c(ztp[i],ztp[i]),tvm)) 
for (i in 1:length(ztp)) lines(trans3d(c(-5,-5),c(0,+0.15),c(ztp[i],ztp[i]),tvm)) 
for (i in 1:length(ztp)) lines(trans3d(c(0,+0.15),c(0,0),c(ztp[i],ztp[i]),tvm))
for (i in 1:length(ztl)) text(trans3d(-5,+0.15,ztp[i],tvm),ztl[i],adj=c(-0.1,0.5),cex=0.8)
text(trans3d(-5,0,0.9,tvm),"Overall attack\nrate",cex=1.0)
text(trans3d(0,-5,0.9,tvm),"B",cex=1.2,font=2)

# Plot for Figure 3
windows(width=12/cm(1),height=16/cm(1))
par(	mai=c(0,0,0,0),		 	# c(bottom, left, top, right)
		xaxs="i",				# Precise axis lengths - no margins on the plot area
		yaxs="i",
		xpd=NA,	
		mgp=c(3,0,0),			# Line gaps on axes, title, labels, axis
		par(cex=1)	)
		
f3cols <- heat.colors(100)
f3xlim <- c(0.1,0.8)
xat <- 1:8/10 
xlab <- c("","20%","","40%","","60%","","80%")
f3ylim <- c(.55,.75)
yat <- 0.55 + 0:4*0.05
ylab <- c("","60% ","","70% ","")
f3col <- heat.colors(100)
subfiglabcex <- 1
subfiglabfont <- 2
subfiglabline <- -1.25
f3tl <- 0.02
lgap <- 0.175
bgap <- 0.1
tgap <- 0.05
rgap <- 0.05
xw <- (1-lgap-rgap)/1
yw <- (1-tgap-bgap)/4
yaxlab_x <- lgap/2
yaxlab_y <- bgap+2*yw
xaxlab_x <- lgap+xw/2
xaxlab_y <- bgap/2
x1 <- lgap
x2 <- 1 - rgap
y1 <- bgap
y2 <- bgap + 1*yw
y3 <- bgap + 2*yw
y4 <- bgap + 3*yw
y5 <- bgap + 4*yw

chtpos <- list(	c(x1,x2,y4,y5),
				c(x1,x2,y3,y4),
				c(x1,x2,y2,y3),
				c(x1,x2,y1,y2)	)

par(fig=chtpos[[1]])
image(	df3_sa_p1[[1]],df3_sa_p1[[2]],t(df3_sa_p1[[3]]),
		ylim=f3ylim,xlim=f3xlim,col=f3cols,axes=FALSE,
		xlab="",ylab="")
axis(1,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(3,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(2,at=yat,labels=ylab,tck=f3tl,las=1)
axis(4,at=yat,labels=rep("",length(ylab)),tck=f3tl)
mtext(	"A",side=2,at=f3ylim[2]-0.025,cex=subfiglabcex,
		font=subfiglabfont,las=1,line=subfiglabline)

par(fig=chtpos[[2]],new=TRUE)
image(	df3_sa_p2[[1]],df3_sa_p2[[2]],t(df3_sa_p2[[3]]),
		ylim=f3ylim,xlim=f3xlim,col=f3cols,axes=FALSE,
		xlab="",ylab="")
axis(1,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(3,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(2,at=yat,labels=ylab,tck=f3tl,las=1)
axis(4,at=yat,labels=rep("",length(ylab)),tck=f3tl)
mtext(	"B",side=2,at=f3ylim[2]-0.025,cex=subfiglabcex,
		font=subfiglabfont,las=1,line=subfiglabline)

par(fig=chtpos[[3]],new=TRUE)
image(	df3_sb_p1[[1]],df3_sb_p1[[2]],t(df3_sb_p1[[3]]),
		ylim=f3ylim,xlim=f3xlim,col=f3cols,axes=FALSE,
		xlab="",ylab="")
axis(1,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(3,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(2,at=yat,labels=ylab,tck=f3tl,las=1)
axis(4,at=yat,labels=rep("",length(ylab)),tck=f3tl)
mtext(	"C",side=2,at=f3ylim[2]-0.025,cex=subfiglabcex,
		font=subfiglabfont,las=1,line=subfiglabline)

par(fig=chtpos[[4]],new=TRUE)
image(df3_sb_p2[[1]],df3_sb_p2[[2]],t(df3_sb_p2[[3]]),
		ylim=f3ylim,xlim=f3xlim,col=f3cols,axes=FALSE,
		xlab="",ylab="")
axis(1,at=xat,labels=xlab,tck=f3tl)
axis(3,at=xat,labels=rep("",length(xlab)),tck=f3tl)
axis(2,at=yat,labels=ylab,tck=f3tl,las=1)
axis(4,at=yat,labels=rep("",length(ylab)),tck=f3tl)
mtext(	"D",side=2,at=f3ylim[2]-0.025,cex=subfiglabcex,
		font=subfiglabfont,las=1,line=subfiglabline)

mtext("Resistant attack rate",side=1,font=2,cex=1.2,line=1.5)
mtext("Overall attack rate",side=2,at=0.95,font=2,cex=1.2,line=2.5)
