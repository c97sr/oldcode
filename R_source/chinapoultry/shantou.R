# Next 	- make up a production quality table for each year of the study for the types of interest
# 		- need to be able to handle sets of strings for data, do by a preprocessing
# 	 	- construct a simple but efficient and generalized 2 type two host differential equation model
# 		- model needs to be able to be solved by odeint, but shouldn't necessarily be slow

# if needed setwd("/Users/sriley/Dropbox/svneclipse/R Source/chinapoultry/")
rm(list=ls(all=TRUE))
options(error=recover)
source("./shantou_funcs.r")
source("/Users/sriley/Dropbox/svneclipse/R Source/SR_misc.r")


# fInt 	<- "/influenza/shantou/intermediate_files/dummy_intermediate.csv"
fInt 	<- "/Volumes/NO\ NAME/data/influenza/shantou/intermediate_files/070829.csv"
fTab1 	<- "./tab1.csv"

mdt 		<- read.csv(fInt)
col 		<- checkPairs(mdt)
cleantype 	<- cleanTypes(mdt) 
mdt 		<- cbind(mdt,Pairs=col,CleanType=cleantype)

subtypes <- c("H5","H6","H9")
hosttypes <- c("Chicken","Waterfowl","Quail","Other")
notypes <- length(subtypes)

binbounds <- c()
for (year in c("2001","2002","2003","2004","2005","2006","2007"))
	for (month in c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")) 
		for (day in c("1")) 
			binbounds <- c(binbounds,paste(day,month,year,sep=""))
binbounds <- as.date(binbounds)
sdate <- binbounds[1]
edate <- binbounds[length(binbounds)]

inc <- genSummary(mdt,hterms=subtypes,tterms=hosttypes,vecBinBounds=binbounds)
write.csv(inc,file=fTab1)

# Plotting stuff
noobs <- dim(inc)[1]

dlabs <- c(	"1Jan1999","1July1999",
			"1Jan2000","1July2000",
			"1Jan2001","1July2001",
			"1Jan2002","1July2002",
			"1Jan2003","1July2003",
			"1Jan2004","1July2004",
			"1Jan2005","1July2005",
			"1Jan2006","1July2006",
			"1Jan2007","1July2007")

nodlabs <- length(dlabs)
dlabsspace <- dlabs
for (i in 1:nodlabs) dlabsspace[i] <- paste(dlabsspace[i]," ",sep="")

minprev <- 0.0001

wwidth <- 20.5/cm(1)
wheight <- 10/cm(1)
pdf(file="./sr_R_plot_1.pdf",width=wwidth,height=wheight)
# X11(width=wwidth,height=wheight)
# par(mai=c(0,0,0,0),xaxs="i",yaxs="i",xpd=FALSE)
# par(mgp=c(0,0,0))
par(cex=1)

par(cex.axis=0.8)
plot(	1:2,type="n",axes=FALSE,xlab="",ylab="",
		xlim=c(sdate,edate),ylim=c(minprev,1),log="y")
axis(	2,pos=sdate,tck=-0.01,las=1,at=10^(-4:0),
		labels=c("0.00001 ","","0.001 ","","1.0 "),line=2)
# mtext(paste(type,host),side=2,at=1,cex=1.0,font=2,las=-1,line=2.0)
axis(	1,pos=minprev,tck=-0.01,at=as.date(dlabs),
		labels=dlabs,outer=TRUE,padj=0.5,las=2,line=1)
# mtext("Time (days)",side=1,line=1)	
chttype <- c("H9","H9")
chthost <- c("Waterfowl","Chicken")
chtcol 	<- c("red","blue")
chtoff	<- c(-5,+5)

for (i in 1:length(chttype)) {
	for (j in 1:noobs) 
		lines(	c((inc$StartBin[j]+inc$EndBin[j])/2+chtoff[i],(inc$StartBin[j]+inc$EndBin[j])/2+chtoff[i]),
				c(	binCI(inc[j,paste(chthost[i],"N",sep=".")],inc[j,paste(chthost[i],chttype[i],sep=".")],0.025),
					binCI(inc[j,paste(chthost[i],"N",sep=".")],inc[j,paste(chthost[i],chttype[i],sep=".")],0.975)),
				col=chtcol[i])
}

# savePlot(filename="~/files/tmp/tmp.png",type="png")

dev.off()
