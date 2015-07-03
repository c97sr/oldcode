# To Dos here
# - Include 1D and 2D histogram plots in the standard report

# Clear the workspace and source mics functions
rm(list=ls())
source("SR_misc.r")

strFileDir <- "~/tmpresults/090417/"		# Working directory
# strFileDir <- "/Volumes/ide01.dynalias.com/results/090305/"
vf <- dir(strFileDir,pattern='trace.out') 	# Generate a list of files
vf <- paste(strFileDir,vf,sep="")			# Add the path back to the filenames

# Load up the data from the file with some thining if required
od_max <- sr_load_mcmc(vf[1],10)

# Set the proportion of the sample to be plotted
samp_prop <- 0.15
maxsample <- dim(od_max)[1]
od <- od_max[max(1,(1-samp_prop)*maxsample):maxsample,]

# Clear windows and plot listed parameters
# Index lnl noInf noSus noImm pAcT pAcP pAcI q_base q_int_1 t_1 i_1
pnames <- c("lnl","q_base","q_int_1","noInf")
nodev <- length(dev.list()); if (nodev > 0) for (i in 1:nodev) dev.off()
for (p in pnames) {
	X11(); plot(od$index,od[,p],type="l",main=p)
}

x <- table(od$noInf)
plot(od$q_int_1,od$q_int_1)
x <- sr1DHist(od$noInf,bins=350-309,min=309.5,max=350.5,login=FALSE)
srPlot1DHist(x,logarg="x")

dim(od)
names(od)

h1 <- sr2DHist(	od_dash$noInf,od_dash$q_base,
				min_x=0,max_x=40,bins_x=40,min_y=0,max_y=1,bins=40)
		
