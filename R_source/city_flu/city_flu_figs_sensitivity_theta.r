rm(list=ls(all=TRUE))
source("city_flu_figs_header.r")

root = "D:\\share\\files\\projects\\flu2005\\results\\060330_theta_sens\\"
pfile = "D:\\share\\files\\projects\\flu2005\\results\\060330_theta_sens\\ParameterList.txt"
params <- read.table(pfile,header=TRUE)
root_gam = "D:\\share\\files\\projects\\flu2005\\results\\060330_gamma_sens\\"
pfile_gam = "D:\\share\\files\\projects\\flu2005\\results\\060330_gamma_sens\\ParameterList.txt"
params_gam <- read.table(pfile,header=TRUE)
stem = "set"
leaf = ".txt"

sind <- 1
eind <- 100
max_i = eind-sind+1
ps = 1

res_None <- flu_processor_dash(root,stem,leaf,"None",sind,eind)
res_Q <- flu_processor_dash(root,stem,leaf,"Q",sind,eind)
res_QI <- flu_processor_dash(root,stem,leaf,"QI",sind,eind)
res_QIA <- flu_processor_dash(root,stem,leaf,"QIA",sind,eind)
res_QIAT <- flu_processor_dash(root,stem,leaf,"QIAC",sind,eind)
res_QIATC <- flu_processor_dash(root,stem,leaf,"QIAC",sind,eind)

res_None_gam <- flu_processor_dash(root_gam,stem,leaf,"None",sind,eind)
res_Q_gam <- flu_processor_dash(root_gam,stem,leaf,"Q",sind,eind)
res_QI_gam <- flu_processor_dash(root_gam,stem,leaf,"QI",sind,eind)
res_QIA_gam <- flu_processor_dash(root_gam,stem,leaf,"QIA",sind,eind)
res_QIAT_gam <- flu_processor_dash(root_gam,stem,leaf,"QIAC",sind,eind)
res_QIATC_gam <- flu_processor_dash(root_gam,stem,leaf,"QIAC",sind,eind)

# Plot routines for figure S5
windows(width=6,height=4)
popsize=250000
par(xpd=FALSE)

top_mar = 0.03
x_gap = 0.15
y_mar = 0.4
chart_height = (1-top_mar-y_mar)
x_mar = 0.1
lhs_mar = 0.06
chart_width = (1-x_mar-lhs_mar-x_gap)/2

left_chart_pos = c(x_mar,x_mar+chart_width,y_mar,y_mar+chart_height)
right_chart_pos = c(x_mar+chart_width+x_gap,1-lhs_mar,y_mar,y_mar+chart_height)

x_limits = c(0.2,0.5)

par(
	fig=left_chart_pos,
	mai=c(0,0,0,0),
	xpd=FALSE,
	mgp=c(0,0.25,0),
	tcl=0.25
)

plot(1:2,xlim=c(0.2,0.5),ylim=c(0,1),type="n",axes=FALSE)
axis(2,pos=0.2,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
mtext("a",side=2,las=1,line=2,at=1,font=2)
mtext("Infection attack rate",side=2,line=2)
axis(1,pos=0)
mtext("Theta",side=1,line=0.75)
points(params$Theta[sind:eind],1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$Theta[sind:eind],1-res_Q$S_365[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$Theta[sind:eind],1-res_QI$S_365[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$Theta[sind:eind],1-res_QIA$S_365[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$Theta[sind:eind],1-res_QIAT$S_365[1:max_i]/popsize,cex=ps,col="orange",pch=20)

par(xpd=NA)
legend(	0.3,y=-0.2,c("None","Quarantine","Quarantine and isolation","Quarantine, isolation and antivirals","Quarantine, isolation, anitvirals and contact tracing"),
		col=c("black","green","blue","cyan","orange"),
		pch=20,pt.cex=2,y.intersp=1.1,yjust=1,ncol=1,bty="n")
par(xpd=FALSE)

par(
	fig=right_chart_pos,
	new=TRUE
)

plot(1:2,xlim=c(0.5,0.75),ylim=c(0,1),type="n",axes=FALSE)
axis(2,pos=0.5,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
mtext("b",side=2,las=1,line=2,at=1,font=2)
mtext("Infection attack rate",side=2,line=2)
axis(1,pos=0)
mtext("Gamma",side=1,line=0.75)
points(params_gam$Gamma[sind:eind],1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params_gam$Gamma[sind:eind],1-res_Q$S_365[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params_gam$Gamma[sind:eind],1-res_QI$S_365[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params_gam$Gamma[sind:eind],1-res_QIA$S_365[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params_gam$Gamma[sind:eind],1-res_QIAT$S_365[1:max_i]/popsize,cex=ps,col="orange",pch=20)

