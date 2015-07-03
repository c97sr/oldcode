source("city_flu_figs_header.r")

root = "C:\\notefiles\\files\\projects\\flu2005\\results\\060630_sens\\"
stem = "set"
leaf = ".txt"
pfile = "C:\\notefiles\\files\\projects\\flu2005\\results\\060630_sens\\ParameterListMidGam.in"
params <- read.table(pfile)

sind <- 1
eind <- 200
max_i = eind-sind+1
ps = 0.5

res_None <- flu_processor_dash(root,stem,leaf,"None",sind,eind)
res_Q <- flu_processor_dash(root,stem,leaf,"Q",sind,eind)
res_QA <- flu_processor_dash(root,stem,leaf,"QA",sind,eind)
res_QI <- flu_processor_dash(root,stem,leaf,"QI",sind,eind)
res_QIA <- flu_processor_dash(root,stem,leaf,"QIA",sind,eind)
res_QIAC <- flu_processor_dash(root,stem,leaf,"QIAC",sind,eind)

chart_pos <- function(xindex,yindex) {

	x_left_mar = 0.1
	x_right_mar = 0.05
	x_gap = 0.15
	x_n_charts = 2
	x_width = (1-x_left_mar-x_right_mar-x_gap*(x_n_charts-1))/x_n_charts
	x_charts_l=rep(0,x_n_charts)
	for (i in 1:x_n_charts) x_charts_l[i]=x_left_mar+(i-1)*x_width+(i-1)*x_gap
	x_charts_r=rep(0,x_n_charts)
	for (i in 1:x_n_charts) x_charts_r[i]=x_charts_l[i]+x_width

	y_top_mar = 0.1
	y_bottom_mar = 0.1
	y_gap = 0.02
	y_n_charts = 5
	y_height = (1-y_top_mar-y_bottom_mar-y_gap*(y_n_charts-1))/y_n_charts
	y_charts_b = rep(0,y_n_charts)
	for (i in 1:y_n_charts) y_charts_b[i] = y_bottom_mar+(i-1)*y_height+(i-1)*y_gap
	y_charts_t = rep(0,y_n_charts)
	for (i in 1:y_n_charts) y_charts_t[i] = y_charts_b[i]+y_height

	c(x_charts_l[xindex],x_charts_r[xindex],y_charts_b[yindex],y_charts_t[yindex])

}

# Plot routines for figure 4
f3width=7
f3height=6*1.5
windows(width=f3width,height=f3height,rescale="R")
popsize=500000

par(xpd=FALSE)
par(cex=0.8)
par(mgp=c(3,0.25,0))
par(tcl=0.25)
par(mai=(c(0,0,0,0)))

x_limits = c(1,3)
r0_factor=1.83/1.8
labelline = 2.5

par(fig=sr_chart_pos(1,1,1,2))
plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("C",side=2,las=1,line=labelline,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_QA$S_365[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

par(fig=chart_pos(1,1),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("E",side=2,las=1,line=labelline,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
axis(1,pos=-1/10)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_QIAC$S_365[1:max_i]/popsize,cex=ps,col="orange",pch=20)

par(fig=chart_pos(1,2),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("D",side=2,las=1,line=labelline,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_QIA$S_365[1:max_i]/popsize,cex=ps,col="cyan",pch=20)

par(fig=chart_pos(1,4),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("B",side=2,las=1,line=labelline,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_QI$S_365[1:max_i]/popsize,cex=ps,col="blue",pch=20)

par(fig=chart_pos(1,5),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,1),type="n",axes=FALSE,ylab="Infection\nattack rate",xlab="R0")
mtext("A",side=2,las=1,line=labelline,at=1,font=2)
axis(2,pos=1,at=(0:4)/4, labels=c("00.0","0.25","0.50","0.75","1.00"),las=1)
points(params$R0[sind:eind]/r0_factor,1-res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,1-res_Q$S_365[1:max_i]/popsize,cex=ps,col="green",pch=20)

par(xpd=NA)
legend(	2,y=1.2,c("None","Q","QI","QA","QIA","QIAC"),
		col=c("black","green","blue","magenta","cyan","orange"),
		pch=20,pt.cex=2,y.intersp=1.1,horiz=TRUE,yjust=0,bty="n")
par(xpd=FALSE)


ymax = 12
y_labels = c("0","","4","","8","","12")
no_y_labels = length(y_labels)
par(fig=chart_pos(2,1),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,ymax),type="n",axes=FALSE,ylab="Total daily\ndoses per person",xlab="R0")
mtext("J",side=2,las=1,line=labelline,at=ymax,font=2)
axis(1,pos=-ymax/10)
axis(2,pos=1,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax, labels=y_labels,las=1)
points(params$R0[sind:eind]/r0_factor,res_Q$Sum_A[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QI$Sum_A[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIA$Sum_A[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIAC$Sum_A[1:max_i]/popsize,cex=ps,col="orange",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QA$Sum_A[1:max_i]/popsize,cex=ps,col="magenta",pch=20)
points(params$R0[sind:eind]/r0_factor,res_None$Sum_A[1:max_i]/popsize,cex=ps,col="black",pch=20)

ymax = 0.012
y_labels = c("0.0%","","0.4%","","0.8%","","1.2%")
no_y_labels = length(y_labels)
par(fig=chart_pos(2,2),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,ymax),type="n",axes=FALSE,ylab="Maximum\nincidence isolation")
mtext("I",side=2,las=1,line=labelline,at=ymax,font=2)
axis(2,pos=1,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax, labels=y_labels,las=1)
points(params$R0[sind:eind]/r0_factor,res_None$Max_I_Inc[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QI$Max_I_Inc[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIA$Max_I_Inc[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QA$Max_I_Inc[1:max_i]/popsize,cex=ps,col="magenta",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIAC$Max_I_Inc[1:max_i]/popsize,cex=ps,col="orange",pch=20)

ymax = 0.08
y_labels = c("0%","","2%","","4%","","6%","","8%")
no_y_labels = length(y_labels)
par(fig=chart_pos(2,3),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,ymax),type="n",axes=FALSE,ylab="Maximum\nprevalence isolation")
mtext("H",side=2,las=1,line=labelline,at=ymax,font=2)
axis(2,pos=1,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax, labels=y_labels,las=1)
points(params$R0[sind:eind]/r0_factor,res_QI$Max_I[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIA$Max_I[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIAC$Max_I[1:max_i]/popsize,cex=ps,col="orange",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QA$Max_I[1:max_i]/popsize,cex=ps,col="magenta",pch=20)
points(params$R0[sind:eind]/r0_factor,res_None$Max_I[1:max_i]/popsize,cex=ps,col="black",pch=20)

ymax = 0.1
y_labels = c("0.0%","2.5%","5.0%","7.5%","10.0%")
no_y_labels = length(y_labels)
par(fig=chart_pos(2,4),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,ymax),type="n",axes=FALSE,ylab="Maximum\nincidence quarantine")
mtext("G",side=2,las=1,line=labelline,at=ymax,font=2)
axis(2,pos=1,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax, labels=y_labels,las=1)
points(params$R0[sind:eind]/r0_factor,res_Q$Max_Q_Inc[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QI$Max_Q_Inc[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QA$Max_Q_Inc[1:max_i]/popsize,cex=ps,col="magenta",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIA$Max_Q_Inc[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIAC$Max_Q_Inc[1:max_i]/popsize,cex=ps,col="orange",pch=20)

ymax = 0.6
y_labels = c("0%","","20%","","40%","","60%")
no_y_labels = length(y_labels)
par(fig=chart_pos(2,5),new=TRUE)
plot(1:2,xlim=x_limits,ylim=c(0,ymax),type="n",axes=FALSE,ylab="Maximum\nprevalence quarantine")
mtext("F",side=2,las=1,line=labelline,at=ymax,font=2)
axis(2,pos=1,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax, labels=y_labels,las=1)
points(params$R0[sind:eind]/r0_factor,res_Q$Max_Q[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QI$Max_Q[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIA$Max_Q[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QA$Max_Q[1:max_i]/popsize,cex=ps,col="magenta",pch=20)
points(params$R0[sind:eind]/r0_factor,res_QIAC$Max_Q[1:max_i]/popsize,cex=ps,col="orange",pch=20)
