
root = "D:\\share\\files\\projects\\flu2005\\results\\051016b\\"
stem = "set"
leaf = ".txt"
pfile = "D:\\share\\files\\projects\\flu2005\\results\\051016b\\parameterList.in"
params <- read.table(pfile)

sind <- 1
eind <- 9
max_i = eind-sind+1

res_None <- flu_processor(root,stem,leaf,"None",sind,eind)
res_Q <- flu_processor(root,stem,leaf,"Q",sind,eind)
res_QI <- flu_processor(root,stem,leaf,"QI",sind,eind)
res_QIA <- flu_processor(root,stem,leaf,"QIA",sind,eind)
res_QIAT <- flu_processor(root,stem,leaf,"QIAT",sind,eind)
res_QIATC <- flu_processor(root,stem,leaf,"QIATC",sind,eind)

dcol=7
scol=3.5

# Plot routines for figure 1
windows(width=dcol,height=dcol*2/3)

par(cex=0.8)
par(mai=c(0.75,0.75,0.1,0.1))
ps = 1
popsize=100000

par(fig=c(0,0.5,0,1))
plot(1:2,xlim=c(1,4),ylim=c(0,1),type="n",xlab="R0",ylab="s_365")
points(params$R0[1:max_i],res_None$S_365[1:max_i]/popsize,cex=ps,col="black",pch=20)
points(params$R0[1:max_i],res_Q$S_365[1:max_i]/popsize,cex=ps,col="green",pch=20)
legend(3,y=1,c("None","Q"),col=c("black","green"),pch=20,pt.cex=2,y.intersp=1.1)

par(fig=c(0.5,1,0.5,1),new=TRUE)
plot(1:2,xlim=c(1,4),ylim=c(0,0.3),type="n",xlab="R0",ylab="max(q)")
points(params$R0[1:max_i],res_Q$Max_Q[1:max_i]/popsize,cex=ps,col="black",pch=20)

par(fig=c(0.5,1,0,0.5),new=TRUE)
plot(1:2,xlim=c(1,4),ylim=c(0,0.3),type="n",xlab="R0",ylab="max(i)")
points(params$R0[1:max_i],res_Q$Max_I[1:max_i]/popsize,cex=ps,col="black",pch=20)

# Plot routines for figure 2
windows(width=dcol,height=dcol)

par(cex=0.8)
par(mai=c(0.75,0.75,0.1,0.1))
ps = 1
popsize=100000

plot(1:2,xlim=c(1,4),ylim=c(0,1),type="n",xlab="R0",ylab="delta(s_365) vs None")
points(params$R0[1:max_i],(res_Q$S_365[1:max_i]-res_None$S_365[1:max_i])/popsize,cex=ps,col="green",pch=20)
points(params$R0[1:max_i],(res_QI$S_365[1:max_i]-res_None$S_365[1:max_i])/popsize,cex=ps,col="blue",pch=20)
points(params$R0[1:max_i],(res_QIA$S_365[1:max_i]-res_None$S_365[1:max_i])/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[1:max_i],(res_QIAT$S_365[1:max_i]-res_None$S_365[1:max_i])/popsize,cex=ps,col="red",pch=20)
points(params$R0[1:max_i],(res_QIATC$S_365[1:max_i]-res_None$S_365[1:max_i])/popsize,cex=ps,col="magenta",pch=20)
for (i in 1:max_i) 
	lines( c(params$R0[i],params$R0[i]), c((res_Q$S_365[i]-res_None$S_365[i])/popsize,(res_QIATC$S_365[i]-res_None$S_365[i])/popsize),lty=1)
legend(1,y=1,c("Q","QA","QT","QTA","QTAC"),col=c("green","blue","cyan","red","magenta"),pch=20,pt.cex=2,y.intersp=1.1,ncol=1)
text(2.96,0.1,"A",pos=1)
text(3.04,0.14,"B",pos=1)

# Plot routines for figure 3
f3width=scol
f3height=scol*3
windows(width=f3width,height=f3height)

par(cex=0.8)

b_mar_in = 0.7
l_mar_in = 0.7
delta=0.015

b_mar_prop = b_mar_in/f3height

fw = c(1,1,1,1,1)
fw = fw/sum(fw)*(1-b_mar_prop)

yps = array(0,c(5))
yps[1] = fw[1]+b_mar_prop
for (i in 2:5) yps[i] = yps[i-1]+fw[i]  

par(fig=c(0,0.95,0,yps[1]-delta))
par(mai=(c(b_mar_in,l_mar_in,0,0)))
plot(1:2,xlim=c(1,4),ylim=c(0,0.1),type="n",axes=FALSE, xlab="R0",ylab="max(CT)")
axis(1,pos=-0.01)
axis(2,pos=1,at=(0:4)/4*0.1, labels=c("0.00","","0.05","","0.10"),las=1)
points(params$R0[1:max_i],res_QIATC$Max_C[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

par(fig=c(0,0.95,yps[1]+delta,yps[2]-delta),new=TRUE)
par(mai=(c(0,l_mar_in,0,0)))
plot(1:2,xlim=c(1,4),ylim=c(0,0.4),type="n",axes=FALSE,ylab="max(Testing)")
axis(2,pos=1,at=(0:4)/4*0.4, labels=c("0.0","0.1","0.2","0.3","0.4"),las=1)
points(params$R0[1:max_i],res_QIA$Max_T[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[1:max_i],res_QIAT$Max_T[1:max_i]/popsize,cex=ps,col="red",pch=20)
points(params$R0[1:max_i],res_QIATC$Max_T[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

par(fig=c(0,0.95,yps[2]+delta,yps[3]-delta),new=TRUE)
par(mai=(c(0,l_mar_in,0,0)))
plot(1:2,xlim=c(1,4),ylim=c(0,15),type="n",axes=FALSE,ylab="total(Treatments)")
axis(2,pos=1,at=(0:6)*10/4, labels=c("0","","5","","10","","15"),las=1)
points(params$R0[1:max_i],res_QI$Sum_A[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[1:max_i],res_QIAT$Sum_A[1:max_i]/popsize,cex=ps,col="red",pch=20)
points(params$R0[1:max_i],res_QIATC$Sum_A[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

par(fig=c(0,0.95,yps[3]+delta,yps[4]-delta),new=TRUE)
par(mai=(c(0,l_mar_in,0,0)))
plot(1:2,xlim=c(1,4),ylim=c(0,0.6),type="n",axes=FALSE,ylab="max(Quarantine)")
axis(2,pos=1,at=(0:3)/3*0.6, labels=c("0.0","0.2","0.4","0.6"),las=1)
points(params$R0[1:max_i],res_Q$Max_Q[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$R0[1:max_i],res_QI$Max_Q[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[1:max_i],res_QIA$Max_Q[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[1:max_i],res_QIAT$Max_Q[1:max_i]/popsize,cex=ps,col="red",pch=20)
points(params$R0[1:max_i],res_QIATC$Max_Q[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

par(fig=c(0,0.95,yps[4]+delta,yps[5]-delta),new=TRUE)
par(mai=(c(0,l_mar_in,0,0)))
plot(1:2,xlim=c(1,4),ylim=c(0,0.2),type="n",axes=FALSE,ylab="max(Isolation)")
axis(2,pos=1,at=(0:4)/5*0.2, labels=c("0.0","","0.1","","0.2"),las=1)
points(params$R0[1:max_i],res_Q$Max_I[1:max_i]/popsize,cex=ps,col="green",pch=20)
points(params$R0[1:max_i],res_QI$Max_I[1:max_i]/popsize,cex=ps,col="blue",pch=20)
points(params$R0[1:max_i],res_QIA$Max_I[1:max_i]/popsize,cex=ps,col="cyan",pch=20)
points(params$R0[1:max_i],res_QIAT$Max_I[1:max_i]/popsize,cex=ps,col="red",pch=20)
points(params$R0[1:max_i],res_QIATC$Max_I[1:max_i]/popsize,cex=ps,col="magenta",pch=20)

## Fig 5
filestem = "D:\\share\\files\\projects\\flu2005\\results\\050913\\"
no_int <- read.table(file=paste(filestem,"DailyIncidence_None.txt",sep=""))
qta_uncon <- read.table(file=paste(filestem,"DailyIncidence_QTA_Inf.txt",sep=""))
qta_100_T <- read.table(file=paste(filestem,"DailyIncidence_QTA_1.txt",sep=""))
qta_tests <- read.table(file=paste(filestem,"numActiveTestServers_QTA_Inf.txt",sep=""))
qta_100_T_tests <- read.table(file=paste(filestem,"numActiveTestServers_QTA_100.txt",sep=""))
no_days = dim(no_int)[1]
no_reals = dim(no_int)[2]
baseline = array(0,c(no_days))
baseline_bounds = array(0,c(2,no_days))
constrained = array(0,c(no_days))
constrained_bounds = array(0,c(2,no_days))
unconstrained = array(0,c(no_days))
unconstrained_bounds = array(0,c(2,no_days))
constrained_tests = array(0,c(no_days))
unconstrained_tests = array(0,c(no_days))
for (i in 1:no_days) {
	baseline[i]=mean(as.numeric(no_int[i,]))
	quant = quantile(as.numeric(no_int[i,]),probs=c(0.025,0.975))
	for (j in 1:2) baseline_bounds[j,i] <- quant[j]
	unconstrained[i]=mean(as.numeric(qta_uncon[i,]))
	quant = quantile(as.numeric(qta_100_T[i,]),probs=c(0.025,0.975))
	for (j in 1:2) unconstrained_bounds[j,i] <- quant[j]
	constrained[i]=mean(as.numeric(qta_100_T[i,]))
	quant = quantile(as.numeric(qta_uncon[i,]),probs=c(0.025,0.975))
	for (j in 1:2) constrained_bounds[j,i] <- quant[j]
	constrained_tests[i]=mean(as.numeric(qta_100_T_tests[i,]))
	unconstrained_tests[i]=mean(as.numeric(qta_tests[i,]))
}

f3width=scol
f3height=2*scol
windows(width=f3width,height=f3height)
par(fig=c(0.1,0.95,0.6,0.99))
par(mai=(c(0.5,0.5,0,0)))
par(mgp=c(1.5,0.15,0))
par(tcl=0.25)

plot(baseline,type="n",col="black",xlim=c(0,200),ylim=c(0,400),lwd=2,axes=FALSE,xlab="Day",ylab="Incidence")
axis(1,pos=-10)
axis(2,pos=0,at=0:4/4*400,labels=c("0","100","200","300","400"),las=1)
mtext("A",side=2,las=1,line=3,at=400,font=2)
## for (i in 1:10) points(no_int[,i],type="l",col="grey",lwd=0.5)
for (i in 1:2) points(baseline_bounds[i,],col="black",lwd=1,type="l")
for (i in 1:2) points(unconstrained_bounds[i,],col="cyan",lwd=1,type="l")
for (i in 1:2) points(constrained_bounds[i,],col="red",lwd=1,type="l")
points(baseline,type="l",col="black",lwd=3)
points(unconstrained,type="l",col="cyan",lwd=3)
points(constrained,type="l",col="red",lwd=3)

par(xpd=NA)

legend(0,y=-120,c("Baseline","Unconstrained","Constrained"),col=c("black","cyan","red"),lwd=2,ncol=1,cex=0.8)

par(xpd=FALSE)

par(fig=c(0.1,0.95,0.1,0.49),new=TRUE)

plot(baseline/10000,type="n",col="black",xlim=c(0,200),ylim=c(0,400),lwd=2,axes=FALSE,xlab="Day",ylab="Incidence")
axis(1,pos=0)
axis(2,pos=0,at=0:4/4*400,labels=c("0","100","200","300","400"),las=1)
mtext("B",side=2,las=1,line=3,at=400,font=2)
for (i in 1:5) points(qta_tests[,i],col="cyan",lwd=1,type="l")
for (i in 1:5) points(qta_100_T_tests[,i],col="red",lwd=1,type="l")


## Fig 6
dcol=7
scol=3.5

filestem = "D:\\share\\files\\projects\\flu2005\\results\\051014b\\"

f3width=scol
f3height=2*scol
windows(width=f3width,height=f3height)
par(mai=(c(0.5,0.5,0,0)))
par(mgp=c(1,0.15,0))
par(tcl=0.25)
par(xpd=TRUE)

## set up plot
ymax = 0.04
no_y_labels = 5
y_labels = c("0","1","2","3","4")
par(fig=c(0.1,0.95,0.66,0.99))
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="Day",ylab="Incidence (%)")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("A",side=2,las=1,line=3,at=ymax,font=2)

list_ints = c("None","Q","QI","QIA","QIAT","QIATC")
list_cols = c("",c(length(list_ints)))
for (i in 1:length(list_ints)) list_cols[i] = c_ass(list_ints[i])
no_ints = length(list_ints)
legend(100,y=ymax,list_ints,col=list_cols,pch=20,pt.cex=2,y.intersp=1.1,horiz=FALSE)


plot_int(filestem,"DailyIncidence",list_ints[1],list_cols[1],TRUE)
for (i in 2:no_ints) plot_int(filestem,"DailyIncidence",list_ints[i],list_cols[1],FALSE)

base_table = table_int(filestem,list_ints)
write.table(base_table,file=paste(filestem,"interventions_table.out",sep=""),sep="\t",row.names=TRUE)

ymax = 0.30
no_y_labels = 4
y_labels = c("0","10","20","30")
par(fig=c(0.1,0.95,0.33,0.65),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="Day",ylab="Prevalence (%)")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("B",side=2,las=1,line=3,at=ymax,font=2)

list_intsB = c("Q","QI","QIA","QIAT","QIATC")
no_intsB = length(list_intsB)

for (i in 1:no_intsB) plot_int(filestem,"numIndividualQuarantined",list_intsB[i],list_cols[1],FALSE)

ymax = 0.1
no_y_labels = 5
y_labels = c("0.0","2.5","5.0","7.5","10.0")
par(fig=c(0.1,0.95,0,0.32),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="Day",ylab="Prevalence (%)")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("C",side=2,las=1,line=3,at=ymax,font=2)

list_intsC = c("None","Q","QI","QIA","QIAT","QIATC")
no_intsC = length(list_intsC)

for (i in 1:no_intsC) plot_int(filestem,"numIndividualIsolated",list_intsC[i],list_cols[1],FALSE)

