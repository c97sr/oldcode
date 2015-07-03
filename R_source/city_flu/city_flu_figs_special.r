source("city_flu_figs_header.r")

# Figure 7
dcol=7
scol=3.5

list_ints = c("None","Q","QI","QIA","QIAT","QIATC")
list_cols = c("",c(length(list_ints)))
for (i in 1:length(list_ints)) list_cols[i] = c_ass(list_ints[i])
no_ints = length(list_ints)

filestem = "D:\\share\\files\\projects\\flu2005\\results\\"

f3width=scol
f3height=2*scol
windows(width=f3width,height=f3height)
par(mai=(c(0.6,0.7,0,0)))
par(mgp=c(2,0.15,0))
par(tcl=0.25)
par(xpd=TRUE)

## set up plot
ymax = 0.04
y_labels = c("0%","","1%","","2%","","3%","","4%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0.67,0.99))
plot(c(0,1),type="n",col="black",xlim=c(0,250),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="",ylab="Daily incidence\nof infection")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("a",side=2,las=1,line=3,at=ymax,font=2)
popsize = 1000000


plot_int(filestem,"051109_base\\BaseCase-11-7\\DailyIncidence","None","black",FALSE,popsize)
plot_int(filestem,"051109_base\\BaseCase-11-7\\DailyIncidence","Q","green",FALSE,popsize)
plot_int(filestem,"051109_comp\\pComp25\\DailyIncidence","Q","red",FALSE,popsize)
plot_int(filestem,"051109_comp\\pComp75\\DailyIncidence","Q","blue",FALSE,popsize)
plot_int(filestem,"051109_base\\BaseCase-11-7\\DailyIncidence","QI","grey",FALSE,popsize)
plot_int(filestem,"051109_comp\\pComp25\\DailyIncidence","QI","cyan",FALSE,popsize)
plot_int(filestem,"051109_comp\\pComp75\\DailyIncidence","QI","magenta",FALSE,popsize)
legend(60,y=ymax,
	c("None","Q p=0.25","Q p=0.50","Q p=0.75","QI p=0.25","QI p=0.50","QI p=0.75"),
	col=c("black","red","green","blue","cyan","grey","magenta"),
	lwd=2,y.intersp=1.0,horiz=FALSE,bty="n",cex=0.75)

list_ints = c("None","Q","QI","QIA","QIAC")
base_table = table_int(paste(filestem,"051109_comp\\pComp25\\",sep=""),list_ints)
write.table(base_table,file=paste(filestem,"interventions_table_25.out",sep=""),sep="\t",row.names=TRUE)
base_table = table_int(paste(filestem,"051109_comp\\pComp75\\",sep=""),list_ints)
write.table(base_table,file=paste(filestem,"interventions_table_75.out",sep=""),sep="\t",row.names=TRUE)

ymax = 0.015
y_labels = c("0%","","0.5%","","1.0%","","1.5%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0.34,0.66),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,250),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="",ylab="Daily incidence\nof infection")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("b",side=2,las=1,line=3,at=ymax,font=2)

plot_int(filestem,"051109_base\\BaseCase-11-7\\DailyIncidence","QIA","red",FALSE,popsize)
plot_int(filestem,"051109_test\\sens60_CapInf\\DailyIncidence","QIAT","green",FALSE,popsize)
plot_int(filestem,"051109_test\\sens90_CapInf\\DailyIncidence","QIAT","blue",FALSE,popsize)
legend(60,y=ymax,c("Constrained 60%","Unconstrained 60%","Unconstrained 90%"),
	col=c("red","green","blue"),lwd=2,y.intersp=1.1,horiz=FALSE,bty="n",cex=0.75)

ymax = 0.12
y_labels = c("0%","","4%","","8%","","12%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0.00,0.33),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,250),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="Day",ylab="Prevalence of\nquarantine")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("c",side=2,las=1,line=3,at=ymax,font=2)

plot_int(filestem,"051109_base\\BaseCase-11-7\\numIndividualQuarantined","QIA","red",FALSE,popsize)
plot_int(filestem,"051109_test\\sens60_CapInf\\numIndividualQuarantined","QIAT","green",FALSE,popsize)
plot_int(filestem,"051109_test\\sens90_CapInf\\numIndividualQuarantined","QIAT","blue",FALSE,popsize)


