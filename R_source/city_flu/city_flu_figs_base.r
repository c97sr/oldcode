source("city_flu_figs_header.r")

## Fig 6
dcol=7
scol=3.5

popsize=1
linethickness = 2

#filestem = "D:\\share\\files\\projects\\flu2005\\results\\060525_base\\"
#filestem = "D:\\share\\files\\projects\\flu2005\\results\\060526_delay_start\\"
filestem = "C:\\notefiles\\files\\projects\\flu2005\\results\\060630_base\\"

f3width=scol
f3height=2*scol
windows(width=f3width,height=f3height)
par(mai=(c(0.5,0.5,0,0)))
par(mgp=c(1.5,0.15,0))
par(tcl=0.25)
par(xpd=TRUE)

## set up plot
ymax = 0.04
y_labels = c("0%","1%","2%","3%","4%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0.66,0.99))
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="",ylab="")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("",side=2,las=1,line=3,at=ymax,font=2)

list_ints = c("None","Q","QI","QA","QIA","QIAC")
list_cols = c("",c(length(list_ints)))
for (i in 1:length(list_ints)) list_cols[i] = c_ass(list_ints[i])
no_ints = length(list_ints)
legend(100,y=ymax,list_ints,col=list_cols,lwd=2,pt.cex=2,y.intersp=1.1,horiz=FALSE,bty="n")


plot_int(filestem,"DailyIncidence",list_ints[1],list_cols[1],FALSE,popsize,thickness=linethickness)
for (i in 2:no_ints) plot_int(filestem,"DailyIncidence",list_ints[i],list_cols[i],FALSE,popsize,thickness=linethickness)

base_table = table_int(filestem,list_ints)
write.table(base_table,file=paste(filestem,"interventions_table.out",sep=""),sep="\t",row.names=TRUE)

ymax = 0.3
y_labels = c("0%","","10%","","20%","","30%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0.33,0.65),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="",ylab="")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("",side=2,las=1,line=3,at=ymax,font=2)

list_intsB = c("Q","QI","QA","QIA","QIAC")
no_intsB = length(list_intsB)

for (i in 1:no_intsB) plot_int(filestem,"numIndividualQuarantined",list_intsB[i],list_cols[i+1],FALSE,popsize,thickness=linethickness)

ymax = 0.02
y_labels = c("0.0%","0.5%","1.0%","1.5%","2.0%")
no_y_labels = length(y_labels)
par(fig=c(0.1,0.95,0,0.32),new=TRUE)
plot(c(0,1),type="n",col="black",xlim=c(0,150),ylim=c(0,ymax),lwd=2,axes=FALSE,xlab="",ylab="")
axis(1,pos=-0.001)
axis(2,pos=0,at=(0:(no_y_labels-1))/(no_y_labels-1)*ymax,labels=y_labels,las=1)
mtext("",side=2,las=1,line=3,at=ymax,font=2)

list_intsC = c("None","Q","QI","QA","QIA","QIAC")
no_intsC = length(list_intsC)

for (i in 1:no_intsC) plot_int(filestem,"numIndividualIsolated",list_intsC[i],list_cols[i],FALSE,popsize,thickness=linethickness)

