pdf("c:\\alpha.pdf", version = "1.4", width = 4, height = 4)
# postscript("c:\\alpha.ps",width = 4, height = 4)
par(mar = c(5, 4, 2, 2))
plot(rnorm(500), rnorm(500), col = rgb(0, 0, 1, 0.2), pch = 16)
dev.off()

#Examples which might not actually work because of missing data sructures, but HAVE worked
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