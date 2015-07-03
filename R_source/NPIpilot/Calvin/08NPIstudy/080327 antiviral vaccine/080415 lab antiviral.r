#080414 add clinical and lab confirmed 2nd case (not finished)
#######
#check the viral load of index
pcr <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_pcr0708_o.csv", header=TRUE) #this year
pcr <- pcr[pcr$member ==0,] #index only

#function for plotting graphs
plot_1g <-  function(plot1,plot2,plot3,name1,name2,name3,title){

par(mar=c(3,5,3,1))
plot( plot1$visitday, plot1$vload, axes=FALSE, type="b", main = title,
xlim=c(0,7), ylim=c(0,8), xlab=" ", ylab="log viral load ", cex.lab=1.2)
axis(1, pos=0, at=0:7, labels=rep(NA, 7), cex.axis=1.2)
axis(2, pos=0, at=c(0,2,4,6,8), labels=c("0","2","4","6","8"), las=1, cex.axis=1.2)

lines(plot2$visitday, plot2$vload, type= "b", col = "red")
lines( plot3$visitday, plot3$vload, type="b",col = "blue")

legend(4, 6, cex=1,
  legend=c(name1,name2,name3), col=c("black","red","blue"),  lty=c(1,1,1)
  #pch = c(17,16,15) ,)
)
#x axis
mtext("Day  0", side=1, at=-0.3, line=1, cex=1.2)
mtext("1", side=1, at=1, line=1, cex=1.2)
mtext("2", side=1, at=2, line=1, cex=1.2)
mtext("3", side=1, at=3, line=1, cex=1.2)
mtext("4", side=1, at=4, line=1, cex=1.2)
mtext("5", side=1, at=5, line=1, cex=1.2)
mtext("6", side=1, at=6, line=1, cex=1.2)
mtext("7", side=1, at=7, line=1, cex=1.2)


}
############

temp <- pcr[pcr$member ==0 & pcr$visit ==1,] 
flutype <- temp[c(1,5)]	#for getting the influenza type
visit1load <- pcr[pcr$member ==0 & pcr$visit ==1,]
visit1load <- visit1load[c(1,6)]
visit2load <- pcr[pcr$member ==0 & pcr$visit ==2,]
visit2load <- visit2load[c(1,6)]
visit3load <- pcr[pcr$member ==0 & pcr$visit ==3,]
visit3load <- visit3load[c(1,6)]

vload123 <-  merge(visit1load, visit2load,  by= "hhID",  all.x= TRUE) #for viral load in 3 visits
names(vload123) <- c("hhID","v1qPCR","v2qPCR")
vload123 <-  merge(vload123, visit3load,  by= "hhID",  all.x= TRUE)
names(vload123) <- c("hhID","v1qPCR","v2qPCR","v3qPCR")

pcr$qPCR <- as.numeric(as.character(pcr$qPCR))
#
hh <- as.character(unique(pcr$hhID))

dat <- data.frame(hhID =0, meanload =0)
for (i in 1:length(hh)){ 
  
  dat[i,1] <- hh[i]
  dat[i,2] <- mean(na.exclude(pcr$qPCR[as.character(pcr$hhID) == hh[i] ]))
}



####################
#for antivirals
av <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_av0708_m.csv", header=TRUE)
av <- av[av$member ==0,]
av <- av[c(1,4)]

dat <- merge(dat, av,  by= "hhID",  all.x= TRUE)
dat$av <- as.character(dat$av)
dat$av[is.na(dat$av)] <- 0
dat$av[dat$av == "amantadine"] <- 1
dat$av[dat$av == "tamiflu"] <- 2
dat$av[dat$av == "relenza"] <- 3
dat$av[dat$av == "ribavirin"] <- 3
dat <- dat[dat$av !=3,]

dat <- merge(dat, flutype,  by= "hhID",  all.x= TRUE)
####################
# check mean viral load with different antiviral

mean(na.exclude(dat$meanload[dat$av == 0 & dat$flu.type=="A"]))
mean(na.exclude(dat$meanload[dat$av == 1& dat$flu.type=="A"]))
mean(na.exclude(dat$meanload[dat$av == 2& dat$flu.type=="A"]))

mean(na.exclude(dat$meanload[dat$av == 0 & dat$flu.type=="B"]))
mean(na.exclude(dat$meanload[dat$av == 1& dat$flu.type=="B"]))
mean(na.exclude(dat$meanload[dat$av == 2& dat$flu.type=="B"]))
#####################

#check 3 visits' averge viral load 
dat <- merge(dat, vload123,  by= "hhID",  all.x= TRUE)
dat$v1qPCR <- log(as.numeric(as.character(dat$v1qPCR)),10)
dat$v2qPCR <- log(as.numeric(as.character(dat$v2qPCR)),10)
dat$v3qPCR <- log(as.numeric(as.character(dat$v3qPCR)),10)

dat$v1qPCR[dat$v1qPCR == -Inf] <-0
dat$v2qPCR[dat$v2qPCR == -Inf] <-0
dat$v3qPCR[dat$v3qPCR == -Inf] <-0


#plots.none <- data.frame( av=0, visit=1:3, vload =0 ) 
#plots.aman <- data.frame( av=1, visit=1:3, vload =0 ) 
#plots.tami <- data.frame( av=2, visit=1:3, vload =0 ) 

#plots.none[1,3] <- mean(na.exclude(dat$v1qPCR[dat$av == 0 & dat$flu.type=="A"]))
#plots.none[2,3] <- mean(na.exclude(dat$v2qPCR[dat$av == 0 & dat$flu.type=="A"]))
#plots.none[3,3] <- mean(na.exclude(dat$v3qPCR[dat$av == 0 & dat$flu.type=="A"]))

#plots.aman[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 1 & dat$flu.type=="A"]))
#plots.aman[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 1 & dat$flu.type=="A"]))
#plots.aman[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 1 & dat$flu.type=="A"]))

#plots.tami[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 2 & dat$flu.type=="A"]))
#plots.tami[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 2 & dat$flu.type=="A"]))
#plots.tami[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 2 & dat$flu.type=="A"]))

#plot_1g(plots.none$vload,plots.aman$vload,plots.tami$vload,"None","Amantadine","Tamiflu", 
#"Mean viral load among different visits (flu A)")


#flu B
#plots.none[1,3] <- mean(na.exclude(dat$v1qPCR[dat$av == 0 & dat$flu.type=="B"]))
#plots.none[2,3] <- mean(na.exclude(dat$v2qPCR[dat$av == 0 & dat$flu.type=="B"]))
#plots.none[3,3] <- mean(na.exclude(dat$v3qPCR[dat$av == 0 & dat$flu.type=="B"]))

#plots.aman[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 1 & dat$flu.type=="B"]))
#plots.aman[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 1 & dat$flu.type=="B"]))
#plots.aman[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 1 & dat$flu.type=="B"]))

#plots.tami[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 2 & dat$flu.type=="B"]))
#plots.tami[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 2 & dat$flu.type=="B"]))
#plots.tami[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 2 & dat$flu.type=="B"]))

#plot_1g(plots.none$vload,plots.aman$vload,plots.tami$vload,"None","Amantadine","Tamiflu", 
#"Mean viral load among different visits (flu B)")

#add housechar data (visit day)
housechar <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_housechar0708_h.csv", header=TRUE)
visitday <- housechar[c(1,11:14)] 
visitday$v1_day <- visitday$v1_day - visitday$clinic_day
visitday$v2_day <- visitday$v2_day - visitday$clinic_day
visitday$v3_day <- visitday$v3_day - visitday$clinic_day

dat <- merge(dat, visitday,  by= "hhID",  all.x= TRUE)

#plot for visit day instead of which visit for influenza A


plots.none <- data.frame( av=rep(0,3), visitday=0, vload =0 ) 
plots.aman <- data.frame( av=rep(1,3), visitday=0, vload =0 ) 
plots.tami <- data.frame( av=rep(2,3), visitday=0, vload =0 ) 

plots.none[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 0 & dat$flu.type=="A"]))
plots.none[1,3] <- mean(na.exclude(dat$v1qPCR[dat$av == 0 & dat$flu.type=="A"]))
plots.none[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 0 & dat$flu.type=="A"]))
plots.none[2,3] <- mean(na.exclude(dat$v2qPCR[dat$av == 0 & dat$flu.type=="A"]))
plots.none[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 0 & dat$flu.type=="A"]))
plots.none[3,3] <- mean(na.exclude(dat$v3qPCR[dat$av == 0 & dat$flu.type=="A"]))

plots.aman[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 1 & dat$flu.type=="A"]))
plots.aman[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 1 & dat$flu.type=="A"]))
plots.aman[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 1 & dat$flu.type=="A"]))
plots.aman[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 1 & dat$flu.type=="A"]))
plots.aman[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 1 & dat$flu.type=="A"]))
plots.aman[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 1 & dat$flu.type=="A"]))

plots.tami[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 2 & dat$flu.type=="A"]))
plots.tami[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 2 & dat$flu.type=="A"]))
plots.tami[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 2 & dat$flu.type=="A"]))
plots.tami[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 2 & dat$flu.type=="A"]))
plots.tami[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 2 & dat$flu.type=="A"]))
plots.tami[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 2 & dat$flu.type=="A"]))

plot_1g(plots.none,plots.aman,plots.tami,"None","Amantadine","Tamiflu", 
"Mean viral load after antiviral prescription (flu A)")

#for influenza B

plots.none <- data.frame( av=rep(0,3), visitday=0, vload =0 ) 
plots.aman <- data.frame( av=rep(1,3), visitday=0, vload =0 ) 
plots.tami <- data.frame( av=rep(2,3), visitday=0, vload =0 ) 

plots.none[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 0 & dat$flu.type=="B"]))
plots.none[1,3] <- mean(na.exclude(dat$v1qPCR[dat$av == 0 & dat$flu.type=="B"]))
plots.none[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 0 & dat$flu.type=="B"]))
plots.none[2,3] <- mean(na.exclude(dat$v2qPCR[dat$av == 0 & dat$flu.type=="B"]))
plots.none[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 0 & dat$flu.type=="B"]))
plots.none[3,3] <- mean(na.exclude(dat$v3qPCR[dat$av == 0 & dat$flu.type=="B"]))

plots.aman[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 1 & dat$flu.type=="B"]))
plots.aman[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 1 & dat$flu.type=="B"]))
plots.aman[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 1 & dat$flu.type=="B"]))
plots.aman[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 1 & dat$flu.type=="B"]))
plots.aman[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 1 & dat$flu.type=="B"]))
plots.aman[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 1 & dat$flu.type=="B"]))

plots.tami[1,2] <- mean(na.exclude(dat$v1_day[dat$av == 2 & dat$flu.type=="B"]))
plots.tami[1,3] <-mean(na.exclude(dat$v1qPCR[dat$av == 2 & dat$flu.type=="B"]))
plots.tami[2,2] <- mean(na.exclude(dat$v2_day[dat$av == 2 & dat$flu.type=="B"]))
plots.tami[2,3] <-mean(na.exclude(dat$v2qPCR[dat$av == 2 & dat$flu.type=="B"]))
plots.tami[3,2] <- mean(na.exclude(dat$v3_day[dat$av == 2 & dat$flu.type=="B"]))
plots.tami[3,3] <-mean(na.exclude(dat$v3qPCR[dat$av == 2 & dat$flu.type=="B"]))

plot_1g(plots.none,plots.aman,plots.tami,"None","Amantadine","Tamiflu", 
"Mean viral load after antiviral prescription (flu B)")



#################
clin2nd <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_clin2nd0708_m.csv", header=TRUE)

clin2nd <- as.data.frame(table(clin2nd$hhID))
names(clin2nd) <-c("hhID","clin2nd")

dat <- merge( dat, clin2nd, by= "hhID" ,  all.x=TRUE) #onsettime, if non-index, onsettime is follow index
dat$clin2nd[is.na(dat$clin2nd)] <- 0

############ labresult
lab2nd <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_lab2nd0708_m.csv", header=TRUE)

lab2nd <- as.data.frame(table(lab2nd$hhID))
names(lab2nd) <-c("hhID","lab2nd")
dat <- merge(  dat, lab2nd, by= "hhID" ,  all.x=TRUE) #onsettime, if non-index, onsettime is follow index

#######################################
#for those have done lab results (as this year lab result is lagging behind
 lab2nd08 <- read.csv(paste(input_dir08,"lab2nd2008.csv", sep=""), header=TRUE) #this year
 donelab <-   lab2nd08[ !(is.na(lab2nd08$V1) & is.na(lab2nd08$V2) & is.na(lab2nd08$V3)),]
 doneno <- c(unique(donelab$hhID))
 doneno <- paste("h080",doneno,sep="")
 doneno <- as.data.frame(doneno)
 names(doneno) <- "hhID"
 dat$lab2nd[is.na(dat$lab2nd) & (dat$hhID %in% doneno$hhID | substr(dat$hhID,1,3) == "h07" )] <- 0
#

#check the viral load of index
pcr <- read.csv("P:/GROUP/NPIstudy/data/combined/2008_4_11_pcr0708_o.csv", header=TRUE) #this year
pcr <- pcr[pcr$member ==0,] 
