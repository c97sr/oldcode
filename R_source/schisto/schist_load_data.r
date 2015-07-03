if (SR_version=="windows") {
	wdir = "D:\\files\\projects\\sars\\amoy\\eclipse\\R Source"
	testfile = "D:\\files\\projects\\schisto\\data\\fromPatrick\\060719\\"
	datDemogFile = "D:\\files\\projects\\schisto\\data\\actual_demog_SR_not_220.txt"
} else if (SR_version=="linux") {
	wdir = "/root/workspace/R Source/"
	testfile = "/mnt/windesk/projects/schisto/data/fromPatrick/060719/"
	datDemogFile = "/mnt/windesk/projects/schisto/data/actual_demog_SR_not_220.txt"
}

datDemog = read.table(datDemogFile,header=TRUE)
global_villpv = rep(0.01,dim(datDemog)[1])

datHumanDemog = read.table(paste(testfile,"HumanDemographics_excl_220.txt",sep=""),header=TRUE,sep=",")
datHumanInf = read.table(paste(testfile,"null4Dynamic-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datHumanExcerpt = as.data.frame(cbind(datHumanDemog[,"Brgy"],datHumanDemog[,"pattern"]))
names(datHumanExcerpt)<-c("Village","Pattern")
pattHuman <- table(datHumanExcerpt)

datCatDemog = read.table(paste(testfile,"CatsDemographics.txt",sep=""),header=TRUE,sep=",")
datCatInf = read.table(paste(testfile,"null4Dynamic-Cats-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datCatExcerpt = as.data.frame(cbind(datCatDemog[,"brgy"],datCatDemog[,"pattern"]))
names(datCatExcerpt)<-c("Village","Pattern")
pattCat <- table(datCatExcerpt)

datDogDemog = read.table(paste(testfile,"DogsDemographics.txt",sep=""),header=TRUE,sep=",")
datDogInf = read.table(paste(testfile,"null4Dynamic-Dogs-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datDogExcerpt = as.data.frame(cbind(datDogDemog[,"brgy"],datDogDemog[,"pattern"]))
names(datDogExcerpt)<-c("Village","Pattern")
pattDog <- table(datDogExcerpt)

datPigDemog = read.table(paste(testfile,"PigsDemographics.txt",sep=""),header=TRUE,sep=",")
datPigInf = read.table(paste(testfile,"null4Dynamic-Pigs-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datPigExcerpt = as.data.frame(cbind(datPigDemog[,"brgy"],datPigDemog[,"pattern"]))
names(datPigExcerpt)<-c("Village","Pattern")
pattPig <- table(datPigExcerpt)

datWtrbuffDemog = read.table(paste(testfile,"WtrbuffsDemographics.txt",sep=""),header=TRUE,sep=",")
datWtrbuffInf = read.table(paste(testfile,"null4Dynamic-Wtrbuffs-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datWtrbuffExcerpt = as.data.frame(cbind(datWtrbuffDemog[,"brgy"],datWtrbuffDemog[,"pattern"]))
names(datWtrbuffExcerpt)<-c("Village","Pattern")
pattWtrbuff <- table(datWtrbuffExcerpt)

datRatDemog = read.table(paste(testfile,"RatsDemographics.txt",sep=""),header=TRUE,sep=",")
datRatInf = read.table(paste(testfile,"null4Dynamic-Rats-postprobStats.txt",sep=""),header=TRUE,sep="\t")
datRatExcerpt = as.data.frame(cbind(datRatDemog[,"brgy"],datRatDemog[,"pattern"]))
names(datRatExcerpt)<-c("Village","Pattern")
pattRat <- table(datRatExcerpt)

list_patterns = list(pattHuman,pattCat,pattDog,pattPig,pattWtrbuff,pattRat)
list_probs = list(datHumanInf,datCatInf,datDogInf,datPigInf,datWtrbuffInf,datRatInf)

lookupPattern(111,3,list_probs[[1]])
listPTS <- list(c(1,2,3),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2))
dataSum <- summeriseCorrectedResults(listPTS,list_patterns,list_probs)
dataSumCIs <- calcConfidenceIntervals(dataSum)
