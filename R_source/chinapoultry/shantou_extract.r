rm(list=ls(all=TRUE))
source("./shantou_funcs.r")

fnFName <- function(s) {paste("/Volumes/NO\ NAME/data/influenza/shantou/data/Shantou 20",s," (Steven).csv",sep="")}
fIntermediateFile <- "./intermediatev2.csv"

yearsinclude <- c("00","01","02","03","04","05","06")
# yearsinclude <- c("00")

matchedtf <- c(	"Sample.Type",
				"Sample.Type",
				"Sample.Type",
				"Animal.and.SampleType",
				"Animal.and.SampleType",
				"Animal.and.SampleType",
				"Animal.and.SampleType") # from "Sample.Type", "Animal.and.SampleType" 

matchedsub <- c("Subtype",
				"Subtype",
				"Subtype",
				"Subtype",
				"HA.subtype",
				"HA.subtype",
				"HA.subtype") # One from "Subtype", "HA.subtype"

mainDat <- extractShantouData(yearsinclude,matchedtf,matchedsub,blocksize=10000)

write.csv(mainDat,file=fIntermediateFile,row.names=FALSE)
