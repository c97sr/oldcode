# To Dos
# - Open and append
# - Include output on the values of the phi variables

rm(list=ls(all=TRUE))
options(error=recover)
require("odesolve")
require("date")

vecArgs <- commandArgs(trailingOnly=TRUE)

cat(vecArgs)

# interpret command line arguments
if (length(vecArgs)>1) stop("Only one argument required") else 
if (length(vecArgs)==0) version <- "local" else
version <- vecArgs[1]

if (version=="local")  dirRH1N1 <- "../" else
if (version=="server") dirRH1N1 <- "../../src/H1N1pdm/"

source(paste(dirRH1N1,"funcs.R",sep=""))
source("default_QLD_params.R")

prevData <- read.csv("Pub_Aus_Data.csv")
prevData$Date <- as.date(as.character(prevData$Date),order="mdy")
incData <- read.csv("NEJM_ICU_QLD_rounded.csv")
incData$date_start <- as.date(as.character(incData$date_start),order="dmy")	
incData$date_finish <- as.date(as.character(incData$date_finish),order="dmy")	
