source("SR_misc.r")
source("resistance_sens/funcs.r")

fileSensData <- "D:/files/projects/influenza/Resistance/sens_data/sensitivity_v1.csv"

pTBounds <- c(0.2,0.4,0.6,0.8)
pTIndices <- c(1,999,1997,2995,3992)

sens_p4 <- read.csv(fileSensData,skip=pTIndices[1]-1,nrows=pTIndices[2]-pTIndices[1])


