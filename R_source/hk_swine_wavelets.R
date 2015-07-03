rm(list=ls(all=TRUE))
options(recover=NULL)
require("wmtsa")

hksd <- read.csv("/Volumes/NO NAME/data/influenza/hk_swine/swine_prevalence.csv")
hksd$Date <- as.date(as.character(hksd$Date),order="dmy")

# First, plot the overall timeseries
plot(hksd$Date,hksd$Total)

# Test that there is a significant decreasing trend
# The p-value of 1.73e-12 for the hksd$Date coefficient is
# good evidence that the incidence is falling. Its not a perfect model
# for the data, but a slope through the points for the mean of the incidence
# does a lot better than a flat line
model1 <- glm(hksd$Total~hksd$Date,poisson)
summary(model1)
#plot(model1)

# Decompose the incidence into wavelets to see if there is a signature for annual periodicity
plot(hksd$Date,hksd$Total,type="l")
swine.ts <- as.ts(as.numeric(hksd$Total),frequency=12)
cwt <- wavCWT(swine.ts)

# Not sure how to interpret the following data
# but I should be able to match up the continuous and the discrete
# dwt <- wavDWT(swine.ts)
# summary(dwt)

pdf(file="./wavelet_power_chart.pdf",width=15/cm(1),height=15/cm(1))
plot(cwt)
dev.off()

cat("On a log2 scale 12 months is ",log2(12)," and 36 moths is ", log2(36))

