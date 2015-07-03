#
# Script to generate two sequences of randomised assignments of 1500+ households to three groups
# group 1 = control group;  group 2 = hand hygiene;   group 3 = hand hygiene plus face masks
# Software used: R Version 2.3.1 
# Date: 4 January 2008
#

set.seed(1234)
tmp <- sample(c(18, 24, 30), 300, replace=TRUE)
y <- NULL
for(i in 1:70) y <- c(y, sample(rep(1:3, tmp[i]/3), tmp[i], replace=FALSE))
length(y)
z <- factor(y, levels=1:3, labels=c("control", "hand hygiene", "masks and HH"))
write.csv(data.frame(hhID=1:length(y), group=y, intervention=z), "C:\\codesQuickVue.csv", row.names=FALSE)

set.seed(5678)
tmp <- sample(c(18, 24, 30), 300, replace=TRUE)
y <- NULL
for(i in 1:70) y <- c(y, sample(rep(1:3, tmp[i]/3), tmp[i], replace=FALSE))
length(y)
z <- factor(y, levels=1:3, labels=c("control", "hand hygiene", "masks and HH"))
write.csv(data.frame(hhID=1:length(y), group=y, intervention=z), "C:\\codesNoQuickVue.csv", row.names=FALSE)


#
# the end
#
#
