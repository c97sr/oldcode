#
# Program to evaluate secondary attack rate and factors affecting the SAR
#

# extract age/sex for contacts and their corresponding index cases,
# the intervention arm, and vaccination status in 2006 OR 2007

demo <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_15_demographic_m.csv")
demo$hhID <- tolower(as.character(demo$hhID))
demo$male <- 1*(demo$sex=="m")

index <- demo[demo$member==0,c("hhID", "age", "male")]
names(index)[2:3] <- c("index.age", "index.male")

contacts <- demo[demo$member!=0,c("hhID", "member", "age", "male")]

housechar <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_15_housechar_h.csv")
housechar$hhID <- tolower(as.character(housechar$hhID))
housechar <- housechar[,c("hhID", "intervention")]
housechar$intervention <- factor(housechar$intervention, levels=1:3, labels=c("control", "mask", "hand"))

vaccine <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_15_baseflu_m.csv")
vaccine$hhID <- tolower(as.character(vaccine$hhID))
vaccine$vaccine <- 1*(vaccine$vaccine06==1 | vaccine$vaccine07==1)
vaccine <- vaccine[vaccine$member==0,c("hhID", "vaccine")]

dim(contacts)
new.demo <- merge(housechar, contacts)
dim(new.demo)
new.demo <- merge(new.demo, vaccine, all=TRUE)
new.demo <- merge(new.demo, index, all=TRUE)
dim(new.demo)
length(unique(new.demo$hhID))

# get the excluded households and the laboratory secondary cases

lab <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_17_home_culture.csv")
lab$hhID <- tolower(as.character(lab$hhID))

index.excl0 <- lab[lab$member==0 & lab$visit==0 & !is.na(lab$culture) & lab$culture=="0",c("hhID", "culture")]
names(index.excl0)[2] <- "culture0"
index.excl1 <- lab[lab$member==0 & lab$visit==1 & !is.na(lab$culture) & lab$culture=="0",c("hhID", "culture")]
names(index.excl1)[2] <- "culture1"

excluded <- merge(index.excl0, index.excl1)
excluded$excluded.index <- 1
excluded <- excluded[,c(1,4)]

second2 <- lab[lab$member!=0 & lab$visit==2 & !is.na(lab$culture),c(1,2,6)]
names(second2)[3] <- "culture2"
second3 <- lab[lab$member!=0 & lab$visit==3 & !is.na(lab$culture),c(1,2,6)]
names(second3)[3] <- "culture3"
second4 <- lab[lab$member!=0 & lab$visit==4 & !is.na(lab$culture),c(1,2,6)]
names(second4)[3] <- "culture4"

second <- merge(second2, second3, all=TRUE)
second <- merge(second, second4, all=TRUE)

second$culture2[is.na(second$culture2)] <- "0"
second$culture3[is.na(second$culture3)] <- "0"
second$culture4[is.na(second$culture4)] <- "0"

second$lab.second <- 1*(second$culture2!="0" |  second$culture3!="0" | second$culture4!="0")

second <- second[,c("hhID", "member", "lab.second")]

dim(new.demo)
new.data <- merge(new.demo, excluded, all=TRUE)
dim(new.data)
new.data <- merge(new.data, second, all=TRUE)
dim(new.data)
length(unique(new.data$hhID))

# !! why does this increase -- are some member IDs in the lab data miscoded??

# get the clinical secondary cases

symptoms <- read.csv("P:\\GROUP\\NPIpilot\\data\\2007_10_15_symptomday_d.csv")
symptoms$hhID <- tolower(as.character(symptoms$hhID))

symptoms$vague <- 1*(symptoms$fever==1 | sum(symptoms$cough + symptoms$sthroat + symptoms$rnose +
  symptoms$tired + symptoms$headache + symptoms$sjoint)>=2)

symptoms$fcs <- 1*(symptoms$fever==1 & (symptoms$cough==1 | symptoms$sthroat==1))

new.data$vague.second <- 0
new.data$fcs.second <- 0
for(i in 1:nrow(new.data)){
  new.data$vague.second[i] <- 1*(sum(symptoms$vague[symptoms$hhID==new.data$hhID[i] &
    symptoms$member==new.data$member[i]])>0)
  new.data$fcs.second[i] <- 1*(sum(symptoms$fcs[symptoms$hhID==new.data$hhID[i] &
    symptoms$member==new.data$member[i]])>0)
}
table(new.data$vague.second, exclude=NULL)
table(new.data$fcs.second, exclude=NULL)


# calculate SARs 

new.data$lab.second[is.na(new.data$lab.second)] <- 0
new.data$vague.second[is.na(new.data$vague.second)] <- 0
new.data$fcs.second[is.na(new.data$fcs.second)] <- 0

dim(new.data)
new.data <- new.data[is.na(new.data$excluded),]   # exclude the index subjects who were -ve at baseline
dim(new.data)

length(unique(new.data$hhID[new.data$intervention=="control"]))
length(unique(new.data$hhID[new.data$intervention=="mask"]))
length(unique(new.data$hhID[new.data$intervention=="hand"]))

lab.table <- table(new.data$lab.second, new.data$intervention)
fcs.table <- table(new.data$fcs.second, new.data$intervention)

lab.table
t(round(t(lab.table)/colSums(lab.table), 2))

fcs.table
t(round(t(fcs.table)/colSums(fcs.table), 2))

# fit regression models

require(gee)

dim(new.data)
newer.data <- new.data[!is.na(new.data$intervention) & !is.na(new.data$index.age) & 
  !is.na(new.data$index.male) & !is.na(new.data$vaccine),]
dim(newer.data)

#!! why do we lose household contacts?

newer.data$index.adult <- cut(newer.data$index.age, c(0,15,100))
newer.data$adult <- cut(newer.data$age, c(0,15,100))

table(newer.data$intervention)
table(newer.data$index.adult)
table(newer.data$index.male)
table(newer.data$vaccine)

gee.l <- gee(lab.second~intervention + index.adult + index.male + vaccine,
  id=factor(hhID), data=newer.data, family="binomial")
#summary(gee.l)
results.l <- data.frame(beta=gee.l$coef, se=sqrt(diag(gee.l[[20]])),
  row.names=names(gee.l$coef))
results.l$OR <- exp(results.l$beta)
results.l$lower.CI <- exp(results.l$beta - 1.96*results.l$se)
results.l$upper.CI <- exp(results.l$beta + 1.96*results.l$se)
results.l$p.value <- 2*pnorm(-1*abs(results.l$beta/results.l$se))
round(results.l, 2)[,-(1:2)]

gee.c <- gee(fcs.second~intervention + index.adult + index.male,
  id=factor(hhID), data=newer.data, family="binomial")
#summary(gee.l)
results.c <- data.frame(beta=gee.c$coef, se=sqrt(diag(gee.c[[20]])),
  row.names=names(gee.c$coef))
results.c$OR <- exp(results.c$beta)
results.c$lower.CI <- exp(results.c$beta - 1.96*results.c$se)
results.c$upper.CI <- exp(results.c$beta + 1.96*results.c$se)
results.c$p.value <- 2*pnorm(-1*abs(results.c$beta/results.c$se))
round(results.c, 2)[,-(1:2)]


#
# the end
#
#
