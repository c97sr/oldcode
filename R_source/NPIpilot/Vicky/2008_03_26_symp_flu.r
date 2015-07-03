
#                                                                   #
# Logistic regression for flu(lab-confirmed)~symptoms, and draw ROC #
#                                                                   #

rawdata <- read.csv("P:\\GROUP\\NPIpilot\\data\\qv paper data\\cleandata\\qvdat.csv",header=T)

data <- rawdata[,c(1,2,5,6,8,9:29,37)]  # scrID, culture, sex, age, self_fever, headache-nsweat, onsettime, QVres, bodytemp, qPCR
data$fever <- 1*(data$bodytemp>=37.8)
for (i in 1:nrow(data)){
   if(is.na(data$fever[i])) data$fever[i] <- data$self_fever[i]
}
data$labconfirm <- 0
data$labconfirm[ ((data$culture=="A"|data$culture=="B"))
         | (!is.na(data$qPCR)&data$qPCR>0) ] <- 1
data <- data[,c(1,29,3,4,28,6:23,5,26,24,25,2,27)]
names(data)[c(24,26)] <- c("feverish","delay")

## fit logistic regression
library(gee)

inf.glm <- glm(labconfirm~fever+headache+lnode+rnose+hvoice+sthroat+cough+phlegm+sbreath+cpain+apain
                   +vomit+diarrhoea+sjoint+pmuscle+dizzy+tired+chill+nsweat,
               data=data, family="binomial")

inf.glm <- glm(labconfirm~fever+cough+rnose,
               data=data, family="binomial")

summary(inf.glm)

results <- data.frame(beta=inf.glm$coef, se=sqrt(diag(vcov(inf.glm))),
  row.names=names(inf.glm$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
round(results, 2) 

data$pred <- predict(inf.glm, data=data, type="response")


library(ROCR)
p1 <- prediction(data$pred, data$labconfirm)
ss <- performance(p1, "sens", "spec")
round(ss@y.values[[1]][3],2)  # sensitivity
round(ss@x.values[[1]][3],2)  # specificity
p2 <- performance(p1, "auc")
round(p2@y.values[[1]],2)     # area under roc

