source ("D:/Work/SVNrepository/NPIpilot/Calvin\\080125 qv paper/070801 qv function.r")

cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")
source ("D:/Work/SVNrepository/NPIpilot/Calvin/080125 qv paper/080211_add_groups.r")

#31/5/07 change QV positive to 1 and negative and na to 0 
QVava <- cdc[(cdc$QVres ==1 |cdc$QVres ==2|cdc$QVres ==3) & !is.na(cdc$QVres),]
cultava <- cdc[!is.na(cdc$cultpos),]
cultpos <-QVava[(QVava$culture=="A"|QVava$culture=="B")& !is.na(QVava$culture),]

#unknown onsettime to be 36 hours and onsettime >48hrs to be 48hrs
cultpos$onsettime[cultpos$onsettime==9] <- NA
cultpos$agegroup[cultpos$agegroup==9] <- NA
#cultpos$onsettime[cultpos$onsettime==5] <- 4


#check for logistic regression

#original lm
glm.fit <- glm( QVpos ~ factor(agegroup)+male+factor(onsettime)+factor(qvdonegp)+fcs
,family = binomial, data=cultpos, 
  subset=(onsettime!=5)
 #subset=(onsettime!=5 & onsettime!=9 & agegroup!=9)
  ) 


results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 2)

#write.csv(round(results.glm, 3), "d:\\r\\temp\\xxx.csv")

# 071119 testing for spline functions
require(Design)

glm.fit <- glm( QVpos ~ rcs(age, 5),
  family = binomial, data=cultpos, subset=(onsettime!=5 & onsettime!=9)
  )  


results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 3)


summary(z)
exp(coef(z))


# for plotting
newdata <- data.frame(age <- 5:75)
plot(predict(glm.fit, newdata))


