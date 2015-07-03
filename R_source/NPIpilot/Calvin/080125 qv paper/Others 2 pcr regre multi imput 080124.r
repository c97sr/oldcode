#run clean clinic data first

library(Hmisc)

#080121 get those with PCR results 
cdc <- read.csv("P:/GROUP/NPIpilot/data/qv paper data/cleandata/qvdat.csv")

# 080122 regroup the data for regression analysis
#cdcPCR$agegroup2[cdcPCR$agegroup2==1 ] <-2

#08/01/31 add new entry agegroup

cdc$agegroup2[cdc$age <=15 &!is.na(cdc$age)]<- 1
cdc$agegroup2[cdc$age >15 & cdc$age <=50 &!is.na(cdc$age)]<-2
cdc$agegroup2[cdc$age >50 &!is.na(cdc$age)]<-3
#cdc$agegroup[is.na(cdc$age)] <-9



cdc$onsettime[cdc$onsettime==9] <- NA
cdc$agegroup2[cdc$agegroup2==9] <- NA

dim(cdc)

#############################

QVava <- cdc[(cdc$QVres ==1 | cdc$QVres ==2 | cdc$QVres ==3) & !is.na(cdc$QVres),]
cultava <- cdc[!is.na(cdc$cultpos),]
cultpos <-QVava[(QVava$culture=="A" | QVava$culture=="B") & !is.na(QVava$culture),]

cultpos$agegroup2 <- factor(cultpos$agegroup2, levels=c(2,1,3))
cultpos$onsettime[cultpos$onsettime==5] <- 4
cultpos$onsettime <- factor(cultpos$onsettime)
cultpos$qvdonegp <- factor(cultpos$qvdonegp)

dim(cultpos)
cultposA <- cultpos[cultpos$cultposA==1,]
dim(cultposA)

#cultposPCR.nomiss <- cultposPCR[cultposPCR$onsettime!=5, ]

glm.fit <- glm( QVposA ~ onsettime + agegroup2 + qvdonegp + male, family = binomial, data=cultposA) 

results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 2)

###########
set.seed(15)
cultposA.i <- aregImpute( ~ onsettime + agegroup2 + qvdonegp + male + QVposA, data=cultposA, n.impute=10)

cultposA.nomiss <- list(cultposA, cultposA, cultposA, cultposA, cultposA,
  cultposA, cultposA, cultposA, cultposA, cultposA)

for(i in 1:10){
	cultposA.nomiss[[i]]$onsettime[is.na(cultposA.nomiss[[i]]$onsettime)] <-
		cultposA.i$imputed$onsettime[,i]
	cultposA.nomiss[[i]]$agegroup2[is.na(cultposA.nomiss[[i]]$agegroup2)] <-
		cultposA.i$imputed$agegroup2[,i] 

}
set.seed(15)
model1 <- list(NA)
for(i in 1:10){
	model1[[i]] <- glm( QVpos ~ onsettime + agegroup2 + qvdonegp + male,
	  family = binomial, data=cultposA.nomiss[[i]]) 
}

#
# Combine imputed results and summarise parameter estimates
#

combine.mi <- function(model, n.impute=10){
	betas <- matrix(c(model[[1]]$coef, model[[2]]$coef, 
		model[[3]]$coef, model[[4]]$coef, model[[5]]$coef, 
		model[[6]]$coef, model[[7]]$coef, model[[8]]$coef, 
		model[[9]]$coef, model[[10]]$coef), byrow=FALSE, ncol=10)
	vars <- matrix(c(diag(vcov(model[[1]])), diag(vcov(model[[2]])), diag(vcov(model[[3]])), 
	diag(vcov(model[[4]])), diag(vcov(model[[5]])), diag(vcov(model[[6]])), diag(vcov(model[[7]])), 
	diag(vcov(model[[8]])), diag(vcov(model[[9]])), diag(vcov(model[[10]]))), byrow=FALSE, ncol=10)
	coef.names <- names(model[[1]]$coef)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	data.frame(beta = mean.coefs,
		lowerCI = mean.coefs - qt(0.975, df=degf)*sqrt(T),
		upperCI = mean.coefs + qt(0.975, df=degf)*sqrt(T),
		p.value = 2*(1 - pt(abs(mean.coefs)/sqrt(T), df=degf)),
		row.names=coef.names)
}

round(combine.mi(model1, n.impute=10), 3)

#
# now add the qPCR for A +ve
#

cultposA$qPCRgp <-cut(cultposA$logqPCR, breaks=c(0,5,6,10))
#cultposA$qPCRgp <-cut(cultposA$logqPCR, breaks=c(0,5,6,7,10))
#breakPCR <-quantile(na.omit(cultposA$logqPCR), probs = seq(0, 1, 0.20)) 
#cultposA$qPCRgp <- cut(cultposA$logqPCR, breaks=c(0,breakPCR[2:5],10)) 
cultposA$qPCRgp <- as.numeric(cultposA$qPCRgp)
cultposA$qPCRgp <- factor(cultposA$qPCRgp,levels=c(2,1,3) )



#########
glm.fit <- glm( QVposA ~ onsettime + agegroup2 + qvdonegp + male + qPCRgp, family = binomial, data=cultposA) 

results.glm <- data.frame(beta=glm.fit$coef, se=sqrt(diag(summary(glm.fit)[[17]])),
  row.names=names(glm.fit$coef))

results.glm$OR <- exp(results.glm$beta)
results.glm$lower.CI <- exp(results.glm$beta - 1.96*results.glm$se)
results.glm$upper.CI <- exp(results.glm$beta + 1.96*results.glm$se)
results.glm$p.value <- 2*pnorm(-1*abs(results.glm$beta/results.glm$se))

round(results.glm, 2)

###

set.seed(15)
cultposA.ipcr <- aregImpute( ~  onsettime + agegroup2 + qvdonegp 
+ male + QVposA + qPCRgp , data=cultposA, n.impute=10)

cultposA.nomisspcr <- list(cultposA, cultposA, cultposA, cultposA, cultposA,
  cultposA, cultposA, cultposA, cultposA, cultposA)

for(i in 1:10){
	cultposA.nomisspcr[[i]]$onsettime[is.na(cultposA.nomisspcr[[i]]$onsettime)] <-
		cultposA.ipcr$imputed$onsettime[,i]
	cultposA.nomisspcr[[i]]$agegroup2[is.na(cultposA.nomisspcr[[i]]$agegroup2)] <-
		cultposA.ipcr$imputed$agegroup2[,i] 
	cultposA.nomisspcr[[i]]$qPCRgp[is.na(cultposA.nomisspcr[[i]]$qPCRgp)] <-
		cultposA.ipcr$imputed$qPCRgp[,i] 

}

set.seed(15)
model2 <- list(NA)
for(i in 1:10){
	model2[[i]] <- glm( QVposA ~ onsettime + agegroup2 + qvdonegp + male + qPCRgp,
	  family = binomial, data=cultposA.nomisspcr[[i]]) 
}

#
# Combine imputed results and summarise parameter estimates
#

combine.mi <- function(model, n.impute=10){
	betas <- matrix(c(model[[1]]$coef, model[[2]]$coef, 
		model[[3]]$coef, model[[4]]$coef, model[[5]]$coef, 
		model[[6]]$coef, model[[7]]$coef, model[[8]]$coef, 
		model[[9]]$coef, model[[10]]$coef), byrow=FALSE, ncol=10)
	vars <- matrix(c(diag(vcov(model[[1]])), diag(vcov(model[[2]])), 
		diag(vcov(model[[3]])), diag(vcov(model[[4]])), diag(vcov(model[[5]])), 
		diag(vcov(model[[6]])), diag(vcov(model[[7]])), diag(vcov(model[[8]])), 
		diag(vcov(model[[9]])), diag(vcov(model[[10]]))), byrow=FALSE, ncol=10)
	coef.names <- names(model[[1]]$coef)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	data.frame(beta = mean.coefs,
		lowerCI = mean.coefs - qt(0.975, df=degf)*sqrt(T),
		upperCI = mean.coefs + qt(0.975, df=degf)*sqrt(T),
		p.value = 2*(1 - pt(abs(mean.coefs)/sqrt(T), df=degf)),
		row.names=coef.names)
}

tmp1 <-combine.mi(model1, n.impute=10)
round(cbind(exp(tmp1[,1:3]), tmp1[,4]), 2)

tmp2 <-combine.mi(model2, n.impute=10)
round(cbind(exp(tmp2[,1:3]), tmp2[,4]), 2)


#round(combine.mi(model1, n.impute=10), 3)
#round(combine.mi(model2, n.impute=10), 3)

table( cultposA$QVposA ,cultposA$qPCRgp)