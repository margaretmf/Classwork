

require(reshape2)
require(plyr)
require(ROCR)
require(Hmisc)
require(FSelector)


# Lab model ---------------------------------

# data l24
data.l24 = read.csv("data.1.labs.csv", header=TRUE)

data.l24$SEPSIS = as.factor(data.l24$SEPSIS)
data.l24$SEPSIS

# cfs - uses best.first.search
features.bestfirst.lab <- cfs(SEPSIS~., data.l24) 
data.bestfirst.lab <- data.l24[,c(features.bestfirst.lab, "X50383.MEDIAN","X50386.MEDIAN","SEPSIS")]
data.bestfirst.lab$SEPSIS = as.factor(data.bestfirst.lab$SEPSIS)


# Modelling

# Create test and train data
# We need to draw samples from control group and target group respectively

# First, get control and target group ids.

control.id.lab <- data.l24[data.l24$SEPSIS==0,]$SUBJECT_ID
target.id.lab <- data.l24[data.l24$SEPSIS==1,]$SUBJECT_ID

require(caret)

set.seed(12345)
control.folds.lab <- createFolds(control.id.lab, k=10)
target.folds.lab <- createFolds(target.id.lab, k=10)

# Initialization for storing metrics
# quantile evaluation
sp.quantile.lab.all <- c()
sp.pct.lab.all <- c()
sp.pop.lab.all <- c()
sp.sp.lab.all <- c()
# auc
auc.lab.all <- c()
# classification
class.metrics.lab <- c()
# coefficients
coeffs.lab.all <- c()

# Cross Validation
for (i in c(1:10)) {
  test.control.id.lab <- control.id.lab[unlist(control.folds.lab[i])]
  test.target.id.lab <- target.id.lab[unlist(target.folds.lab[i])]
  train.control.id.lab <- control.id.lab[unlist(control.folds.lab[!c(1:10) %in% i])] 
  train.target.id.lab <- target.id.lab[unlist(target.folds.lab[!c(1:10) %in% i])] 
  
  train.id.lab <- c(train.control.id.lab, train.target.id.lab)
  test.id.lab <- c(test.control.id.lab, test.target.id.lab)
  length(train.id.lab)+length(test.id.lab)
  nrow(data.l24)
  
  train.ind.lab <- which(data.l24$SUBJECT_ID %in% train.id.lab)
  test.ind.lab <- which(data.l24$SUBJECT_ID %in% test.id.lab)
  
  # Build models using bestfirst selection of features
  data.train.lab <- data.bestfirst.lab[train.ind.lab,]
  data.test.lab <- data.bestfirst.lab[test.ind.lab,]
  
  # Regression Trees
  library(rpart)
  library(ROCR)
  library(RWeka)
  
  # Logistic Model Trees
  lmt.lab <- LMT(SEPSIS~., data.train.lab, control="-A")  
  summary(lmt.lab)
  lmt.lab
  
  # Prediction
  lmt.predict.lab <- predict(lmt.lab, data.test.lab, type="probability")
  lmt.predict.lab <- lmt.predict.lab[,2]
  lmt.pred.lab <- prediction(lmt.predict.lab, data.test.lab$SEPSIS)
  lmt.perf.lab <- performance(lmt.pred.lab, "tpr", "fpr")
  lmt.auc.lab <- performance(lmt.pred.lab, "auc")
  plot(lmt.perf.lab)
  lmt.auc.lab
  
  # Evaluation
  # record proportion of sepsis to whole population
  prop.lab <- length(target.id.lab)/ (length(target.id.lab) + length(control.id.lab))
  #length(train.target.id)/length(train.control.id)
  # pick a threshold
  
  hist(lmt.predict.lab)
  
  th.lab <- 1
  error <- 0.005
  while (th.lab > 0) {
    sp <- lmt.predict.lab[lmt.predict.lab >= th.lab]
    nonsp <- lmt.predict.lab[lmt.predict.lab < th.lab]
    prop.predict <- length(sp) / (length(sp)+length(nonsp))
    if (abs(prop.predict - prop.lab) < error){
      break
    }else {
      th.lab <- th.lab - 0.005
    }
  }

  # assign class
  lmt.predict.class.lab <- lmt.predict.lab
  lmt.predict.class.lab[lmt.predict.lab >= th.lab] <- 1
  lmt.predict.class.lab[lmt.predict.lab < th.lab] <- 0
  lmt.predict.class.lab <- as.factor(lmt.predict.class.lab)
  class.result.lab <- confusionMatrix(lmt.predict.class.lab, data.test.lab$SEPSIS, positive="1")
  
  # Quantile evaluation
  # create a new column that label which quantile each observation is in
  lmt.predict.lab.quantiles <- cut(lmt.predict.lab, 
                                    breaks=quantile(lmt.predict.lab, probs=seq(0,1, by=0.1)), 
                                    include.lowest=TRUE)
  lmt.predict.lab.quantiles.mean <- cut2(lmt.predict.lab, 
                                          g=10, 
                                          levels.mean=T)
  lmt.predict.lab.df <- data.frame(data.test.lab, lmt.predict.lab)
  lmt.predict.lab.df$quantiles <- lmt.predict.lab.quantiles
  lmt.predict.lab.df$quantiles.mean <- as.numeric(as.character(lmt.predict.lab.quantiles.mean))
  require(ggplot2)
  lmt.predict.sp.pct.lab <- ddply(lmt.predict.lab.df, .(quantiles), 
                                   function(x){return(c(quantiles.mean=x$quantiles.mean[1],
                                                        percent.in.quantile=sum(x$SEPSIS==1)/nrow(x),
                                                        percent.in.pop=sum(x$SEPSIS==1)/nrow(data.test.lab),
                                                        percent.in.sepsis=sum(x$SEPSIS==1)/sum(data.test.lab$SEPSIS==1)))})
  
  # Save metrics
  # quantile
  sp.quantile.lab.all <- cbind(sp.quantile.lab.all, lmt.predict.sp.pct.lab[,2])
  sp.pct.lab.all <- cbind(sp.pct.lab.all, lmt.predict.sp.pct.lab[,3])
  sp.pop.lab.all<- cbind(sp.pop.lab.all, lmt.predict.sp.pct.lab[,4])
  sp.sp.lab.all <- cbind(sp.sp.lab.all, lmt.predict.sp.pct.lab[,5])
  # auc
  auc.lab.all <- rbind(auc.lab.all, lmt.auc.lab@y.values[[1]])
  # class metrics
  class.metrics.lab <- rbind(class.metrics.lab, class.result.lab[4][[1]])
  
}

sp.quantile.lab.all <- apply(sp.quantile.lab.all, 1, mean)
sp.pct.lab.all <- apply(sp.pct.lab.all, 1, mean)
sp.pop.lab.all <- apply(sp.pop.lab.all, 1, mean)
sp.sp.lab.all <- apply(sp.sp.lab.all, 1, mean)

sp.df <- data.frame(quantiles=c(1:10), quantiles.mean=sp.quantile.lab.all, percent.in.quantile=sp.pct.lab.all,
                    percent.in.pop=sp.pop.lab.all, percent.in.target=sp.sp.lab.all)

# save results
# quantile evaluation
write.csv(sp.df, "quantile_final.24lmt.lab.csv", row.names=F)
# auc
write.csv(auc.lab.all, "auc_final.24lmt.lab.csv", row.names=F)
# class evaluation
write.csv(class.metrics.lab,"class_final.24lmt.lab.csv", row.names=F)


predict.plot.lab <- ggplot(sp.df, aes(x=quantiles.mean, y=percent.in.quantile))
predict.plot.lab <- predict.plot.lab + geom_point()
predict.plot.lab <- predict.plot.lab + labs(x="Mean of Quantile", y="Percent in Quantile") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=70, vjust=0.5)) 
predict.plot.lab


# Vital model ---------------------------------

# data v24
data.v24 = read.csv("data.1.vitals.csv", header=TRUE)

data.v24$SEPSIS = as.factor(data.v24$SEPSIS)
data.v24$SEPSIS

# cfs - uses best.first.search
features.bestfirst.vital <- cfs(SEPSIS~., data.v24) 
data.bestfirst.vital <- data.v24[,c(features.bestfirst.vital, "SEPSIS")]
data.bestfirst.vital$SEPSIS = as.factor(data.bestfirst.vital$SEPSIS)


# Modelling

# Create test and train data
# We need to draw samples from control group and target group respectively

# First, get control and target group ids.

control.id.vital <- data.v24[data.v24$SEPSIS==0,]$SUBJECT_ID
target.id.vital <- data.v24[data.v24$SEPSIS==1,]$SUBJECT_ID

require(caret)

set.seed(12345)
control.folds.vital <- createFolds(control.id.vital, k=10)
target.folds.vital <- createFolds(target.id.vital, k=10)

# Initialization for storing metrics
# quantile evaluation
sp.quantile.vital.all <- c()
sp.pct.vital.all <- c()
sp.pop.vital.all <- c()
sp.sp.vital.all <- c()
# auc
auc.vital.all <- c()
# classification
class.metrics.vital <- c()
# coefficients
coeffs.vital.all <- c()

# Cross Validation
for (i in c(1:10)) {
  test.control.id.vital <- control.id.vital[unlist(control.folds.vital[i])]
  test.target.id.vital <- target.id.vital[unlist(target.folds.vital[i])]
  train.control.id.vital <- control.id.vital[unlist(control.folds.vital[!c(1:10) %in% i])] 
  train.target.id.vital <- target.id.vital[unlist(target.folds.vital[!c(1:10) %in% i])] 
  
  train.id.vital <- c(train.control.id.vital, train.target.id.vital)
  test.id.vital <- c(test.control.id.vital, test.target.id.vital)
  length(train.id.vital)+length(test.id.vital)
  nrow(data.v24)
  
  train.ind.vital <- which(data.v24$SUBJECT_ID %in% train.id.vital)
  test.ind.vital <- which(data.v24$SUBJECT_ID %in% test.id.vital)
  
  # Build models using bestfirst selection of features
  data.train.vital <- data.bestfirst.vital[train.ind.vital,]
  data.test.vital <- data.bestfirst.vital[test.ind.vital,]
  
  # Regression Trees
  library(rpart)
  library(ROCR)
  library(RWeka)
  
  # Logistic Model Trees
  lmt.vital <- LMT(SEPSIS~., data.train.vital, control="-A")  
  summary(lmt.vital)
  lmt.vital
  
  # Prediction
  lmt.predict.vital <- predict(lmt.vital, data.test.vital, type="probability")
  lmt.predict.vital <- lmt.predict.vital[,2]
  lmt.pred.vital <- prediction(lmt.predict.vital, data.test.vital$SEPSIS)
  lmt.perf.vital <- performance(lmt.pred.vital, "tpr", "fpr")
  lmt.auc.vital <- performance(lmt.pred.vital, "auc")
  plot(lmt.perf.vital)
  lmt.auc.vital
  
  # Evaluation
  # record proportion of sepsis to whole population
  prop.vital <- length(target.id.vital)/ (length(target.id.vital) + length(control.id.vital))
  #length(train.target.id)/length(train.control.id)
  # pick a threshold
  
  hist(lmt.predict.vital)
  
  th.vital <- 1
  error <- 0.005
  while (th.vital > 0) {
    sp <- lmt.predict.vital[lmt.predict.vital >= th.vital]
    nonsp <- lmt.predict.vital[lmt.predict.vital < th.vital]
    prop.predict <- length(sp) / (length(sp)+length(nonsp))
    if (abs(prop.predict - prop.vital) < error){
      break
    }else {
      th.vital <- th.vital - 0.005
    }
  }
  
  # assign class
  lmt.predict.class.vital <- lmt.predict.vital
  lmt.predict.class.vital[lmt.predict.vital >= th.vital] <- 1
  lmt.predict.class.vital[lmt.predict.vital < th.vital] <- 0
  lmt.predict.class.vital <- as.factor(lmt.predict.class.vital)
  class.result.vital <- confusionMatrix(lmt.predict.class.vital, data.test.vital$SEPSIS, positive="1")
  
  # Quantile evaluation
  # create a new column that label which quantile each observation is in
  lmt.predict.vital.quantiles <- cut(lmt.predict.vital, 
                                   breaks=quantile(lmt.predict.vital, probs=seq(0,1, by=0.1)), 
                                   include.lowest=TRUE)
  lmt.predict.vital.quantiles.mean <- cut2(lmt.predict.vital, 
                                         g=10, 
                                         levels.mean=T)
  lmt.predict.vital.df <- data.frame(data.test.vital, lmt.predict.vital)
  lmt.predict.vital.df$quantiles <- lmt.predict.vital.quantiles
  lmt.predict.vital.df$quantiles.mean <- as.numeric(as.character(lmt.predict.vital.quantiles.mean))
  require(ggplot2)
  lmt.predict.sp.pct.vital <- ddply(lmt.predict.vital.df, .(quantiles), 
                                  function(x){return(c(quantiles.mean=x$quantiles.mean[1],
                                                       percent.in.quantile=sum(x$SEPSIS==1)/nrow(x),
                                                       percent.in.pop=sum(x$SEPSIS==1)/nrow(data.test.vital),
                                                       percent.in.sepsis=sum(x$SEPSIS==1)/sum(data.test.vital$SEPSIS==1)))})
  
  # Save metrics
  # quantile
  sp.quantile.vital.all <- cbind(sp.quantile.vital.all, lmt.predict.sp.pct.vital[,2])
  sp.pct.vital.all <- cbind(sp.pct.vital.all, lmt.predict.sp.pct.vital[,3])
  sp.pop.vital.all<- cbind(sp.pop.vital.all, lmt.predict.sp.pct.vital[,4])
  sp.sp.vital.all <- cbind(sp.sp.vital.all, lmt.predict.sp.pct.vital[,5])
  # auc
  auc.vital.all <- rbind(auc.vital.all, lmt.auc.vital@y.values[[1]])
  # class metrics
  class.metrics.vital <- rbind(class.metrics.vital, class.result.vital[4][[1]])
  
}

sp.quantile.vital.all <- apply(sp.quantile.vital.all, 1, mean)
sp.pct.vital.all <- apply(sp.pct.vital.all, 1, mean)
sp.pop.vital.all <- apply(sp.pop.vital.all, 1, mean)
sp.sp.vital.all <- apply(sp.sp.vital.all, 1, mean)

sp.df <- data.frame(quantiles=c(1:10), quantiles.mean=sp.quantile.vital.all, percent.in.quantile=sp.pct.vital.all,
                    percent.in.pop=sp.pop.vital.all, percent.in.target=sp.sp.vital.all)

# save results
# quantile evaluation
write.csv(sp.df, "quantile_final.24lmt.vital.csv", row.names=F)
# auc
write.csv(auc.vital.all, "auc_all_final.24lmt.vital.csv", row.names=F)
# class evaluation
write.csv(class.metrics.vital,"class_all_final.24lmt.vital.csv", row.names=F)


predict.plot.vital <- ggplot(sp.df, aes(x=quantiles.mean, y=percent.in.quantile))
predict.plot.vital <- predict.plot.vital + geom_point()
predict.plot.vital <- predict.plot.vital + labs(x="Mean of Quantile", y="Percent in Quantile") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=70, vjust=0.5)) 
predict.plot.vital


# All model ---------------------------------
library(digest)

# FOR LACTATE CLEARANCE
# # data a24
# cbind(data.l24, data.v24)[1:5,]
# #data.a24 = read.csv("data.2.all.csv", header=TRUE)
# data.a24 = cbind(data.l24, data.v24)
# data.a24 <- data.a24[!duplicated(lapply(data.a24, summary))]
# data.a24$SEPSIS = as.factor(data.a24$SEPSIS)
# table(data.a24$SEPSIS)
# # cfs - uses best.first.search
# features.bestfirst.all <- cfs(SEPSIS~., data.a24) 
# data.bestfirst.all <- data.a24[,c(features.bestfirst.all, "X50383.MEDIAN","X50386.MEDIAN", "SEPSIS")]
# data.bestfirst.all$SEPSIS = as.factor(data.bestfirst.all$SEPSIS)


# data v24
data.a24 = read.csv("data.1.all.csv", header=TRUE)

data.a24$SEPSIS = as.factor(data.a24$SEPSIS)
data.a24$SEPSIS

# cfs - uses best.first.search
features.bestfirst.all <- cfs(SEPSIS~., data.a24) 
data.bestfirst.all <- data.a24[,c(features.bestfirst.all, "SEPSIS")]
data.bestfirst.all$SEPSIS = as.factor(data.bestfirst.all$SEPSIS)

# Modelling

# Create test and train data
# We need to draw samples from control group and target group respectively

# First, get control and target group ids.

control.id.all <- data.a24[data.a24$SEPSIS==0,]$SUBJECT_ID
target.id.all <- data.a24[data.a24$SEPSIS==1,]$SUBJECT_ID

require(caret)

set.seed(12345)
control.folds.all <- createFolds(control.id.all, k=10)
target.folds.all <- createFolds(target.id.all, k=10)

# Initialization for storing metrics
# quantile evaluation
sp.quantile.all.all <- c()
sp.pct.all.all <- c()
sp.pop.all.all <- c()
sp.sp.all.all <- c()
# auc
auc.all.all <- c()
# classification
class.metrics.all <- c()
# coefficients
coeffs.all.all <- c()

# Cross Validation
for (i in c(1:10)) {
  test.control.id.all <- control.id.all[unlist(control.folds.all[i])]
  test.target.id.all <- target.id.all[unlist(target.folds.all[i])]
  train.control.id.all <- control.id.all[unlist(control.folds.all[!c(1:10) %in% i])] 
  train.target.id.all <- target.id.all[unlist(target.folds.all[!c(1:10) %in% i])] 
  
  train.id.all <- c(train.control.id.all, train.target.id.all)
  test.id.all <- c(test.control.id.all, test.target.id.all)
  length(train.id.all)+length(test.id.all)
  nrow(data.a24)
  
  train.ind.all <- which(data.a24$SUBJECT_ID %in% train.id.all)
  test.ind.all <- which(data.a24$SUBJECT_ID %in% test.id.all)
  
  # Build models using bestfirst selection of features
  data.train.all <- data.bestfirst.all[train.ind.all,]
  data.test.all <- data.bestfirst.all[test.ind.all,]
  
  # Regression Trees
  library(rpart)
  library(ROCR)
  library(RWeka)
  
  # Logistic Model Trees
  lmt.all <- LMT(SEPSIS~., data.train.all, control="-A")  
  summary(lmt.all)
  lmt.all
  
  # Prediction
  lmt.predict.all <- predict(lmt.all, data.test.all, type="probability")
  lmt.predict.all <- lmt.predict.all[,2]
  lmt.pred.all <- prediction(lmt.predict.all, data.test.all$SEPSIS)
  lmt.perf.all <- performance(lmt.pred.all, "tpr", "fpr")
  lmt.auc.all <- performance(lmt.pred.all, "auc")
  plot(lmt.perf.all)
  lmt.auc.all
  
  # Evaluation
  # record proportion of sepsis to whole population
  prop.all <- length(target.id.all)/ (length(target.id.all) + length(control.id.all))
  #length(train.target.id)/length(train.control.id)
  # pick a threshold
  
  hist(lmt.predict.all)
  
  th.all <- 1
  error <- 0.005
  while (th.all > 0) {
    sp <- lmt.predict.all[lmt.predict.all >= th.all]
    nonsp <- lmt.predict.all[lmt.predict.all < th.all]
    prop.predict <- length(sp) / (length(sp)+length(nonsp))
    if (abs(prop.predict - prop.all) < error){
      break
    }else {
      th.all <- th.all - 0.005
    }
  }
  
  # assign class
  lmt.predict.class.all <- lmt.predict.all
  lmt.predict.class.all[lmt.predict.all >= th.all] <- 1
  lmt.predict.class.all[lmt.predict.all < th.all] <- 0
  lmt.predict.class.all <- as.factor(lmt.predict.class.all)
  class.result.all <- confusionMatrix(lmt.predict.class.all, data.test.all$SEPSIS, positive="1")
  
  # Quantile evaluation
  # create a new column that label which quantile each observation is in
  lmt.predict.all.quantiles <- cut(lmt.predict.all, 
                                     breaks=quantile(lmt.predict.all, probs=seq(0,1, by=0.1)), 
                                     include.lowest=TRUE)
  lmt.predict.all.quantiles.mean <- cut2(lmt.predict.all, 
                                           g=10, 
                                           levels.mean=T)
  lmt.predict.all.df <- data.frame(data.test.all, lmt.predict.all)
  lmt.predict.all.df$quantiles <- lmt.predict.all.quantiles
  lmt.predict.all.df$quantiles.mean <- as.numeric(as.character(lmt.predict.all.quantiles.mean))
  require(ggplot2)
  lmt.predict.sp.pct.all <- ddply(lmt.predict.all.df, .(quantiles), 
                                    function(x){return(c(quantiles.mean=x$quantiles.mean[1],
                                                         percent.in.quantile=sum(x$SEPSIS==1)/nrow(x),
                                                         percent.in.pop=sum(x$SEPSIS==1)/nrow(data.test.all),
                                                         percent.in.sepsis=sum(x$SEPSIS==1)/sum(data.test.all$SEPSIS==1)))})
  
  # Save metrics
  # quantile
  sp.quantile.all.all <- cbind(sp.quantile.all.all, lmt.predict.sp.pct.all[,2])
  sp.pct.all.all <- cbind(sp.pct.all.all, lmt.predict.sp.pct.all[,3])
  sp.pop.all.all<- cbind(sp.pop.all.all, lmt.predict.sp.pct.all[,4])
  sp.sp.all.all <- cbind(sp.sp.all.all, lmt.predict.sp.pct.all[,5])
  # auc
  auc.all.all <- rbind(auc.all.all, lmt.auc.all@y.values[[1]])
  # class metrics
  class.metrics.all <- rbind(class.metrics.all, class.result.all[4][[1]])
  
}

sp.quantile.all.all <- apply(sp.quantile.all.all, 1, mean)
sp.pct.all.all <- apply(sp.pct.all.all, 1, mean)
sp.pop.all.all <- apply(sp.pop.all.all, 1, mean)
sp.sp.all.all <- apply(sp.sp.all.all, 1, mean)

sp.df <- data.frame(quantiles=c(1:10), quantiles.mean=sp.quantile.all.all, percent.in.quantile=sp.pct.all.all,
                    percent.in.pop=sp.pop.all.all, percent.in.target=sp.sp.all.all)

# save results
# quantile evaluation
write.csv(sp.df, "quantile_final.24lmt.all.csv", row.names=F)
# auc
write.csv(auc.all.all, "auc_final.24lmt.all.csv", row.names=F)
# class evaluation
write.csv(class.metrics.all,"class_final.24lmt.all.csv", row.names=F)


predict.plot.all <- ggplot(sp.df, aes(x=quantiles.mean, y=percent.in.quantile))
predict.plot.all <- predict.plot.all + geom_point()
predict.plot.all <- predict.plot.all + labs(x="Mean of Quantile", y="Percent in Quantile") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=70, vjust=0.5)) 
predict.plot.all


