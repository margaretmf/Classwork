


# Group members: Margaret, Joseph, Young Gee 
# Data: US News College Rankings and Reviews - Fall 2013 Acceptance Rate
# Mini-project


####################  Reading in and cleaning the dataset ########################################

## Reading in and creating a data frame

#school <- read.csv("SchoolData.csv", header = T)
school <- read.table("SchoolData.txt", header=T, sep='\t')

## Converting setting variable into categorical.

school$settingIndex <- c(0)

school$setting <- as.character(school$setting)

a <- school$setting
school$settingIndex <- ifelse(a == "suburban", 3,
                              ifelse(a == "urban" | a == "city", 2,
                                     ifelse(a == "rural", 1, 0)))

## Converting academicCal variable into categorical.

school$academicCalIndex <- c(0)

school$academicCal <- as.character(school$academicCal)
summary(as.factor(school$academicCal))

b <- school$academicCal
school$academicCalIndex <- ifelse(b == "semester", 1,
                                  ifelse(b == "quarter", 2, 3))

## Checking the variable classes
sapply(school,class)

## Changing to categorical variables 
school$region <- as.factor(school$region)
school$settingIndex <- as.factor(school$settingIndex)
school$academicCalIndex <- as.factor(school$academicCalIndex)
school$type <- as.factor(school$type)

## Delete unneeded variables: state, setting, academicCal
school <- school[,-c(2,5,6)]

## Divide the list into 4 groups: national, liberal, national.testing, liberal.testing
## Combine into two data frames: colleges (testing data) and colleges.testing (testing data)
national <- school[1:40,]
liberal <- school[53:82,]
national.testing <- school[41:45,]
liberal.testing <- school[83:87,]
colleges <- rbind(national, liberal)
colleges.testing <- rbind(national.testing,liberal.testing)



####################  Method 1: Making a linear model with AIC and BIC method ########################################

colleges.lm <- lm(acceptRate~+region+type+rank+tuitionIS+tuitionOS+enrollment+
                    yearFounded+retentionRate+gradRate+FacperStu+classLT20+classGT50+counselorScore+endowment+
                    perFemale+AvgDebt+settingIndex+academicCalIndex, data = na.omit(colleges))

# testing multicollinearity
vif(colleges.lm) > 5

# summary of linear model
summary(colleges.lm)

# testing stepAIC 
library(MASS)
step1 <- stepAIC(colleges.lm, direction="both", k=2)

# Step:  AIC=-443.51
# acceptRate ~ region + type + tuitionOS + enrollment + retentionRate + 
# gradRate + classLT20 + counselorScore + perFemale + academicCalIndex

# testing stepBIC
step2 <- stepAIC(colleges.lm, direction="both", k=log(20))

# Step:  AIC=-429.57
# acceptRate ~ region + type + tuitionOS + enrollment + retentionRate + 
# gradRate + classLT20 + counselorScore + perFemale + academicCalIndex

# creating a reduced model
colleges.reduced.lm <- lm(acceptRate ~ region + type + tuitionOS + enrollment + retentionRate + 
                            gradRate + classLT20 + counselorScore + perFemale + academicCalIndex, data = colleges)
summary(colleges.reduced.lm)

# prediction
predict(colleges.reduced.lm, colleges.testing, interval="conf")
colleges.testing[,9]
# missed: observations 44, 83, 84



#################### Method 2: making a linear model with regsubsets() function method ####################################

library(leaps)
c.sub <- regsubsets(acceptRate~+region+type+rank+tuitionIS+tuitionOS+enrollment+
                      yearFounded+retentionRate+gradRate+FacperStu+classLT20+classGT50+counselorScore+endowment+
                      perFemale+AvgDebt+settingIndex+academicCalIndex, data = na.omit(colleges), nbest=2, nvmax = 20)

rsq <- summary(c.sub)$rsq
adjr2 <- summary(c.sub)$adjr2
cp <- summary(c.sub)$cp
bic <- summary(c.sub)$bic
cbind(summary(c.sub)$which,rsq, adjr2,cp,bic)

regsubsets.lm <- lm(acceptRate ~ region+enrollment+retentionRate+gradRate+FacperStu+classLT20+counselorScore+perFemale+academicCalIndex, data = colleges)
summary(regsubsets.lm)
vif(regsubsets.lm)

## Linear Model with Interaction Variable 'type:tuitionOS'
regsubsets2.lm <- lm(acceptRate~region+type:tuitionOS+enrollment+retentionRate+gradRate+FacperStu+classLT20+counselorScore+perFemale+academicCalIndex, data=colleges)
summary(regsubsets2.lm)

# Linear Model without "FacperStu"
regsubsets3.lm <- lm(acceptRate~region+type:tuitionOS+enrollment+retentionRate+gradRate+classLT20+counselorScore+perFemale+academicCalIndex, data=colleges)
summary(regsubsets3.lm)

# Comparing the BIC values of these two models (regsubsets3.lm is lower)
library(nlme)
BIC(regsubsets2.lm)
BIC(regsubsets3.lm)

# prediction
predict(regsubsets3.lm, colleges.testing, interval="predict")
colleges.testing[,9]



########################### Regression Diagnostics using Graphics  ###################################

# Analysis of variance
anova(regsubsets3.lm)
# This shows all the variables except enrollment have some correlation with the response variable.
# The enrollment variable, however, should not be deleted because it can have some impact on 
# other variables. 

# Scatterplots for each variables against acceptRate
layout(matrix(1:4, nrow = 2))
plot(acceptRate~region+type:tuitionOS+enrollment+retentionRate+gradRate+FacperStu+classLT20+counselorScore+perFemale+academicCalIndex, data=colleges)
# These 11 scatterplots show the relationship between the response variable and 
# each explanatory variable. Some variables such as couselorScore, gradRate, and 
# retentionRate show strong correlation (either positive or negative). In contrast, 
# variables such as tuitionOS, enrollment, and perFemale don’t have strong 
# relationships.  Dropping those weakly related variables is invalid because it is 
# not always possible to see the relationship from a scatterplot, especially when 
# the effect of other predictors has not been accounted for.

#  Plot of residuals against fitted values
dev.off()
plot(regsubsets3.lm,pch=20,cex=.2,which=1)
# If the points in a residual plot are randomly dispersed around the horizontal axis, a 
# linear regression model is appropriate for the data; otherwise, a non-linear model is 
# more appropriate. This plot indicates a fairly random pattern; this data is fit to a
# linear model. Observations 53, 55, and 57 give rather large residual values.

# Normal probability plot of residuals
plot(regsubsets3.lm,pch=20,cex=.2,which=2)
# This plot is a graphical technique for assessing whether or not a data set is
# approximately normally distributed. The points in middle on this plot form a nearly 
# linear pattern, which indicates that the normal distribution is a good model for this 
# data set. However few data at the both tails have tendency to be far away from the 
# straight line. It shows observation 53,55, and 57 are extreme again. This indicates 
# that one of major assumptions is violated such as it is not a random drawings, or 
# from a fixed distribution.

# Index plot of Cook’s distances
plot(regsubsets3.lm,pch=20,cex=.2,which=4)
# Cook’s distance can have different thresholds based on the data. In this case, since 
# there are 70 observations, the threshold can be calculated by 4/N = 0.057. It 
# suggests that observation 73, 79, and 80 have undue influence on the estimated 
# regression coefficients, but the three outliers (53,55,57) identified previously do not. 

###########################################################################################################################

# Additional Commentary

# A regression model with the interaction variable "type:tuitionIS" instead of "type:tuitionOS"
# Both models show the positive correlation between acceptance rate and tuition, particularly at public universities
regsubsets3.lm <- lm(acceptRate~region+type:tuitionIS+enrollment+retentionRate+gradRate+classLT20+counselorScore+perFemale+academicCalIndex, data=colleges)
summary(regsubsets3.lm)


# Reexamining "perFemale" variable without the military academies and the all-female colleges shows that
# there is not a strong correlation between the two variables
x <- colleges[colleges$perFemale>0.25 & colleges$perFemale<0.8,]
x.lm <- lm(acceptRate~perFemale, data=x)
summary(x.lm)

# Plots of the "perFemale" with the training data acceptance rate and testing data acceptance rate
# indicate a possible explanation for the incorrect prediction of Barnard College's acceptance rate.
plot(colleges$acceptRate~colleges$perFemale)
plot(colleges.testing$acceptRate~colleges.testing$perFemale)



##################################################################################################


## Linear model with subset of data
colleges3.lm <- lm(acceptRate~type+tuitionOS+
                     gradRate+retentionRate+
                     academicCalIndex, data = colleges.no.mil)
summary(colleges3.lm)


predict(colleges3.lm, national.testing, interval="conf",na.action=na.pass, expand.na = F)
national.testing[,12]

predict(colleges3.lm, liberal.testing, interval="conf",na.action=na.pass, expand.na = F)
liberal.testing[,12]

## Testing variables with just national universities
n.lm <- lm(acceptRate~type+tuitionOS+
             gradRate+retentionRate+
             academicCalIndex, data=national)
summary(n.lm)

predict(n.lm, national.testing, interval="conf",na.action=na.pass, expand.na = F)
national.testing[,12]

## Testing variables with just liberal arts colleges
l.lm <- lm(acceptRate~type+
             gradRate+retentionRate+
             academicCalIndex+tuitionOS, data=liberal)
summary(l.lm)

# Eliminate tuitionOS variable because it is not significant to the model
l.lm <- lm(acceptRate~type+
             gradRate+retentionRate+
             academicCalIndex, data=liberal)
summary(l.lm)

predict(l.lm, liberal.testing, interval="conf",na.action=na.pass, expand.na = F)
liberal.testing[,12]
