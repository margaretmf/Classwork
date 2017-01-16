
# Group Statistical Modeling Class
# Group members: Youngwoon, Mike, Xuandong, Margaret
# Data: http://www.healthdata.org/mcpa


getwd()

# Linear regression part

# Import the data set
zambia <- read.csv("Zambia.csv", header=T, stringsAsFactors = F, 
                   quote = "\"", sep=",", na.strings=c("null","N/A", "-1"))
zambia

# Get only observations of 2010
zambia.2010 <- zambia[which(zambia$Year == 2010),]

# Remove the national part
zambia.2010 <- zambia.2010[-73,]


unique(zambia.2010$Province)

# Separate training data
zambia.2010.training <- zambia.2010[1:60,]

# data for prediction
zambia.predict <- zambia.2010[61:72,]


# Exclude Household_size, IRS, Pentavalent for setting up full model due to the missing value
# that will cause another missing values for prediction
zambia.lm <- lm(U5_mortality ~ Antenatal1 + Antenatal4 + Schooling + Schooling_women +
                  BCG + Underweight + DPT + Breastfeeding + Female_households +
                  Electricity + Fuel + Sanitation + Wall + IPT + ITN_owner + ITN_owner_or_IRS + 
                  ITN_use + ITN_use_or_IRS + Measles + In_school + In_school_women + 
                  Polio + Skilled_birth, data=zambia.2010.training)

# P-values are very large. does not seem appropriate
summary(zambia.lm)


# Check the multicollinearity
library(car)
vif(zambia.lm)
# multicollinearity exist in the full model

# Perform model selection by using stepwise
library(MASS)
stepAIC(zambia.lm, direction="both",k=2)  
# Antenatal1, Antenatal4, Schooling, BCG, Underweight, Breastfeeding, Female_households, Electricity, ITN_owner, ITN_use_or_IRS, Skilled_birth
stepAIC(zambia.lm, direction="both",k=log(60))
# Antenatal1, Schooling, Underweight, Breastfeeding, Female_households, Electricity


new.zambia1.lm <- lm(U5_mortality ~ Antenatal1 + Schooling + BCG + Underweight + Breastfeeding + Female_households + Electricity + ITN_owner + ITN_use_or_IRS + Skilled_birth, data = zambia.2010.training)
new.zambia2.lm <- lm(U5_mortality ~ Antenatal1 + Schooling + Underweight + Breastfeeding + Female_households + Electricity, data = zambia.2010.training)

# Confirm the new models are better
summary(new.zambia1.lm)  # higher adj R^2 = 0.7079 ; Skilled_birth is not significant and neither is ITN_use_or_IRS
summary(new.zambia2.lm)  # P-values show that the variables are significant; adj R^2 = 0.6709 (removes Skilled_birth, BCG, ITN_owner, ITN_use_or_IRS)

# Remove Skilled_birth and ITN_use_or_IRS from new.zambia1.lm to improve the first model
new.zambia3.lm <- lm(U5_mortality ~ Antenatal1 + Schooling + BCG + Underweight + Breastfeeding + Female_households + Electricity + ITN_owner + ITN_use_or_IRS, data = zambia.2010.training)
summary(new.zambia3.lm)  # adj R^2 = 0.7042; all p-values are significant except for ITN_use_or_IRS whose p-value is > 0.05 but < 0.1

# Durbin-Watson test
dwtest(new.zambia3.lm)  # stepAIC with Skilled_birth and ITN_use_or_IRS removed
dwtest(new.zambia2.lm)  # stepBIC
# DW value is closed to 2 for new.zambia2.lm which means the variables are fairly independent. Also p-value is 0.1613 which
# is greater than 0.05 thus, we can conclude that correlation doesn't seem exist

# Check outliers
plot(new.zambia3.lm)
plot(new.zambia2.lm)

# Perform partial F-test to confirm the new model is better
anova(zambia.lm, new.zambia3.lm)  # better
anova(zambia.lm,new.zambia2.lm)  # better

# Predict the motality rate
predict <- predict(new.zambia2.lm, zambia.predict, interval="prediction")
predict
write.csv(predict, "mortalityrate_prediction.csv")



#####################################################################################

# Time series part

# Get only observations of 2010
zambia.national <- zambia[which(zambia$Level == "National"),]

zambia.national

# Before we fitting ARIMA model, we need to check our time series is stationary
# By definition of time series, Time series is stationary if its mean level and variance stay steady over time

install.packages("fpp")
install.packages("forecast")

library(fpp)
library(forecast)

attach(zambia.national)

# Drop the variables that include missing values
zambia.national$Household_size <- NULL
zambia.national$IRS <- NULL
zambia.national$Pentavalent<- NULL

# create the time series data
zambia.ts.motal <- ts(U5_mortality, start=1990, end=2010)
zambia.ts.motal

plot(zambia.ts.motal, ylab = "Under 5 Mortality Rate", xlab = "Time (1990-2010)", main = "Time Series Plot \n for Under 5 Mortality Rate")


# Check autocorrelation function
Acf(zambia.ts.motal, main = "Autocorrelation in Zambia Time Series \n for Under 5 Mortality Rate")
# each 1 lag means it shows the autocorrelation between time Yt and Yt-1
# But here, we got 1 at the time lag1 means the autocorrelation between Y1990, and Y1990 is not suprisingly 1.

# Check partial autocorrelation function
# If we assume we have Y = b0 + b1*x1 +b2*x2 in this case we check autocorrelation while we remove another time
# effect(X2), so we only compare X1 vs X2, X2 vs X3 and so on.
Pacf(zambia.ts.motal, "Partial Autocorrelation in Zambia \n Time Series for Under 5 Mortality Rate")

Box.test(zambia.ts.motal, lag=20, type="Ljung-Box")
adf.test(zambia.ts.motal, alternative="stationary")
kpss.test(zambia.ts.motal)

# we can use below to fit the ARIMA model
# p is the number of autoregressive terms,
# d is the number of nonseasonal differences, and
# q is the number of lagged forecast errors in the prediction equation.
fit <- arima(zambia.ts.motal, order=c(p, d, q)
             
             # However, this helps us to find appropriate ARIMA model automatically similar as automatic 
             # reduced variable, step-wise function for regression model 
             # Automated forecasting using an ARIMA model
             fit <- auto.arima(zambia.ts.motal)
             fit
             
             library(forecast)
             
             # Forecast the next 20 years
             forecast(fit, 20)
             write.csv(forecast(fit,20), "forecasted.csv")
             
             # As we expected, the motality rate keeps decreasing
             plot(forecast(fit, 20), xlab = "Year", ylab = "Forecasted Under-5 Mortality Rate")
             
             
             # ARIMA (2 2 1)
             
             # AR(p)
             # AR(2): Yt is explained by its past data
             # so we got p=2 which means the any current Yt is explained by
             # Yt-1 and Yt-2 (Y2005 <- affected by Y2004, Y2003 mostly)
             # Yt = alpha1*Yt-1 + alpha2*Yt-2
             
             # MA(q)
             # MA(2): the time series data is affected by error terms of past data
             # so we got q=2 which means any current Yt is explained by
             # Yt = e(t) - beta1*e(t-1) - beta2*e(t-2)
             
             
             # ARIMA means mixed version of AR and MA
             
             # In order to fit ARIMA model, the data must be stationary first.
             # If the model is not stationary, any result of time series analysis will be distorted.
             # Thus, if it is non-stationary we take log or difference (ex. Yt-Yt-1)
             # For example, once we have seasonal, monthly, or yearly variation,
             # then we should take log to the variable
             
             
             
             ######################## the multivariate time series model(failed)
             
             
             
             
