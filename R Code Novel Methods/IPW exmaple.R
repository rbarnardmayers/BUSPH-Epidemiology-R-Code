# Simple data exmple
setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 12 - MSMs II")


library(dplyr)
library(geepack)
library(haven)

data <- read_sas("nhefs.sas7bdat")

# manipulate data
data <- data %>% 
  mutate(match = 1, 
         cens = case_when(is.na(wt82) ~ 1, 
                          TRUE ~ 0))

data1 <- subset(data, !is.na(wt82))

# Explore data ----
log.out1 <- lm(wt82_71 ~ qsmk, data = data)

summary(data[data$qsmk == 1, ]$age)
summary(data[data$qsmk == 0, ]$age)

summary(data[data$qsmk == 1, ]$wt71)
summary(data[data$qsmk == 0, ]$wt71)

summary(data[data$qsmk == 1, ]$smokeintensity)
summary(data[data$qsmk == 0, ]$smokeintensity)

summary(data[data$qsmk == 1, ]$smokeyrs)
summary(data[data$qsmk == 0, ]$smokeyrs)

summary(data[data$qsmk == 1, ]$sex)
summary(data[data$qsmk == 0, ]$sex)

summary(data[data$qsmk == 1, ]$race)
summary(data[data$qsmk == 0, ]$race)

summary(data[data$qsmk == 1, ]$education)
summary(data[data$qsmk == 0, ]$education)


# Treatment Weights (unstabalized)----
# Regression model
log.out2 <- glm(qsmk ~ sex + race + age + age^2 + factor(education) + smokeintensity + 
                  smokeintensity^2 + smokeyrs^2 +
                  smokeyrs + factor(exercise) + factor(active) + wt71 + wt71^2, 
                data = data1, family = binomial(link = "logit"))

# saving the predicted probabilities
data1$p_qsmk <- log.out2$fitted.values

# generating the weights based off the predicted values
# nonstabalized weights, if had exposure, weight is 1 / predicted probability of having exposure
# # if did not have exposure, weight is 1 / predicted probability of NOT having exposure
data1$w <- ifelse(data1$qsmk == 1, 1/data1$p_qsmk, 
                 ifelse(data1$qsmk == 0, 1/(1-data1$p_qsmk), NA))

mean(data1$w) # should be 2, is 1.997331
sd(data1$w) # should be 1.5, is 1.379866
range(data1$w) # range should be 1.05 to 16.7, is  1.047523 to 13.405788

# Computing IPTW association 
log.out3 <- geeglm(wt82_71 ~ qsmk, data=data1, weights = w,
                    id = factor(seqn))

# force out sandwich estimate of the variance
tidy(log.out3) 

# Diagnositcs 
library(survey)

wtd.table(data1$qsmk, weights = data1$w)
wtd.table(data1$sex, weights = data1$w)

table(data1[data1$sex == 1 & data1$race == 0, ]$age, data1[data1$sex == 1 & data1$race == 0, ]$qsmk)


# Treatment Weights (stabalized)----
# Regression model
log.out_den <- glm(qsmk ~ sex + race + age + age^2 + factor(education) + smokeintensity + 
                  smokeintensity^2 + smokeyrs^2 +
                  smokeyrs + factor(exercise) + factor(active) + wt71 + wt71^2, 
                data = data, family = binomial(link = "logit"))

log.out_num <- glm(qsmk ~ 1, 
                data = data, family = binomial(link = "logit"))

# saving the predicted probabilities
data$pn_qsmk <- log.out_num$fitted.values
data$pd_qsmk <- log.out_den$fitted.values

# generating the weights based off the predicted values
# nonstabalized weights, if had exposure, weight is 1 / predicted probability of having exposure
# # if did not have exposure, weight is 1 / predicted probability of NOT having exposure
data$sw <- ifelse(data$qsmk == 1, data$pn_qsmk/data$pd_qsmk, 
                 ifelse(data$qsmk == 0, (1-data$pn_qsmk)/(1-data$pd_qsmk), NA))

mean(data$sw) # 0.999
sd(data$sw) # 0.259
range(data$sw) # 0.336, 3.52

# Computing IPTW association 
log.out6 <- geeglm(wt82_71 ~ qsmk, data=data, weights = sw,
                   id = factor(seqn))

# force out sandwich estimate of the variance
tidy(log.out6)

# Diagnositcs 
library(survey)

wtd.table(data$qsmk, weights = data$sw)
wtd.table(data$sex, weights = data$sw)

table(data[data$sex == 1 & data$race == 0, ]$age, data[data$sex == 1 & data$race == 0, ]$qsmk)

# Continuation weights ----
data2 <- data
data2$wt82_71 <- ifelse(data2$death == 1, NA, data2$wt82_71)
data2$age <- ifelse(data2$death == 1, data2$age - 7, data2$age)
data2$wt71 <- ifelse(data2$death == 1, data2$wt71 + 4, data2$wt71)
data2$cens <- ifelse(data2$death ==1, 1, data2$cens)

prop.table(table(data2$qsmk, data2$cens),1)
# 0.1981682 with qsmk = 0, 0.2710280 with qsmk = 1

mean(data2[data2$cens == 1, ]$wt71) # 76.4274
mean(data2[data2$cens == 0, ]$wt71) # 70.55735

# censoring/continuation model 
library(DescTools) 
# logistic regression for denominator of continuation weights
log.out7 <- glm(cens ~ factor(qsmk) + factor(sex) + factor(race) + age + age^2 + 
                  factor(education) + smokeintensity + smokeintensity^2 + 
                  smokeyrs + smokeyrs^2 + factor(exercise) + factor(active) + wt71 + wt71^2, 
                data = data2, family = binomial(link = "logit"))

Cstat(log.out7) # 0.8383339

# adding predicted values from regression to dataset
data2$pd_cens <- log.out7$fitted.values

# get numerator for continuation weights
log.out8 <- glm(cens ~ qsmk, data = data2, family = binomial(link="logit"))
data2$pn_cens <- log.out8$fitted.values

# need to re run the treatment weights  
log.out_den8 <- glm(qsmk ~ factor(sex) + factor(race) + age + age^2 + factor(education) + smokeintensity + 
                     smokeintensity^2 + smokeyrs^2 +
                     smokeyrs + factor(exercise) + factor(active) + wt71 + wt71^2, 
                   data = data2, family = binomial(link = "logit"))

log.out_num8 <- glm(qsmk ~ 1, 
                   data = data2, family = binomial(link = "logit"))

data2$pd_qsmk <- log.out_den8$fitted.values
data2$pn_qsmk <- log.out_num8$fitted.values

data3 <- subset(data2, cens = 0)

data3$sw_a <- ifelse(data3$qsmk == 1, data3$pn_qsmk / data3$pd_qsmk, 
                     ifelse(data3$qsmk == 0, (1-data3$pn_qsmk) / (1-data3$pd_qsmk), NA))

data3$sw_c <- data3$pn_cens / data3$pd_cens
data3$sw = data3$sw_a * data3$sw_c

mean(data3$sw) # should be 1 and is not
sd(data3$sw) # should be 1 and is not

# Compute IPTCW 

log.out9 <- geeglm(wt82_71 ~ qsmk, data=data3, weights = sw,
                              id = factor(seqn))



