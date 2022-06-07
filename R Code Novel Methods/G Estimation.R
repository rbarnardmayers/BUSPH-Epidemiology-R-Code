library(dplyr)
library(haven)

setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 13 - G Estimation")

nhefs <- read_sas("nhefs.sas7bdat")

nhefs <- nhefs %>% mutate(older = ifelse(age > 50, 1, ifelse(age <= 50, 0, NA)), 
                          qsmkolder = qsmk * older)

log.out1 <- lm(wt82_71 ~ qsmk + older + qsmkolder, data = nhefs)

summary(log.out1)$coefficients
  
# Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)  2.9914484  0.2609718 11.4627257 2.895592e-29
# qsmk         3.0649598  0.5502326  5.5702985 2.989613e-08
# older       -3.7534739  0.5038557 -7.4495017 1.542974e-13
# qsmkolder   -0.2024598  0.9273881 -0.2183118 8.272147e-01

# E[Y | A=0, L = 0] = 
summary(log.out1)$coefficients[1,1]
# E[Y | A=0, L = 1] = 
summary(log.out1)$coefficients[1,1] + summary(log.out1)$coefficients[3,1]
# E[Y | A=1, L = 0] = 
summary(log.out1)$coefficients[1,1] + summary(log.out1)$coefficients[2,1]
# E[Y | A=1, L = 1] = 
summary(log.out1)$coefficients[1,1] + summary(log.out1)$coefficients[2,1] + 
summary(log.out1)$coefficients[3,1] + summary(log.out1)$coefficients[4,1]










