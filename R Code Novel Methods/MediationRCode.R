setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 5 - Mediation II")
library(haven)
library(mediation)
library(medflex)

allvitd2 <- UPBdata

# setting the mediator model and creating and expanded counterfactual dataset ---- 
medFit <- glm(negaff ~ factor(attbin) + gender + educ + age,
              family = gaussian, data = UPBdata)
expData <- neWeight(medFit)

# Fitting the natural effects model ----
# need to use the expanded dataset and the counterfactual exposure variables
neMod1 <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age,
                  family = binomial("logit"), expData = expData)

# Exponentiating the model parameter estimates provides estimates that 
# can be interpreted as odds ratios.
summary(neMod1)

# beta 1 is NDE
# beta 2 is NIE

# IMPUTATION METHOD
impFit <- glm(UPB ~ factor(attbin) + negaff + gender + educ + age,
              family = binomial("logit"), data = UPBdata)

expData <- neImpute(impFit)

neMod1 <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age,
                  family = binomial("logit"), expData = expData, se = "robust")
# beta 1 is NDE
# beta 2 is NIE


# In class ----
allvitd2 <- read.csv("allvitd2.csv")
#allvitd2 <- read_sas("~/Downloads/allvitd2.sas7bdat", NULL)
#write.csv(allvitd2, "~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 5/allvitd2.csv")

# Predicting Outcome from Exposure and Mediator (Total Effect)
out.fit <- glm(lethal10 ~ qv25*chronbin + agebld_yr + bmi94 + gleas_cat, data = allvitd2)

# Predicting Mediator from Exposure
med.fit <- glm(chronbin ~ qv25 + agebld_yr + bmi94 + gleas_cat, data = allvitd2)

# Mediation analysis
med.out <- mediate(med.fit, out.fit, treat="qv25", mediator="chronbin")
summary(med.out)$d1
summary(med.out)$z1

# ACME = indirect
# ADE = direct
# Total Effect = total effect

# HW 2 ----
data1 <- read_dta(file.choose())
data1$educ_L <- ifelse(data1$educ == 1, 1, 0)
data1$educ_M <- ifelse(data1$educ == 2, 1, 0)
# Predicting Outcome from Exposure and Mediator (Total Effect)
out.fit <- glm(UPB ~ attbin + negaff + gender + educ_L + educ_M + age, data = data1)

# Predicting Mediator from Exposure
med.fit <- lm(negaff ~ attbin, data = data1)

# Mediation analysis
med.out <- mediate(med.fit, out.fit, treat="attbin", mediator="negaff")
summary(med.out)
summary(med.out)$z1