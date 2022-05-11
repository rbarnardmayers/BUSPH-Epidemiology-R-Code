# Set up ----
library(purrr)
library(dplyr)
library(naniar)
library(reshape2)
library(broom)
library(table1)
library(medflex)

setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Final Project/")

ltm <- read.csv("LTM_FinalData.csv")

# Subsetting ----
# Only look at nulliparous Women without diabetes
ltm <- ltm %>% subset(Q610_1 == 1 & Diabetes == 2 & !is.na(PregnancyWeight_Gain)  & !is.na(Race_Ethnicity))

# Mediation Analysis (imputation approach) ----

# setting the mediator model and creating and expanded counterfactual dataset ---- 
impData <- neImpute(Cesarean ~ GestDiab + Babyweight_Ounces + PregnancyWeight_Gain + Age + Prepregnancy_Weight + Race_Ethnicity,
              family = binomial("logit"), nMed = 2, data = ltm)

# Fitting the natural effects model ----
# need to use the expanded dataset and the counterfactual exposure variables
neMod1 <- neModel(Cesarean ~ GestDiab0 + GestDiab1 +  Age + Prepregnancy_Weight + Race_Ethnicity,
                  family = binomial("logit"), expData = impData)

# Exponentiating the model parameter estimates provides estimates that 
# can be interpreted as odds ratios.
summary(neMod1)

# Joint Mediation 
# NDE = 1.052624 using Told Large 95% CI:  0.648 1.721
# NIE = 1.091257 using Told Large 95% CI:  0.954 1.261

# NDE = 1.170689 using babyweight  95% CI: 0.7633 1.883
# NIE = 0.9824385 using babyweight  95% CI: 0.918 1.047

impData1 <- neImpute(Cesarean ~ GestDiab + ToldLarge + PregnancyWeight_Gain + Age + Prepregnancy_Weight + Race_Ethnicity,
                    family = binomial("logit"), nMed = 2, data = ltm)

# Fitting the natural effects model ----
# need to use the expanded dataset and the counterfactual exposure variables
neMod2 <- neModel(Cesarean ~ GestDiab0 + GestDiab1 +  Age + Prepregnancy_Weight + Race_Ethnicity,
                  family = binomial("logit"), expData = impData1)



