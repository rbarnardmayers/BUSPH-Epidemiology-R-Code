# Final Project 
# Set up ----
library(purrr)
library(dplyr)
library(naniar)
library(reshape2)
library(broom)
library(table1)
library(medflex)
library(invgamma)

setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Final Project/")

ltm <- read.csv("LTM_FinalData.csv")
ltm_m <- ltm[, c('Cesarean', 'GestDiab','Babyweight_Ounces', 'Induce', 'WasLarge' , 'Age', 'Prepregnancy_Weight', "PregnancyWeight_Gain",'Race_Ethnicity')]
ltm_m$Cesarean <- as.factor(ltm_m$Cesarean)
ltm_m$GestDiab <- as.factor(ltm_m$GestDiab)


# Function ----
getCI.int <- function(model, n1, n2){
  coefs <- summary(model)$coefficients[,1]
  X <- model.matrix(model)
  dof <- nrow(X) - ncol(X)
  coefs_var <- vcov(model)
  halfCI <- qt(0.975, dof) * sqrt(coefs_var[n1,n1]+coefs_var[n2,n2]+2*coefs_var[n1,n2])
  return(as.vector(c(coefs[n1]+coefs[n2]-halfCI, coefs[n1]+coefs[n2]+halfCI)))
}

get.SE <- function(list){
  lower <- list[1]
  upper <- list[2]
  se <- (upper - lower)/3.92
  return(se)
}

# Table 1 ----
table1(~ Age_f + Race_Ethnicity_f + Insurance_f  + 
         Income_f + Cesarean_f + Babyweight_Ounces + 
         Babyweight_Categorical_f +  WasLarge_f + 
         ToldLarge_f + Induce_f +
         PregnancyWeight_Gain + Prepregnancy_Weight |
         GestDiab_f, data=ltm)

# Conventional Analysis ----
ltm$Cesarean <- as.factor(ltm$Cesarean)
ltm$WasLarge <- as.factor(ltm$WasLarge)
ltm$ToldLarge <- as.factor(ltm$ToldLarge)

GC.crude.out <- glm(Cesarean ~ GestDiab, family = binomial(link = "logit"), 
                 data = ltm)

GC.adj.out <- glm(Cesarean ~ GestDiab + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain, 
                  family = binomial(link = "logit"),
                 data = ltm) 

GM.crude.out <- lm(Babyweight_Ounces ~ GestDiab, 
                   data = ltm)

GM.adj.out <- lm(Babyweight_Ounces ~ GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight,
                 data = ltm) 

MC.crude.out <- glm(Cesarean ~ Babyweight_Ounces, family = binomial(link = "logit"), 
                           data = ltm)

MC.adj.out <- glm(Cesarean ~ Babyweight_Ounces * GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight, 
                        family = binomial(link = "logit"),
                         data = ltm) 

# Prepare the columns
OR_1 <- round(exp(coef(GC.crude.out)), 2)
OR_2 <- round(exp(coef(GC.adj.out)), 2)
OR_3 <- round(coef(GM.crude.out), 2)
OR_4 <- round(coef(GM.adj.out), 2)
OR_5 <- round(exp(coef(MC.crude.out)), 2)
OR_6 <- round(exp(coef(MC.adj.out)), 2)
OR <- c(OR_1,OR_2,OR_3,OR_4,OR_5,OR_6)

CI_1 <- round(exp(confint(GC.crude.out)), 2)
CI_2 <- round(exp(confint(GC.adj.out)), 2)
CI_3 <- round(confint(GM.crude.out), 2)
CI_4 <- round(confint(GM.adj.out), 2)
CI_5 <- round(exp(confint(MC.crude.out)), 2)
CI_6 <- round(exp(confint(MC.adj.out)), 2)

CI <- as.data.frame(rbind(CI_1,CI_2,CI_3,CI_4,CI_5,CI_6))
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")

export <- cbind(OR, CI)
export$Model <- c(rep("GC Crude", length(OR_1)), 
                  rep("GC Adjusted", length(OR_2)), 
                  rep("GM Crude", length(OR_3)), 
                  rep("GM Adjusted", length(OR_4)), 
                  rep("MC Crude", length(OR_5)), 
                  rep("MC Adjusted", length(OR_6)))
write.csv(export, "Regression_Results.csv")

# Mediation Analyses ----
# Mediation Analysis (imputation approach) 
# setting the mediator model and creating and expanded counterfactual dataset 
impData <- neImpute(Cesarean ~ GestDiab * Babyweight_Ounces + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                    family = binomial("logit"), nMed = 1, data = ltm_m)

# Fitting the natural effects model
# need to use the expanded dataset and the counterfactual exposure variables
neMod1 <- neModel(Cesarean ~ GestDiab0 + GestDiab1 + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                  family = binomial("logit"), expData = impData)
summary(neMod1)

Results_Mediation <- data.frame(Analysis = c("NDE", "NIE"), 
                                Estimate = c(exp(summary(neMod1)$coefficients[2,1]), exp(summary(neMod1)$coefficients[3,1])),
                                Lower = c(exp(confint(neMod1)[2,1]), exp(confint(neMod1)[3,1])),
                                Upper = c(exp(confint(neMod1)[2,2]), exp(confint(neMod1)[3,2])), 
                                logOR = c(summary(neMod1)$coefficients[2,1], summary(neMod1)$coefficients[3,1]),
                                selogOR = c(summary(neMod1)$coefficients[2,2], summary(neMod1)$coefficients[3,2]))

write.csv(Results_Mediation, "Results_Mediation.csv")

# Sensitivity Analyses Table 1s ----
ltm <- read.csv("LTM_SA12.csv")
table1(~ Age_f + Race_Ethnicity_f + Insurance_f  + 
         Income_f + Cesarean_f + Babyweight_Ounces + 
         Babyweight_Categorical_f +  WasLarge_f + 
         ToldLarge_f + Induce_f + 
         PregnancyWeight_Gain + Prepregnancy_Weight |
         GestDiab_f, data=ltm)

# Sensitivity Analysis Was Large DONE ----
ltm <- read.csv("LTM_FinalData.csv")

MC2.crude.out <- glm(Cesarean ~ WasLarge, family = binomial(link = "logit"), 
                     data = ltm)

MC2.adj.out <- glm(Cesarean ~ WasLarge + GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain, 
                   family = binomial(link = "logit"),
                   data = ltm) 

GM2.crude.out <- glm(WasLarge ~ GestDiab, family = binomial(link = "logit"),
                     data = ltm)

GM2.adj.out <- glm(WasLarge ~ GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight,
                   family = binomial(link = "logit"), data = ltm) 

# Prepare the columns
OR_1 <- round(exp(coef(GC.crude.out)), 2)
OR_2 <- round(exp(coef(GC.adj.out)), 2)
OR_3 <- round(exp(coef(GM2.crude.out)), 2)
OR_4 <- round(exp(coef(GM2.adj.out)), 2)
OR_5 <- round(exp(coef(MC2.crude.out)), 2)
OR_6 <- round(exp(coef(MC2.adj.out)), 2)
OR <- c(OR_1,OR_2,OR_3,OR_4,OR_5,OR_6)

CI_1 <- round(exp(confint(GC.crude.out)), 2)
CI_2 <- round(exp(confint(GC.adj.out)), 2)
CI_3 <- round(exp(confint(GM2.crude.out)), 2)
CI_4 <- round(exp(confint(GM2.adj.out)), 2)
CI_5 <- round(exp(confint(MC2.crude.out)), 2)
CI_6 <- round(exp(confint(MC2.adj.out)), 2)

CI <- as.data.frame(rbind(CI_1,CI_2,CI_3,CI_4,CI_5,CI_6))
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")

export1 <- cbind(OR, CI)
export1$Model <- c(rep("GC Crude", length(OR_1)), 
                  rep("GC Adjusted", length(OR_2)), 
                  rep("GM Crude", length(OR_3)), 
                  rep("GM Adjusted", length(OR_4)), 
                  rep("MC Crude", length(OR_5)), 
                  rep("MC Adjusted", length(OR_6)))
write.csv(export1, "Regression_Results_WasLarge.csv")

# Mediation Was Large DONE ----
impData2 <- neImpute(Cesarean ~ GestDiab * WasLarge + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                     family = binomial("logit"), nMed = 1, data = ltm_m)

# Fitting the natural effects model
# need to use the expanded dataset and the counterfactual exposure variables
neMod3 <- neModel(Cesarean ~ GestDiab0 + GestDiab1 + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                  family = binomial("logit"), expData = impData2)

Results_Mediation_2 <- data.frame(Analysis = c("SA NDE", "SA NIE"), 
                                  Estimate = c(exp(summary(neMod3)$coefficients[2,1]), 
                                               exp(summary(neMod3)$coefficients[3,1])),
                                  Lower = c(exp(confint(neMod3)[2,1]), exp(confint(neMod3)[3,1])),
                                  Upper = c(exp(confint(neMod3)[2,2]), exp(confint(neMod3)[3,2])), 
                                  logOR = c(summary(neMod3)$coefficients[2,1], 
                                            summary(neMod3)$coefficients[3,1]),
                                  selogOR = c(summary(neMod3)$coefficients[2,2], 
                                              summary(neMod3)$coefficients[3,2]))
write.csv(Results_Mediation_2, "Results_Mediation_WasLarge.csv")


# Sensitivity Analysis Told Large  ----
crude.out <- glm(Cesarean ~ ToldLarge, family = binomial(link = "logit"), 
                     data = ltm)

adj.out <- glm(Cesarean ~ ToldLarge * GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight + Race_Ethnicity, 
               family = binomial(link = "logit"),
                   data = ltm) 

# Prepare the columns
OR_1 <- round(exp(coef(crude.out)), 2)
OR_2 <- round(exp(coef(adj.out)), 2)
OR_3 <- c(exp(adj.out$coefficients[2] + adj.out$coefficients[8]))
OR <- c(OR_1,OR_2, OR_3)

CI_1 <- round(exp(confint(crude.out)), 2)
CI_2 <- round(exp(confint(adj.out)), 2)
CI_3 <- round(exp(getCI.int(adj.out, 2, 8)), 2)

CI <- as.data.frame(rbind(CI_1,CI_2, CI_3))

# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")

export1 <- cbind(OR, CI)
export1$Model <- c(rep("Crude", length(OR_1)), 
                   rep("Adjusted", length(OR_2)), 
                   "interaction")

write.csv(export1, "Regression_Results_ToldLarge.csv")


# Mediation with Induction ----
impData3 <- neImpute(Cesarean ~ GestDiab * Induce + Babyweight_Ounces + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                     family = binomial("logit"), nMed = 2, data = ltm_m)

# Fitting the natural effects model
# need to use the expanded dataset and the counterfactual exposure variables
neMod4 <- neModel(Cesarean ~ GestDiab0 + GestDiab1 + Age + Prepregnancy_Weight + Race_Ethnicity + PregnancyWeight_Gain,
                  family = binomial("logit"), expData = impData3)

Results_Mediation_Induce <- data.frame(Analysis = c("Original NDE", "Original NIE"), 
                                  Estimate = c(exp(summary(neMod4)$coefficients[2,1]), 
                                               exp(summary(neMod4)$coefficients[3,1])),
                                  Lower = c(exp(confint(neMod4)[2,1]), exp(confint(neMod4)[3,1])),
                                  Upper = c(exp(confint(neMod4)[2,2]), exp(confint(neMod4)[3,2])), 
                                  logOR = c(summary(neMod4)$coefficients[2,1], 
                                            summary(neMod4)$coefficients[3,1]),
                                  selogOR = c(summary(neMod4)$coefficients[2,2], 
                                              summary(neMod4)$coefficients[3,2]))
write.csv(Results_Mediation_Induce, "Results_Mediation_Induce.csv")

# Validation ----
table(ltm$ToldLarge_f, ltm$WasLarge) 
