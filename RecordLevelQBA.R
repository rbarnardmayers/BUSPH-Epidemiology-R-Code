# Set up ----
library(purrr)
library(dplyr)
library(naniar)
library(reshape2)
library(broom)
library(table1)
library(mediation)
library(medflex)

setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Final Project/")

ltm <- read.csv("LTM_QBA.csv")
# ltm <- read.csv("LTM_QBA_SA12.csv")

# QBA parameters ----
j = 10000  # number of iterations to run "observed" data through
k = j     # number of iterations to run misclassification process (always same as j)

misc <- ltm %>% slice(rep(1:n(), each = j))
misc$iteration <- rep(1:j, nrow(ltm))

# internal validation  ----

se = 33/(33 + 21) # 0.6111111
sp = 0.8 # 649/(649+193)  # 0.7707838

# Freq counts that mirror 2x2 table (1 record per combination of e and d)---
test <- as.data.frame(table(misc[misc$iteration == 1,]$d, misc[misc$iteration == 1,]$e))
colnames(test) <- c("d", "e_obs", "Freq")

# 1 dataset with all four combinations (1 row, 4 columns)
test2 <- data.frame(ed_o = test[test$d == 1 & test$e_obs == 1,]$Freq,
                    eud_o = test[test$d == 0 & test$e_obs == 1,]$Freq,
                    ued_o = test[test$d == 1 & test$e_obs == 0,]$Freq,
                    ueud_o = test[test$d == 0 & test$e_obs == 0,]$Freq)

# Calculate se and sp, and expected truths for each iteration of data

test3 <- data.frame(iteration = 1:k, 
                    ed_o = rep(test2$ed_o, k), 
                    eud_o = rep(test2$eud_o, k),
                    ued_o = rep(test2$ued_o, k), 
                    ueud_o = rep(test2$ueud_o, k))

# se is 603*0.884 = alpha (TPs) and 603-(603*0.884) = beta (FNs)

# Classification Parameters  ----

set.seed(13548); test3$se_UD <- rbeta(k, se*226, 226-(226*se)) 
set.seed(8463132);test3$sp_UD <- rbeta(k, sp*670, 670-(670*sp)) 

set.seed(82309); test3$se_D <- rbeta(k, se*226, 226-(226*se)) 
set.seed(22934); test3$sp_D <- rbeta(k, sp*670, 670-(670*sp)) 

# Calculate expected "truth" based on se and sp----
test3$ed_t <- (test3$ed_o-(1-test3$sp_D)*(test3$ed_o+test3$ued_o))/(test3$se_D-(1-test3$sp_D))
test3$ued_t <-(test3$ed_o+test3$ued_o)-test3$ed_t
test3$eud_t <-(test3$eud_o-(1-test3$sp_UD)*(test3$eud_o+test3$ueud_o))/(test3$se_UD-(1-test3$sp_UD))
test3$ueud_t <-(test3$eud_o+test3$ueud_o)-test3$eud_t

#expected misclassification adjusted----
test3$syst_exp_adj <- (test3$ed_t /test3$eud_t) / (test3$ued_t /test3$ueud_t)

# sample prevalence
set.seed(64521); test3$prev_e_d <- rbeta(k, test3$ed_t, test3$eud_t)
set.seed(894352132); test3$prev_e_ud <- rbeta(k, test3$eud_t, test3$ueud_t)

# Calculate PPV and NPV ----
# want to make negative values NA for PPV and NPV
test3$PPV_d <- ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                      (test3$se_D* test3$prev_e_d)/((test3$se_D* test3$prev_e_d)+((1-test3$sp_D)*(1-test3$prev_e_d))))
test3$NPV_d <-  ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                       (test3$sp_D * (1-test3$prev_e_d))/(((1-test3$se_D) * test3$prev_e_d)+(test3$sp_D*(1-test3$prev_e_d))))
test3$PPV_ud <- ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                       (test3$se_UD * test3$prev_e_ud)/((test3$se_UD* test3$prev_e_ud)+((1-test3$sp_UD)*(1-test3$prev_e_ud))))
test3$NPV_ud <- ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                       (test3$sp_UD * (1-test3$prev_e_ud))/(((1-test3$se_UD) * test3$prev_e_ud)+(test3$sp_UD*(1-test3$prev_e_ud)))) 

# Create new dataset without "true" cells 
misc2 <- merge(test3, misc, by="iteration")

m = nrow(misc2)

set.seed(8465); misc2$e_exp_t_PPVd <- as.numeric(rbernoulli(m, misc2$PPV_d))
set.seed(4357); misc2$e_exp_t_NPVd <- as.numeric(rbernoulli(m, 1-misc2$NPV_d))
set.seed(6467); misc2$e_exp_t_PPVud <- as.numeric(rbernoulli(m, misc2$PPV_ud))
set.seed(684312); misc2$e_exp_t_NPVud <- as.numeric(rbernoulli(m, 1-misc2$NPV_ud))

misc2$e_exp_t <- ifelse(misc2$e == 1 & misc2$d == 1, misc2$e_exp_t_PPVd, 
                        ifelse(misc2$e == 0 & misc2$d == 1, misc2$e_exp_t_NPVd, 
                               ifelse(misc2$e == 1 & misc2$d == 0, misc2$e_exp_t_PPVud, 
                                      ifelse(misc2$e == 0 & misc2$d == 0, misc2$e_exp_t_NPVud, NA))))
misc2$d <- as.factor(misc2$d)
# misc2 <- misc2 %>% select(-c(misc2$e_exp_t_NPVd, misc2$e_exp_t_NPVud, misc2$e_exp_t_PPVd, misc2$e_exp_t_PPVud))

# Crude RR ----
output_obs_c <- glm(d~e, data=misc2[misc2$iteration ==1,], family=binomial(link='logit')) # Just do this for iteration 1

# Save the output
estim0 <- data.frame(Parameter="e_exp_t", e_conv_c = output_obs_c$coefficients[2],
                     conv_syst_c = summary(output_obs_c)$coefficients[2,2])
# Adjusted RR ----
output_obs_a <- glm(d~e + GestDiab + Age + PregnancyWeight_Gain + Prepregnancy_Weight, 
                    data=misc2[misc2$iteration ==1,], 
                    family=binomial(link='logit')) # Just do this for iteration 1

# Save the output
estim1 <- data.frame(Parameter="e_exp_t", e_conv_a = output_obs_a$coefficients[2],
                     conv_syst_a = summary(output_obs_a)$coefficients[2,2])

# Obtain RR estimates for each iteration----
library(tidyr)
estim2 <- misc2 %>% 
  group_by(iteration) %>% 
  do(model = tidy(glm(d~e_exp_t + GestDiab + Age + PregnancyWeight_Gain + 
                        Prepregnancy_Weight , data = ., 
                      family = binomial(link = "logit")))[2,2:3]) %>% 
  unnest(cols = model) 
colnames(estim2) <- c("iteration", "e_syst", "stderr_syst")
estim2$Parameter <- rep("e_exp_t", k)

# Dataset with only se, sp, and syst_exp_adj 
estim3 <- distinct(misc2[, names(misc2) %in% c('iteration', 
                                               'syst_exp_adj', 'se_D', 'sp_D', 'se_UD', 'sp_UD')])

# merge all data together to get information for each iteration
estim4 <- merge(merge(merge(estim0, estim1, by="Parameter"), estim2, by="Parameter"), estim3, by="iteration")

# Calculate OR intervals for random, systematic, and total error 
set.seed(5684); 
estim4$or_rand_c <- exp(estim4$e_conv_c - estim4$conv_syst_c*rnorm(j))
set.seed(215456846); 
estim4$or_rand_a <- exp(estim4$e_conv_a - estim4$conv_syst_a*rnorm(j))
estim4$or_syst <- estim4$syst_exp_adj
set.seed(586); 
estim4$or_tot <- exp(estim4$e_syst - rnorm(j)*estim4$stderr_syst)

# Create a neat output dataset with median as well as 2.5th and 97.5th percentiles 
output_data <- data.frame(Analysis = c("Random Error Crude","Random Error Adjusted" , "Systematic Error", "Total Error"), 
                          Estimate = c(quantile(estim4$or_rand_c, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_rand_a, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_syst, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_tot, 0.5, na.rm = TRUE)),
                          Percentile_lower = c(quantile(estim4$or_rand_c, 0.025, na.rm = TRUE),
                                               quantile(estim4$or_rand_a, 0.025, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.025, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.025, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(estim4$or_rand_c, 0.975, na.rm = TRUE),
                                               quantile(estim4$or_rand_a, 0.975, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.975, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.975, na.rm = TRUE)))

output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower
output_data

write.csv(output_data, "QBA_TOLDLARGE.csv")
# write.csv(output_data, "QBA_TOLDLARGE_SA12.csv")

# View Output
hist(estim4$se_D, breaks = 30)
hist(estim4$sp_D, breaks = 30)
hist(estim4$se_UD, breaks = 30)
hist(estim4$sp_UD, breaks = 30)
hist(estim4$or_syst, breaks = 30)
hist(estim4$or_tot, breaks = 30)
