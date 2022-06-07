# Questions ----
# Have PPV, can we use that instead of sp and se? 
# How to pick distribution for se and sp given mean and 95% confidence intervals?

# Final Project QBA ----
library(purrr)
library(dplyr)
library(naniar)
library(reshape2)
library(broom)
library(trapezoid)

j = 100  # number of iterations to run "observed" data through

a = 22  # corresponds to the a cell in a 2x2 contingency table
b = 15  # corresponds to the b cell in a 2x2 contingency table
c = 23  # corresponds to the c cell in a 2x2 contingency table
d= 55  # corresponds to the d cell in a 2x2 contingency table
n = sum(a,b,c,d)  # total number of people in dataset

# Create your simulated data frame (one row per observations per iteration = j*n rows)
misc <- data.frame(Obs = rep(1:j, n), iteration=rep(1:j, each=n), c = rep(as.numeric(rbernoulli(1)), n*j), e_obs = rep(c(rep(1,a+c), rep(0,b+d)),j),
                   d = rep(c(rep(1, a), rep(0, c), rep(1, b), rep(0,d)),j))

# Freq counts that mirror 2x2 table (1 record per combination of e and d)
test <- as.data.frame(table(misc[misc$iteration == 1,]$d, misc[misc$iteration == 1,]$e_obs))
colnames(test) <- c("d", "e_obs", "Freq")

# 1 dataset with all four combinations (1 row, 4 columns)
test2 <- data.frame(ed_o = test[test$d == 1 & test$e_obs == 1,]$Freq,
                    eud_o = test[test$d == 0 & test$e_obs == 1,]$Freq,
                    ued_o = test[test$d == 1 & test$e_obs == 0,]$Freq,
                    ueud_o = test[test$d == 0 & test$e_obs == 0,]$Freq)

# Calculate se and sp, and expected truths for each iteration of data
k = j     # number of iterations to run misclassification process (always same as j)

test3 <- data.frame(iteration = 1:k, 
                    ed_o = rep(test2$ed_o, k), 
                    eud_o = rep(test2$eud_o, k),
                    ued_o = rep(test2$ued_o, k), 
                    ueud_o = rep(test2$ueud_o, k))

# Data from validation paper (NY NDEM)----

test3$se_D <- rbeta(k, 85, 3)
test3$sp_D <- rbeta(k, 95, 7)

test3$se_UD <- rbeta(k, 95, 5)
test3$sp_UD <- rbeta(k, 90, 10)

# Calculate expected "truth" based on se and sp----
test3$ed_t <- (test3$ed_o-(1-test3$sp_D)*(test3$ed_o+test3$ued_o))/(test3$se_D-(1-test3$sp_D))
test3$ued_t <-(test3$ed_o+test3$ued_o)-test3$ed_t
test3$eud_t <-(test3$eud_o-(1-test3$sp_UD)*(test3$eud_o+test3$ueud_o))/(test3$se_UD-(1-test3$sp_UD))
test3$ueud_t <-(test3$eud_o+test3$ueud_o)-test3$eud_t

#expected misclassification adjusted----
test3$syst_exp_adj <- (test3$ed_t /test3$eud_t) / (test3$ued_t /test3$ueud_t)

# sample prevalence
test3$prev_e_d <- rbeta(k, test3$ed_t, test3$eud_t)
test3$prev_e_ud <- rbeta(k, test3$eud_t, test3$ueud_t)

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
misc2 <- merge(select(test3, -c(ed_t, ued_t, eud_t, ueud_t)), misc, by="iteration")

# Calculate expected cells based off PPV and NPV
misc2 <- misc2[order(-misc2$d, -misc2$e_obs),]

m = nrow(misc2)

misc2$e_exp_t_PPVd <- as.numeric(rbernoulli(m, misc2$PPV_d))
misc2$e_exp_t_NPVd <- as.numeric(rbernoulli(m, 1-misc2$NPV_d))
misc2$e_exp_t_PPVud <- as.numeric(rbernoulli(m, misc2$PPV_ud))
misc2$e_exp_t_NPVud <- as.numeric(rbernoulli(m, 1-misc2$NPV_ud))

misc2$e_exp_t <- ifelse(misc2$e_obs == 1 & misc2$d == 1, misc2$e_exp_t_PPVd, 
                        ifelse(misc2$e_obs == 0 & misc2$d == 1, misc2$e_exp_t_NPVd, 
                               ifelse(misc2$e_obs == 1 & misc2$d == 0, misc2$e_exp_t_PPVud, 
                                      ifelse(misc2$e_obs == 0 & misc2$d == 0, misc2$e_exp_t_NPVud, NA))))
misc2$d <- as.factor(misc2$d)

# Crude RR ----
output_obs <- glm(d~e_obs, data=misc2[misc2$iteration ==1,], family=binomial(link='logit')) # Just do this for iteration 1

# Save the output
estim <- data.frame(Parameter="e_exp_t", e_conv = output_obs$coefficients[2],
                     conv_syst = summary(output_obs)$coefficients[2,2])

# Obtain RR estimates for each iteration----
library(tidyr)
estim2 <- misc2 %>% 
  group_by(iteration) %>% 
  do(model = tidy(glm(d~e_exp_t, data = ., family = binomial(link = "logit")))[2,2:3]) %>% 
  unnest(cols = model) 
colnames(estim2) <- c("iteration", "e_syst", "stderr_syst")
estim2$Parameter <- rep("e_exp_t", k)

# Dataset with only se, sp, and syst_exp_adj 
estim3 <- distinct(select(misc2, c(iteration, syst_exp_adj, se_D, sp_D, se_UD, sp_UD)))

# merge all data together to get information for each iteration
estim4 <- merge(merge(estim, estim2, by="Parameter"), estim3, by="iteration")

# Calculate OR intervals for random, systematic, and total error 
estim4$or_rand <- exp(estim4$e_conv - estim4$conv_syst*rnorm(j))
estim4$or_syst <- estim4$syst_exp_adj
estim4$or_tot <- exp(estim4$e_syst - rnorm(j)*estim4$stderr_syst)

# Create a neat output dataset with median as well as 2.5th and 75th percentiles 
output_data <- data.frame(Analysis = c("Random Error Crude","Systematic Error", "Total Error"), 
                          Estimate = c(quantile(estim4$or_rand, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_syst, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_tot, 0.5, na.rm = TRUE)),
                          Percentile_lower = c(quantile(estim4$or_rand, 0.025, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.025, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.025, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(estim4$or_rand, 0.975, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.975, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.975, na.rm = TRUE)))

output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower
output_data

write.csv(output_data, "Results/Individual/NonDiffMis.csv")

# View Output
hist(estim4$se_D, breaks = 30)
hist(estim4$sp_D, breaks = 30)
hist(estim4$se_UD, breaks = 30)
hist(estim4$sp_UD, breaks = 30)
hist(estim4$or_syst, breaks = 30)
hist(estim4$or_tot, breaks = 30)
