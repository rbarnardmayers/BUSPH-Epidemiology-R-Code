###### R CODE FOR RECORD LEVEL BIAS ANALYSIS: EXPOSURE MISCLASSIFICATION #############
# Created by: 		Ruby Barnard-Mayers
# Date created: 	Jan 12, 2022
# Purpose: 			Simple R code for individual record level PBA for exposure misclassification
#########################################################################################

library(purrr)
library(dplyr)
j = 100  # number of iterations to run "observed" data through
k = j     # number of iterations to run misclassification process (always same as j)

# Data information ----
# IF YOU HAVE INDIVIDUAL LEVEL DATA ALREADY, uncomment the lines below and skip lines 15-27
# data <- read.csv(file.choose) # upload as csv
# data <- rename(data, e_obs = exposure, d = outcome)  # rename exposure and outcome columns to e_obs and d
# misc <- rep(cbind(data),j)

# if you need to simulate record level data from a 2x2 table

a = 215  # corresponds to the a cell in a 2x2 contingency table
b = 1449  # corresponds to the b cell in a 2x2 contingency table
c = 668  # corresponds to the c cell in a 2x2 contingency table
d= 4296  # corresponds to the d cell in a 2x2 contingency table
n = sum(a,b,c,d)  # total number of people in dataset

# Create your simulated data frame (one row per observations per iteration = j*n rows)
misc <- data.frame(Obs = rep(1:j, n), iteration=rep(1:j, each=n), c = rep(as.numeric(rbernoulli(1)), n*j), e_obs = rep(c(rep(1,a+c), rep(0,b+d)),j),
                        d = rep(c(rep(1, a), rep(0, c), rep(1, b), rep(0,d)),j))

# Freq counts that mirror 2x2 table (1 record per combination of e and d)
test <- as.data.frame(table(misc[misc$iteration == 1,]$d, misc[misc$iteration == 1,]$e_obs))
test <- rename(test, d = Var1)
test <- rename(test, e_obs = Var2)

# 1 dataset with all four combinations (1 row, 4 columns)
test2 <- data.frame(ed_o = test[test$d == 1 & test$e_obs == 1,]$Freq,
                    eud_o = test[test$d == 0 & test$e_obs == 1,]$Freq,
                    ued_o = test[test$d == 1 & test$e_obs == 0,]$Freq,
                    ueud_o = test[test$d == 0 & test$e_obs == 0,]$Freq)

# Calculate se and sp, and expected truths for each iteration of data


test3 <- data.frame(iteration = 1:k, ed_o = rep(test2$ed_o, k), eud_o = rep(test2$eud_o, k), ued_o = rep(test2$ued_o, k), ueud_o = rep(test2$ueud_o, k))

# Correlated data ----
# This is for differential
# For non-differential, run test3$se_UD <- test3$se_D and test3$sp_UD <- test3$sp_D 
# Sensitivity
# test3$rho <- 0.8  # Correlation coefficient
# test3$Z1 <- rnorm(1)
# test3$U1 <- pnorm(test3$Z1)
# test3$Z2 <- rnorm(1)
# test3$U2 <- pnorm(test3$Z2)
# test3$se_D <- pbeta(test3$U1, 50.6, 14.3)  # beta distribution (centered on 50.6/(50.6+14.3))
# test3$se_UD <- pbeta(test3$U2, 50.6, 14.3)  # beta distribution (centered on 50.6/(50.6+14.3)) 

# Specificity
# test3$Z3 <- rnorm(1)
# test3$U3 <- pnorm(test3$Z3)
# test3$Z4 <- rnorm(1)
# test3$U4 <- pnorm(test3$Z4)
# test3$sp_D <- pbeta(test3$U3, 70, 1)  # beta distribution  (centered on 70/(70+1)) 
# test3$sp_UD <- pbeta(test3$U4, 70, 1)  # beta distribution (centered on 70/(70+1)) 
# ----

# This is for non-differential non-correlated. For differential, change lines 70 and 71 
# These data would come from a valdiation dataset 
test3$se_D <- rbeta(k, 50.6,14.3) #  (centered on 0.8/(0.8+0.1))
test3$sp_D <- rbeta(k, 70,1) #  (centered on 32/(32+8)) 

test3$se_UD <- test3$se_D
test3$sp_UD <- test3$sp_D
# test3$se_UD <- rbeta(k, 95,5)
# test3$sp_UD <- rbeta(k, 90,10)

# Calculate expected "truth" based on se and sp
test3$ed_t <- (test3$ed_o-(1-test3$sp_D)*(test3$ed_o+test3$ued_o))/(test3$se_D-(1-test3$sp_D))
test3$ued_t <-(test3$ed_o+test3$ued_o)-test3$ed_t
test3$eud_t <-(test3$eud_o-(1-test3$sp_UD)*(test3$eud_o+test3$ueud_o))/(test3$se_UD-(1-test3$sp_UD))
test3$ueud_t <-(test3$eud_o+test3$ueud_o)-test3$eud_t

#expected misclassification adjusted
test3$syst_exp_adj <- (test3$ed_t /test3$eud_t) / (test3$ued_t /test3$ueud_t)

# sample prevalence
test3$prev_e_d <- ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                         rbeta(k, test3$ed_t, test3$eud_t))
test3$prev_e_ud <- ifelse(test3$ed_t <= 0 | test3$eud_t <= 0 | test3$ued_t <= 0 | test3$ueud_t <= 0, NA,
                         rbeta(k, test3$eud_t, test3$ueud_t))

# Calculate PPV and NPV
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
misc2$e_exp_t <- ifelse(misc2$e_obs == 1 & misc2$d == 1, rbernoulli(misc2$PPV_d, p=misc2$PPV_d), 
                        ifelse(misc2$e_obs == 0 & misc2$d == 1, rbernoulli(1-misc2$NPV_d,p=1-misc2$NPV_d), 
                              ifelse(misc2$e_obs == 1 & misc2$d == 0, rbernoulli(misc2$PPV_ud, p=misc2$PPV_ud), 
                                      ifelse(misc2$e_obs == 0 & misc2$d == 0, rbernoulli(1-misc2$NPV_ud,p=1-misc2$NPV_ud), 0))))

# Obtain RR estimate of crude analysis
output_obs <- glm(d~e_obs, data=misc2[misc2$iteration ==1,], family=binomial(link='logit')) # Just do this for iteration 1
# Save the output
estim1 <- data.frame(Parameter="e_exp_t", e_conv = output_obs$coefficients[2],
                          conv_syst = summary(output_obs)$coefficients[2,2])

# Obtain RR estimates for each iteration
e_syst_c <- c()
stderr_syst_c <- c()
for(i in 1:k){
  output_exp <- glm(d~e_exp_t, data=misc2[misc2$iteration == i,], family=binomial(link='logit')) # Do this for every iteration
  e_syst_c <- c(e_syst_c, summary(output_exp)$coefficients[2,1])
  stderr_syst_c <- c(stderr_syst_c, summary(output_exp)$coefficients[2,2])
  # print(i)  # Optional if you want to track how many iterations have been completed
} 

# Save the output
estim2<- data.frame(iteration = 1:k, Parameter = rep("e_exp_t", k), e_syst = e_syst_c, 
                          stderr_syst= stderr_syst_c)

# Dataset with only se, sp, and syst_exp_adj 
estim3 <- distinct(select(misc2[misc2$Obs == 1, ], c(iteration, syst_exp_adj, se_D, sp_D, se_UD, sp_UD)))
estim3$Parameter <- rep("e_exp_t", nrow(estim3))

# merge all data together to get information for each iteration
estim_4 <- merge(estim1, estim2, by="Parameter")
estim4 <- merge(estim_4, estim3, by="iteration")

# Calculate OR intervals for random, systematic, and total error 
estim4$or_rand <- exp(estim4$e_conv - estim4$conv_syst*rnorm(1))
estim4$or_syst <- estim4$syst_exp_adj
estim4$or_tot <- exp(estim4$e_syst - rnorm(1)*estim4$stderr_syst)

# Create a neat output dataset with median as well as 2.5th and 75th percentiles 
output_data <- data.frame(Analysis = c("Random Error", "Systematic Error", "Total Error"), 
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
# View Output
summary(estim4$se_D)
summary(estim4$sp_D)
summary(estim4$sp_D)
par(mfrow=c(1,2))
hist(estim4$se_D)
hist(estim4$sp_D)
hist(estim4$se_UD)
hist(estim4$sp_UD)
par(mfrow=c(1,1))
output_data
