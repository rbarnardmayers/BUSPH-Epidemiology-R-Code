# Class 2 R Code
library(dplyr)
library(survey)
library(purrr)
library(EnvStats)

# Warm up exercise ----
# N = 300
# pr_C = 0.50
# Pr_E = 0.2 + 0.5*c 
# Pr_D = 0.1 + 0.3*c + 0.2*e
# se = rbeta(1, 90, 10)
# sp = rbeta(1, 95,5)

set.seed(124543)
k=300
c=as.numeric(rbernoulli(k, 0.5))
e=as.numeric(rbernoulli(k, 0.2+0.5*c))
d=as.numeric(rbernoulli(k, 0.01+0.30*c+0.2*e))

dat <- as.data.frame(cbind(id=1:k, e, d, c))

dat$c_misc <- ifelse(dat$c == 1, as.numeric(rbernoulli(k, p=0.9)), 
                     rbernoulli(k, p=1-0.95))
set.seed(NULL)

# Linear Regression ---
summary(lm(dat$e~dat$c))$coefficient[2,1]
summary(lm(dat$d~dat$c))$coefficient[2,1]
summary(lm(dat$d~dat$e))$coefficient[2,1]
summary(lm(dat$d~dat$e + dat$c))$coefficient[2,1]
summary(lm(dat$d~dat$e+dat$c_misc))$coefficient[2,1]

# Confounder Misclassification ----
j = 1000
misc <- dat %>% slice(rep(1:n(), j))
misc$iteration <- rep(1:j, each=nrow(dat))

test <- as.data.frame(table(misc[misc$iteration == 1,]$d, 
                            misc[misc$iteration == 1,]$e, 
                            misc[misc$iteration == 1,]$c))

test <- rename(test, d = Var1)
test <- rename(test, e = Var2)
test <- rename(test, c = Var3)


# 1 dataset with all four combinations (1 row, 8 columns)
test2 <- data.frame(cp_ed_o = test[test$d == 1 & test$e == 1 & test$c == 1,]$Freq,
                    cp_eud_o = test[test$d == 0 & test$e == 1 & test$c == 1,]$Freq,
                    cp_ued_o = test[test$d == 1 & test$e == 0 & test$c == 1,]$Freq,
                    cp_ueud_o = test[test$d == 0 & test$e == 0 & test$c == 1,]$Freq,
                    
                    cn_ed_o = test[test$d == 1 & test$e == 1 & test$c == 0,]$Freq,
                    cn_eud_o = test[test$d == 0 & test$e == 1 & test$c == 0,]$Freq,
                    cn_ued_o = test[test$d == 1 & test$e == 0 & test$c == 0,]$Freq,
                    cn_ueud_o = test[test$d == 0 & test$e == 0 & test$c == 0,]$Freq)

test3 <- data.frame(iteration = 1:j, cp_ed_o = rep(test2$cp_ed_o, j), 
                    cp_eud_o = rep(test2$cp_eud_o, j), 
                    cp_ued_o = rep(test2$cp_ued_o, j), 
                    cp_ueud_o = rep(test2$cp_ueud_o, j),
                    cn_ed_o = rep(test2$cn_ed_o, j), 
                    cn_eud_o = rep(test2$cn_eud_o, j), 
                    cn_ued_o = rep(test2$cn_ued_o, j), 
                    cn_ueud_o = rep(test2$cn_ueud_o, j))


# This is for non-differential non-correlated. For differential, change lines 70 and 71 
# These data would come from a valdiation dataset 
test3$se <- rbeta(j, 90, 10)  
test3$sp <- rbeta(j, 95, 5) 

test3$cp_ed_t=(test3$cp_ed_o-(1-test3$sp)*(test3$cp_ed_o+test3$cn_ed_o))/(test3$se-(1-test3$sp))
test3$cp_ued_t=(test3$cp_ued_o-(1-test3$sp)*(test3$cp_ued_o+test3$cn_ued_o))/(test3$se-(1-test3$sp))
test3$cp_eud_t=(test3$cp_eud_o-(1-test3$sp)*(test3$cp_eud_o+test3$cn_eud_o))/(test3$se-(1-test3$sp))
test3$cp_ueud_t=(test3$cp_ueud_o-(1-test3$sp)*(test3$cp_ueud_o+test3$cn_ueud_o))/(test3$se-(1-test3$sp))

test3$cn_ed_t=(test3$cp_ed_o+ test3$cn_ed_o)-test3$cp_ed_t
test3$cn_ued_t=(test3$cp_ued_o+ test3$cn_ued_o)-test3$cp_ued_t 
test3$cn_eud_t=(test3$cp_eud_o+ test3$cn_eud_o)-test3$cp_eud_t
test3$cn_ueud_t=(test3$cp_ueud_o+ test3$cn_ueud_o)-test3$cp_ueud_t

test3$cn = test3$cn_ed_t+test3$cn_ued_t+test3$cn_eud_t+test3$cn_ueud_t
test3$cp = test3$cp_ed_t+test3$cp_ued_t+test3$cp_eud_t+test3$cp_ueud_t

#expected misclassification adjusted
test3$syst_exp_adj_smr = (test3$cp_ed_t+ test3$cn_ed_t) / 
  ((test3$cp_ued_t * test3$cp_eud_t / test3$cp_ueud_t) + (test3$cn_ued_t * test3$cn_eud_t / test3$cn_ueud_t))

test3$syst_exp_adj =(test3$cp_ed_t * test3$cp_ueud_t / test3$cp + test3$cn_ed_t * test3$cn_ueud_t / test3$cp) / 
  (test3$cp_ued_t * test3$cp_eud_t / test3$cn + test3$cn_ued_t * test3$cn_eud_t / test3$cn)

# sample prevalence
test3$prev_c_e_d <- ifelse(test3$cp_ed_t > 0 & test3$cn_ed_t > 0,
                         rbeta(j, test3$cp_ed_t, test3$cn_ed_t), NA)

test3$prev_c_ue_d <- ifelse( test3$cp_ued_t > 0 & test3$cn_ued_t > 0,
                            rbeta(k, test3$cp_ued_t, test3$cn_ued_t), NA)

test3$prev_c_e_ud <- ifelse( test3$cp_eud_t > 0 & test3$cn_eud_t > 0,
                            rbeta(k, test3$cp_eud_t, test3$cn_eud_t), NA)

test3$prev_c_ue_ud <- ifelse( test3$cp_ueud_t > 0 & test3$cn_ueud_t > 0,
                             rbeta(k, test3$cp_ueud_t, test3$cn_ueud_t), NA)

# Calculate PPV and NPV
# want to make negative values NA for PPV and NPV
test3$PPV_ed <- ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                         test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                       (test3$se*test3$prev_c_e_d)/((test3$se* test3$prev_c_e_d)+((1-test3$sp)*(1-test3$prev_c_e_d))), NA)
test3$NPV_ed <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                          test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                        (test3$sp*(1-test3$prev_c_e_d))/(((1-test3$se) * test3$prev_c_e_d)+(test3$sp*(1-test3$prev_c_e_d))), NA)

test3$PPV_ued <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                           test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                         (test3$se*test3$prev_c_ue_d)/((test3$se* test3$prev_c_ue_d)+((1-test3$sp)*(1-test3$prev_c_ue_d))), NA)

test3$NPV_ued <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                           test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                         (test3$sp*(1-test3$prev_c_ue_d))/(((1-test3$se) * test3$prev_c_ue_d)+(test3$sp*(1-test3$prev_c_ue_d))), NA)

test3$PPV_eud <- ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                          test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                        (test3$se*test3$prev_c_e_ud)/((test3$se* test3$prev_c_e_ud)+((1-test3$sp)*(1-test3$prev_c_e_ud))), NA)

test3$NPV_eud <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                           test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                         (test3$sp*(1-test3$prev_c_e_ud))/(((1-test3$se) * test3$prev_c_e_ud)+(test3$sp*(1-test3$prev_c_e_ud))), NA)

test3$PPV_ueud <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                            test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                          (test3$se*test3$prev_c_ue_ud)/((test3$se*test3$prev_c_ue_ud)+((1-test3$sp)*(1-test3$prev_c_ue_ud))), NA)

test3$NPV_ueud <-  ifelse(test3$cp_ed_t > 0 & test3$cp_ued_t > 0 & test3$cp_eud_t > 0 & test3$cp_ueud_t > 0 & 
                            test3$cn_ed_t > 0 & test3$cn_ued_t > 0 & test3$cn_eud_t > 0 & test3$cn_ueud_t > 0,
                          (test3$sp*(1-test3$prev_c_ue_ud))/(((1-test3$se) * test3$prev_c_ue_ud)+(test3$sp*(1-test3$prev_c_ue_ud))), NA)

# Create new dataset without "true" cells 
misc2 <- merge(test3, misc, by="iteration")

# Calculate expected cells based off PPV and NPV
misc2$c_exp_t <- ifelse(misc2$e == 1 & misc2$d == 1 & misc2$c_misc == 1, rbernoulli(misc2$PPV_ed), 
                        ifelse(misc2$e == 1 & misc2$d == 1 & misc2$c_misc == 0, rbernoulli(1-misc2$NPV_ed), 
                               ifelse(misc2$e == 1 & misc2$d == 0 & misc2$c_misc == 1, rbernoulli(misc2$PPV_eud), 
                                      ifelse(misc2$e == 1 & misc2$d == 0 & misc2$c_misc == 0, rbernoulli(1-misc2$NPV_eud), 
                                             ifelse(misc2$e == 0 & misc2$d == 1 & misc2$c == 1, rbernoulli(misc2$PPV_ued),
                                                ifelse(misc2$e == 0 & misc2$d == 1 & misc2$c == 0, rbernoulli(misc2$PPV_ued),
                                                       ifelse(misc2$e == 0 & misc2$d == 0 & misc2$c == 1, rbernoulli(misc2$PPV_ueud),
                                                              ifelse(misc2$e == 0 & misc2$d == 0 & misc2$c == 0, rbernoulli(misc2$PPV_ueud), NA))))))))

# Obtain RR estimate of crude analysis
# Just do this for iteration 1
output_obs <- glm(d~e+c_misc, data=misc2[misc2$iteration ==1,], family=binomial(link='logit')) 

# Save the output
estim1 <- data.frame(Parameter="c_exp_t", c_conv = output_obs$coefficients[2],
                     conv_syst = summary(output_obs)$coefficients[2,2])

# Obtain RR estimates for each iteration
c_syst_c <- c()
stderr_syst_c <- c()
for(i in 1:k){
  output_exp <- glm(d~e+c_exp_t, data=misc2[misc2$iteration == i,], family=binomial(link='logit')) # Do this for every iteration
  c_syst_c <- c(c_syst_c, summary(output_exp)$coefficients[2,1])
  stderr_syst_c <- c(stderr_syst_c, summary(output_exp)$coefficients[2,2])
  # print(i)  # Optional if you want to track how many iterations have been completed
} 

# Save the output
estim2<- data.frame(iteration = 1:k, Parameter = rep("c_exp_t", k), c_syst = c_syst_c, 
                    stderr_syst= stderr_syst_c)

# Dataset with only se, sp, and syst_exp_adj 
estim3 <- select(misc2, c(iteration, syst_exp_adj, se, sp))
estim3$Parameter <- rep("c_exp_t", nrow(estim3))

# merge all data together to get information for each iteration
estim_4 <- merge(estim1, estim2, by="Parameter")
estim4 <- merge(estim_4, estim3, by="iteration")

# Calculate OR intervals for random, systematic, and total error 
estim4$or_rand <- exp(estim4$c_conv - estim4$conv_syst*rnorm(1))
estim4$or_syst <- estim4$syst_exp_adj
estim4$or_tot <- exp(estim4$c_syst - rnorm(1)*estim4$stderr_syst)

# Create a neat output dataset with median as well as 2.5th and 75th percentiles 
output_data <- data.frame(Analysis = c("Random Error", "Systematic Error", "Total Error"), 
                          Estimate = c(quantile(estim4$or_rand, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_syst, 0.5, na.rm = TRUE),
                                       quantile(estim4$or_tot, 0.5, na.rm = TRUE)),
                          Percentile_lower = c(quantile(estim4$or_rand, 0.25, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.25, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.25, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(estim4$or_rand, 0.75, na.rm = TRUE),
                                               quantile(estim4$or_syst, 0.75, na.rm = TRUE),
                                               quantile(estim4$or_tot, 0.75, na.rm = TRUE)))

output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower

# View Output
output_data


# Checking work by comparing to SAS Data ----
setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/R Code/SAS Data 2")

estim1
estim2
estim3
estim4
estim5
misc
misc2
test3 


