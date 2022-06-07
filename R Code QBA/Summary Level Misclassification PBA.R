library(dplyr)
library(trapezoid)

i = 1000000 # number of iterations

a = 215  # corresponds to the a cell in a 2x2 contingency table
b = 668 # corresponds to the b cell in a 2x2 contingency table
c = 1449  # corresponds to the c cell in a 2x2 contingency table
d= 4296  # corresponds to the d cell in a 2x2 contingency table

alpha_1 = 50.6 # alpha for a beta distribution of the sensitivitiiy 
beta_1 = 14.3  # beta for a beta distribution of the sensitivitiiy 
alpha_0 = 70 # alpha for a beta distribution of the specificity 
beta_0 = 1  # beta for a beta distribution of the specificity 

miscdata <- data.frame(id = 1:i, ed_o = rep(a, i), eud_o = rep(b, i),
                       ued_o = rep(c, i), ueud_o = rep(d, i),
                       se_D = rbeta(i, alpha_1, beta_1), sp_D = rbeta(i, alpha_0, beta_0))
miscdata$se_UD <- miscdata$se_D
miscdata$sp_UD <- miscdata$sp_D

# Calcaulate the crude OR and SE(LN(OR))
miscdata$conv_or <- (miscdata$ed_o*miscdata$ueud_o)/(miscdata$eud_o*miscdata$ued_o) 
miscdata$conv_std <- sqrt(1/miscdata$ed_o + 1/miscdata$eud_o + 1/miscdata$ued_o + 1/miscdata$ueud_o)

#expected based on sens and spec
miscdata$ed_t=(miscdata$ed_o-(1-miscdata$sp_D)*(miscdata$ed_o+miscdata$ued_o))/(miscdata$se_D-(1-miscdata$sp_D))
miscdata$ued_t=(miscdata$ed_o+miscdata$ued_o)-miscdata$ed_t
miscdata$eud_t=(miscdata$eud_o-(1-miscdata$sp_UD)*(miscdata$eud_o+miscdata$ueud_o))/(miscdata$se_UD-(1-miscdata$sp_UD))
miscdata$ueud_t=(miscdata$eud_o+miscdata$ueud_o)-miscdata$eud_t

#sample prevalence data
miscdata$prev_e_d = rbeta(i,miscdata$ed_t, miscdata$ued_t)
miscdata$prev_e_ud = rbeta(i,miscdata$eud_t, miscdata$ueud_t)

# PPV and NPV
miscdata$prev_e_d <- ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, miscdata$prev_e_d)
miscdata$prev_e_ud <- ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, miscdata$prev_e_ud)
miscdata$PPV_d <- ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, 
                         (miscdata$sp_D * (1-miscdata$prev_e_d))/(((1-miscdata$se_D) * miscdata$prev_e_d)+(miscdata$sp_D*(1-miscdata$prev_e_d))))
miscdata$NPV_d <- ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, 
                      (miscdata$se_D*miscdata$prev_e_d)/((miscdata$se_D* miscdata$prev_e_d)+((1-miscdata$sp_D)*(1-miscdata$prev_e_d))))

miscdata$PPV_ud <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, 
                          (miscdata$se_UD* miscdata$prev_e_ud)/((miscdata$se_UD* miscdata$prev_e_ud)+((1-miscdata$sp_UD)*(1-miscdata$prev_e_ud))))

miscdata$NPV_ud <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA, 
                           (miscdata$sp_UD * (1-miscdata$prev_e_ud))/(((1-miscdata$se_UD) * miscdata$prev_e_ud)+(miscdata$sp_UD*(1-miscdata$prev_e_ud))))
# sample cell counts
miscdata$x1 <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,rbinom(i, prob = miscdata$PPV_d, size = miscdata$ed_o))
miscdata$x2 <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,rbinom(i, prob = miscdata$PPV_ud, size = miscdata$eud_o))
miscdata$x3 <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,rbinom(i, prob = miscdata$NPV_d, size = miscdata$ued_o))
miscdata$x4 <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,rbinom(i, prob = miscdata$NPV_ud, size = miscdata$ueud_o))

# sampled misclassfied adjusted OR
miscdata$ed_t_s <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,(miscdata$x1 + miscdata$ued_o-miscdata$x3))
miscdata$eud_t_s <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,(miscdata$x2 + miscdata$ueud_o-miscdata$x4))
miscdata$ued_t_s <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,(miscdata$x3 + miscdata$ed_o-miscdata$x1))
miscdata$ueud_t_s <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,(miscdata$x4 + miscdata$eud_o-miscdata$x2))

miscdata$e_syst <-  ifelse(miscdata$ed_t <= 0 | miscdata$eud_t <= 0 | miscdata$ued_t <= 0 | miscdata$ueud_t <= 0, NA,
                           (miscdata$ed_t_s/ miscdata$eud_t_s) / (miscdata$ued_t_s/miscdata$ueud_t_s))
# conventional random error OR
miscdata$or_rand <- ifelse(miscdata$ed_t_s <= 0 | miscdata$eud_t_s <= 0 | miscdata$ued_t_s <= 0 | miscdata$ueud_t_s <= 0, NA,
                           exp(log(miscdata$conv_or - miscdata$conv_std*rnorm(i))))
# expected misclassification adjusted OR
miscdata$or_syst <- ifelse(miscdata$ed_t_s <= 0 | miscdata$eud_t_s <= 0 | miscdata$ued_t_s <= 0 | miscdata$ueud_t_s <= 0, NA,
                           (miscdata$ed_t /miscdata$eud_t) / (miscdata$ued_t /miscdata$ueud_t))
# standard error systematic adjusted
miscdata$stderr_syst <- ifelse(miscdata$ed_t_s <= 0 | miscdata$eud_t_s <= 0 | miscdata$ued_t_s <= 0 | miscdata$ueud_t_s <= 0, NA,
                               sqrt(1/miscdata$ed_t_s+1/miscdata$eud_t_s+1/miscdata$ued_t_s+1/miscdata$ueud_t_s))
# total error estimates
miscdata$or_tot <- ifelse(miscdata$ed_t_s <= 0 | miscdata$eud_t_s <= 0 | miscdata$ued_t_s <= 0 | miscdata$ueud_t_s <= 0, NA,
                               exp(log(miscdata$e_syst) - rnorm(i)*miscdata$stderr_syst))

output_data <- data.frame(Analysis = c("Random Error", "Systematic Error", "Total Error"), 
                          Percentile_lower = c(quantile(miscdata$or_rand, 0.25, na.rm = TRUE),
                                          quantile(miscdata$or_syst, 0.25, na.rm = TRUE),
                                          quantile(miscdata$or_tot, 0.25, na.rm = TRUE)),
                          Estimate = c(quantile(miscdata$or_rand, 0.5, na.rm = TRUE),
                                               quantile(miscdata$or_syst, 0.5, na.rm = TRUE),
                                               quantile(miscdata$or_tot, 0.5, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(miscdata$or_rand, 0.75, na.rm = TRUE),
                                               quantile(miscdata$or_syst, 0.75, na.rm = TRUE),
                                               quantile(miscdata$or_tot, 0.75, na.rm = TRUE)))
output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower

# View Output
summary(miscdata$se_D)
summary(miscdata$sp_D)
par(mfrow=c(1,2))
hist(miscdata$se_D)
hist(miscdata$sp_D)
par(mfrow=c(1,1))
output_data