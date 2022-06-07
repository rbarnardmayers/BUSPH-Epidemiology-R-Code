###### SAS CODE FOR SUMMARY LEVEL BIAS ANALYSIS: EXPOSURE MISCLASSIFICATION #############
# Created by: 		Matt Fox
# Date created: 	Sept 19, 2018
# Date modified: 	January 5, 2022
# Purpose: 			SAS code for summary level PBA for exposure misclassification
##########################################################################################
# NOTES: code to be edited is listed below between lines that say 
#        "#####EDIT#####" and "####END EDIT####"
#        code calculates odds ratios, modify for risk ratios or risk differences
#########################################################################################
library(dplyr)
library(purrr)
library(trapezoid)

j = 1000000
a = 105
b = 527 
c = 85
d = 93 
m = a + b
n = c + d

conf <- data.frame(iteration = 1:j,
                   p_c1_e1 = rtrapezoid(j, 0.7, 0.75,0.85, 0.9), 
                   p_c1_e0 = rtrapezoid(j, 0.03, 0.04,0.07,0.10), 
                   p_RRc = rtrapezoid(j, 0.5, 0.6,0.7, 0.8))

conf$m1 <- conf$p_c1_e1*m
conf$n1 <- conf$p_c1_e0*n
conf$n0 <- n - conf$n1
conf$m0 <- m - conf$m1
conf$a1 <- conf$p_RRc * conf$m1 * a / (conf$p_RRc * conf$m1 + m-conf$m1)
conf$c1 <- conf$p_RRc * conf$n1 * c / (conf$p_RRc * conf$n1 + n-conf$n1)  
conf$b1 <- conf$m1 - conf$a1
conf$d1 <- conf$n1 - conf$c1
conf$a0 <- a-conf$a1
conf$d0 <- d-conf$d1
conf$b0 <- b - conf$b1
conf$c0 <- c - conf$c1

conf$a_conf_prob = conf$a1/a
conf$b_conf_prob = conf$b1/b
conf$c_conf_prob = conf$c1/c
conf$d_conf_prob = conf$d1/d

conf$a1_exp_s <- rbinom(j, a, conf$a_conf_prob)
conf$b1_exp_s <- rbinom(j, b, conf$b_conf_prob)
conf$c1_exp_s <- rbinom(j, c, conf$c_conf_prob)
conf$d1_exp_s <- rbinom(j, d, conf$d_conf_prob)

conf$a0_exp_s <- a - conf$a1_exp_s
conf$b0_exp_s <- b - conf$b1_exp_s
conf$c0_exp_s <- c - conf$c1_exp_s
conf$d0_exp_s <- d - conf$d1_exp_s

conf$m1_exp_s <- conf$a1_exp_s + conf$b1_exp_s
conf$m0_exp_s <- conf$a0_exp_s + conf$b0_exp_s

conf$n1_exp_s <- conf$c1_exp_s + conf$d1_exp_s
conf$n0_exp_s <- conf$c0_exp_s + conf$d0_exp_s

conf$check <- ifelse(conf$a1_exp_s > 0 & conf$b1_exp_s > 0 & 
                       conf$c1_exp_s > 0 & conf$d1_exp_s > 0 & 
                       conf$a0_exp_s > 0 & conf$b0_exp_s > 0 & 
                       conf$c0_exp_s > 0 & conf$d0_exp_s > 0, 1, 0)

# to see the rows with 0 check 
nrow(conf[conf$check ==0,])
conf <- conf[conf$check == 1,]
# conventional random error OR and SD
conf$rr_conv <- (a/m)/(c/n) 
conf$conv_std <- sqrt (1/a - 1/m + 1/b - 1/n)

# Expected confounder adjusted, and adjusted OR (Using MH not SMR)
conf$rr_rand <- exp(log(conf$rr_conv) - conf$conv_std*rnorm(1))
conf$rr_syst_exp <- (conf$a0*conf$n0/(conf$m0+conf$n0)+conf$a1*conf$n1/(conf$m1+conf$n1))/(conf$c0*conf$m0/(conf$m0+conf$n0)+conf$c1*conf$m1/(conf$m1+conf$n1))
# Probabilistic confounder adjusted (simulated);
conf$rr_syst <- (conf$a0_exp_s*conf$n0_exp_s/(conf$m0_exp_s+conf$n0_exp_s)+conf$a1_exp_s*conf$n1_exp_s/(conf$m1_exp_s+conf$n1_exp_s))/
  (conf$c0_exp_s*conf$m0_exp_s/(conf$m0_exp_s+conf$n0_exp_s)+conf$c1_exp_s*conf$m1_exp_s/(conf$m1_exp_s+conf$n1_exp_s))

# Standard error systematic adjusted
conf$stderr_syst <- sqrt(((((conf$a0+conf$c0)*conf$m0*conf$n0/(conf$m0+conf$n0)**2-conf$a0*conf$c0/(conf$m0+conf$n0))) + 
                            (((conf$a1+conf$c1)*conf$m1*conf$n1/(conf$m1+conf$n1)**2-conf$a1*conf$c1/(conf$m1+conf$n1))))/
                   (((conf$a0*conf$n0/(conf$m0+conf$n0)+conf$a1*conf$n1/(conf$m1+conf$n1)))*(conf$c0*conf$m0/(conf$m0+conf$n0)+conf$c1*conf$m1/(conf$m1+conf$n1))))

# Total Error Estimate
conf$rr_tot <- exp(log(conf$rr_syst) - (rnorm(1)*conf$stderr_syst))

# Output

output_data <- data.frame(Analysis = c("Random Error", "Systematic Error", "Total Error"), 
                          Percentile_lower = c(quantile(conf$rr_rand, 0.25, na.rm = TRUE),
                                               quantile(conf$rr_syst_exp, 0.25, na.rm = TRUE),
                                               quantile(conf$rr_tot, 0.25, na.rm = TRUE)),
                          Esitmate = c(quantile(conf$rr_rand, 0.5, na.rm = TRUE),
                                       quantile(conf$rr_syst_exp, 0.5, na.rm = TRUE),
                                       quantile(conf$rr_tot, 0.5, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(conf$rr_rand, 0.75, na.rm = TRUE),
                                               quantile(conf$rr_syst_exp, 0.75, na.rm = TRUE),
                                               quantile(conf$rr_tot, 0.75, na.rm = TRUE)))

output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower

output_data
par(mfrow=c(1,1))
hist(conf$p_RRc)
hist(conf$p_c1_e0)
hist(conf$p_c1_e1)

summary(conf$p_RRc)
summary(conf$p_c1_e0)
summary(conf$p_c1_e1)
