###### R CODE FOR RECORD LEVEL BIAS ANALYSIS: UNCONTROLLED CONFOUNDING #############
# Created by: 		Ruby Barnard-Mayers
# Date created: 	Jan 12, 2022
# Purpose: 			Simple R code for individual record level PBA for uncontrolled confounding
#########################################################################################

library(purrr)
library(dplyr)
j = 100  # number of iterations to run

# Data information ----
# IF YOU HAVE INDIVIDUAL LEVEL DATA ALREADY, uncomment the lines below and skip lines 15-27
# conf <- read.csv(file.choose) # upload as csv
# conf <- rename(conf, e_obs = exposure, d = outcome)  # rename exposure and outcome columns to e_obs and d
# conf <- cbind(conf. j)

# if you need to simulate record level data from a 2x2 table

a = 105  # corresponds to the a cell in a 2x2 contingency table
b = 85  # corresponds to the b cell in a 2x2 contingency table
c = 527  # corresponds to the c cell in a 2x2 contingency table
d= 93  # corresponds to the d cell in a 2x2 contingency table
n = b + d  # total number of people in unexposed
m = a + c # total number of people in exposed

# These values (22,23,15,55) correspond to (a,b,c,d) in a 2x2 table where 115 is the total N
conf <- data.frame(Obs = rep(1:j, n+m), iteration=rep(1:j, each=n+m), e = rep(c(rep(1,a+c), rep(0,b+d)),j),
                   y = rep(c(rep(1, a), rep(0, c), rep(1, b), rep(0,d)),j))


# Freq counts that mirror 2x2 table (1 record per combination)
test <- as.data.frame(table(conf[conf$iteration == 1,]$y, conf[conf$iteration == 1,]$e))
test <- rename(test, y = Var1)
test <- rename(test, e = Var2)

# 1 dataset with all four combinations
test2 <- data.frame(a = test[test$y == 1 & test$e == 1,]$Freq,
                    b = test[test$y == 0 & test$e == 1,]$Freq,
                    c = test[test$y == 1 & test$e == 0,]$Freq,
                    d = test[test$y == 0 & test$e == 0,]$Freq)

k = j 
test3 <- data.frame(iteration = 1:k)

test3$p_c1_e1 = rtrapezoid(j, 0.7, 0.75,0.85, 0.9)
test3$p_c1_e0 = rtrapezoid(j, 0.03, 0.04,0.07,0.10)
test3$p_RRc = rtrapezoid(j, 0.5, 0.6,0.7, 0.8)


test3$m1 <- test3$p_c1_e1*m
test3$n1 <- test3$p_c1_e0*n
test3$n0 <- n - test3$n1
test3$m0 <- m - test3$m1
test3$a1 <- test3$p_RRc * test3$m1 * a / (test3$p_RRc * test3$m1 + m-test3$m1)
test3$b1 <- test3$p_RRc * test3$n1 * b / (test3$p_RRc * test3$n1 + n-test3$n1)  
test3$c1 <- test3$m1 - test3$a1
test3$d1 <- test3$n1 - test3$b1
test3$a0 <- a-test3$a1
test3$d0 <- d-test3$d1
test3$b0 <- b - test3$b1
test3$c0 <- c - test3$c1

test3$syst_exp_adj = a / ((test3$m1 * test3$b1 /test3$n1) +(test3$m0 * test3$b0 /test3$n0))


test3$stderr_mh=sqrt(((((test3$a0+test3$b0)*test3$m0*test3$n0/(test3$m0+test3$n0)^2-test3$a0*test3$b0/(test3$m0+test3$n0))) + 
                        (((test3$a1+test3$b1)*test3$m1*test3$n1/(test3$m1+test3$n1)^2-test3$a1*test3$b1/(test3$m1+test3$n1))))/ 
                       (((test3$a0*test3$n0/(test3$m0+test3$n0)+test3$a1*test3$n1/(test3$m1+test3$n1)))*(test3$b0*test3$m0/(test3$m0+test3$n0)+test3$b1*test3$m1/(test3$m1+test3$n1))))

test3$a_conf_prob = test3$a1/a
test3$b_conf_prob = test3$b1/b
test3$c_conf_prob = test3$c1/c
test3$d_conf_prob = test3$d1/d

test3$check <- ifelse(test3$a_conf_prob > 0 & test3$b_conf_prob > 0 & 
                        test3$c_conf_prob > 0 & test3$d_conf_prob > 0, 1, 0)

test3 <- test3[test3$check == 1,]

conf2 <- merge(test3, conf, by="iteration")

conf2$conf <- ifelse(conf2$e == 1 & conf2$y == 1, rbernoulli(conf2$a_conf_prob), 
                        ifelse(conf2$e == 0 & conf2$y == 1, rbernoulli(1-conf2$b_conf_prob), 
                               ifelse(conf2$e == 1 & conf2$y == 0, rbernoulli(conf2$c_conf_prob), 
                                      ifelse(conf2$e == 0 & conf2$y == 0, rbernoulli(1-conf2$d_conf_prob), 0))))


# Obtain RR estimate of crude analysis
output_obs <- glm(y~e, data=conf2[conf2$iteration ==1,], family=binomial(link='logit')) # Just do this for iteration 1

# Save the output
estim1 <- data.frame(Parameter="e", e_conv = output_obs$coefficients[2],
                     stderr_conv = summary(output_obs)$coefficients[2,2])

# Obtain RR estimates for each iteration
e_syst_c <- c()
stderr_syst_c <- c()
for(i in 1:k){
  output_exp <- glm(y~e+conf, data=conf2[conf2$iteration == i,], family=binomial(link='logit')) # Do this for every iteration
  e_syst_c <- c(e_syst_c, summary(output_exp)$coefficients[2,1])
  stderr_syst_c <- c(stderr_syst_c, summary(output_exp)$coefficients[2,2])
} 
# Save the output
estim2<- data.frame(iteration = 1:k, Parameter = rep("e", k), e_syst = e_syst_c, 
                    stderr_syst= stderr_syst_c)

# Dataset with only se, sp, and syst_exp_adj 
estim3 <- distinct(select(conf2[conf2$Obs == 1, ], c(iteration, syst_exp_adj, stderr_mh, p_c1_e1, p_c1_e0, p_RRc)))
estim3$Parameter <- rep("e", nrow(estim3))

# merge all data together to get information for each iteration
test4 <- merge(estim1, estim2, by="Parameter")
estim4 <- merge(test4, estim3, by="iteration")

# Calculate OR intervals for random, systematic, and total error 
estim4$rr_rand <- exp(estim4$e_conv - estim4$stderr_conv*rnorm(1))
estim4$rr_syst <- estim4$syst_exp_adj
estim4$rr_tot <- exp(estim4$e_syst - rnorm(1)*estim4$stderr_mh)

# Create a neat output dataset with 2.5, .5 and .75 percentiles 
output_data <- data.frame(Analysis = c("Random Error", "Systematic Error", "Total Error"), 
                          Percentile_lower = c(quantile(estim4$rr_rand, 0.25, na.rm = TRUE),
                                               quantile(estim4$rr_syst, 0.25, na.rm = TRUE),
                                               quantile(estim4$rr_tot, 0.25, na.rm = TRUE)),
                          Estimate = c(quantile(estim4$rr_rand, 0.5, na.rm = TRUE),
                                       quantile(estim4$rr_syst, 0.5, na.rm = TRUE),
                                       quantile(estim4$rr_tot, 0.5, na.rm = TRUE)),
                          Percentile_Upper = c(quantile(estim4$rr_rand, 0.75, na.rm = TRUE),
                                               quantile(estim4$rr_syst, 0.75, na.rm = TRUE),
                                               quantile(estim4$rr_tot, 0.75, na.rm = TRUE)))

output_data$Width <- output_data$Percentile_Upper/output_data$Percentile_lower

# View Output
summary(estim4$se_D)
summary(estim4$sp_D)
par(mfrow=c(1,2))
hist(estim4$se_D, breaks = 50)
hist(estim4$sp_D, breaks = 50)

output_data

