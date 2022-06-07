# Random Error 
# Set up ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 3")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")
library(infer)

# p curves ----
j = 10000
data2 <- data.frame(replicate = 1:j) 
s = 100
pr_E = 0.3
pr_D = 0.2
RD_ED = 0.5

E = rbern(pr_E, n=1000)
D = rbern(pr_D + RD_ED * E, n=1000)

dag1 <- data.frame(E, D) 
data1 <- rep_sample_n(dag1, s, replace = FALSE, reps = j, prob = NULL)

data2[,paste0("pval", s, pr_E, pr_D, RD_ED)] <- 0

for(i in data2$replicate){
  reg <- summary(glm(D~E, data = data1[data1$replicate == i, ]))
  #reg <- chisq.test(table(data1[data1$replicate == i, ]$E, data1[data1$replicate == i, ]$D))
  #data2[data2$replicate == i, ]$pval <- reg$p.value
  data2[data2$replicate == i, paste0("pval", s, pr_E, pr_D, RD_ED) ] <- reg$coefficients[2,4]
  #data2[data2$replicate == i, ]$coef <- reg$coefficients[2]
}

hist(data2[,paste0("pval", s, pr_E, pr_D, RD_ED)], breaks = 30, 
     main = paste0("pval", s,"_", pr_E,"_", pr_D,"_", RD_ED))

# NDEM ----
j = 1000
sp = 0.7
se = 0.8
data2 <- data.frame(replicate = 1:j) 

s = 100
pr_E = 0.3
pr_D = 0.2
RD_ED = 0.5

E = rbern(pr_E, n=1000)
D = rbern(pr_D + RD_ED * E, n=1000)
E_misc = rbern((1-sp) + (se - (1-sp)) * E)

dag1 <- data.frame(E, D, E_misc) 
data1 <- rep_sample_n(dag1, s, replace = FALSE, reps = j, prob = NULL)

data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "misc")] <- 0
data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "true")] <- 0

for(i in data2$replicate){
  reg_misc <- summary(glm(D~E_misc, data = data1[data1$replicate == i, ]))
  reg_true <- summary(glm(D~E, data = data1[data1$replicate == i, ]))
  #reg <- chisq.test(table(data1[data1$replicate == i, ]$E, data1[data1$replicate == i, ]$D))
  #data2[data2$replicate == i, ]$pval <- reg$p.value
  data2[data2$replicate == i, paste0("pval", s, pr_E, pr_D, RD_ED, "misc") ] <- reg_misc$coefficients[2,4]
  data2[data2$replicate == i, paste0("pval", s, pr_E, pr_D, RD_ED, "true") ] <- reg_true$coefficients[2,4]
  #data2[data2$replicate == i, ]$coef <- reg$coefficients[2]
}

hist(data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "misc")], breaks = 30, 
     main = paste0("Misclassified"))
hist(data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "true")], breaks = 30, 
     main = paste0("Truth"))

# DEM ----
j = 1000
sp_D = 0.7
se_D = 0.8
sp_UD = 0.8
se_UD = 0.95
data2 <- data.frame(replicate = 1:j) 
s = 100
pr_E = 0.3
pr_D = 0.2
RD_ED = 0.5

E = rbern(pr_E, n=1000)
D = rbern(pr_D + RD_ED * E, n=1000)
E_misc = rbern((1-se_UD) + ((sp_D - (1-sp_UD)) * D) + ((se_UD - (1-sp_UD)) * E) + 
                ((se_D - (se_UD - (1-sp_UD)) -  (1-sp_UD)) * E * D))

dag1 <- data.frame(E, D, E_misc) 
data1 <- rep_sample_n(dag1, s, replace = FALSE, reps = j, prob = NULL)

data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "misc")] <- 0
data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "true")] <- 0

for(i in data2$replicate){
  reg_misc <- summary(glm(D~E_misc, data = data1[data1$replicate == i, ]))
  reg_true <- summary(glm(D~E, data = data1[data1$replicate == i, ]))
  #reg <- chisq.test(table(data1[data1$replicate == i, ]$E, data1[data1$replicate == i, ]$D))
  #data2[data2$replicate == i, ]$pval <- reg$p.value
  data2[data2$replicate == i, paste0("pval", s, pr_E, pr_D, RD_ED, "misc") ] <- reg_misc$coefficients[2,4]
  data2[data2$replicate == i, paste0("pval", s, pr_E, pr_D, RD_ED, "true") ] <- reg_true$coefficients[2,4]
  #data2[data2$replicate == i, ]$coef <- reg$coefficients[2]
}

hist(data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "misc")], breaks = 30, 
     main = paste0("Misclassified"))
hist(data2[,paste0("pval", s, pr_E, pr_D, RD_ED, "true")], breaks = 30, 
     main = paste0("Truth"))

# Sample size and power ----
pr_D  = 0.2
RD_ED = 0.1
alpha = 0.05
beta = 0.2
n1 = 680
n0 = 680
j = 1000

population <- data.frame(E = c(rep(1, n1), rep(0, n0)))
population$D <- rbern(pr_D + RD_ED * population$E, 680*2)
twobytwo.strat(population, risk = FALSE)
chisq.test(table(population$E, population$D)) 

population2 <- data.frame(E = c(rep(1, n1), rep(0, n0))) %>% slice(rep(1:n(), j)) 
population2$iteration <- rep(1:j, each = n0 + n1) 

D <- list()
for(i in 1:j){
D <- c(D, rbern(pr_D + RD_ED * population2[population2$iteration == i, ]$E, n1 + n0))
}

population2$D <- as.numeric(D)

populationref <- data.frame(iteration = 1:j, pval = NA, OR = NA)

for(i in populationref$iteration){
  reg <- chisq.test(table(population2[population2$iteration == i, ]$E,
                          population2[population2$iteration == i, ]$D))
  populationref[populationref$iteration == i, ]$pval <- reg$p.value
  populationref[populationref$iteration == i, ]$OR <- calc.OR(subset(population2, 
                                                                     population2$iteration == i))
}

populationref$include <- ifelse(populationref$pval <= 0.05, 1, 0)
populationref$OR_f_g <- ifelse(populationref$OR >= 1.71, 1, 0)
populationref$OR_f <- ifelse(populationref$OR_f_g == 1 & populationref$include == 1, 1, 
                             ifelse(populationref$OR_f_g == 0 & populationref$include == 1, 0, NA))

prop.table(table(populationref$include))
prop.table(table(populationref$OR_f))
prop.table(table(populationref$OR_f_g))

# Bootstrapping ----
j = 10000
n = 50000
k = 100

pop <- data.frame(height= rnorm(n, mean = 65, sd = 5))
pop <- rep_sample_n(pop, k, replace = FALSE, reps = j, prob = NULL)

pop1 <- pop %>% group_by(replicate) %>%
  mutate(mean = mean(height), 
         sd = sd(height), 
         lower = mean - 1.96 * sd/sqrt(20-1),
         upper = mean + 1.96 * sd/sqrt(20-1),
         quantile(mean(pop$height), 0.025),
         quantile(mean(pop$height), 0.75))

hist(pop1$mean)

# boostrapping examples ----
library(broom)
n = 100
j = 1000
C = rbern(0.3, n)
E = rbern(0.1 + 0.25 * C, n)
D = rbern(0.07 + C * 0.15, n)

data1 <- data.frame(C, D, E)
summary(glm(D~E, data = data1))

data2 <- rep_sample_n(data1, n , replace = TRUE, reps = j, prob = NULL)

data3 <- data2 %>% 
  group_by(replicate) %>% 
    do(model = tidy(glm(D~E, data=.))[2,2:3]) %>% 
  unnest(cols = model) 

mean(data3$estimate)

