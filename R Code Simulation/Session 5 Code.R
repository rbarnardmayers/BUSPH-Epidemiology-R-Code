# Random Error
# Set up ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 3")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")
library(infer)

#  Challenge (not discussed) ----
U1 = rbern(0.2)
U2 = rbern(0.4)
W = rbern(0.3 + 0.1 * U1 + 0.2 * U2)
Z = rbern(0.2 + 0.4 * W)
X = rbern(0.1 + 0.2 * Z + 0.3 * U2)
Y = rbern(0.1 * 0.2 * X + 0.3 * U1)

dag1 <- data.frame(U1, U2, W, Z, X, Y)
lm(Z~W, data = dag1)$coefficients[2] 
lm(Y~Z, data = dag1)$coefficients[2] 
lm(Y~X, data = dag1)$coefficients[2] 

# Random error (small n) ----
C1 = rbern(0.3, n = 500)
C1_alt = rnorm(n = 500)
C2 = rbern(0.4 + 0.1 * C1, n = 500)
C2_alt = rbern(0.4 + 0.1* C1_alt, n = 500)
E = rbern(0.2 + 0.3 * C1, n = 500)
E_alt = rbern(0.2 + 0.3 * C1_alt, n = 500)
D = rbern(0.1 + 0.3 * E + 0.5 * C2, n = 500)
D_alt = rbern(0.1 + 0.3 * E_alt + 0.5 * C2_alt, n = 500)

dag1 <- data.frame(C1, C2, E, D, C1_alt, C2_alt, E_alt, D_alt)
lm(D~E + C2, data = dag1)$coefficients[2]
lm(D~E + C1, data = dag1)$coefficients[2]

lm(D_alt~E_alt + C2_alt, data = dag1)$coefficients[2]
lm(D_alt~E_alt + C1_alt, data = dag1)$coefficients[2]

# Random Error Set up----
a = 30
b = 20
c = 470 
d = 480
n1 = a + c 
n0 = b + d
n = n1 + n0
m1 = a + b
m2 = c + d

data1 <- data.frame(E = c(rep(1, n1), rep(0, n0)),
                    D = c(rep(1, a), rep(0, c), rep(1, b), rep(0, d)))

twobytwo.strat(data1, risk=FALSE)

# Hypergeometric distribution
p.val <- rhyper(1000, 50, 500, 1000)

(choose(m1, a) * choose(n1 - m1, n1 - a)) / choose(n, n1)

# Dataset 2 ----
a = 20
b = 10
c = 80 
d = 90
n1 = a + c 
n0 = b + d
n = n1 + n0
m1 = a + b
m2 = c + d
j = 100000

data2 <- data.frame(E = c(rep(1, n1), rep(0, n0)),
                    D = c(rep(1, a), rep(0, c), rep(1, b), rep(0, d))) # p-value = 0.07471

data3 <- data2 %>% slice(rep(1:n(), each = j)) 
data3$iteration <- rep(1:j, nrow(data2))
data3$k <- runif(n*j)

twobytwo.strat(data3[data3$iteration == 2, ], risk=FALSE)

data3 <- data3[order(data3$k),]
data3 <- data3[order(data3$iteration),]

data4 <- data3 

data4$E_new <- rep(c(rep(1, n1), rep(0, n0)), j)
data4$Count <- rep(1, j*n)

data5 <- data4[data4$E_new == 1 & data4$D == 1, ] %>%
  group_by(iteration) %>%
    summarise(Frequency = sum(Count))

data5$Bool <- data5$Bool <- ifelse(data5$Frequency >= 20, 1, 0)

p.upper <- prop.table(table(data5$Bool))[2]
p.twoside <- 2 * p.upper

# Confidence Intervals ----
n = 1000
j = 20

pop <- data.frame(height= rnorm(n, mean = 65, sd = 5))
pop <- rep_sample_n(pop, 20, replace = FALSE, reps = 10000, prob = NULL)

pop1 <- pop %>% group_by(replicate) %>%
  mutate(mean = mean(height), 
         sd = sd(height), 
         lower = mean - 1.96 * sd/sqrt(20-1),
         upper = mean + 1.96 * sd/sqrt(20-1))

pop1$include <- ifelse(pop1$lower <= 65 & pop1$upper >= 65, 1, 0)
prop.table(table(pop1$include))




