# Set up ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 1")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")

# Simulation from DAGs

# Example 1 (G -> Z -> A -> X, G -> X, Z -> X)----
G = rbern(p = 0.25)
Z = rbern(p = 0.5 + 0.1*G)
A = rbern(p = 0.15 + 0.08*Z)
X = rbern(p = 0.1 + 0.02*A + 0.07*Z + 0.05*G)

dat <- data.frame(G, Z, A, X)

# Crude RDs from linear regression
lm(Z~G, data = dat)$coefficients[2]
lm(A~Z, data = dat)$coefficients[2]
lm(X~Z, data = dat)$coefficients[2]
lm(X~G, data = dat)$coefficients[2]
lm(X~A, data = dat)$coefficients[2]

# Adjusted

# Can't interpret beta of G as causal, since Z is a mediator
# Table 2 fallacy
lm(X~Z + G, data = dat)$coefficients[2]

# Can't interpret beta of Z as causal, since A is a mediator
# Table 2 fallacy
lm(X~A + Z, data = dat)$coefficients[2]

# Can't interpret effects of G and Z as causal, since controlling for mediators in this model 
# Table 2 fallacy
lm(X~A + G + Z, data = dat)$coefficients[2]

# Negative values ----
G = rbern(p = 0.25)
Z = rbern(p = 0.5 + 0.1*G)
A = rbern(p = 0.15 + 0.08*Z)
X = rbern(p = 0.1 - 0.08*A + 0.07*Z + 0.05*G)

dat <- data.frame(G, Z, A, X)

lm(X~Z+G, data = dat)$coefficients[2]

# Dag 1 (C -> E -> M -> Y, C -> M ) ----
C = rbern(10000)
E = rbern(10000, 0.1 + 0.2*C) # RD_CE = 0.2
M = rbern(10000, 0.15 + 0.3*C + 0.2*E) # RD_EM = 0.2, RD_CM = 0.3
Y = rbern(10000, 0.1 + 0.4*M) # RD_MY = 0.4
# Expect E on Y to be 0.2*0.4 = 0.08

dag1 <- data.frame(C, E, M, Y)
lm(Y ~ E, data = dag1)$coefficients[2]
lm(Y ~ E + C, data = dag1)$coefficients[2]
lm(Y ~ E + M, data = dag1)$coefficients[2] # Adjusting for M removes effect of E on Y

# Dag 2 (E -> M -> Y, C -> M, C -> Y) -----
C = rbern(10000)
E = rbern(10000, 0.4)
M = rbern(10000, 0.15 + 0.3*C + 0.2*E)
Y = rbern(10000, 0.1 + 0.4*M + 0.3*C)

dag2 <- data.frame(C, E, M, Y)
lm(Y ~ E, data = dag2)$coefficients[2]
lm(Y ~ E + C, data = dag2)$coefficients[2]
lm(Y ~ E + C + M, data = dag2)$coefficients[2] # Adjusting for M removes effect of E on Y

# DAG 3 Data to DAG to Simulation (not touched on in class)----
E = rbern(10000, 0.2)
Y = rbern(10000, 0.1 + 0.3*E)

dag3 <- data.frame(E, Y)
lm(Y~E, data = dag3)$coefficients[2]

# DAG 4 (C -> E, C -> Y)  ----
C = rbern(0.4)
E = rbern(0.2 + 0.6 * C)
Y = rbern(0.2 + 0.5 * C)

dag4 <- data.frame(C, E, Y)
lm(Y~E, data = dag4)$coefficients[2] # Confounded RD is 0.3
lm(Y~E+C, data = dag4)$coefficients[2] # Unconfounded RD is 0.01

# DAG 5 (C -> E -> Y, C -> Y) ----
C = rbern(0.4)
E = rbern(0.2 + 0.6 * C)
Y = rbern(0.2 + 0.5 * C + 0.3 * E)

dag5 <- data.frame(C, E, Y)
lm(Y~E, data = dag5)$coefficients[2] 
lm(Y~E+C, data = dag5)$coefficients[2] 

# DAG 6 (F -> E -> Y, F -> Y, C -> E, C -> Y) ----
C = rbern(p = 0.4)
F = rbern(p = 0.2)
E = rbern(p = 0.2 + 0.5 * C + 0.2 * F)
Y = rbern( p = 0.2 + 0.5 * C + 0.1 * E + 0.1 * F)

dag6 <- data.frame(C, E, Y, F)
crude <- lm(Y ~ E, data = dag6)$coefficients[2] 
adjusted <- lm(Y ~ E + C + F, data = dag6) $coefficients[2]

# Risk Difference Due to Confounder = 
crude - adjusted

# Optional HW exercise for more complex confounding ----
C = rbern(0.4)
F = rbern(0.2)
E = rbern(0.2 + 0.5 * C + 0.2 * F)
Y = rbern(0.2 + 0.5 * C + 0.1 * E + 0.1 * F)

dag6 <- data.frame(C, E, Y, F)
lm(Y ~ E, data = dag6)$coefficients[2]
lm(Y ~ E + C + F, data = dag6)$coefficients[2] 
lm(Y ~ E + C, data = dag6)$coefficients[2] 
lm(Y ~ E + F, data = dag6)$coefficients[2] 

L = rbern(100000)
C = rbern(0.4 + 0.2 * L)
F = rbern(0.2 + 0.3 * L )
E = rbern(0.2 + 0.5 * C + 0.2 * F)
Y = rbern(0.2 + 0.5 * C + 0.1 * E + 0.1 * F)

dag7 <- data.frame(C, E, Y, F, L)
lm(Y ~ E, data = dag7)$coefficients[2] 
lm(Y ~ E + C, data = dag7)$coefficients[2] 
lm(Y ~ E + C + F, data = dag7)$coefficients[2] 
lm(C ~ F, data = dag7)$coefficients[2]

# Option HW to mess around with strength of relationships in DAG 5 ----
C = rbern(0.04)
E = rbern(0.2 + 0.4 * C)
Y = rbern(0.3 + 0.4 * C + 0.03 * E)

dag5 <- data.frame(C, E, Y)
lm(Y~E, data = dag5)$coefficients[2] 
lm(Y~E+C, data = dag5)$coefficients[2] 

# Mediation Analysis ----
library(mediation)
# !!! This mediation package only works on the linear scale !!! 
# DAG ( A -> M -> Y, A -> Y)
A = rbern(0.3, n=1000)
M = rbern(0.1 + 0.5 * A, n=1000)
Y = rbern(0.1 + 0.3 * A + 0.2 * M, n=1000)

dat <- data.frame(exp = A, med = M, out = Y)

out.fit <- lm(out~exp + med, data = dat)
med.fit <- lm(med~exp, data = dat)

summary(out.fit)
summary(med.fit)

summary(mediate(med.fit, out.fit, treat="exp", mediator="med"))






