# Set up ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")
library(table1)

# twobytwo.strat(dag1$E, dag1$Y, dag1$C)

# controlling for parent/child of confounder ----
M = rbern(0.4)
C = rbern(0.01+0.95*M)
E = rbern( 0.1 + 0.7*C)
Y = rbern( 0.1 + 0.7*C)
Z = rbern(0.1 + 0.7*C)

dag1 <- data.frame(M, C, E, Y, Z)
lm(Y~E, data=dag1) # biased (0.5)
lm(Y~E + M, data=dag1) # does a little to mitigate the bias (0.5 -> 0.34)
# Weaker RD_MC means less reduction in bias correction (0.5 -> 0.4)
lm(Y~E + Z, data=dag1) # does a little to mitigate the bias (0.5 -> 0.33)
# Weaker RD_ZC means less reduction in bias correction (0.5 -> 0.4), opposite true as well
lm(Y~E + C, data=dag1) # does a lot to mitigate the bias (0.5 -> 0)

# controlling for parent/child of collider ----
M = rbern( 0.4)
E = rbern( 0.4)
Y = rbern( 0.4)
C = rbern( 0.1 + 0.3 * M + 0.3 * E + 0.3 * Y)
Z = rbern( 0.1 + 0.7 * C)

dag2 <- data.frame(M, C, E, Y, Z)
lm(Y~E, data=dag2) # Unbiased (0)
lm(Y~E + C, data=dag2) # biased (0 -> -0.1) 
lm(Y~E + M, data=dag2) # no effect
lm(Y~E + Z, data=dag2) # biased (0 -> -0.05) weakening RDZC = less bias (0 -> -0.002), vice versa

twobytwo.strat(dag2, out = "Y", strat = "C")

# Instrumental variables ----

C = rbern(0.25)
I = rbern(0.25)
M = rbern(0.1 + 0.5 * I + 0.3 * C)
Y = rbern(0.1 + 0.4 * M + 0.4 * C)

dag3 <- data.frame(C, I , M, Y)
lm(Y~M, data = dag3) # Crude = 0.509
lm(Y~M + C, data = dag3) # Adjusted = 0.40130 
lm(Y~I, data = dag3) # Estimate = 0.2011
lm(M~I, data = dag3) # Estimate = 0.5008
# Instrumental Estimate =  0.2011/0.5008 = 0.4

I = rbern(0.25)
X = rbern(0.1 + 0.6 * I)
M = rbern(0.1 + 0.5 * I)
Y = rbern(0.1 + 0.4 * M)

dag4 <- data.frame(I, X, M, Y)
lm(Y~M, data = dag4) # Crude = 0.3976

lm(Y~I, data = dag4) # Estimate = 0.1986
lm(M~I, data = dag4) # Estimate = 0.5003
# Instrumental Estimate with I =  0.1986/0.5003 = 0.4

lm(X~I, data = dag4) # Estimate = 0.60096
lm(M~I, data = dag4) # Estimate = 0.49865
lm(Y~X, data = dag4) # Estimate = 0.1188
# Instrumental Estimate with X = 0.1188/(0.60096*0.49865)  = 0.3964378

# Instrumental with violations ----

I = rbern(0.25)
M = rbern(0.1 + 0.1 * I) 
Y = rbern(0.1 + 0.04 * I + 0.4 * M)
Y2 = rbern(0.1 + 0.1 * I + 0.4 * M)

dag5 <- data.frame(I, M, Y, Y2)
lm(Y~I, data=dag5) # with RD_IM = 0.1 = 0.4391, with RD_IM = 0.7 = 0.4391
lm(M~I, data=dag5) # with RD_IM = 0.1 = 0.09975 , with RD_IM = 0.7 = 0.7011
lm(Y~M, data=dag5) # with RD_IM = 0.1 = 0.4679 , with RD_IM = 0.7 = 0.2888
lm(Y ~ M + I, data=dag5) # with RD_IM = 0.1 = 0.4679 , with RD_IM = 0.7 = 0.2888
# Want to recover M on Y with instrumental I 
# with RD_IM low: 0.4391/0.09975 = 4.4
# with RD_IM high: 0.4391/0.7011 = 0.972258 

#Reduce RD_IY
lm(Y2~I, data=dag5) # with RD_IM = 0.1 = 0.000544, with RD_IM = 0.7 = 0.3794
lm(M~I, data=dag5) # with RD_IM = 0.1 = 0.10009 , with RD_IM = 0.7 = 0.7001
# Instrumental Analysis = 

# low RD_IM = 0.000544/0.10009 =0.005435108
# high RD_IM = 0.3794/0.7001 = 0.5419226

# Instrumental with common between exposure and Y  -----

C = rbern(0.25)
I = rbern(0.5)
M = rbern(0.3 + 0.5 * I + 0.2 * C)
Y = rbern(0.1 + 0.4 * M + 0.3 * C)

dag6 <- data.frame(C, I , M, Y)
lm(Y~M, data = dag6) 
lm(Y~M + I, data = dag6) 
lm(Y~M + I + C, data = dag6) 
lm(Y~I, data = dag6) 
lm(M~I, data = dag6) 
# better off with crude , adjusted for I, or I as IV when C is unmeasured? 

# RD(CM) - 
# Crude M-Y = 0.3554
# Adjusted for I, more bias = 0.34135
# Instrumental Estimate =  0.1992/0.5011 = 0.3975254

# RD(CM) +
# Crude = 0.4467
# Adjusted for I, more bias  = 0.46324
# IV Analysis =  0.1998/0.4982 = 0.4010438

# RD(IM) weak
# Crude = 0.4688
# IV = 0.3970588

# Homework (common cause of I and Y w/o direct effect) ----

X = rbern(0.6)
I = rbern(0.45 + 0.6 * X)
C = rbern(0.2)
M = rbern(0.2 + 0.3 * I + 0.3 * C) 
Y = rbern(0.05 + 0.2 * X + 0.3 * M + 0.4*C)

dagHW <- data.frame(I, M, Y, X, C)
a <- lm(Y~I, data=dagHW)$coefficients[2] #biased
b <- lm(Y~I + X, data=dagHW)$coefficients[2]# unbiased
c <- lm(M~I, data=dagHW)$coefficients[2] # biased
#How much bias is present in crude 
a /  c 

#change assumptions

# Optional Homework (Front door approach) ----
# Multiply ZY|X and XZ to get X on Y 
U = rbern(0.3)
X = rbern(0.25 + 0.3 * U)
Z = rbern(0.1 + 0.5 * X)
Y = rbern(0.1 + 0.3 * U + 0.1 * Z)

dagHW2 <- data.frame(U, X, Z, Y)
lm(Y~Z + X,  data=dagHW2) # 0.10150
lm(Z~X, data=dagHW2)  #  0.5008
lm(Y~X + U, data = dagHW2) #  0.05177











