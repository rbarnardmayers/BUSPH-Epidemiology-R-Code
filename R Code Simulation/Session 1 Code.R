# Set working dictionary and load libraries ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 1")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")

# Setting Seeds ----
# set.seed resets after the subsequent line
set.seed(1348798);
a = rnorm(10) 
set.seed(1348798);
b = rnorm(10) 
ifelse(a == b,"yes","no") # a and b will be the same

set.seed(1348798);
a = rnorm(10)
b = rnorm(10) 
ifelse(a == b,"yes","no") # a and b will NOT be the same

# Random Distributions with Histograms (Default options)----
# Uniform distribution
X <- runif(n= 100000) # min = 1, max = 1

# Triangular distribution
Y <- rtri(n=100000) # min = 0, max = 1, mode = 0.5

# Normal distribution
Z <- rnorm(n=100000) # mean = 0, sd = 1

# Bernoulli distribution
A <- as.numeric(rbernoulli(n=100000)) # p = 0.5

# Set plot margins so that all four histograms show in one window
par(mfrow=c(2,2))
hist(X)
hist(Y)
hist(Z)
hist(A)
par(mfrow=c(1,1))

# Random Distributions with Histograms and changing parameters (Option 1)----
X <- 25 + 6*runif(10000)  # min = 25, max = 31
Y <- 66 + 22*rtri(100000, mode = 0.2) # min = 66, max = 88, mode = 0.2
Z <- 800 + 50*rnorm(10000) # mean is 800, sd is 50 

# Set plot margins so that all three histograms show in one window
par(mfrow=c(2,2))
hist(X)
hist(Y)
hist(Z)
par(mfrow=c(1,1))

# Random Distributions with Histograms and changing parameters (Option 2)----
X <- runif(10000, min = 25, max = 31)  # min = 25, max = 31
Y <- rtri(100000, min = 66, max = 88, mode = 0.2) # min = 66, max = 88, mode = 0.2
Z <- rnorm(10000, mean = 800, sd = 50) # mean is 800, sd is 50 

# Set plot margins so that all three histograms show in one window
par(mfrow=c(2,2))
hist(X)
hist(Y)
hist(Z)
par(mfrow=c(1,1))

# Simulating dataset (no confounder)----
E = as.numeric(rbernoulli(10000))
D = as.numeric(rbernoulli(10000, 0.25+0.25*Exp))

first <- data.frame(E, D)
twobytwo.strat(first)
# Simulating dataset with confounder (no EMM on difference scale)----
Conf = as.numeric(rbernoulli(10000))
Exp = as.numeric(rbernoulli(10000, 0.25 + 0.25*Exp))
Dis = as.numeric(rbernoulli(10000, 0.05 + 0.2*Exp + 0.05*Conf)) 

first <- data.frame(Conf, Exp, Dis)

# Probabilities to Regression ----
A = rbernoulli(100000)
B = rbernoulli(100000)

X = rbernoulli(100000, p = exp(log(0.1) + log(0.3/0.1)*A + log(0.4/0.1)*B + log(8/12)*A*B))
Y = rbernoulli(100000, 0.1 + 0.2*A+0.3*B+0.2*A*B)

# Regressions ----
# can add this to the end of any lines below to get the Exp coefficient only $coefficients[2]
lm(Dis~Exp, data=first)
lm(Dis~Exp + Conf, data=first)

# can add this to the end of any lines below to get the Exp coefficient only $coefficients[2]
# need to exponentiate to translate to OR or RR - exp()
glm(Dis~Exp, data=first, family=poisson(link = "log")) # For risk ratios
glm(Dis~Exp, data=first, family=binomial(link = "logit")) # For risk ratios

# Create EMM dataset from Excel and Simulate ----
# no EMM in RD, EMM in RR, confounding is present
Conf = as.numeric(rbernoulli(20000))
Exp = as.numeric(rbernoulli(20000, p = 0.3 + 0.4*Conf))
Dis1 = as.numeric(rbernoulli(20000, p = 0.2 + 0.2*Exp + 0.1*Conf))
Dis2 = as.numeric(rbernoulli(20000, p = exp(log(0.2) + log(2)*Exp + log(1.5)*Conf + log(2.5/(1.5*2))*Conf*Exp)))

# Create data frame 
test <- data.frame(Conf, Exp, Dis1, Dis2)

# run the models to get RD and RR estimates
lm(Dis1~Exp+Conf, data=test)$coefficients[2] # linear
logistic <- glm(Dis2~Exp*Conf, data=test,family=poisson(link = "log")) # log link 
exp(summary(logistic)$coefficients[2])

#no EMM in RR, EMM in RD, confounding is present
Conf = as.numeric(rbernoulli(20000))
Exp = as.numeric(rbernoulli(20000, p = 0.3 + 0.4*Conf))
Dis1 = as.numeric(rbernoulli(20000, p = 0.2 + 0.2*Exp + 0.1*Conf))
Dis2 = as.numeric(rbernoulli(20000, p = exp(log(0.2) + log(2)*Exp + log(1.5)*Conf)))

test <- data.frame(Conf, Exp, Dis1, Dis2)

lm(Dis1~Exp+Conf, data=test)$coefficients[2]
logistic <- glm(Dis2~Exp*Conf, data=test,family=poisson(link = "log"))
exp(summary(logistic)$coefficients[2])

#EMM in RD and EMM in RR, confounding is present
Conf = as.numeric(rbernoulli(20000))
Exp = as.numeric(rbernoulli(20000, p = 0.3 + 0.4*Conf))
Dis1 = as.numeric(rbernoulli(20000, p = 0.2 + 0.3*Exp + 0.1*Conf + 0.1*Exp*Conf))
Dis2 = as.numeric(rbernoulli(20000, p = exp(log(0.2) + log(2.5)*Exp + log(1.5)*Conf + log(3.5/(1.5*2.5))*Exp*Conf)))

test <- data.frame(Conf, Exp, Dis1, Dis2)

lm(Dis1~Exp*Conf, data=test)$coefficients[2]
logistic <- glm(Dis2~Exp*Conf, data=test,family=poisson(link = "log"))
exp(summary(logistic)$coefficients[2])



