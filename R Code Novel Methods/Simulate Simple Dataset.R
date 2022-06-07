library(dplyr)
library(survey)
library(purrr)
library(EnvStats)

# First go round (not relevant) ----
first <- data.frame(ID = c(1,2,3,4), Exp = c(1,1,0,0), Out = c(1,0,1,0), Count=c(25,75,50,50))


svytable(Out ~ Exp, design=wtdfirst)
wtdfirst <- svydesign(id= ID, weights=Count, data=first)

k=100
a=25
b=75
c=50
d=50

second <- data.frame(id=1:k, Exp=c(rep(1, a+c), rep(0, b+d)), Out = c(rep(1,a), rep(0,c), rep(1,b), rep(0,d)))

# Probabilistic Simulation (Not relevant)----

set.seed(100); seed1 <- data.frame(id = 1:k, fixed_n=runif(k))
set.seed(NULL); seed2 <- data.frame(id = 1:k, fixed_n=runif(k))
 
test <- data.frame(id = 1:k, X=runif(k), Y = rtri(k), Z=rnorm(k), A=rbinom(k, 1,0.8))

test <- data.frame(id = 1:k, X=25+6*runif(k), Y =66+22*rtri(k), Z=800+50*rnorm(k))


first <- data.frame(id = 1:k, exp=rbinom(k,1,0.5))
first$dis <- as.numeric(ifelse(first$exp == 1, rbernoulli(1, p=0.25), rbernoulli(1)))
                            
# Corresponds to DAG with C -> E -> D and C -> D ----
# pr(C) = 0.5 
# pr(E|C-)=0.25
# pr(D|E-, C-) = 0.05
# RD_EC = 0.25
# RD_DC = 0.15
# RD_ED = 0.20

k=100000
c=as.numeric(rbernoulli(k, 0.5))
e=as.numeric(rbernoulli(k, 0.25+0.25*c))
d=as.numeric(rbernoulli(k, 0.05+0.20*e+0.15*c))

dat <- as.data.frame(cbind(id=1:k, e, d, c))

# Linear Regression ---
summary(lm(dat$d~dat$e)) 
summary(lm(dat$d~dat$e+dat$c)) 


# Simulate a dataset with detail from slide 49 (in class exercise) ----
# Corresponds to DAG with M-bias
# pr(a) = 0.3
# pr(b) = 0.3
# pr(e | a-) = 0.8, RD_EA = -0.60
# pr(d | b-) = 0.10, RD_DB = 0.65
# pr(m | a-, b-) = 0.20, RD_AM = 0.6, RD_BM = -0.15

k=10000
a=as.numeric(rbernoulli(k, 0.3))
b=as.numeric(rbernoulli(k, 0.3))
e=as.numeric(rbernoulli(k, 0.80 - 0.6*a))
d=as.numeric(rbernoulli(k, .10+0.65*b))
m=as.numeric(rbernoulli(k, .20+0.60*a-0.15*b))

dat <- as.data.frame(cbind(id=1:k, a, b, e, d, m))

# Linear Regression ---
summary(lm(dat$d ~ dat$e)) 
summary(lm(dat$d ~ dat$e + dat$m))


# As a function ----

sim_mbias <- function(
  # Corresponds to DAG with M-bias
  pr_a = 0.3,
  pr_b = 0.3,
  pr_e_a0 = 0.8,
  RD_EA = -0.60,
  pr_d_b0  = 0.10,
  RD_DB = 0.65,
  pr_m_a0b0 = 0.20,
  RD_AM = 0.6,
  RD_BM = -0.15){
  
  k=10000
  a=as.numeric(rbernoulli(k, pr_a))
  b=as.numeric(rbernoulli(k, pr_b))
  e=as.numeric(rbernoulli(k, pr_e_a0 + RD_EA*a))
  d=as.numeric(rbernoulli(k, pr_d_b0+RD_DB*b))
  m=as.numeric(rbernoulli(k, pr_m_a0b0+RD_AM*a+RD_BM*b))
  
  dat <- as.data.frame(cbind(id=1:k, a, b, e, d, m))
  
  # Linear Regression ---
  crude_est <- summary(lm(dat$d ~ dat$e))$coefficients[2,1]
  adj_est <- summary(lm(dat$d ~ dat$e + dat$m))$coefficients[2,1]

  print(paste0("Crude Estimate: ", crude_est))
  print(paste0("Adjusted Estimate: ", adj_est))
  print(paste0("RRc: ", crude_est/adj_est))
  
  
}

# Results ----
# > sim_mbias(RD_AM = -0.5, RD_BM = 0.5, RD_EA = -0.6, RD_DB = -0.6)
# [1] "Crude Estimate: 0.0107231098380173"
# [1] "Adjusted Estimate: 0.0191658143678722"
# [1] "RRc: 0.559491479579003"



