# IPW in R ----
setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 11 - MSMs I")
source("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2/Functions.R")

# load libraries ----
library(dplyr)
library(geepack)

# read in data ----
# fram.data <- read.csv("frmgham2.csv")
C0 = rbern(0.2, n = 4)
E0 = rbern(0.1 + 0.2 * C0, n = 4)
C1 = rbern(0.2 + 0.7 * C0 + 0.3 * E0, n = 4)
E1 = rbern(0.1 + 0.2 * C1 + 0.4 * E0, n = 4)
C2 = rbern(0.2 + 0.4 * C1 + 0.3 * E1 , n = 4) 
E2= rbern(0.1 + 0.2 * C2 + 0.3 * E1, n = 4)

Y = rbern(0.05 + 0.08 * C0 + 0.07 * E0 + 0.07 * C1 + 0.1 * C2 + 0.2* E1 + 0.2 * E2, n = 4)
  
conf = c(C0, C1, C2)
exp = c(E0, E1, E2)
out = rep(Y, 3)

data <- data.frame(ID = rep(1:4, 3), Visit = rep(0:2, each = 4), Exp = exp, Out = out, C = conf)

