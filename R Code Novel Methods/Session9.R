# Set up ----
library(dplyr)
library(ivpack)
library(haven)

i = 10000

# Create Data Frame ----
U = rnorm(i)
Z = rnorm(i)
W = rnorm(i, mean = 0, sd = 10)
E= rnorm(i, mean = 0, sd = 10)
X = 2 - 13 * U + 10 * Z + W
Y = 5 + 4 * X - .5 * U + E

one <- data.frame(U, Z, X, Y)
Sarah <- data.frame(Z, X, Y)

write.csv(Sarah, "Sarah.csv")
# Run Regression ----
 
log.out1 <- lm(Y~X, data = one) # 0.8381 
log.out2 <- lm(Y~X + U, data = one) # 0.4873 

# Do IV Analysis----

log.out3 <- lm(x~z, data = data)

data$fitted <- log.out3$fitted.values

log.out4 <- lm(y ~ fitted, data = data) # 8.589

# ITT 
 itt.out <- lm(Y~Z, data = one)
 
# First Stage
 
 first.out <- lm(X~Z, data = one)
 
 beta = summary(itt.out)$coefficients[2,1] / summary(first.out)$coefficients[2,1]








