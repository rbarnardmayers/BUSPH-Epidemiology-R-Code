library(dplyr)

i = 10000
dist <- data.frame(Obs=1:i,
                   X=rbinom(n=i,size=100,prob=0.5),
                   Y=rbeta(i,10,90),
                   Z=runif(i, min=0.8, max=1))
par(mfrow=c(2,2))
hist(dist$Z)
hist(dist$Y)
hist(dist$X)


j = 10000
test <- data.frame(Obs = 1:j,
                   X=runif(i),
                   Y=runif(i, min=0.5, max=0.7))
par(mfrow=c(1,2))
hist(test$X, breaks=seq(0,1,by=0.02))
hist(test$Y, breaks=20)

k = 10000
a = 45 
b = 94
c = 257
d = 945
randomerr <- data.frame(Obs=1:k,
                        stderr = rep(sqrt(1/a-1/c+1/b-1/d), k),
                        Y = rnorm(k, mean=sqrt(1/a-1/c+1/b-1/d)))

hist(randomerr$Y)
summary(randomerr$Y)

par(mfrow=c(1,1))




