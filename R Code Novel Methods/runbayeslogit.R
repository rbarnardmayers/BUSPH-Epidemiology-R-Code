#Load JAGS library
library(rjags)

#Create Data
case=c(rep(1,48),rep(1,17),rep(0,32),rep(0,37))
exposed=c(rep(1,48),rep(0,17),rep(1,32),rep(0,37))
freq.model<-glm(case~exposed,family=binomial)

summary(freq.model)
confint(freq.model)


#Run Bayes Model

# tell JAGS what variables to monitor
parameters<-c("b0","b1")

# get all the data in one place
data.for.jags<-list('x'=exposed,'y'=case,'N'=length('x'))

#compile the model
jags.compiled<- jags.model(data = data.for.jags,
                           file = "/Volumes/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 3 - Bayes/logistic.r")

jags.output<-jags.samples(model=jags.compiled,variable.names = parameters,n.iter=10000)
#get pt est and CI
exp(quantile(jags.output$b1,c(.025,.5,.975)))
#trace plot
plot(jags.output$b1)


