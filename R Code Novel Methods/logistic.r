model {
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    p[i]<-1/(1+exp(-(b0 + b1*x[i])))
  }
  b0 ~ dnorm( 0 , 1.0E-6 )
  b1 ~ dnorm(log(1.133) , 1/(((log(1.173) - log(1.095))/(2*1.96))^2))
}

