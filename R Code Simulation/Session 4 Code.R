# Set up ----- 
setwd("~/Google Drive/My Drive/PhD/Sem 4/Simulation/Day 2")
source("Functions.R")

# Case Control Sampling ----
m = 40000
simdat <- data.frame(C = c(rep(1, m/2), rep(0, m/2)), 
                   E = c(rep(1, m/4), rep(0, m/4), rep(1, m/4), rep(0,m/4)),
                   D = c(rep(1, 4000), rep(0, 6000), rep(1, 2000), rep(0, 8000), 
                         rep(1, 400), rep(0, 9600), rep(1, 200), rep(0, 9800)))

twobytwo.strat(simdat, strat = "C")

# logit link
glm(D ~ E * C, data=simdat, family = binomial(link = "logit"))

# log link
glm(D ~ E * C, data=simdat, family = binomial(link="log"))

# Select only cases from simdat
simdat_case <- simdat[simdat$D == 1,]
simdat_case$Case <- 1 

# Scenario 1: Case-control sampling with EMM on OR scale  ----
simdat_control <- simdat[simdat$D == 0,]
simdat_control$Case <- 0
simdat_control <- simdat_control %>% sample_frac(.1)

simdat_new_1 <- rbind(simdat_case, simdat_control)
glm(Case ~ E * C, data=simdat_new_1, family = binomial(link = "logit"))

# Scenario 2: Case Cohort sampling with EMM on OR scale ----
simdat_control <- simdat %>% sample_frac(.1)
simdat_control$Case <- 0

simdat_new_2 <- rbind(simdat_case, simdat_control)
glm(Case ~ E * C, data=simdat_new_2, family = binomial(link = "logit"))

# Looping ----

all_data <- data.frame()

for(p_A in seq(0.1, 0.8, by=0.1)){
  for(e_CE in seq(.1, .4, by=0.1)){
    C = rbern(p_A)
    E = rbern(0.2 + e_CE)
    D = rbern(-.1 + 0.1*E + 0.3*C)
    if(is_empty(all_data)){
    all_data <- data.frame(C, E, D, prev = c(rep(p_A, length(C))), effect = c(rep(e_CE, length(C))))} else{
      all_data <- rbind(all_data, data.frame(C, E, D, prev = c(rep(p_A, length(C))), effect = c(rep(e_CE, length(C)))))
    }
  }
}

# Misclassification (NDEM) ----
se_D = 0.70
sp_D = 0.95
se_UD = se_D
sp_UD = sp_D
E <- rbern(0.3)
D <- rbern(0.1 + 0.5 * E)
# Two options for observed E: 
# Option 1) 
E_obs <- ifelse(E == 1, rbern(se_D), 
                     ifelse(E == 0, rbern(1-sp_D), NA))
# Option 2) <-- better
E_obs_2  = rbern((1-sp_D) + (se_D - (1-sp_D)) * E)
misc <- data.frame(E, D, E_obs, E_obs_2)

twobytwo.strat(misc, out="E_obs_2", risk=FALSE)

glm(D~E, data = misc) # Estimate =  0.4999
glm(D~E_obs, data = misc) # Estimate = 0.3691 
glm(D~E_obs_2, data = misc) # Estimate = 0.3685 

(se_D - (1-sp_D))

# Misclassification (DEM) ----
se_D = 0.70
sp_D = 0.95
se_UD = 0.65
sp_UD = 0.97

E = rbern(0.5)
D = rbern(0.1 + 0.3 * E)
E_obs = rbern((1-se_UD) + ((sp_D - (1-sp_UD)) * D) + ((se_UD - (1-sp_UD)) * E) + 
                ((se_D - (se_UD - (1-sp_UD)) -  (1-sp_UD)) * E * D))

misc2 <- data.frame(E, D, E_obs)

lm(D~E, data = misc2)
lm(D~E_obs, data = misc2)

