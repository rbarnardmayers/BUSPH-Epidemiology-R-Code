library(haven)
library(MatchIt)
library(cobalt)
library(optmatch)
library(dplyr)

setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 6 - Propensity Scores I")
# setwd("G:/My Drive/PhD/Sem 4/Novel Methods/Session 6")
lalonde <- read_sas("lalonde.sas7bdat")


# Crude analysis ----
crude.out <- lm(re78 ~ treated, data = lalonde)
adj.out <- lm(re78 ~ treated + education + black  + married  + 
                hispanic  + u74 +  u75  + age  + re74 + 
                re75, data = lalonde)

# propensity model ----
prop.out <- glm(treated ~ education + black  + married  + 
                hispanic  + u74 +  u75  + age  + re74 + 
                re75, data = lalonde)

lalonde$pscore_1 <- prop.out$fitted.values

# Trimming ----
trim1 <- quantile(lalonde[lalonde$treated == 1, ]$pscore_1, 0.975)
trim2 <- quantile(lalonde[lalonde$treated == 0, ]$pscore_1, 0.025)

lalonde1 <- lalonde %>% mutate(predvals_cat = case_when(pscore_1 > trim1 | pscore_1 < trim2  ~ 1, 
                                               TRUE ~ 0))

lalonde2 <- subset(lalonde1, lalonde1$predvals_cat == 0 )

# Deciles ----

lalonde2$decile <- ntile(lalonde2$pscore_1, 10)

decile.out <- lm(re78 ~ treated + decile, data = lalonde2)

# PS Matched ----

m.out <- matchit(treated ~  education + black + married  + 
                   hispanic  + u74 +  u75  + age  + re74 + 
                   re75, data=lalonde, method = "exact")

m.out <- matchit(treated ~  education + black + married  + 
                   hispanic  + u74 +  u75  + age  + re74 + 
                   re75, data=lalonde, method = "nearest", caliper = 0.1)

# summary(m.out)
# plot(m.out, type = "jitter")
# plot(m.out, type = "hist")

lalonde$pscore_2 <- m.out$distance

second.out <- lm(re78~ treated, data = m.data)

m.data <- match.data(m.out, distance = 'prop.score')

bal.plot(m.out, var.name = "re75", 
         which = 'both', grid = TRUE)

View(m.data)

# Variable Ratio Matching ----
m.out2 <- matchit(treated ~  education + black + married  + 
                    hispanic  + u74 +  u75  + age  + re74 + 
                    re75, data=lalonde, method = "nearest", ratio = 4)

m.data2 <- match.data(m.out2, distance = 'prop.score')
merge.2 <- data.frame(table(m.data2$subclass))
colnames(merge.2) <- c("subclass", "Freq")

m.data2 <- merge(m.data2, merge.2, by = "subclass")

matched.out2 <- lm(re78~treated + Freq, data = m.data2)

exp(summary(matched.out2)$coefficients)
exp(confint(matched.out2))



