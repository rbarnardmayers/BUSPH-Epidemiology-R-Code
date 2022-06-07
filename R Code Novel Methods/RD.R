# Regression Discontinuity Exercise
setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Session 8")

library(haven)
library(dplyr)
library(ggplot2)
library(rdd)

sard <- read_sas("south_africa_rd_updated.sas7bdat")

# Step 1: Dataset setup ----

sard2 <- sard %>% 
  
  rename(Y1 = single_sub_firstline_nrti, 
         Y2 = death24, 
         Y3 = comp24, 
         Y4 = cd4_change24, 
         Y5 = fail18) %>%
  
  mutate(TRT_rec = case_when(firstnrti == "TDF" ~ 1, TRUE~0), 
         T = case_when(cz > 0 ~ 1,
                       TRUE ~ 0), 
         Z_bin = round(cz/7) * 7) %>% # days / 7 makes it weeks 
  subset(cz > -365.25 & cz < 365.25) %>% 
  mutate(cz_x_t = cz * T) %>%
  
  # IKbandwidth(sard2$cz, sard2$TRT_rec)
  
  # negative cz represents old drug era, positive is new drug era
  group_by(Z_bin) %>%
  
  mutate(trt_mean = mean(TRT_rec), 
         y1_mean = mean(Y1, na.rm = TRUE), 
         y2_mean = mean(Y2, na.rm = TRUE),
         y3_mean = mean(Y3, na.rm = TRUE),
         y4_mean = mean(Y4, na.rm = TRUE),
         y5_mean = mean(Y5, na.rm = TRUE)) 

sard2 <- sard2[order(sard2$Z_bin), ]

#Step 2: Check RD assumptions ----

# histogram 
p <- ggplot(sard2, aes(x = cz)) + 
  geom_histogram(color="darkblue", bins = (365.25*2/7)) +
  geom_vline(xintercept = 0, color = "red") + 
  labs(title = "Histogram displaying number of ART initiations in South Africa (n=6,418) 
before and after the introduction of tenofovir in first-line therapy")
p

# could suggest could be treatment delay in this sample becaues increase in 
# histogram post change in guidelines (bunching is present)

# regressions og covariates ----
age.out <- lm(age1 ~ cz + T + cz_x_t , data = sard2)

# (Intercept)           cz            T       cz_x_t  
# 35.857355    -0.001884    -0.628869     0.004866  

# above threshold average age is 35.857355
# below threshold average age is 35.857355 - 0.628869 = 35.22849 
# p-value is  summary(age.out)$coefficients[3,4]

wght.out <- lm(base_wgt1 ~ cz + T + cz_x_t , data = sard2)

# (Intercept)           cz            T       cz_x_t  
# 64.235230    -0.004697     0.331597     0.009133  

# above threshold average weight is 64.235230
# below threshold average weight is 64.235230 + 0.331597 = 64.56683
# p-value is  summary(wght.out)$coefficients[3,4]

cd4.out <- lm(base_cd4 ~ cz + T + cz_x_t , data = sard2)

# (Intercept)           cz            T       cz_x_t  
# 156.35137      0.04532     12.13061     -0.06049 

# above threshold average cd4 is 156.35137
# below threshold average cd4 is 156.35137 + 12.13061 = 168.482
# p-value is  summary(cd4.out)$coefficients[3,4]

hb.out <- lm(base_hb ~ cz + T + cz_x_t , data = sard2)

# (Intercept)           cz            T       cz_x_t  
# 11.001301    -0.000696     0.089101     0.003034  

# above threshold average hb is 11.001301
# below threshold average hb is 11.001301 + 0.089101 = 11.0904
# p-value is  summary(hb.out)$coefficients[3,4]

male.out <- lm(male ~ cz + T + cz_x_t , data = sard2)

# (Intercept)           cz            T       cz_x_t  
# 0.237474    -0.002409     0.038863     0.002578  

# above threshold proportion of males is 0.237474
# below threshold proportion of males is 0.237474 + 0.038863 = 0.276337
# p-value is  summary(male.out)$coefficients[3,4]

# Regression discontinuity ----

# TRT_rec

linreg_tr.out <- lm(TRT_rec ~ cz + T + cz_x_t, data = sard2)
linreg_y1.out <- lm(Y1 ~ cz + T + cz_x_t, data = sard2)
linreg_y2.out <- lm(Y2 ~ cz + T + cz_x_t, data = sard2)
linreg_y3.out <- lm(Y3 ~ cz + T + cz_x_t, data = sard2)
linreg_y4.out <- lm(Y4 ~ cz + T + cz_x_t, data = sard2)
linreg_y5.out <- lm(Y5 ~ cz + T + cz_x_t, data = sard2)

# summary(linreg2.out)

sard2$pred_tr <- linreg_tr.out$fitted.values
sard2$pred_y1 <- linreg_y1.out$fitted.values
sard2$pred_y2 <- linreg_y2.out$fitted.values
sard2$pred_y3 <- linreg_y3.out$fitted.values

sard2$pred_y4 <- ifelse(is.na(sard2$Y4), NA, linreg_y4.out$fitted.values)

sard2$pred_y5 <- linreg_y5.out$fitted.values

DCdensity(sard2$cz, 0)

# Plotting ----
# TRT_rec

cols1 <- c("#1170AA", "#55AD89")#, "#EF6F6A")

 trtrect.plot <- sard2 %>% ggplot() +
   geom_point(size = 0.2, aes(x = cz, y = pred_tr, color = factor(T), fill = factor(T))) + 
   geom_point(size = 0.5, aes(x = Z_bin, y = trt_mean, color = factor(T), fill = factor(T)))+
   scale_color_manual( values = cols1) +
   geom_vline(xintercept = 0, color = "red") +
   labs(title = "Regression discontinuity showing the probability of receiving tenofovir in South Africa (n=6,418)",
        x = "Days since WHO guidelines change", 
        y = "Proportion initiatied on TDF")
 trtrect.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Y1 - Y5
y.plot <- sard2 %>% ggplot() +
  geom_point(size = 0.2, aes(x = cz, y = pred_y2, color = factor(T), fill = factor(T))) + 
  geom_point(size = 0.5, aes(x = Z_bin, y = y2_mean, color = factor(T), fill = factor(T)))+
  scale_color_manual( values = cols1) +
  geom_vline(xintercept = 0, color = "red") +
  labs(title = "Regression discontinuity showing the probability of y3 in South Africa (n=6,418)",
       x = "Days since WHO guidelines change", 
       y = "Number who initiatied on TDF")

y.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Y4
y.plot <- sard2 %>% subset(!is.na(pred_y4)) %>% ggplot() +
  geom_point(size = 0.2, aes(x = cz, y = pred_y4, color = factor(T), fill = factor(T))) + 
  geom_point(size = 0.5, aes(x = Z_bin, y = y4_mean, color = factor(T), fill = factor(T)))+
  scale_color_manual(values = cols1) +
  geom_vline(xintercept = 0, color = "red") +
  labs(title = "Regression discontinuity showing the probability of y4 in South Africa (n=6,418)",
       x = "Days since WHO guidelines change", 
       y = "Number who initiatied on TDF")

y.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))


# CACE 
 # y2 = -2.725e-02/0.70 
# 0.7 is coming from uptake at the threshold 
# attirition = 






