library(ggplot2)
setwd("~/Google Drive/My Drive/PhD/Sem 4/Novel Methods/Final Project/")

# QBA Results ----
results <- read.csv("QBA_TOLDLARGE.csv")
results <- rename(results, Lower = Percentile_lower, Upper = Percentile_Upper)
results$label <- paste0("OR = ", round(results$OR, 2))

#create forest plot
ggplot(data=results, aes(y=Model, x=log(OR),
                             xmin=log(Lower), 
                             xmax=log(Upper))) +
  geom_point(colour = "magenta") + 
  geom_text(label = results$label, nudge_y = 0.25) + 
  geom_errorbarh(height=.1) + labs(y = "", x = "logOR") + 
  #scale_y_continuous(labels = results$Analysis) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 11))

# QBA Sens ----
results <- read.csv("QBA_Sens.csv")
results <- rename(results, Lower = Percentile_lower, Upper = Percentile_Upper)
results$label <- paste0("OR = ", round(results$Estimate, 2))
results[5,] <- c(5, "Conventional Analysis", 1.06, 0.69, 1.61, 1.61/0.69, "OR = 1.06")

#create forest plot
ggplot(data=results, aes(y=Analysis, x=log(Estimate),
                         xmin=log(Lower), 
                         xmax=log(Upper))) +
  geom_point(colour = "magenta") + 
  geom_text(label = results$label, nudge_y = 0.25) + 
  geom_errorbarh(height=.1) + labs(y = "", x = "logOR") + 
  #scale_y_continuous(labels = results$Analysis) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 11))