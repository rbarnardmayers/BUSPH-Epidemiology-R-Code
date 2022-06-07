library(dplyr)
library(survey)
library(purrr)
library(EnvStats)
library(epiR)


rbern <- function(p = 0.5, n = 1000000){
  require(purrr)
  # Function that returns a vector of size n (default is 1,000,000) 
  # each number corresponds to a bernoulli trial using probability p  (defualty is 0.5)
  return(as.numeric(rbernoulli(n, p)))
}

twobytwo.strat <- function(info, exposure = "E", out = "D", strat = NULL, risk = TRUE){
  # Creates a two by two table of exposure and out (out as rows, exposure as cols)
  # assumse dichotomous variables only
  # exposure is the name of your exposure column (default is E)
  # out is the name of your outcome column (default is D)
  # strat is the name of your stratifying variable (default is none)
  # if strat is given (another column name in info dataset), stratified tables will print as well
  # if risk = TRUE, risks, RRs, RDs, and ORs will be printed by strata
  # if risk = FALSE, no risk calculations will be shown
  # variables must be binary for all variables
  options(scipen = 100)
  
  if(is.null(strat)){ # if no stratification variable 
    dat = data.frame(Exposure = info[[exposure]], Outcome = info[[out]], 
                     Count = rep(1, length(info[[out]]))) 
    dat2 <- data.frame("a" = sum(dat[dat$Exposure == 1 & dat$Outcome == 1, ]$Count),
                       "b" = sum(dat[dat$Exposure == 0 & dat$Outcome == 1, ]$Count),
                       "c" = sum(dat[dat$Exposure == 1 & dat$Outcome == 0, ]$Count),
                       "d" = sum(dat[dat$Exposure == 0 & dat$Outcome == 0, ]$Count))
    dat2$Total_e <- dat2$a + dat2$c
    dat2$Total_ue <- dat2$b + dat2$d
    dat2$Total_d <- dat2$a + dat2$b
    dat2$Total_ud <- dat2$c + dat2$d
    dat2$Total <- sum(dat2$Total_ud +  dat2$Total_d)
    
    if(risk){
    dat2$Risk_e <- dat2$a / dat2$Total_e
    dat2$Risk_ue <- dat2$b / dat2$Total_ue
    dat2$RD <- dat2$Risk_e - dat2$Risk_ue
    dat2$RR <- dat2$Risk_e/dat2$Risk_ue
    dat2$OR <- (dat2$a * dat2$d)/(dat2$b * dat2$c)
    
    dimnames = list(c("D+", "D-", "Total", "Risk", "RD", "RR", "OR"), c("E+", "E-", "Total"))
    
    print("Overall")
    print(matrix(c(dat2$a, dat2$c, dat2$Total_e, dat2$Risk_e, dat2$RD, 
                   dat2$RR, round(dat2$OR, 2), dat2$b, dat2$d, dat2$Total_ue, dat2$Risk_ue, 
                   "ref", "ref", "ref", dat2$Total_d, dat2$Total_ud, dat2$Total, "", "", "", ""),
                 nrow = 7, dimnames = dimnames), quote=FALSE)
    
    } else{
      dimnames = list(c("D+", "D-", "Total"), c("E+", "E-", "Total"))
      
      print("Overall")
      print(matrix(c(dat2$a, dat2$c, dat2$Total_e, dat2$b, dat2$d, dat2$Total_ue, dat2$Total_d, dat2$Total_ud, dat2$Total), nrow = 3, dimnames = dimnames), quote=FALSE)
    }
    
  } else { # if there is stratification by a variable
    dat = data.frame(Exposure = info[[exposure]], Outcome = info[[out]], 
                     Stratify = info[[strat]], Count = rep(1, length(info[[out]])))
    dat2 <- data.frame("a" = sum(dat[dat$Exposure == 1 & dat$Outcome == 1, ]$Count),
                       "b" = sum(dat[dat$Exposure == 0 & dat$Outcome == 1, ]$Count),
                       "c" = sum(dat[dat$Exposure == 1 & dat$Outcome == 0, ]$Count),
                       "d" = sum(dat[dat$Exposure == 0 & dat$Outcome == 0, ]$Count),
                       "a_1" = sum(dat[dat$Exposure == 1 & dat$Outcome == 1 & dat$Stratify == 1, ]$Count),
                       "b_1" = sum(dat[dat$Exposure == 0 & dat$Outcome == 1 & dat$Stratify == 1, ]$Count),
                       "c_1" = sum(dat[dat$Exposure == 1 & dat$Outcome == 0 & dat$Stratify == 1, ]$Count),
                       "d_1" = sum(dat[dat$Exposure == 0 & dat$Outcome == 0 & dat$Stratify == 1, ]$Count),
                       "a_0" = sum(dat[dat$Exposure == 1 & dat$Outcome == 1 & dat$Stratify == 0, ]$Count),
                       "b_0" = sum(dat[dat$Exposure == 0 & dat$Outcome == 1 & dat$Stratify == 0, ]$Count),
                       "c_0" = sum(dat[dat$Exposure == 1 & dat$Outcome == 0 & dat$Stratify == 0, ]$Count),
                       "d_0" = sum(dat[dat$Exposure == 0 & dat$Outcome == 0 & dat$Stratify == 0, ]$Count))
    dat2$Total_e <- dat2$a + dat2$c
    dat2$Total_ue <- dat2$b + dat2$d
    dat2$Total_d <- dat2$a + dat2$b
    dat2$Total_ud <- dat2$c + dat2$d
    dat2$Total <- sum(dat2$Total_ud +  dat2$Total_d)
    
    dat2$Total_1e <- dat2$a_1 + dat2$c_1
    dat2$Total_1ue <- dat2$b_1 + dat2$d_1
    dat2$Total_0e <- dat2$a_0 + dat2$c_0
    dat2$Total_0ue <- dat2$b_0 + dat2$d_0
    
    dat2$Total_1d <- dat2$a_1 + dat2$b_1
    dat2$Total_1ud <- dat2$c_1 + dat2$d_1
    dat2$Total_1 <- sum(dat2$Total_1ud +  dat2$Total_1d)
    dat2$Total_0d <- dat2$a_0 + dat2$b_0
    dat2$Total_0ud <- dat2$c_0 + dat2$d_0
    dat2$Total_0 <- sum(dat2$Total_0ud +  dat2$Total_0d)
    
    if(risk){
    dat2$Risk_e <- dat2$a / dat2$Total_e
    dat2$Risk_ue <- dat2$b / dat2$Total_ue
    dat2$Risk <- dat2$Total_d / dat2$Total
    dat2$RD <- dat2$Risk_e - dat2$Risk_ue
    dat2$RR <- dat2$Risk_e/dat2$Risk_ue
    dat2$OR <- (dat2$a * dat2$d)/(dat2$b * dat2$c)
    
    dat2$Risk_1e <- dat2$a_1 / dat2$Total_1e
    dat2$Risk_1ue <- dat2$b_1 / dat2$Total_1ue
    dat2$Risk1 <- dat2$Total_1d / dat2$Total_1ud
    
    dat2$Risk_0e <- dat2$a_0 / dat2$Total_0e
    dat2$Risk_0ue <- dat2$b_0 / dat2$Total_0ue
    dat2$Risk0 <- dat2$Total_0d / dat2$Total_0ud
    
    
    dat2$RD_1 <- dat2$Risk_1e - dat2$Risk_1ue
    dat2$RD_0 <- dat2$Risk_0e - dat2$Risk_0ue
    dat2$RR_1 <- dat2$Risk_1e/dat2$Risk_1ue
    dat2$RR_0 <- dat2$Risk_0e/dat2$Risk_0ue
    dat2$OR_1 <- (dat2$a_1 * dat2$d_1)/(dat2$b_1 * dat2$c_1)
    dat2$OR_0 <- (dat2$a_0 * dat2$d_0)/(dat2$b_0 * dat2$c_0)
    
    dimnames = list(c("D+", "D-", "Total", "Risk", "RD", "RR", "OR"), c("E+", "E-", "Total"))
    
    print("Overall")
    print(matrix(c(dat2$a, dat2$c, dat2$Total_e, dat2$Risk_e, dat2$RD, 
                   dat2$RR, round(dat2$OR, 2), dat2$b, dat2$d, dat2$Total_ue, dat2$Risk_ue, 
                   "ref", "ref", "ref", dat2$Total_d, dat2$Total_ud, dat2$Total, dat2$Risk,
                   "", "", ""),
                 nrow = 7, dimnames = dimnames), quote=FALSE)
    
    print("Stratify = 1")
    print(matrix(c(dat2$a_1, dat2$c_1, dat2$Total_1e, dat2$Risk_1e, dat2$RD_1, 
                   dat2$RR_1, round(dat2$OR_1, 2), dat2$b_1, dat2$d_1, dat2$Total_1ue, dat2$Risk_1ue, 
                   "ref", "ref", "ref", dat2$Total_1d, dat2$Total_1ud, dat2$Total_1, dat2$Risk1,
                   "", "", ""), 
                 nrow = 7, dimnames = dimnames), quote=FALSE)
    
    print("Stratify = 0")
    print(matrix(c(dat2$a_0, dat2$c_0, dat2$Total_0e, dat2$Risk_0e, dat2$RD_0, 
                   dat2$RR_0, round(dat2$OR_0, 2), dat2$b_0, dat2$d_0, dat2$Total_0ue, dat2$Risk_0ue, 
                   "ref", "ref", "ref", dat2$Total_0d, dat2$Total_0ud, dat2$Total_0,dat2$Risk0,
                  "", "", ""), 
                 nrow = 7, dimnames = dimnames), quote=FALSE)
    } else {
      dimnames = list(c("D+", "D-", "Total"), c("E+", "E-"))
    
    print("Overall")
    print(matrix(c(dat2$a, dat2$c, dat2$Total_e, dat2$b, dat2$d, dat2$Total_ue), nrow = 3, dimnames = dimnames), quote=FALSE)
    
    print("Stratify = 1")
    print(matrix(c(dat2$a_1, dat2$c_1, dat2$Total_1e, dat2$b_1, dat2$d_1, dat2$Total_1ue), nrow = 3, dimnames = dimnames), quote=FALSE)
    
    print("Stratify = 0")
    print(matrix(c(dat2$a_0, dat2$c_0, dat2$Total_0e, dat2$b_0, dat2$d_0, dat2$Total_0ue), nrow = 3, dimnames = dimnames), quote=FALSE)
  }
}} 


calc.OR <- function(info, exposure = "E", out = "D"){
  # Calculates the crude odds ratio for a two by two table 
  info = data.frame(E = info[[exposure]], D = info[[out]], 
                   Count = rep(1, length(info[[out]])))
   info$Count <- rep(1, nrow(info))
  a <- sum(info[info$E == 1 & info$D == 1, ]$Count)
  b <- sum(info[info$E == 0 & info$D == 1, ]$Count)
  c <- sum(info[info$E == 1 & info$D == 0, ]$Count)
  d <- sum(info[info$E == 0 & info$D == 0, ]$Count)
  OR <- (a*d)/(b*c)
  return(OR)
  
}
