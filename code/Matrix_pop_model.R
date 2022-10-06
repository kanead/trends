######################################################
# Vulture Trends Matrix Model
# Survival rates are for African white-backed vultures
# 6 scenarios exploring the worst and best case for the
# 3 parks, Katavi, Ruaha and Nyerere 
######################################################
# clean everything first
rm(list = ls())

# load required packages
library(popbio)

#### Matrix values ----------------------------------------
# these are equal across the six matrices 
f1 <- 0 # n females fledged per female aged 1 to 4
f2 <- 0.6 / 2 # n females fledged per female aged >=5
s0 <- 0.42 # survival prob. from fledging until the following year

#### Katavi Worse Case ----
s1 <- 0.613 # survival prob. from age 1-4
s2 <- 0.613 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_KW <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_KW) <- c("1","2","3","4","5+")
A_pre_KW

#### Katavi Best Case ----
s1 <- 0.834 # survival prob. from age 1-5
s2 <- 0.834 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_KB <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_KB) <- c("1","2","3","4","5+")
A_pre_KB


#### Ruaha worst case ---- 
s1 <- 0.894 # survival prob. from age 1-5
s2 <- 0.894 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_RW <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_RW) <- c("1","2","3","4","5+")
A_pre_RW

#### Ruaha best case ---- 
s1 <- 0.894 # survival prob. from age 1-5
s2 <- 0.894 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_RB <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_RB) <- c("1","2","3","4","5+")
A_pre_RB

#### Nyerere worst case ----
s1 <- 0.801 # survival prob. from age 1-5
s2 <- 0.801 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_NW <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_NW) <- c("1","2","3","4","5+")
A_pre_NW

#### Nyerere best case ----
s1 <- 0.853 # survival prob. from age 1-5
s2 <- 0.853 # survival prob. for older birds

# Pre-breeding pulse matrix model
A_pre_NB <-matrix(c(
  s0 * f1, s0 * f1, s0 * f1, s0 * f1, s0 * f2,
  s1, 0, 0, 0, 0,
  0, s1, 0, 0, 0, 
  0, 0, s1, 0, 0,  
  0, 0, 0, s1, s2),nrow=5,byrow=TRUE)

colnames(A_pre_NB) <- c("1","2","3","4","5+" )
A_pre_NB

#' combine the matrices
vul_mats <- list(A_pre_KW, A_pre_KB, A_pre_RW, A_pre_RB, A_pre_NW, A_pre_NB)
vul_mats

#' run the eigen analysis on each 
eigen_results <- lapply(vul_mats, eigen.analysis)

#' pull out the lambda values 
get_lambda <- function(x){
 lambdas <- eigen_results[[x]]$lambda1
}
runs <- c(1:6) # matrix indexed by number
lambda_values <- lapply(runs, get_lambda)
names(lambda_values) <- c("Katavi_Worst", "Katavi_Best", "Ruaha_Worst", 
                          "Ruaha_Best", "Nyerere_Worst", "Nyerere_Best")
data.frame(lambda_values)
