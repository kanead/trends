######################################################
# Vulture Vulture Survival Analysis
######################################################

##clean##
rm(list=ls())
################
##load libraries##
################

library(survival)
library(survminer)

#load and prep data----
d<-read.csv("E:/Vulture/data/vulture4.csv")


#make end year a factor
d$year<-factor(d$dyear)


#relevel location
d$Location<-factor(d$Mort.Loc)
d$relLoc <-relevel(d$Location, ref = "Ruaha")



#perform kapaln-meier ----

km_min <- survfit(Surv(Time, min.censor) ~ 1, data = d)
summary(km_min)
print(km_min)
summary(km_min, times = c(1,30,60,180*(1:16)))
summary(km_min, times = c(1,30*(1:48)))


km_max <- survfit(Surv(Time, max.censor) ~ 1, data = d)
summary(km_max)
print(km_max)
summary(km_max, times = c(1,30,60,180*(1:16)))
summary(km_max, times = c(1,365.25*(1:4)))



#get KM survival table

summary(survfit(Surv(Time, max.censor) ~1, data = d), times = 365.25)
summary(survfit(Surv(Time, min.censor) ~ 1, data = d), times = 365.25)


#get median survival

survfit(Surv(Time, max.censor) ~ 1, data = d)
survfit(Surv(Time, min.censor) ~ 1, data = d)


#cox proportional hazards modeling----

#remove the bird that flew to Botswana
d1<-filter(d, !Location=="Other")

#test location as a fixed term

cox.max <- coxph(Surv(Time, max.censor) ~ relLoc, data=d1)
summary(cox.max)

cox.min <- coxph(Surv(Time, min.censor) ~ relLoc, data=d1)
summary(cox.min)

#test year as a fixed term

cox.max <- coxph(Surv(Time, max.censor) ~ year, data=d)
summary(cox.max)

cox.min <- coxph(Surv(Time, min.censor) ~ year, data=d)
summary(cox.min)

#get probabilities per year
surv_adjustedcurves(cox.max, data = d1, method = "average", variable = "year")
surv_adjustedcurves(cox.min, data = d1, method = "average", variable = "year")

###survival probabilities across locations using KM-----


#create datasets for each location
ruaha<-filter(d, Mort.Loc=="Ruaha")
k<-filter(d, Mort.Loc=="Katavi")
s<-filter(d, Mort.Loc=="Selous")



#survival to 1 year
summary(survfit(Surv(Time, min.censor) ~ 1, data = ruaha), times = 365.25)
summary(survfit(Surv(Time, max.censor) ~ 1, data = ruaha), times = 365.25)


summary(survfit(Surv(Time, min.censor) ~ 1, data = s), times = 365.25)
summary(survfit(Surv(Time, max.censor) ~ 1, data = s), times = 365.25)


summary(survfit(Surv(Time, min.censor) ~ 1, data = k), times = 365.25)
summary(survfit(Surv(Time, max.censor) ~ 1, data = k), times = 365.25)




#survival to 2 years
summary(survfit(Surv(Time, min.censor) ~ 1, data = ruaha), times = 730)
summary(survfit(Surv(Time, max.censor) ~ 1, data = ruaha), times =730)


summary(survfit(Surv(Time, min.censor) ~ 1, data = s), times = 730)
summary(survfit(Surv(Time, max.censor) ~ 1, data = s), times = 730)


summary(survfit(Surv(Time, min.censor) ~ 1, data = k), times = 730)
summary(survfit(Surv(Time, max.censor) ~ 1, data = k), times = 730)



##perform log-rank test using time event data set for calendar years----

dt<-read.csv("E:/Vulture/data/time event.csv")

dt$min.surv <- Surv(dt$Time, dt$min.censor)
dt$max.surv <- Surv(dt$Time, dt$max.censor)

survdiff(formula=min.surv ~ year, data = dt)
survdiff(formula=max.surv ~ year, data = dt)




