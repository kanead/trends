#' Bayesian regression on scavenger survey data
rm(list = ls())
graphics.off()

#' load packages
library(tidyverse)
library(brms)
library(readxl)
library(tidybayes)
library(interactions)

#' load data
mydata <-
  read_xlsx("data/SurveyDataNov2021_amended.xlsx", sheet = "DataAmendR")
mydata %>% print(n = Inf)

tail(mydata$HV)
mydata[is.na(mydata)] <- 0
tail(mydata$HV)

#' observations per park?
mydata %>% group_by(NP) %>% count()

#' change the date to a continuous variable
mydata$DateFinal <-
  as.POSIXct(mydata$DateFinal, format = "%Y-%m-%d")

mydata <- transform(
  mydata,
  ndate = as.numeric(DateFinal),
  nyear  = as.numeric(format(DateFinal, '%Y')),
  nmonth = as.numeric(format(DateFinal, '%m')),
  doy    = as.numeric(format(DateFinal, '%j'))
)


#' make the park, carcass and transect factors
mydata$NP <- as.factor(mydata$NP)
levels(mydata$NP)
mydata$StandardTransect  <- as.factor(mydata$StandardTransect)
levels(mydata$StandardTransect)
mydata$CarcassPres <- as.factor(mydata$CarcassPres)
mydata$nmonth <- as.factor(mydata$nmonth)

#' arrange by date
mydata <- arrange(mydata, DateFinal)

#' convert the dates into months such that the start month is 0 in 2013
#' and the end month is 97 in 2021
#' first write the function
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}
#' then apply it 
mydata$ndate <- c(0, elapsed_months(mydata[-1, 1], mydata[1, 1]))
tail(mydata)

#' take a look at the data by comparing one species across the different 
#' types of time format 
ggplot(mydata, aes(x = ndate, y = AWB)) + geom_point(aes(colour = CarcassPres)) +
  theme_bw() + geom_smooth(method = "lm") + facet_wrap( ~ NP)

ggplot(mydata, aes(x = ndate, y = AWB)) + geom_point(aes(colour = CarcassPres)) +
  theme_bw() + geom_smooth(method = "lm") + facet_wrap( ~ NP)

#' select just what we need
mydata <-
  mydata %>% select(
    WHV,
    TE,
    LFV,
    HV,
    Bateleur,
    AWB,
    ndate,
    nyear,
    NP,
    Season,
    CarcassPres,
    StandardTransect,
    Tlength
  )

####' Fit models with more informative priors ----

newprior <-
  c(
    set_prior("normal(0,.8)", class = "b", coef = "Intercept"),
    set_prior("normal(0,.8)", class = "b", coef = "NPRuaha"),
    set_prior("normal(0,.8)", class = "b", coef = "NPSelous")#,
    # set_prior("normal(0,.05)", class = "b", coef = "ndate"),
    # set_prior("normal(0,.05)", class = "b", coef = "ndate:NPRuaha"),
    # set_prior("normal(0,.05)", class = "b", coef = "ndate:NPSelous")
  )

####' African white-backed vulture Population Trends ----

nbglm_awb_default <-
  brm(
    AWB ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' plotting conditional effects
plotsAWB <- plot(conditional_effects(nbglm_awb_default, prob = 0),
                 ask = FALSE,
                 points = T,
                 offset = T)

plotsAWB[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("African white-backed vulture") + 
  theme_classic()

#' extract the random effects
ranef(nbglm_awb_default)

#' extract the fixed effects
fixef(nbglm_awb_default) %>% round(3)

#' bayesian R^2
bayes_R2(nbglm_awb_default) %>% round(3)

#' export them
awbv_results_fixed <- fixef(nbglm_awb_default) %>% round(3)
awbv_results_fixed <- as_tibble(awbv_results_fixed,rownames=NA) %>% 
  rownames_to_column()
awbv_results_fixed
class(awbv_results_fixed)
write.csv(file = "results/awbv_results1.csv", x = awbv_results_fixed, row.names = F)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_awb_default, ndraws = 100) # posterior predictive check
# 
# # plot posterior intervals
# mcmc_plot(nbglm_awb_default)

#' combine baseline posterior draws with differences
#' combine baseline posterior draws with differences
#' need to do this because the baseline of Katavi has to be added to the
#' other park coefficients for n_date because these are just the differences 
#' from the baseline rather than the absolute value in their own right 
draws_awb_def <- as_draws_df(nbglm_awb_default)

ggplot(draws_awb_def, aes(b_ndate)) + geom_histogram()
ggplot(draws_awb_def, aes(b_ndate + `b_ndate:NPRuaha`)) + geom_histogram()
ggplot(draws_awb_def, aes(b_ndate + `b_ndate:NPSelous`)) + geom_histogram()

#' median of the combined distributions
summary(draws_awb_def$b_ndate) %>% round(3)
quantile(draws_awb_def$b_ndate, c(.025, .975))
median(draws_awb_def$b_ndate) %>% round(3)
median(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPRuaha`) %>% round(3)
median(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPSelous`) %>% round(3)

median_change <- c(median(draws_awb_def$b_ndate),
                   median(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPRuaha`),
                   median(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_awb_def$b_ndate < 0) %>% round(3)
mean(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPRuaha` < 0) %>% round(3)
mean(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPSelous` < 0) %>% round(3)

neg_draw <- c(mean(draws_awb_def$b_ndate < 0),
              mean(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPRuaha` < 0),
              mean(draws_awb_def$b_ndate + draws_awb_def$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
awbv_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
awbv_results
rownames(awbv_results) <- c("Katavi", "Ruaha", "Selous")
colnames(awbv_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
awbv_results
write.csv("results/awbv_results2.csv",x = awbv_results, row.names = F)

#' export results using shinyapp
# launch_shinystan(nbglm_awb_default)
#' I saved the fixed effects for the sake of space in the Latex table 


####' Bateleur Population Trends ----

#' fit the model
nbglm_bat_default <-
  brm(
    Bateleur ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' check the models
# plot(nbglm_bat_default)


#' plotting conditional effects
plotsBat <- plot(conditional_effects(nbglm_bat_default, prob = 0),
                 ask = FALSE,
                 points = T,
                 offset = T)

plotsBat[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("Bateleur") + 
  theme_classic()

#' extract the random effects
ranef(nbglm_bat_default)

#' extract the fixed effects
fixef(nbglm_bat_default) %>% round(3)

#' bayesian R^2
bayes_R2(nbglm_bat_default) %>% round(3)

#' export them
bat_results_fixed <- fixef(nbglm_bat_default) %>% round(3)
bat_results_fixed <- as_tibble(bat_results_fixed,rownames=NA) %>% 
  rownames_to_column()
bat_results_fixed
class(bat_results_fixed)
write.csv(file = "results/bat_results1.csv", x = bat_results_fixed, row.names = F)

#' model fit
bayes_R2(nbglm_bat_default) %>% round(3)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_bat_default, ndraws = 100) # posterior predictive check

# plot posterior intervals
# mcmc_plot(nbglm_bat_default)

#' combine baseline posterior draws with differences
draws_bat <- as_draws_df(nbglm_bat_default)

ggplot(draws_bat, aes(b_ndate)) + geom_histogram()
ggplot(draws_bat, aes(b_ndate + `b_ndate:NPRuaha`)) + geom_histogram()
ggplot(draws_bat, aes(b_ndate + `b_ndate:NPSelous`)) + geom_histogram()

#' median of the combined distributions
summary(draws_bat$b_ndate) %>% round(3)
quantile(draws_bat$b_ndate, c(.025, .975))
median(draws_bat$b_ndate) %>% round(3)
median(draws_bat$b_ndate + draws_bat$`b_ndate:NPRuaha`) %>% round(3)
median(draws_bat$b_ndate + draws_bat$`b_ndate:NPSelous`) %>% round(3)

median_change <- c(median(draws_bat$b_ndate),
                   median(draws_bat$b_ndate + draws_bat$`b_ndate:NPRuaha`),
                   median(draws_bat$b_ndate + draws_bat$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_bat$b_ndate < 0) %>% round(3)
mean(draws_bat$b_ndate + draws_bat$`b_ndate:NPRuaha` < 0) %>% round(3)
mean(draws_bat$b_ndate + draws_bat$`b_ndate:NPSelous` < 0) %>% round(3)

neg_draw <- c(mean(draws_bat$b_ndate < 0),
              mean(draws_bat$b_ndate + draws_bat$`b_ndate:NPRuaha` < 0),
              mean(draws_bat$b_ndate + draws_bat$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
bat_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
bat_results
rownames(bat_results) <- c("Katavi", "Ruaha", "Selous")
colnames(bat_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
bat_results
write.csv("results/bat_results2.csv",x = bat_results, row.names = F)

####' Hooded vulture population trends -----
nbglm_hv_default <-
  brm(
    HV ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' check the models
# plot(nbglm_hv_default)

#' plotting conditional effects
plotsHV <- plot(conditional_effects(nbglm_hv_default, prob = 0),
                ask = FALSE,
                points = T,
                offset = T)

plotsHV[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("Hooded vulture") + 
  theme_classic()

#' extract the random effects
ranef(nbglm_hv_default)

#' extract the fixed effects
fixef(nbglm_hv_default) %>% round(3)

#' export them
hv_results_fixed <- fixef(nbglm_hv_default) %>% round(3)
hv_results_fixed <- as_tibble(hv_results_fixed,rownames=NA) %>% 
  rownames_to_column()
hv_results_fixed
class(hv_results_fixed)
write.csv(file = "results/hv_results1.csv", x = hv_results_fixed, row.names = F)

#' model fit
bayes_R2(nbglm_hv_default) %>% round(3)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_hv_default, ndraws = 100) # posterior predictive check

# plot posterior intervals
# mcmc_plot(nbglm_hv_default)

#' combine baseline posterior draws with differences
draws_hv <- as_draws_df(nbglm_hv_default)

hist(draws_hv$b_ndate)
hist(draws_hv$b_ndate + draws_hv$`b_ndate:NPRuaha`)
hist(draws_hv$b_ndate + draws_hv$`b_ndate:NPSelous`)

#' median of the combined distributions
summary(draws_hv$b_ndate) %>% round(3)
quantile(draws_hv$b_ndate, c(.025, .975))
median(draws_hv$b_ndate) %>% round(3)
median(draws_hv$b_ndate + draws_hv$`b_ndate:NPRuaha`) %>% round(3)
median(draws_hv$b_ndate + draws_hv$`b_ndate:NPSelous`) %>% round(3)


median_change <- c(median(draws_hv$b_ndate),
                   median(draws_hv$b_ndate + draws_hv$`b_ndate:NPRuaha`),
                   median(draws_hv$b_ndate + draws_hv$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_hv$b_ndate < 0) %>% round(3)
mean(draws_hv$b_ndate + draws_hv$`b_ndate:NPRuaha` < 0) %>% round(3)
mean(draws_hv$b_ndate + draws_hv$`b_ndate:NPSelous` < 0) %>% round(3)

neg_draw <- c(mean(draws_hv$b_ndate < 0),
              mean(draws_hv$b_ndate + draws_hv$`b_ndate:NPRuaha` < 0),
              mean(draws_hv$b_ndate + draws_hv$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
hv_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
hv_results
rownames(hv_results) <- c("Katavi", "Ruaha", "Selous")
colnames(hv_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
hv_results
write.csv("results/hv_results2.csv",x = hv_results, row.names = F)


####' Lappet faced vulture population trends -----
nbglm_lfv_default <-
  brm(
    LFV ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' check the models
# plot(nbglm_lfv_default)

#' plotting conditional effects
plotsLF <- plot(conditional_effects(nbglm_lfv_default, prob = 0),
                ask = FALSE,
                points = T,
                offset = T)

plotsLF[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("Lappet-faced vulture") + 
  theme_classic()

#' extract the random effects
ranef(nbglm_lfv_default) %>% round(3)

#' extract the fixed effects
fixef(nbglm_lfv_default) %>% round(3)

# bayes R2
bayes_R2(nbglm_lfv_default) %>% round(3)

#' export them
lfv_results_fixed <- fixef(nbglm_lfv_default) %>% round(3)
lfv_results_fixed <- as_tibble(lfv_results_fixed,rownames=NA) %>% 
  rownames_to_column()
lfv_results_fixed
class(lfv_results_fixed)
write.csv(file = "results/lfv_results1.csv", x = lfv_results_fixed, row.names = F)

#' model fit
bayes_R2(nbglm_lfv_default) %>% round(3)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_lfv_default, ndraws = 100) # posterior predictive check

# plot posterior intervals
# mcmc_plot(nbglm_lfv_default)

#' combine baseline posterior draws with differences
draws_lfv <- as_draws_df(nbglm_lfv_default)

hist(draws_lfv$b_ndate)
hist(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPRuaha`)
hist(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPSelous`)

#' median of the combined distributions
median(draws_lfv$b_ndate) %>% round(3)
median(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPRuaha`) %>% round(3)
median(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPSelous`) %>% round(3)

median_change <- c(median(draws_lfv$b_ndate),
                   median(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPRuaha`),
                   median(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_lfv$b_ndate < 0)
mean(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPRuaha` < 0)
mean(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPSelous` < 0)

neg_draw <- c(mean(draws_lfv$b_ndate < 0),
              mean(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPRuaha` < 0),
              mean(draws_lfv$b_ndate + draws_lfv$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
lfv_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
lfv_results
rownames(lfv_results) <- c("Katavi", "Ruaha", "Selous")
colnames(lfv_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
lfv_results
write.csv("results/lfv_results2.csv",x = lfv_results, row.names = F)

####' Tawny Eagle population trends -----

nbglm_te_default <-
  brm(
    TE ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' check the models
# plot(nbglm_te_default)

#' plotting conditional effects
plotsTawny <- plot(conditional_effects(nbglm_te_default, prob = 0),
                   ask = FALSE,
                   points = T,
                   offset = T)

plotsTawny[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("Tawny eagle") + 
  theme_classic()

#' extract the random effects
ranef(nbglm_te_default) %>% round(3)

#' extract the fixed effects
fixef(nbglm_te_default) %>% round(3)

#' bayes R2
bayes_R2(nbglm_te_default) %>% round(3)

#' export them
te_results_fixed <- fixef(nbglm_te_default) %>% round(3)
te_results_fixed <- as_tibble(te_results_fixed,rownames=NA) %>% 
  rownames_to_column()
te_results_fixed
class(te_results_fixed)
write.csv(file = "results/te_results1.csv", x = te_results_fixed, row.names = F)

#' model fit
bayes_R2(nbglm_te_default) %>% round(3)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_te_default, ndraws = 100) # posterior predictive check

# plot posterior intervals
# mcmc_plot(nbglm_te_default)

#' combine baseline posterior draws with differences
draws_te <- as_draws_df(nbglm_te_default)

hist(draws_te$b_ndate)
hist(draws_te$b_ndate + draws_te$`b_ndate:NPRuaha`)
hist(draws_te$b_ndate + draws_te$`b_ndate:NPSelous`)

#' median of the combined distributions
median(draws_te$b_ndate) %>% round(3)
median(draws_te$b_ndate + draws_te$`b_ndate:NPRuaha`) %>% round(3)
median(draws_te$b_ndate + draws_te$`b_ndate:NPSelous`) %>% round(3)

median_change <- c(median(draws_te$b_ndate),
                   median(draws_te$b_ndate + draws_te$`b_ndate:NPRuaha`),
                   median(draws_te$b_ndate + draws_te$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_te$b_ndate < 0)
mean(draws_te$b_ndate + draws_te$`b_ndate:NPRuaha` < 0)
mean(draws_te$b_ndate + draws_te$`b_ndate:NPSelous` < 0)

neg_draw <- c(mean(draws_te$b_ndate < 0),
              mean(draws_te$b_ndate + draws_te$`b_ndate:NPRuaha` < 0),
              mean(draws_te$b_ndate + draws_te$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
te_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
te_results
rownames(te_results) <- c("Katavi", "Ruaha", "Selous")
colnames(te_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
te_results
write.csv("results/te_results2.csv",x = te_results, row.names = F)


####' White-headed Vulture population trends -----

#' fit the model
nbglm_wh_default <-
  brm(
    WHV ~  0 + Intercept + # this allows control of prior on the intercept
      ndate * NP + Season + CarcassPres +
      (1 | NP / StandardTransect) +
      offset(log(Tlength)),
    family = negbinomial,
    data = mydata,
    chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    prior = newprior
  )

#' check the models
# plot(nbglm_wh_default)

#' plotting conditional effects
plotsWH <- plot(conditional_effects(nbglm_wh_default, prob = 0),
                ask = FALSE,
                points = T,
                offset = T)

plotsWH[[5]] +xlab("time (months)") + ylab("counts") + 
  ggtitle("White headed vulture") + 
  theme_classic()


#' extract the random effects
ranef(nbglm_wh_default) %>% round(3)

#' extract the fixed effects
fixef(nbglm_wh_default) %>% round(3)

#' bayes R2
bayes_R2(nbglm_wh_default) %>% round(3)

#' export them
wh_results_fixed <- fixef(nbglm_wh_default) %>% round(3)
wh_results_fixed <- as_tibble(wh_results_fixed,rownames=NA) %>% 
  rownames_to_column()
wh_results_fixed
class(wh_results_fixed)
write.csv(file = "results/wh_results1.csv", x = wh_results_fixed, row.names = F)


#' model fit
bayes_R2(nbglm_wh_default) %>% round(3)

# check whether the distributional assumption of the model is reasonable
# pp_check(nbglm_wh_default, ndraws = 100) # posterior predictive check

# plot posterior intervals
# mcmc_plot(nbglm_wh_default)

#' combine baseline posterior draws with differences
draws_wh <- as_draws_df(nbglm_wh_default)

hist(draws_wh$b_ndate)
hist(draws_wh$b_ndate + draws_wh$`b_ndate:NPRuaha`)
hist(draws_wh$b_ndate + draws_wh$`b_ndate:NPSelous`)

#' median of the combined distributions
median(draws_wh$b_ndate) %>% round(3)
median(draws_wh$b_ndate + draws_wh$`b_ndate:NPRuaha`) %>% round(3)
median(draws_wh$b_ndate + draws_wh$`b_ndate:NPSelous`) %>% round(3)

median_change <- c(median(draws_wh$b_ndate),
                   median(draws_wh$b_ndate + draws_wh$`b_ndate:NPRuaha`),
                   median(draws_wh$b_ndate + draws_wh$`b_ndate:NPSelous`)
) ; median_change %>% round(3)

#' convert to % change
perc_change <- (exp(median_change)-1)*100
#' convert to 12 month % change
yearly_change <- (exp(median_change*12)-1)*100

#' proportion of posterior that's negative
mean(draws_wh$b_ndate < 0)
mean(draws_wh$b_ndate + draws_wh$`b_ndate:NPRuaha` < 0)
mean(draws_wh$b_ndate + draws_wh$`b_ndate:NPSelous` < 0)

neg_draw <- c(mean(draws_wh$b_ndate < 0),
              mean(draws_wh$b_ndate + draws_wh$`b_ndate:NPRuaha` < 0),
              mean(draws_wh$b_ndate + draws_wh$`b_ndate:NPSelous` < 0) 
); neg_draw %>% round(3)

#' collate these derived results and export them 
wh_results <- matrix(c(median_change,perc_change,yearly_change,neg_draw),ncol=4) %>% round(3)
wh_results
rownames(wh_results) <- c("Katavi", "Ruaha", "Selous")
colnames(wh_results) <- c("estimate (combined)", "monthly change (%)", "yearly change (%)", "Prop posterior negative")
wh_results
write.csv("results/wh_results2.csv",x = wh_results, row.names = F)

####' compile the results ----
#' first for the summary of the regressions by species
data_path <- "results/Priors3/"   # path to the data

files1 <- dir(data_path, pattern = "*1") # get file names
length(files1)

results_data1 <- files1  %>% 
  # read in all the files, appending the path before the file name
  map( ~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind)
results_data1

write.csv(x = results_data1, file = "regression_summary.csv", row.names = F)

#' then for the proportion of the parameters on the dates that are negative 
files2 <- dir(data_path, pattern = "*2") # get file names
length(files2)

results_data2 <- files2  %>% 
  # read in all the files, appending the path before the file name
  map( ~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind)
results_data2

write.csv(x = results_data2, file = "prop_negative_summary.csv", row.names = F)

