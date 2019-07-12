library(dplyr)
library(rstan)
library(shinystan)

zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv") %>% 
  mutate(nofish = factor(nofish),
         livebait = factor(livebait),
         camper = factor(camper))

x <- zinb %>% select(nofish:child)

dat <- list(
  N = nrow(zinb),
  K = ncol(x),
  x = x,
  y = zinb$count)


fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)

fit <- sampling(sm, data = dat)
fit

