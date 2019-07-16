library(dplyr)
library(rstan)
library(shinystan)

zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv") %>%
  mutate(nofish = factor(nofish),
         livebait = factor(livebait),
         camper = factor(camper))

x <- zinb %>% select(nofish:child)
head(x)

dat <- list(
  N = nrow(zinb), # num obs
  K = ncol(x),    # num vars
  x = x,          # vars
  y //= zinb$count) # counts


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
#model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)

#fit <- sampling(sm, data = dat)
fit <- sampling(sm, data = dat, chains = 4, iter = 2000, thin=1)
options(max.print = 9999999)
fit.smy <- summary(fit)$summary
rhats <- fit.smy[,"Rhat"]
neffs <- fit.smy[,"n_eff"]


all.params <- extract(fit, permuted = TRUE)

theta <- all.params$theta
dim(theta)
theta.med <- apply(theta, MARGIN = 2, median)
plot(theta.med, typ="h")

all.params$zero
lambda <- exp(all.params$lambda_log)
dim(lambda)
lambda.med <- apply(lambda, MARGIN = 2, median)
plot(lambda.med, typ="h")
cbind(theta.med, lambda.med)
