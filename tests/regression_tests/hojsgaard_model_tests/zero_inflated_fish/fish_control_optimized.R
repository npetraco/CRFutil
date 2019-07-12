library(dplyr)
library(rstan)
library(shinystan)

zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv") %>%
  mutate(nofish = factor(nofish),
         livebait = factor(livebait),
         camper = factor(camper))

x <- zinb %>% select(nofish:child)
head(x)


y         <- zinb$count
zero.idxs <- which(y==0)
y.zero    <- y[zero.idxs]
y.nonzero <- y[-zero.idxs]
x.z       <- x[zero.idxs,]
x.nz      <- x[-zero.idxs,]
N.zero    <- length(y.zero)
N.nonzero <- length(y.nonzero)

dat <- list(
  N_zero = N.zero,       # num obs
  N_nonzero = N.nonzero,
  K = ncol(x),           # num vars
  x_z = x.z,             # vars
  x_nz = x.nz,
  y_zero = y.zero,       # counts
  y_nonzero = y.nonzero
)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
#model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson_optimized3.stan"), model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)

#fit <- sampling(sm, data = dat)
fit <- sampling(sm, data = dat, chains = 8, iter = 20000, thin=10)
options(max.print = 9999999)
fit
# ~ 5.7 sec for 4000 samples

params <- extract(fit, permuted = TRUE)

theta.z <- params$theta_z
dim(theta.z)
theta.z.med <- apply(theta.z, MARGIN = 2, median)
plot(theta.z.med, typ="h")

theta.nz <- params$theta_nz
dim(theta.nz)
theta.nz.med <- apply(theta.nz, MARGIN = 2, median)
plot(theta.nz.med, typ="h")

