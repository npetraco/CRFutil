library(rstan)

set.seed(689934)

N <- 5000
x <- rnorm(N, 10, 1)
X = t(data.matrix(data.frame(x, x * x)))

M <- 2
beta = matrix(c(2.5, -1), nrow=M, ncol=1)
alpha <- -0.275
sigma <- 0.8

mu <- t(X) %*% beta + alpha
y = sapply(1:N, function(n) rnorm(1, mu[n], sigma))
y

stan_rdump(c("N", "M", "X", "y"), file="regr.data.R")
setwd(paste0("/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/","tests"))

writeLines(readLines("regr.stan"))
input_data <- read_rdump("regr.data.R")
input_data
class(input_data)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("stan_utility.R")


fit <- stan(file='regr.stan', data=input_data, seed=483892929)
fit
check_treedepth(fit)
check_treedepth(fit, 1)
print(fit)
check_energy(fit)
check_divergences(fit)
check_div(fit)

sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
stepsizes <- sapply(sampler_params, function(x) x[1,'stepsize__'])
names(stepsizes) <- list("Chain 1", "Chain 2", "Chain 3" ,"Chain 4")
stepsizes

n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients

partition <- partition_div(fit)
params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(params$'beta[1]', params$'beta[2]',
     col=c_dark_trans, pch=16, cex=0.8, xlab="beta[1]", ylab="beta[2]",
     xlim=c(1.5, 3), ylim=c(-1.1, -0.9))
points(beta[1,1], beta[2,1],
       col=c_mid, pch=17, cex=2)

fit@sim$samples
library(shinystan)
launch_shinystan(fit)

?check_treedepth

check_hmc_diagnostics(fit)


setwd(paste0("/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/","tests"))
qr_fit <- stan(file='qr_regr.stan', data=input_data, seed=483892929)
check_treedepth(qr_fit)
check_energy(qr_fit)
check_div(qr_fit)

sampler_params <- get_sampler_params(qr_fit, inc_warmup=FALSE)
qr_stepsizes <- sapply(sampler_params, function(x) x[1,'stepsize__'])
names(qr_stepsizes) <- list("Chain 1", "Chain 2", "Chain 3" ,"Chain 4")
qr_stepsizes

n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients

partition <- partition_div(qr_fit)
params <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(params$'beta_tilde[1]', params$'beta_tilde[2]',
     col=c_dark_trans, pch=16, cex=0.8, xlab="beta_tilde[1]", ylab="beta_tilde[2]")


