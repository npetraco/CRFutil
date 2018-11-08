library(rstan)

# Simulate some fake data
a.tru   <- 3.2
b1.tru  <- 0.76
b2.tru  <- 2.34
sig.tru <- 0.12

alps  <- rnorm(100, mean = a.tru,  sd = sig.tru)
bets1 <- rnorm(100, mean = b1.tru, sd = sig.tru)
bets2 <- rnorm(100, mean = b2.tru, sd = sig.tru)

x1  <- sample(c(-1,1), 100, replace = T)
x2  <- sample(c(-1,1), 100, replace = T)
mus <- alps + bets1*x1 + bets2*x2
y   <- sapply(1:100, function(xx){rnorm(1, mean = mus[xx], sd = sig.tru)})

plot(x1,y)
plot(x2,y)
# Standard MLE linear regression:
mle.fit <- lm(y ~ x1 + x2)
summary(mle.fit)
abline(mle.fit)


x1.new   <- sample(c(-1,1), 100, replace = T)
x2.new   <- sample(c(-1,1), 100, replace = T)
dat <- list(N    = 100,
            x1    = x1,
            x2    = x2,
            y    = y,
            Nnew = 100,
            x1new = x1.new,
            x2new = x2.new
)


# Compile and run the model:
stan.c <-"
data {
  int<lower=0> N;
  vector[N]    x1;
  vector[N]    x2;
  vector[N]    y;
  int<lower=0> Nnew;
  vector[Nnew] x1new;
  vector[Nnew] x2new;
}
parameters {
  real alpha;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}
model {
  alpha ~ cauchy(0,1);
  beta1 ~ cauchy(0,5);
  beta2 ~ cauchy(0,5);
  sigma ~ normal(0,1);
  y ~ normal(alpha + beta1 * x1 + beta2 * x2, sigma);
}
generated quantities {
  vector[Nnew] mu_pred;
  vector[Nnew] y_tilde;

  mu_pred = alpha + beta1 * x1new + beta2 * x2new;

  for (i in 1:Nnew){
    y_tilde[i] = normal_rng(mu_pred[i], sigma);
  }
}"

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# Compile and run the model
model.c <- stanc(model_code = stan.c, model_name = 'model1')
sm <- stan_model(stanc_ret = model.c, verbose = T)
#bayes.fit <- sampling(sm, data = dat, iter=10000, thin = 1, chains = 4)
bayes.fit <- sampling(sm, data = dat, iter=5000, thin = 1, chains = 4, control = list(adapt_delta = 0.9))
bayes.fit

# Examine the sampling output in more detail:
alpha <- extract(bayes.fit,"alpha")[[1]]
plot(alpha, typ="l", main="Trace")
hist(alpha, bre=80, probability = T)
acf(alpha)
mean(alpha)
median(alpha)
a.tru

beta <- extract(bayes.fit,"beta")[[1]]
plot(beta, typ="l", main="Trace")
hist(beta, bre=80, probability = T)
acf(beta)
mean(beta)
median(beta)
b1.tru

sigma <- extract(bayes.fit,"sigma")[[1]]
plot(sigma, typ="l", main="Trace")
hist(sigma, bre=80, probability = T)
acf(sigma)
mean(sigma)
median(sigma)
sig.tru

